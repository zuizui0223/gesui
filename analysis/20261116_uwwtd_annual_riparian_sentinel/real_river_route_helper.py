#!/usr/bin/env python3
"""Generate ±1.5 km river-path locations around original 16 outfalls, source ONLY.

Uses the original signed river-chain ID/order and location, not the original
mid-reach U1/D1 GBIF site centres. Local metric projection per outfall;
actual HydroRIVERS geometries, directed reach adjacency and endpoint gaps
must support the path. NO GBIF / EIV / TP / butterfly outcomes read, and no
actual wastewater-delivered N/P source is inferred.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,math
from pathlib import Path
import geopandas as gpd
import numpy as np
from pyproj import CRS,Transformer
from shapely.geometry import LineString,MultiLineString,Point
from shapely.ops import linemerge,transform

STUDY_SHA="9cf5757b96864b5918d205c0d1a9bd8734c85605522aab3763b1ed457d1af335"
GEO_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
SOURCE_SCHEMA="GESUI_20261116_REAL_LOCAL_1500M_HYDROWASTE_HYDRORIVERS_SOURCE_ONLY_V1"
TARGET_METRES=1500.
MAX_JOIN_GAP_M=100.
MAX_OUTFALL_SNAP_M=100.
LOCAL_SITE_RADIUS_M=250.
MIN_UNIQUE_FAC=10
MIN_COUNTRIES=3

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def one(dir,name):
    files=list(dir.rglob(name))
    if len(files)!=1:raise RuntimeError("Original source unavailable/duplicated: "+name)
    return files[0]
def load_source(source):
    path=one(source,"candidate_site_design.json")
    if sha(path)!=STUDY_SHA:raise RuntimeError("Original outfall/chain source ID coordinate SHA changed")
    j=json.loads(path.read_text())
    sites=j["selected"]
    if len(sites)!=16 or len({str(s["wwtp_id"]) for s in sites})!=16:
        raise RuntimeError("Original source selected exact 16 pilot outfalls altered")
    if collections.Counter(s["country"] for s in sites)!={"AUT":4,"CZE":4,"DEU":4,"DNK":4}:
        raise RuntimeError("Original source-only pilot national selections altered")
    if collections.Counter(s["exposure_stratum"] for s in sites)!={"high":8,"low":8}:
        raise RuntimeError("Original source-only pilot high/low model fraction strata altered")
    for item in sites:
        if item["river_hydraulic_chain"]!="U2>U1>OUTFALL_REACH>D1":
            raise RuntimeError("Original outfall source did not have four linked river units")
    return sites

def line(g):
    if isinstance(g,LineString):ret=g
    elif isinstance(g,MultiLineString):
        ret=linemerge(g)
        if not isinstance(ret,LineString):
            raise ValueError("Real river reach MultiLineString is physically disconnected")
    else:raise ValueError("Original source river geometry is neither LineString nor connected Multiline")
    if ret.is_empty or len(ret.coords)<2:
        raise ValueError("Real source river reach path empty/degenerate")
    return ret

def load_original_rivers(path,ids):
    if sha(path)!=GEO_SHA:raise RuntimeError("Original public HydroRIVERS Europe GeoParquet SHA changed")
    d=gpd.read_parquet(path)
    if d.crs is None:raise RuntimeError("Original river network CRS missing")
    d=d.rename(columns={k:k.upper() for k in d.columns if k!=d.geometry.name})
    for k in ("HYRIV_ID","NEXT_DOWN"):
        if k not in d.columns:raise RuntimeError("Original converted GeoParquet missing original source topology "+k)
    subset=d[d["HYRIV_ID"].astype(int).isin(ids)].to_crs("EPSG:4326")
    if len(subset)!=len(ids):
        raise RuntimeError("Original selected river reach source ID missing/duplicated "+str((len(subset),len(ids))))
    out={}
    for r in subset.itertuples(index=False):
        orig=int(r.HYRIV_ID)
        if orig in out:raise RuntimeError("Original river id repeated")
        geom=getattr(r,d.geometry.name)
        out[orig]=(int(r.NEXT_DOWN),line(geom))
    if set(out)!=set(ids):raise RuntimeError("Original geometric reach source coverage not complete")
    return out

def ends(segment):
    a=list(segment.coords)
    return (Point(a[0]),Point(a[-1]))
def join_closest(a,b):
    endsa=ends(a);endsb=ends(b)
    cand=sorted((float(endsa[i].distance(endsb[j])),i,j)
       for i in (0,1) for j in (0,1))
    return cand[0],cand[1][0]

def oriented(g,reverse):
    return LineString(list(g.coords)[::-1]) if reverse else g

def project_metric(item,river):
    original_ids=[int(item[k]) for k in ("u2_river_id","u1_river_id","river_id","d1_river_id")]
    if len(set(original_ids))!=4:
        return None,{"stop":"source_chain_duplicate_reach"}
    for a,b in zip(original_ids[:-1],original_ids[1:]):
        if river[a][0]!=b:
            return None,{"stop":"real_Next_DOWN_graph_does_not_match_original_source_ids",
                         "invalid_pair":[a,b],"actual_nextdown":river[a][0]}
    lon=float(item["outfall_lon"]);lat=float(item["outfall_lat"])
    crs=CRS.from_proj4(f"+proj=aeqd +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m +no_defs")
    to=Transformer.from_crs("EPSG:4326",crs,always_xy=True)
    back=Transformer.from_crs(crs,"EPSG:4326",always_xy=True)
    mlines=[line(transform(to.transform,river[r][1])) for r in original_ids]
    joininfo=[join_closest(mlines[i],mlines[i+1]) for i in range(3)]
    gaps=[j[0][0] for j in joininfo]
    if any(x>MAX_JOIN_GAP_M for x in gaps):
        return None,{"stop":"source_endpoint_gap_exceeds_100m","junction_gaps_m":gaps}
    # For each ordered pair A -> B the closest endpoint of A must be A's
    # downstream END and the closest endpoint of B its upstream START.
    endidx=[j[0][1] for j in joininfo]
    startidx=[j[0][2] for j in joininfo]
    if startidx[0]==endidx[1] or startidx[1]==endidx[2]:
        return None,{"stop":"interior_river_reach_needs_same_end_for_both_directions",
                     "junction_endpoint_source_indices":{"starts":startidx,"ends":endidx}}
    rev=[endidx[0]==0,startidx[0]==1,startidx[1]==1,startidx[2]==1]
    segments=[oriented(g,flip) for g,flip in zip(mlines,rev)]
    if any(g.length<1 for g in segments):
        return None,{"stop":"source_river_segment_length_lt_1m"}
    point=Point(*to.transform(lon,lat))
    o=segments[2]
    source_snap=float(point.distance(o))
    if source_snap>MAX_OUTFALL_SNAP_M:
        return None,{"stop":"hydrowaste_source_outfall_snap_over_100m",
                     "snap_outfall_m":source_snap}
    routecoords=list(segments[0].coords)
    segment_starts=[0.]
    gap=[]
    for i in range(1,len(segments)):
        connector=float(Point(routecoords[-1]).distance(Point(segments[i].coords[0])))
        gap.append(connector)
        if connector>MAX_JOIN_GAP_M:
            return None,{"stop":"oriented_lines_have_large_source_endpoint_gap_m","oriented_gaps":gap}
        segment_starts.append(LineString(routecoords).length+connector)
        routecoords.extend(segments[i].coords)
    route=LineString(routecoords)
    q=float(o.project(point))
    outlet_pos=segment_starts[2]+q
    if outlet_pos<TARGET_METRES or route.length-outlet_pos<TARGET_METRES:
        return None,{"stop":"insufficient_source_contiguous_mainstem_1p5km_each_side",
                     "available_upstream_m":outlet_pos,
                     "available_downstream_m":route.length-outlet_pos}
    pu=route.interpolate(outlet_pos-TARGET_METRES)
    pd=route.interpolate(outlet_pos+TARGET_METRES)
    # Don't accept coordinates generated on short straight connectors
    # if their spatial position differs from the real hydrologic line.
    upoff=min(float(pu.distance(g)) for g in segments)
    dnoff=min(float(pd.distance(g)) for g in segments)
    if max(upoff,dnoff)>5:
        return None,{"stop":"source_target_falls_on_gap_not_actual_river",
                     "up_off_line_m":upoff,"down_off_line_m":dnoff}
    separation=float(pu.distance(pd))
    if separation<=2*LOCAL_SITE_RADIUS_M:
        return None,{"stop":"local_250m_radius_up_down_buffers_overlap",
                     "up_down_euclidean_m":separation}
    up_lon,up_lat=back.transform(pu.x,pu.y)
    down_lon,down_lat=back.transform(pd.x,pd.y)
    local={
       "original_outfall_id":str(item["wwtp_id"]),
       "country":item["country"],
       "original_source_outfall_lat":lat,"original_source_outfall_lon":lon,
       "original_source_outfall_reach":original_ids[2],
       "original_source_modelled_water_fraction":item.get("effluent_fraction_proxy",None),
       "original_source_facility_status":str(item["source_status"]),
       "actual_reported_2018_or_2022_NP_available_in_original_selected_panel":False,
       "new_UPSTREAM_1500m_lat":float(up_lat),
       "new_UPSTREAM_1500m_lon":float(up_lon),
       "new_DOWNSTREAM_1500m_lat":float(down_lat),
       "new_DOWNSTREAM_1500m_lon":float(down_lon),
       "source_outfall_to_reach_geometric_snap_m":source_snap,
       "source_U2_U1_OUTFALL_D1_endpoint_gaps_m":gaps,
       "source_original_river_true_link_and_orientation_checked":True,
       "original_up_down_river_path_distance_m":2*TARGET_METRES,
       "new_up_down_Euclidean_distance_m":separation,
       "candidate_future_GBIF_buffer_radius_m":LOCAL_SITE_RADIUS_M,
       "GBIF_or_biology_outcomes_examined":False,
       "not_wwtp_nutrient_concentration_or_ecological_effect":True}
    return local,{"stop":None,"outfall_snap_m":source_snap,
                  "up_down_euclidean_m":separation,
                  "oriented_source_endpoint_gap_max_m":max(gaps),
                  "source_route_length_m":route.length}

def main(args):
    args.output.mkdir(parents=True,exist_ok=True)
    selected=load_source(args.source)
    ids={int(q[k]) for q in selected for k in ("u2_river_id","u1_river_id","river_id","d1_river_id")}
    rivers=load_original_rivers(args.rivers,ids)
    accepted=[];failed=[];diagnostics=[]
    for q in selected:
        points,qa=project_metric(q,rivers)
        diagnostics.append({"outfall":str(q["wwtp_id"]),"country":q["country"],**qa})
        if points is None:failed.append({"outfall":str(q["wwtp_id"]),"country":q["country"],**qa})
        else:accepted.append(points)
    countries=sorted({q["country"] for q in accepted})
    passes=len(accepted)>=MIN_UNIQUE_FAC and len(countries)>=MIN_COUNTRIES
    result={"schema":SOURCE_SCHEMA,
       "status":"PASS_OUTFALL_1500M_GEOGRAPHY_ONLY" if passes else "STOP_PHYSICAL_LOCAL_SOURCE_GEOGRAPHY",
       "original_selected_16_units":len(selected),
       "original_16_selected_four_reach_unique_IDS":len(ids),
       "original_source_study_SHA256":STUDY_SHA,
       "original_HydroRIVERS_thirdparty_GeoParquet_SHA256":GEO_SHA,
       "source_true_chain_outfalls_geometry_pass":len(accepted),
       "source_true_chain_outfalls_geometry_fail":len(failed),
       "source_true_chain_passing_country_groups":countries,
       "source_outcome_blind_minimum_facilities":MIN_UNIQUE_FAC,
       "source_outcome_blind_minimum_countries":MIN_COUNTRIES,
       "source_true_outfall_1500m_up_down_both_250m_buffers":passes,
       "source_rejected_details":failed,
       "diagnostics_all_source_outfalls":diagnostics,
       "coordinate_selection_uses_true_directed_HYRIV_and_NEXT_DOWN":True,
       "old_U2_U1_D1_reach_centroid_GBIF_study_not_reused_as_local_sites":True,
       "thirdparty_HydroRIVERS_v1_conversion_NOT_official_byte_topology_verified":True,
       "original_actual_yearly_reported_wwtp_NP_at_receiver_not_qualified":True,
       "new_local_GBIF_and_field_biota_not_downloaded":True,
       "source_no_ecological_effect_estimated":True}
    path=args.output/"ORIGINAL16_REAL_WWTP_LOCAL_1500M_RIVER_SITE_GEOMETRY.json"
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    if accepted:
        # User-facing original public source coordinates only, no third-party
        # full downloaded river network geometry rehosting.
        import csv
        with (args.output/"ORIGINAL16_LOCAL_UPSTREAM_DOWNSTREAM_1500M_WWTP_SOURCE_SITES.csv").open(
             "w",encoding="utf-8",newline="") as f:
             w=csv.DictWriter(f,fieldnames=list(accepted[0]))
             w.writeheader();w.writerows(sorted(accepted,key=lambda a:a["original_outfall_id"]))
    print("GESUI_ORIGINAL16_REAL_LOCAL_1500M_WWTP_GEO="+json.dumps({
       "status":result["status"],
       "original_units":len(selected),
       "source_pass_n":len(accepted),
       "countries":countries,
       "source_fail_reasons":collections.Counter(q["stop"] for q in failed),
       "no_biology_or_pollutant_exposure_fitted":True},default=str,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--rivers",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
