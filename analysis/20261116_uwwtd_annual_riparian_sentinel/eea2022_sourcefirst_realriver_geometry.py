#!/usr/bin/env python3
"""True 2022 EEA measured NP outfall source-first local ±1.5km river geometry.

No NDVI, biology, aquatic EQR numerical response, GBIF flora or butterflies
accessed. Original EEA measured annual N/P source alone selects n1156; actual
HydroRIVERS next-down graph/coordinates define physical reaches and riparian
positions; original measured annual mass is not a local same-day chemical dose.

Source-only ecological study opportunity; DO NOT fit ecological treatment.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,math,sys
from pathlib import Path
import geopandas as gpd
import numpy as np
from shapely.geometry import box,Point
from shapely.ops import transform

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_uwwtd_annual_riparian_sentinel"))
from real_river_route_helper import project_metric,line,GEO_SHA

EEA_BLOB_SHA="1b4fbdfd8516e4040615a6cc1737febb1b6a33f0"
SCHEMA="GESUI_20261116_EEA2022_MEASURED_NP_ONE_FW_OUTFALL_DIRECTED_HYDRORIVERS_SOURCE_GEO_V1"
MAX_SNAP_M=100
MIN_SECOND_GAP_M=25
MIN_SEPARATION_BETWEEN_OUTFALLS_M=5000
PILOT_N_PER_COUNTRY=24
GEOMETRY_GATE_TOTAL=24
GEOMETRY_GATE_COUNTRY_MIN=8

def gitblob(b):
    return hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def distance_m(a,b):
    lat1,lon1=float(a[0]),float(a[1]);lat2,lon2=float(b[0]),float(b[1])
    p=math.pi/180;dlat=(lat2-lat1)*p;dlon=(lon2-lon1)*p
    k=math.sin(dlat/2)**2+math.cos(lat1*p)*math.cos(lat2*p)*math.sin(dlon/2)**2
    return 12742017.6*math.asin(min(1.,math.sqrt(k)))
def original_source(path):
    raw=path.read_bytes()
    if gitblob(raw)!=EEA_BLOB_SHA:
        raise RuntimeError("Original authentic PR23 2022 EEA UWWTD measured outfall manifest Gitblob changed")
    with path.open(newline="",encoding="utf-8-sig") as f:rows=list(csv.DictReader(f))
    if len(rows)!=1156 or collections.Counter(r["countryCode"] for r in rows)!={"DE":630,"AT":526}:
        raise RuntimeError("2022 exact reported measured both-positive source plant cohort changed")
    if len({r["original_EEA_plant_source_key"] for r in rows})!=1156:
        raise RuntimeError("Original source unique facility duplicate")
    out=[]
    for r in rows:
        if (r["source_one_original_FW_outfall"]!="True" or
            r["original_source_cycle"]!="EEA_2022_UWWTD_public"):
            raise RuntimeError("Original source single aquatic FW outlet label invalid")
        lat=float(r["discharge_lat_WGS84"]);lon=float(r["discharge_lon_WGS84"])
        n=float(r["original_reported_measured_N_t_year"])
        p=float(r["original_reported_measured_P_t_year"])
        if not (-90<=lat<=90 and -180<=lon<=180 and all(map(math.isfinite,[lat,lon,n,p])) and n>0 and p>0):
            raise RuntimeError("Source original annual positive measured mass and coordinate invalid")
        out.append({"reporter":r["original_EEA_reporter_code"],
          "id":r["original_EEA_plant_source_key"],
          "discharge_id":r["original_EEA_discharge_ID"],
          "country":r["countryCode"],"lat":lat,"lon":lon,"reportedN2022":n,"reportedP2022":p})
    return out

def geom_original(path):
    if sha(path)!=GEO_SHA:
        raise RuntimeError("Original source SHA of version-locked HydroRIVERS GeoParquet changed")
    gd=gpd.read_parquet(path)
    if gd.crs is None:
        raise RuntimeError("No original HydroRIVERS CRS")
    gd=gd.rename(columns={k:k.upper() for k in gd.columns if k!=gd.geometry.name})
    needed=("HYRIV_ID","NEXT_DOWN","UPLAND_SKM","DIS_AV_CMS")
    if any(k not in gd.columns for k in needed):
        raise RuntimeError("Original HydroRIVERS native drainage rank/graph source missing")
    gd=gd.to_crs("EPSG:4326")
    gd=gd[gd.geometry.notna() & ~gd.geometry.is_empty].reset_index(drop=True)
    ids=gd.HYRIV_ID.astype(np.int64).to_numpy()
    if len(ids)!=len(set(ids)):
        raise RuntimeError("Original HydroRIVERS cannot map unique river IDs")
    rid={int(v):i for i,v in enumerate(ids)}
    downstream={}
    up=collections.defaultdict(list)
    area={}
    riverq={}
    for i,v in enumerate(ids):
        nxt=int(gd.iloc[i].NEXT_DOWN) if math.isfinite(float(gd.iloc[i].NEXT_DOWN)) else -1
        downstream[int(v)]=nxt
        a=gd.iloc[i].UPLAND_SKM
        dis=gd.iloc[i].DIS_AV_CMS
        area[int(v)]=float(a) if math.isfinite(float(a)) else -1
        riverq[int(v)]=float(dis) if math.isfinite(float(dis)) else -1
        if nxt>0 and nxt in rid:up[nxt].append(int(v))
    return gd,rid,downstream,up,area,riverq

def snap_candidate(s,geo,rid,spatial):
    lat=s["lat"];lon=s["lon"]
    # Approximate metre conversion only for finding the nearest geometric source;
    # the final 1.5 km route uses precise locally centered WGS84 AEQD projection.
    dlat=0.006;dlon=0.006/max(.2,math.cos(math.radians(lat)))
    ix=spatial.query(box(lon-dlon,lat-dlat,lon+dlon,lat+dlat))
    if len(ix)==0:return None,{"reason":"no_original_river_reach_near_2022_outfall"}
    mx=111132.;my=111320*math.cos(math.radians(lat))
    p=Point(0,0)
    trials=[]
    for i in ix:
        g=geo.geometry.iloc[int(i)]
        def metric(x,y,z=None):
            xx=(np.asarray(x)-lon)*my
            yy=(np.asarray(y)-lat)*mx
            return (xx,yy) if z is None else (xx,yy,z)
        gmetric=transform(metric,g)
        dist=float(gmetric.distance(p))
        trials.append((dist,int(geo.HYRIV_ID.iloc[int(i)])))
    trials.sort()
    first=trials[0]
    if first[0]>MAX_SNAP_M:
        return None,{"reason":"original_2022_EEA_outlet_farther_than_100m_from_HYDRORIVERS",
                "nearest_river_path_offset_m":round(first[0],1)}
    if len(trials)>1 and trials[1][0]-first[0]<MIN_SECOND_GAP_M:
        return None,{"reason":"ambiguous_two_nearby_HydroRIVERS_reaches",
                "nearest_river_path_offset_m":round(first[0],1),
                "next_best_source_offset_m":round(trials[1][0],1)}
    return first[1],{"nearest_river_path_offset_m":round(first[0],3),
           "next_best_source_offset_m":round(trials[1][0],3) if len(trials)>1 else None,
           "source_candidates_within_geobbox":len(trials)}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    source=original_source(a.eea)
    gd,rid,downstream,up,area,riverq=geom_original(a.rivers)
    spatial=gd.sindex
    accepted=[];failed=[];snapstats=collections.Counter()
    for count,s in enumerate(source,1):
        idriver,detail=snap_candidate(s,gd,rid,spatial)
        if idriver is None:
            failed.append({"id":s["id"],"country":s["country"],**detail})
            snapstats[detail["reason"]]+=1
            continue
        nxt=downstream[idriver]
        if nxt not in rid:
            failed.append({"id":s["id"],"country":s["country"],
                "reason":"source_river_downstream_reach_not_in_original_graph"})
            snapstats["source_river_downstream_reach_not_in_original_graph"]+=1
            continue
        def best(q):
            options=up.get(q,[])
            return max(options,key=lambda x:(area[x],riverq[x],-x)) if options else None
        u1=best(idriver)
        u2=best(u1) if u1 is not None else None
        if u1 is None or u2 is None or len({u2,u1,idriver,nxt})!=4:
            failed.append({"id":s["id"],"country":s["country"],
                 "reason":"source_mainstem_two_upstream_or_one_downstream_reach_missing"})
            snapstats["source_mainstem_two_upstream_or_one_downstream_reach_missing"]+=1
            continue
        chain=[u2,u1,idriver,nxt]
        original={int(v):(downstream[v],line(gd.geometry.iloc[rid[v]])) for v in chain}
        item={"wwtp_id":s["id"],"country":"DEU" if s["country"]=="DE" else "AUT",
            "outfall_lat":s["lat"],"outfall_lon":s["lon"],
            "u2_river_id":u2,"u1_river_id":u1,"river_id":idriver,"d1_river_id":nxt,
            "source_status":"2022_EEA_reported_annual_measured_both_N_and_P",
            "effluent_fraction_proxy":None}
        coords,reason=project_metric(item,original)
        if coords is None:
            q=reason["stop"]
            failed.append({"id":s["id"],"country":s["country"],
                   "reason":q,**detail})
            snapstats[q]+=1
            continue
        coords.update({"original_EEA2022_plant_id":s["id"],
            "original_EEA2022_discharge_id":s["discharge_id"],
            "original_EEA2022_reporter":s["reporter"],
            "source_N_t_year_reported_measured_2022":s["reportedN2022"],
            "source_P_t_year_reported_measured_2022":s["reportedP2022"],
            "reported_source_year":2022,
            "reported_source_measurement_is_yearly_not_monthly":True,
            "source_GIS_original_river_snapped_reach":idriver,
            "source_nearest_original_reach_m":detail["nearest_river_path_offset_m"],
            "source_original_upstream_reach_ids":[u2,u1],
            "source_true_graph_downstream_reach_id":nxt,
            "EEA_original_countryCode":s["country"],
            "landscape_study_no_GBIF_or_Sentinel_NDVI_opened":True})
        accepted.append(coords)
        if count%200==0:print("EEA_2022_NP_OUTFALL_HYDRORIVERS_SOURCE_AUDIT_ROWS="+str(count),flush=True)
    if len(accepted)+len(failed)!=len(source):
        raise RuntimeError("Source original measured annual outfall audit did not exhaust original 1156")
    # Demarcate an independent pilot without selecting on satellite values or
    # signs of discharge mass; geographic cells get first opportunity.
    by_country=collections.defaultdict(list)
    for q in accepted:by_country[q["EEA_original_countryCode"]].append(q)
    def hashid(q):
        return hashlib.sha256(q["original_EEA2022_plant_id"].encode()).hexdigest()
    chosen=[];in_cell=set()
    for country in ("AT","DE"):
        studies=sorted(by_country[country],key=lambda q:(hashid(q),q["original_EEA2022_plant_id"]))
        sel=[]
        for q in studies:
            cell=(country,math.floor(float(q["original_source_outfall_lat"])/2),
                         math.floor(float(q["original_source_outfall_lon"])/2))
            if (cell not in in_cell and all(distance_m(
                (q["original_source_outfall_lat"],q["original_source_outfall_lon"]),
                (v["original_source_outfall_lat"],v["original_source_outfall_lon"])
                )>=MIN_SEPARATION_BETWEEN_OUTFALLS_M for v in sel)):
                in_cell.add(cell);sel.append(q)
                if len(sel)>=PILOT_N_PER_COUNTRY:break
        for q in studies:
            if len(sel)>=PILOT_N_PER_COUNTRY:break
            if q in sel:continue
            if all(distance_m(
                (q["original_source_outfall_lat"],q["original_source_outfall_lon"]),
                (v["original_source_outfall_lat"],v["original_source_outfall_lon"])
                )>=MIN_SEPARATION_BETWEEN_OUTFALLS_M for v in sel):
                sel.append(q)
        chosen.extend(sel)
    counts=collections.Counter(q["EEA_original_countryCode"] for q in chosen)
    passgeo=(len(chosen)>=GEOMETRY_GATE_TOTAL and counts["AT"]>=GEOMETRY_GATE_COUNTRY_MIN
             and counts["DE"]>=GEOMETRY_GATE_COUNTRY_MIN)
    manifest={
       "schema":SCHEMA,
       "source_status":"PASS_SOURCE_GEOMETRY_FOR_STAC_METADATA_ONLY" if passgeo else "STOP_SOURCE_GEOMETRY_UNSUPPORTED",
       "source_EEA2022_original_manifest_Gitblob":EEA_BLOB_SHA,
       "source_HydroRIVERS_thirdparty_converted_SHA256":GEO_SHA,
       "original_source_2022_measuredNP_oneFW_outlets":len(source),
       "original_country_source_outlets":dict(collections.Counter(q["country"] for q in source)),
       "original_source_real_gis_route_pass_n":len(accepted),
       "original_source_real_gis_route_pass_by_country":dict(collections.Counter(q["EEA_original_countryCode"] for q in accepted)),
       "source_nearest_line_or_route_failure_counts":dict(snapstats),
       "original_source_real_gis_route_fail_n":len(failed),
       "strict_nearest_river_source_snap_m":MAX_SNAP_M,
       "strict_reach_identification_second_nearest_gap_m":MIN_SECOND_GAP_M,
       "genuine_upstream_downstream_path_distance_each_m":1500,
       "both_local_buffers_250m_nonoverlapping":True,
       "geometry_gate_total_at_least":GEOMETRY_GATE_TOTAL,
       "geometry_gate_each_country_at_least":GEOMETRY_GATE_COUNTRY_MIN,
       "original_response_blind_selected_pilot_max_per_country":PILOT_N_PER_COUNTRY,
       "independent_pilot_selected_n":len(chosen),
       "independent_pilot_selected_country_n":dict(counts),
       "original_source_geographically_deduplicated_min_source_gap_m":MIN_SEPARATION_BETWEEN_OUTFALLS_M,
       "original_source_pilot_pass":passgeo,
       "observed_year2022_EEA_annual_mass_positive_N_and_P_proves_same_day_delivery":False,
       "actual_receiving_river_TNTP_and_leaf_N_NOT_observed":True,
       "source_WWTP_annual_measured_report_NOT_equivalent_Sentinel_NDVI_effect":True,
       "no_sentinel_reflectance_NDVI_or_biological_outcomes_read":True,
       "thirdparty_HydroRIVERS_conversion_NOT_official_byte_topology_verified":True}
    (a.output/"EEA2022_NP_SOURCEFIRST_LOCAL_DIRECTED_RIVER_GEOMETRY_GATE.json").write_text(
       json.dumps(manifest,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    # Full selected source derived locations, no rehosting original source full geometry.
    import csv
    if chosen:
        p=a.output/"EEA2022_SATELLITE_SOURCE_BLIND_INDEPENDENT_OUTFALL_PILOT.csv"
        with p.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(chosen[0]))
            w.writeheader();w.writerows(chosen)
        manifest["pilot_source_CSV_SHA256"]=sha(p)
        (a.output/"EEA2022_NP_SOURCEFIRST_LOCAL_DIRECTED_RIVER_GEOMETRY_GATE.json").write_text(
          json.dumps(manifest,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    # Audit the FULL unchanged source-qualified 282 receiving outfalls rather
    # than only the 48 geography-balanced satellite pilot. Source-only.
    # Keeping the full eligible denominator avoids falsely concluding
    # that no aquatic monitoring exists based on the pilot sample alone.
    if accepted:
        p_all=a.output/"EEA2022_ALL_SOURCE_QUALIFIED_UWWTD_NP_LOCAL_RIVER_SITES.csv"
        with p_all.open("w",encoding="utf-8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(accepted[0]))
            w.writeheader()
            w.writerows(sorted(accepted,key=lambda z:z["original_EEA2022_plant_id"]))
        manifest["original_all_true_receiving_outfall_geo_rows"]=len(accepted)
        manifest["original_all_source_valid_geography_CSV_SHA256"]=sha(p_all)
        (a.output/"EEA2022_NP_SOURCEFIRST_LOCAL_DIRECTED_RIVER_GEOMETRY_GATE.json").write_text(
            json.dumps(manifest,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_EEA2022_MEASURED_NP_LOCAL_RIVER_SOURCE_GATE="+json.dumps({
         "status":manifest["source_status"],"n_source_annual_measuredNP":len(source),
         "n_genuine_local_geometry_pass":len(accepted),
         "n_by_country":manifest["original_source_real_gis_route_pass_by_country"],
         "failure_counts":dict(snapstats),
         "n_pilot_geography_before_Sentinel":len(chosen),
         "pilot_countries":dict(counts),
         "no_satellite_or_biology_outcomes":True},ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--eea",required=True,type=Path)
    p.add_argument("--rivers",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
