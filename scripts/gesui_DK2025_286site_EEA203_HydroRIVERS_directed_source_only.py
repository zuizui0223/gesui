#!/usr/bin/env python3
"""ORIGINAL source-only Denmark EBMS 2025 WWTP -> river -> survey site feasibility.

Compares 203 original official 2022 EEA single-FW outfall coordinates and
286 original 2025 eBMS public Event-only source sites against the same SHA-
locked third-party HydroRIVERS Europe v1-derived network as prior gesui PR #9.

Reads ZERO butterfly individualCount/taxon/family outcomes. Hydrological
downstream site-point relation is NOT evidence butterflies encounter the water.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

ORIGINAL_SRC="results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
ORIGINAL_NEW399="results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json"
RIVER_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
ORIGINAL_EBMS_SHA="2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b"
DANISH_EXTENDED_BBOX=(6.,53.,17.,59.)
SNAP_METRES=500.
GEODESIC_OUTFALL_RADIUS_KM=5.
ROUTE_LIMIT_KM=50.
DOWNSTREAM_PRIMARY_KM=10.
FAMILY_RESPONSE_MIN_SITE=12
FAMILY_RESPONSE_MIN_WWTP=8
FAMILY_RESPONSE_MIN_COMPONENTS=5
FAMILY_RESPONSE_MIN_REGIONS=3
UPSTREAM_MIN_SITE=8
UPSTREAM_MIN_WWTP=4
UPSTREAM_MIN_REGIONS=3

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for d in iter(lambda:f.read(1024*1024),b""):
            h.update(d)
    return h.hexdigest()

def hav(a,b):
    r=math.pi/180
    a1,a2=a["lat"]*r,b["lat"]*r
    dz=(a["lat"]-b["lat"])*r
    dx=(a["lon"]-b["lon"])*r
    f=math.sin(dz/2)**2+math.cos(a1)*math.cos(a2)*math.sin(dx/2)**2
    return 12742*math.asin(min(1.,math.sqrt(max(0.,f))))

def walk(start,target,graph,limit_km=50,max_hops=750):
    """Approximate river reach-path length; same reach returns unknown, not 0."""
    if start==target:return None
    cur=int(start)
    seen=set()
    km=0.
    for _ in range(max_hops):
        if cur in seen or cur not in graph:return None
        seen.add(cur)
        nxt,length=graph[cur]
        if not(math.isfinite(length) and length>=0.):return None
        km+=length
        if km>limit_km:return None
        if nxt==target:return km
        if nxt is None or nxt<=0:return None
        cur=nxt
    return None

def snap_rows(rows, g, maxdist):
    """Tie-break nearest river reach by distance then immutable HydroRIVERS ID."""
    if not rows:return []
    pts=pd.DataFrame(rows)
    if pts.source_row_id.duplicated().any():
        raise ValueError("Source point identifiers duplicated")
    geoms=gpd.GeoDataFrame(pts,geometry=gpd.points_from_xy(pts.lon,pts.lat),crs="EPSG:4326")
    proj=g[["HYRIV_ID","geometry"]].to_crs("EPSG:3035")
    res=gpd.sjoin_nearest(geoms.to_crs("EPSG:3035"),proj,how="left",
                           max_distance=maxdist,distance_col="river_snap_m")
    res=res.sort_values(["source_row_id","river_snap_m","HYRIV_ID"],
                        na_position="last").drop_duplicates("source_row_id")
    out=[]
    for r in res.itertuples(index=False):
        try:
            reach=int(r.HYRIV_ID) if pd.notna(r.HYRIV_ID) else None
        except (ValueError,TypeError):reach=None
        out.append({"source_row_id":r.source_row_id,
                    "HYRIV_ID":reach,
                    "snap_m":float(r.river_snap_m) if pd.notna(r.river_snap_m) else None})
    if len(out)!=len(rows):raise ValueError("Nearest-river merge cardinality not equal source locations")
    return sorted(out,key=lambda r:r["source_row_id"])

def run(root,riverfile):
    p=root/ORIGINAL_SRC
    q=root/ORIGINAL_NEW399
    original=json.loads(p.read_text(encoding="utf-8"))
    cohort=json.loads(q.read_text(encoding="utf-8"))
    if original["original_eBMS_2026_v1_4_source"]["source_archive_SHA256"]!=ORIGINAL_EBMS_SHA:
        raise ValueError("Public pre-outcome eBMS source version changed")
    if not original["butterfly_species_2023_2024_occurrence_extension_READ"] is False:
        raise ValueError("Original site source was not independent of butterfly results")
    if not cohort["NEW2025_Lycaenidae_or_any_butterfly_family_individualCount_NOT_READ"]:
        raise ValueError("Original new399 source cohort contaminated by butterfly counts")
    if not cohort["frozen_before_2025_NEW399_butterfly_counts_gate"]["2025_NEW399_butterfly_family_count_outcome_permitted"]:
        raise ValueError("Original new399 source and independent effort gate failed")
    if cohort["original_2023_24_Event_source_only_735site_SHA256"]!=sha(p):
        raise ValueError("Original eBMS Event-only site coordinates changed")
    sources=original["original_2022_EEA_DK_qualifying_Fw_outfalls"]
    if len(sources)!=203 or len({s["outfall_id_SHA256"] for s in sources})!=203:
        raise ValueError("Expected exactly 203 unique EEA registered Danish source IDs")
    allsites=original["original_2023_2024_repeated_event_site_source_metadata"]
    loc={s["site_id_SHA256"]:s for s in allsites}
    candidates=cohort["all_original_399_source_only_new_site_2025_documented_visit_metadata"]
    if len(candidates)!=399 or len(loc)!=735:
        raise ValueError("Original source site population changed")
    eligible=[s for s in candidates if s["pre_2025_family_any_occurrence_linked_site_source_gate_pass"]]
    if len(eligible)!=286 or Counter(x["source_band"] for x in eligible)!=Counter(
            {"near5":70,"remote10":92,"intermediate":124}):
        raise ValueError("Same original 286 survey-day support cohort must be retained")
    sites=[]
    for s in eligible:
        key=s["site_id_SHA256"]
        v=loc.get(key)
        if v is None:raise ValueError("Original public source-only point not in 735 event sites")
        check=hashlib.sha256(f"{v['lat']:.5f}|{v['lon']:.5f}".encode()).hexdigest()
        if check!=key:raise ValueError("Original public source coordinate hash mismatch")
        if abs(float(s["source_distance_km"])-float(v["nearest_qualified_EEA2022_DK_FW_outfall_km"]))>1e-7:
            raise ValueError("Original static WWTP proximity changed")
        sites.append({
            "source_row_id":"S"+key,
            "site_id_SHA256":key,
            "lat":float(v["lat"]),"lon":float(v["lon"]),
            "original_country_1degree_region":s["original_1degree_region"],
            "original_reported_2025_any_occurrence_dates":int(s["original_linked_occurrence_dates"]),
            "original_geodesic_near_distance_km":float(s["source_distance_km"]),
            "original_geodesic_source_band":s["source_band"],
            "original_near_registered_WWTP_component_SHA256":s.get("5km_original_WWTP_spatial_component_SHA256")
        })
    outlets=[{"source_row_id":"W"+x["outfall_id_SHA256"],"source_id_SHA256":x["outfall_id_SHA256"],
              "lat":float(x["lat"]),"lon":float(x["lon"])} for x in sources]
    if sha(riverfile)!=RIVER_SHA:
        raise ValueError("HydroRIVERS EU v1-derived third-party GeoParquet original SHA changed")
    g=gpd.read_parquet(riverfile)
    if g.crs is None:
        raise ValueError("River GeoParquet lacks original projection")
    g=g.to_crs("EPSG:4326")
    g=g.rename(columns={k:k.upper() for k in g.columns if k!=g.geometry.name})
    if {"HYRIV_ID","NEXT_DOWN","LENGTH_KM"}-set(g.columns):
        raise ValueError("HydroRIVERS original graph schema missing")
    g=g.cx[DANISH_EXTENDED_BBOX[0]:DANISH_EXTENDED_BBOX[2],
             DANISH_EXTENDED_BBOX[1]:DANISH_EXTENDED_BBOX[3]].copy()
    for name in ("HYRIV_ID","NEXT_DOWN","LENGTH_KM"):
        g[name]=pd.to_numeric(g[name],errors="coerce")
    duplicate=g[g.HYRIV_ID.duplicated(keep=False) & g.HYRIV_ID.notna()]
    if len(duplicate):
        q=duplicate.groupby("HYRIV_ID").agg({"NEXT_DOWN":"nunique","LENGTH_KM":"nunique"})
        if ((q.NEXT_DOWN>1)|(q.LENGTH_KM>1)).any():
            raise ValueError("Duplicated HydroRIVERS reach IDs with conflicting direction")
    g=g[g.HYRIV_ID.notna()].sort_values("HYRIV_ID").drop_duplicates("HYRIV_ID")
    graph={int(z.HYRIV_ID):(int(z.NEXT_DOWN) if pd.notna(z.NEXT_DOWN) else 0,float(z.LENGTH_KM))
           for z in g[["HYRIV_ID","NEXT_DOWN","LENGTH_KM"]].itertuples(index=False)}
    if len(graph)<100:
        raise ValueError("Unexpectedly empty Denmark HydroRIVERS source graph")
    site_snap=snap_rows(sites,g,SNAP_METRES)
    outlet_snap=snap_rows(outlets,g,SNAP_METRES)
    bysite={x["source_row_id"]:x for x in site_snap}
    bysource={x["source_row_id"]:x for x in outlet_snap}
    # No 2025 butterfly observations are touched; each site-source candidate
    # is determined exclusively by original public source points and flow graph.
    pair_rows=[]
    records=defaultdict(list)
    for site in sites:
        sk=bysite[site["source_row_id"]]
        near_count=0
        for src in outlets:
            km=hav(site,src)
            if km>GEODESIC_OUTFALL_RADIUS_KM+1e-12:continue
            near_count+=1
            ss=bysource[src["source_row_id"]]
            relation="NOT_EVALUABLE"
            routed=None
            opposite=None
            if sk["HYRIV_ID"] is None:
                relation="SITE_REACH_NOT_SNAPPED"
            elif ss["HYRIV_ID"] is None:
                relation="SOURCE_REACH_NOT_SNAPPED"
            elif sk["HYRIV_ID"]==ss["HYRIV_ID"]:
                relation="SAME_REACH_ORDER_UNKNOWN"
            else:
                routed=walk(ss["HYRIV_ID"],sk["HYRIV_ID"],graph,ROUTE_LIMIT_KM)
                opposite=walk(sk["HYRIV_ID"],ss["HYRIV_ID"],graph,ROUTE_LIMIT_KM)
                if routed is not None:
                    relation=("DOWNSTREAM_CONNECTED_LE10KM" if routed<=DOWNSTREAM_PRIMARY_KM
                             else "DOWNSTREAM_CONNECTED_10TO50KM")
                elif opposite is not None:
                    relation=("UPSTREAM_CONNECTED_LE10KM" if opposite<=DOWNSTREAM_PRIMARY_KM
                             else "UPSTREAM_CONNECTED_10TO50KM")
                else:
                    relation="BOTH_SNAPPED_NO_VERIFIED_DIRECTED_ROUTE_LE50KM"
            r={
                "original_public_site_SHA256":site["site_id_SHA256"],
                "original_registered_outfall_SHA256":src["source_id_SHA256"],
                "original_geodesic_WWTP_to_site_km":km,
                "registered_outfall_to_HydroRIVERS_snap_m":ss["snap_m"],
                "public_site_to_HydroRIVERS_snap_m":sk["snap_m"],
                "original_source_outfall_HYRIV_ID":ss["HYRIV_ID"],
                "original_butterfly_monitor_centroid_HYRIV_ID":sk["HYRIV_ID"],
                "directed_source_to_site_reach_approx_route_km":routed,
                "directed_site_to_source_reach_approx_route_km":opposite,
                "physical_source_2022_to_site2025_class":relation,
                "source_original_geographic_1degree_region":site["original_country_1degree_region"],
                "original_preoutcome_2025_any_occurrence_date_support":site["original_reported_2025_any_occurrence_dates"],
                "same_reach_order_ambiguous_and_not_downstream":sk["HYRIV_ID"]==ss["HYRIV_ID"] if sk["HYRIV_ID"] is not None and ss["HYRIV_ID"] is not None else False,
                "source_near_WWTP_component_SHA256":site["original_near_registered_WWTP_component_SHA256"],
                "no_butterfly_FAMILY_count_accessed":True
            }
            pair_rows.append(r)
            records[site["site_id_SHA256"]].append(r)
        if site["original_geodesic_source_band"]=="near5" and near_count==0:
            raise ValueError("Original site declared <=5km but no matching original outfall")
    if not pair_rows:
        raise ValueError("No original near-source/outfall pairs; source cohort inconsistent")
    total_sources=Counter(r["physical_source_2022_to_site2025_class"] for r in pair_rows)
    downstream=[r for r in pair_rows if r["physical_source_2022_to_site2025_class"]=="DOWNSTREAM_CONNECTED_LE10KM"]
    upstream=[r for r in pair_rows if r["physical_source_2022_to_site2025_class"]=="UPSTREAM_CONNECTED_LE10KM"]
    near_sites=[s for s in sites if s["original_geodesic_source_band"]=="near5"]
    near_siteids={s["site_id_SHA256"] for s in near_sites}
    down_sites={r["original_public_site_SHA256"] for r in downstream if r["original_public_site_SHA256"] in near_siteids}
    up_sites={r["original_public_site_SHA256"] for r in upstream if r["original_public_site_SHA256"] in near_siteids}
    down_sources={r["original_registered_outfall_SHA256"] for r in downstream if r["original_public_site_SHA256"] in near_siteids}
    up_sources={r["original_registered_outfall_SHA256"] for r in upstream if r["original_public_site_SHA256"] in near_siteids}
    down_clusters={r["source_near_WWTP_component_SHA256"] for r in downstream if r["original_public_site_SHA256"] in near_siteids and r["source_near_WWTP_component_SHA256"] is not None}
    ds_region={r["source_original_geographic_1degree_region"] for r in downstream if r["original_public_site_SHA256"] in near_siteids}
    us_region={r["source_original_geographic_1degree_region"] for r in upstream if r["original_public_site_SHA256"] in near_siteids}
    passed=bool(
        len(down_sites)>=FAMILY_RESPONSE_MIN_SITE and len(down_sources)>=FAMILY_RESPONSE_MIN_WWTP
        and len(down_clusters)>=FAMILY_RESPONSE_MIN_COMPONENTS and len(ds_region)>=FAMILY_RESPONSE_MIN_REGIONS
        and len(up_sites)>=UPSTREAM_MIN_SITE and len(up_sources)>=UPSTREAM_MIN_WWTP
        and len(us_region)>=UPSTREAM_MIN_REGIONS
    )
    categories=defaultdict(Counter)
    for s in sites:
        reached=records[s["site_id_SHA256"]]
        if not reached:
            categories[s["original_geodesic_source_band"]]["NO_REGISTERED_QUALIFYING_OUTFALL_WITHIN5KM"]+=1
            continue
        for relation in set(r["physical_source_2022_to_site2025_class"] for r in reached):
            categories[s["original_geodesic_source_band"]][relation]+=1
    out={
        "schema":"gesui_DK2025_eBMS_new399_source286_EEA2022_registeredWWTP203_HydroRIVERS_v1_directed_RECEPTOR_FEASIBILITY_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_source_sha256":{"original_eBMS_Event_only_EEA2022_outfalls":sha(p),
                                   "original_2025_new399_source_only_site_cohort":sha(q),
                                   "third_party_Europe_HydroRIVERS_v1_geoparquet":sha(riverfile)},
        "source_only_no_2025_butterfly_taxon_or_individualCount_read":True,
        "third_party_HydroRIVERS_not_official_verified":True,
        "registered_outfall_2022_status_not_measured_actual_2025_water_receptor_exposure":True,
        "site_500m_nearest_river_is_NOT_necessarily_riparian_or_butterfly_puddling_habitat":True,
        "same_reach_pair_excluded_as_unpositioned":True,
        "2025_original_source_band_sites":dict(Counter(x["original_geodesic_source_band"] for x in sites)),
        "original_Danish_2022_EEA_source_points_n":len(outlets),
        "original_available_HydroRIVERS_reaches_in_Denmark_border_bbox":len(graph),
        "site_HydroRIVERS_reach_snap_le500m_n":sum(s["HYRIV_ID"] is not None for s in site_snap),
        "WWTP_HydroRIVERS_reach_snap_le500m_n":sum(s["HYRIV_ID"] is not None for s in outlet_snap),
        "near_site_HydroRIVERS_snap_le500m_n":sum(bysite[x["source_row_id"]]["HYRIV_ID"] is not None for x in near_sites),
        "site_outfall_geodesic_within5km_pairs_n":len(pair_rows),
        "site_outfall_geodesic_within5km_pair_relation_counts":dict(total_sources),
        "original_source_band_per_site_any_relation_counts":{key:dict(val) for key,val in categories.items()},
        "source_reach_path_distance_threshold_km":[10,25,50],
        "strict_site10km_DOWNSTREAM_original_WWTP_groups":{
            "distinct_near_source_sites":len(down_sites),
            "distinct_reported_FW_outfall_IDs":len(down_sources),
            "distinct_original_near_WWTP_5km_components":len(down_clusters),
            "independent_original_1degree_regions":len(ds_region),
            "sites_SHA256":sorted(down_sites),
            "source_SHA256":sorted(down_sources),
            "near_source_5km_component_SHA256":sorted(down_clusters),
            "regions":sorted(ds_region)},
        "strict_site10km_UPSTREAM_direction_control":{
            "distinct_near_source_sites":len(up_sites),
            "distinct_reported_FW_outfall_IDs":len(up_sources),
            "independent_original_1degree_regions":len(us_region),
            "sites_SHA256":sorted(up_sites),
            "source_SHA256":sorted(up_sources),
            "regions":sorted(us_region)},
        "fixed_source_status_pre_butterfly_FEASIBILITY_GATE":{
            "min_distinct_downstream_near_eBMS_sites":FAMILY_RESPONSE_MIN_SITE,
            "min_distinct_downstream_registered_WWTP_outfalls":FAMILY_RESPONSE_MIN_WWTP,
            "min_distinct_downstream_5km_WWTP_components":FAMILY_RESPONSE_MIN_COMPONENTS,
            "min_downstream_1degree_regions":FAMILY_RESPONSE_MIN_REGIONS,
            "min_distinct_upstream_near_eBMS_site_placebo":UPSTREAM_MIN_SITE,
            "min_upstream_registered_WWTP_outfalls":UPSTREAM_MIN_WWTP,
            "min_upstream_1degree_regions":UPSTREAM_MIN_REGIONS,
            "new_hydrologically_directed_2025_Lycaenidae_family_comparison_Permitted":passed},
        "source_only_original_286_butterfly_monitor_centroids_nearest_HYRIV":[
            {"site_id_SHA256":s["site_id_SHA256"],"original_source_band":s["original_geodesic_source_band"],
             "site_500m_river_snap_m":bysite[s["source_row_id"]]["snap_m"],
             "river_reach_ID":bysite[s["source_row_id"]]["HYRIV_ID"],
             "original_linked_2025_monitor_dates":s["original_reported_2025_any_occurrence_dates"]}
            for s in sites],
        "source_only_original_203_EEA_outfall_nearest_HYRIV":[
            {"outfall_SHA256":s["source_id_SHA256"],
             "outfall_500m_river_snap_m":bysource[s["source_row_id"]]["snap_m"],
             "river_reach_ID":bysource[s["source_row_id"]]["HYRIV_ID"]}
            for s in outlets],
        "all_outfall_site_directed_relation_records_without_butterfly_responses":pair_rows,
        "no_biological_WWTP_positive_or_negative_effect_evaluated":True
    }
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--river-geoparquet",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    value=run(root,args.river_geoparquet)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(value,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
       "expected_all_source_sites":value["2025_original_source_band_sites"],
       "EEA_sources":value["original_Danish_2022_EEA_source_points_n"],
       "rivers_in_bbox":value["original_available_HydroRIVERS_reaches_in_Denmark_border_bbox"],
       "site_snaps":value["site_HydroRIVERS_reach_snap_le500m_n"],
       "outfall_snaps":value["WWTP_HydroRIVERS_reach_snap_le500m_n"],
       "near_source_site_to_outfall_relations":value["site_outfall_geodesic_within5km_pair_relation_counts"],
       "strict_downstream":{k:v for k,v in value["strict_site10km_DOWNSTREAM_original_WWTP_groups"].items() if not isinstance(v,list)},
       "strict_upstream":{k:v for k,v in value["strict_site10km_UPSTREAM_direction_control"].items() if not isinstance(v,list)},
       "outcome_source_gate":value["fixed_source_status_pre_butterfly_FEASIBILITY_GATE"]
    },indent=2))

if __name__=="__main__":
    main()
