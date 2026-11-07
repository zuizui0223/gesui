#!/usr/bin/env python3
"""Hydrologic validation of Peters 2026 TN/TP stations near 16 *fixed* outfalls.

Input Peters source-metadata candidate table comes from original 42369645480
and HydroRIVERS v1 EU-derived GeoParquet with immutable SHA-256. This is a
source-only gate, NOT evidence of actual measured N/P exposure, because the
reported station time envelope is not sufficient for per-year measured values.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import math
import geopandas as gpd
import numpy as np
import pandas as pd

PARQUET_SHA256="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
ORIGINAL_COVERAGE_RUN=42369645480
MAX_STATION_SNAP_M=500
MAX_NETWORK_KM=(10,25,50)
SOURCE_OUTFALLS=16

def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for x in iter(lambda:f.read(1024*1024),b""):
            h.update(x)
    return h.hexdigest()

def walk_to(start,target,river,max_hops=500):
    if start==target:return 0.
    cur=int(start)
    seen=set()
    km=0.
    for _ in range(max_hops):
        if cur in seen or cur not in river: return None
        seen.add(cur)
        nxt,length=river[cur]
        if not math.isfinite(length) or length<0:return None
        km+=length
        if nxt==target:return km
        if not nxt or nxt<0:return None
        cur=nxt
    return None

def relation(snapped,outfall,river):
    if snapped==outfall:return ("OUTFALL_REACH_UNPOSITIONED",0)
    up=walk_to(snapped,outfall,river)
    if up is not None:return ("UPSTREAM_CONNECTED",up)
    dn=walk_to(outfall,snapped,river)
    if dn is not None:return ("DOWNSTREAM_CONNECTED",dn)
    return ("NOT_CONNECTED_IN_SOURCE_GRAPH",None)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--candidates",required=True,type=Path)
    p.add_argument("--prior-coverage",required=True,type=Path)
    p.add_argument("--design",required=True,type=Path)
    p.add_argument("--river-parquet",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    coverage=json.loads(args.prior_coverage.read_text())
    assert coverage["study_outfalls"]==SOURCE_OUTFALLS and coverage["study_sites"]==48
    assert coverage["all_source_files_md5_verified"] is True
    design=json.loads(args.design.read_text())
    selected=design["selected"]
    assert len(selected)==SOURCE_OUTFALLS
    sha=sha256(args.river_parquet)
    if sha!=PARQUET_SHA256:raise ValueError("HydroRIVERS mirror checksum differs")
    candidates=pd.read_csv(args.candidates,dtype={"source_station_id":str,"site_id":str})
    assert len(candidates)==coverage["station_candidate_pairs_within_10km"]
    assert set(candidates.nutrient)=={"TN","TP"}
    stations=candidates.drop_duplicates(["nutrient","source_station_id"]).copy()
    stations["station_key"]=stations.nutrient+":"+stations.source_station_id
    assert stations.station_key.nunique()==len(stations)
    gdf=gpd.read_parquet(args.river_parquet)
    names={str(c).lower():str(c) for c in gdf.columns}
    for col in ("hyriv_id","next_down","length_km"):
        if col not in names:raise ValueError("HydroRIVERS missing "+col)
    gdf=gdf.rename(columns={names["hyriv_id"]:"HYRIV_ID",names["next_down"]:"NEXT_DOWN",
                            names["length_km"]:"LENGTH_KM"})
    if gdf.crs is None:raise ValueError("HydroRIVERS lacks CRS")
    gdf=gdf.to_crs("EPSG:4326")
    gdf=gdf.cx[-13:27,43:65].copy()
    gdf["HYRIV_ID"]=pd.to_numeric(gdf.HYRIV_ID,errors="raise").astype("int64")
    gdf["NEXT_DOWN"]=pd.to_numeric(gdf.NEXT_DOWN,errors="coerce").fillna(0).astype("int64")
    gdf["LENGTH_KM"]=pd.to_numeric(gdf.LENGTH_KM,errors="coerce")
    gdf=gdf.drop_duplicates("HYRIV_ID")
    river={int(z.HYRIV_ID):(int(z.NEXT_DOWN),float(z.LENGTH_KM))
           for z in gdf[["HYRIV_ID","NEXT_DOWN","LENGTH_KM"]].itertuples(index=False)}
    pts=gpd.GeoDataFrame(stations[["station_key"]].copy(),
                        geometry=gpd.points_from_xy(stations.source_station_lon,
                                                    stations.source_station_lat),
                        crs="EPSG:4326")
    projection="EPSG:3035"
    snapped=gpd.sjoin_nearest(pts.to_crs(projection),
        gdf[["HYRIV_ID","geometry"]].to_crs(projection),
        how="left",max_distance=MAX_STATION_SNAP_M,distance_col="snap_distance_m")
    snapped=snapped.sort_values(["station_key","snap_distance_m","HYRIV_ID"],
                                na_position="last").drop_duplicates("station_key")
    station_snap={}
    for row in snapped.itertuples():
        if pd.isna(row.HYRIV_ID):
            station_snap[row.station_key]={"snapped_hyriv_id":None,"station_snap_m":None}
        else:
            station_snap[row.station_key]={
                "snapped_hyriv_id":int(row.HYRIV_ID),
                "station_snap_m":float(row.snap_distance_m)
            }
    byout={int(w["wwtp_id"]):w for w in selected}
    rows=[]
    for v in candidates.itertuples(index=False):
        sample=station_snap[str(v.nutrient)+":"+str(v.source_station_id)]
        water=byout[int(v.wwtp_id)]
        related,km=relation(sample["snapped_hyriv_id"],int(water["river_id"]),river) \
           if sample["snapped_hyriv_id"] else ("NO_RIVER_SNAP_WITHIN_500M",None)
        rec={
          "nutrient":str(v.nutrient),
          "outfall_id":int(v.wwtp_id),
          "plant_site_id":str(v.site_id),
          "plant_site_role":str(v.site_role),
          "plant_site_river_id":int(v.river_reach_id),
          "water_quality_station_id":str(v.source_station_id),
          "station_lat":float(v.source_station_lat),
          "station_lon":float(v.source_station_lon),
          "distance_to_plant_site_km":float(v.candidate_haversine_distance_km),
          "source_start":str(v.source_first_observation),
          "source_end":str(v.source_last_observation),
          "n_reported_measurements_full_period":int(v.source_observations_reported),
          "candidate_stage":"metadata_date_window_overlap_not_measured_years_verified",
          "source_snapped_hyriv_id":sample["snapped_hyriv_id"],
          "station_to_river_snap_m":sample["station_snap_m"],
          "network_position_relative_to_outfall":related,
          "network_path_distance_km":km,
          "official_HydroRIVERS_bytes_or_topology_independently_verified":False,
        }
        rows.append(rec)
    result=pd.DataFrame(rows)
    args.out.mkdir(parents=True,exist_ok=True)
    result.to_csv(args.out/"source_candidate_station_river_connections.csv",index=False,float_format="%.8g")
    summary={
       "schema":"gesui_Peters2026_TNTP_station_to_HydroRIVERS_source_network_20261107",
       "timestamp_utc":datetime.now(timezone.utc).isoformat(),
       "source_Peters_figshare_v5":"10.25400/lincolnuninz.30546371.v5",
       "source_station_metadata_original_run":ORIGINAL_COVERAGE_RUN,
       "source_metadata_sha256":sha256(args.prior_coverage),
       "source_candidate_table_sha256":sha256(args.candidates),
       "source_outfall_design_sha256":sha256(args.design),
       "HydroRIVERS_derived_europe_sha256":sha,
       "HydroRIVERS_official_original_independently_verified":False,
       "source_sites_48_outfalls_16":True,
       "source_snapping_threshold_m":MAX_STATION_SNAP_M,
       "candidate_station_site_pairs_within_10km":len(result),
       "unique_monitoring_source_stations_by_nutrient":{
          k:int(g["water_quality_station_id"].nunique()) for k,g in result.groupby("nutrient")
       },
       "unique_stations_snapped_within_500m_by_nutrient":{
          k:int(g.drop_duplicates("water_quality_station_id").source_snapped_hyriv_id.notna().sum())
          for k,g in result.groupby("nutrient")
       },
       "position_count_by_nutrient":{
          k:{str(i):int(j) for i,j in g.network_position_relative_to_outfall.value_counts().items()}
          for k,g in result.groupby("nutrient")
       },
       "network_distance_km_ceilings":{},
       "measurement_values_or_dated_rows_loaded":False,
       "plant_response_results_not_examined":True,
       "warning":"These are network-positioned station metadata candidates only; real TN/TP observations in 2011–24 have not been verified. Other WWTP discharges and land-use/effort are uncontrolled."
    }
    for km in MAX_NETWORK_KM:
        counts={}
        for nutrient,g in result.groupby("nutrient"):
            x=g[g.network_position_relative_to_outfall.isin(
                ["UPSTREAM_CONNECTED","DOWNSTREAM_CONNECTED"])&
                g.network_path_distance_km.le(km)]
            outs={}
            for v in x.itertuples(index=False):
                outs.setdefault(v.outfall_id,set()).add(v.network_position_relative_to_outfall)
            counts[nutrient]={
                "source_outfalls_with_connected_any":len(outs),
                "source_outfalls_with_upstream_and_downstream_connected":
                    sum(v=={"UPSTREAM_CONNECTED","DOWNSTREAM_CONNECTED"} for v in outs.values()),
                "source_outfalls_upstream_supported":sum("UPSTREAM_CONNECTED" in v for v in outs.values()),
                "source_outfalls_downstream_supported":sum("DOWNSTREAM_CONNECTED" in v for v in outs.values()),
                "unique_connected_stations":int(x.water_quality_station_id.nunique())
            }
        summary["network_distance_km_ceilings"][str(km)]=counts
    (args.out/"NETWORK_SOURCE_QC.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print("GESUI_RIVER_CONFIRMED_NUTRIENT_SOURCE_COVERAGE",json.dumps(summary,sort_keys=True),flush=True)
if __name__=="__main__":main()
