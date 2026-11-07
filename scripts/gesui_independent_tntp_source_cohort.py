#!/usr/bin/env python3
"""Independent SOURCE-ONLY nutrient-ready outfalls from 3569 real river chains.

Recompute HydroWASTE/derived HydroRIVERS source candidate inventory exactly,
then select monitoring-ready sites using original Peters 2026 v5 measured
TN/TP station inventory (available=1, actual public station file exists).
No GBIF or plant response read; original 16 outfalls excluded.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile
import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import requests
from gesui_river_source_gate import (
    SOURCE_HYDROWASTE_BLOB,BBOX,BANNED_STATUS,COUNTRIES,git_blob_sha,
    read_river_subset,cohort_candidates
)
from gesui_nutrient_station_network_gate import walk_to

DOI="10.25400/lincolnuninz.30546371.v5"
ARTICLE="https://api.figshare.com/v2/articles/30546371"
RIVER_SHA256="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
NEEDED=("TN.2_meta.csv","TP.2_meta.csv","TN.2_years.csv","TP.2_years.csv",
        "TN.2_-_Harmonised.zip","TP.2_-_Harmonised.zip")
PRIMARY_MAX_NETWORK_KM=25
PRIMARY_STATION_SNAP_M=500
MIN_OBSERVATIONS_REPORTED=12
R=6371.0088

def get(url,max_bytes):
    r=requests.get(url,stream=True,timeout=(15,90),
      headers={"User-Agent":"gesui-independent-monitored-outfalls-coverage/0.1"})
    r.raise_for_status();b=bytearray()
    for x in r.iter_content(131072):
        b.extend(x)
        if len(b)>max_bytes:
            r.close();raise ValueError("Source exceeds configured bound")
    r.close();return bytes(b)

def xyz(lat,lon):
    a=np.deg2rad(np.asarray(lat,dtype=float))
    b=np.deg2rad(np.asarray(lon,dtype=float))
    return np.column_stack([np.cos(a)*np.cos(b),
                            np.cos(a)*np.sin(b),np.sin(a)])
def chord(km):
    return 2*math.sin(km/(2*R))

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--hydrowaste",required=True,type=Path)
    p.add_argument("--river-parquet",required=True,type=Path)
    p.add_argument("--old16",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    if git_blob_sha(a.hydrowaste)!=SOURCE_HYDROWASTE_BLOB:
        raise ValueError("Original 58502 HydroWASTE source hash mismatch")
    if hashlib.sha256(a.river_parquet.read_bytes()).hexdigest()!=RIVER_SHA256:
        raise ValueError("Original source-verified European river GeoParquet changed")
    old=json.loads(a.old16.read_text())
    excluded={int(v["wwtp_id"]) for v in old["selected"]}
    assert len(excluded)==16
    ww=pd.read_csv(a.hydrowaste,encoding="latin1",low_memory=False)
    if len(ww)!=58502:raise ValueError("Original facility denominator changed")
    for k in ("HYRIV_ID","QUAL_LOC","QUAL_WASTE","WASTE_DIS","RIVER_DIS","COAST_10KM"):
        ww[k]=pd.to_numeric(ww[k],errors="coerce")
    eligible=ww[
        ww.CNTRY_ISO.isin(COUNTRIES)
        & (~ww.STATUS.astype(str).str.lower().isin(BANNED_STATUS))
        & (ww.QUAL_LOC<=2)&(ww.QUAL_WASTE<=3)
        & (ww.WASTE_DIS>0)&(ww.RIVER_DIS>=0.5)
        & ww.HYRIV_ID.notna()
        & (ww.COAST_10KM==0)
        & ww.LAT_OUT.between(BBOX[1],BBOX[3])
        & ww.LON_OUT.between(BBOX[0],BBOX[2])
    ].copy()
    river,upstream,n_rivers=read_river_subset(a.river_parquet)
    river_graph={int(k):(int(v.NEXT_DOWN) if pd.notna(v.NEXT_DOWN) else 0,
                         float(v.LENGTH_KM) if pd.notna(v.LENGTH_KM) else float("nan"))
                 for k,v in river.items()}
    source_candidates,chain_valid=cohort_candidates(eligible,river,upstream)
    assert len(source_candidates)==3569, "Original fixed 3569 source-only non-overlap geography changed"
    candid=[x for x in source_candidates if int(x["wwtp_id"]) not in excluded]
    a.out.mkdir(parents=True,exist_ok=True)
    article=json.loads(get(ARTICLE,500000))
    if article["doi"]!=DOI:raise ValueError("Global nutrient v5 version shifted")
    fidx={x["name"]:x for x in article["files"]}
    frames={}
    files=[]
    for name in NEEDED:
        f=fidx[name]
        raw=get(f["download_url"],32_000_000)
        digest=hashlib.md5(raw).hexdigest()
        if digest!=f["supplied_md5"]:raise ValueError("Original Peters source MD5 mismatch")
        if name.endswith(".csv"):
            frames[name]=pd.read_csv(io.BytesIO(raw),low_memory=False,dtype={"WQ.ID":str})
        else:
            prefix=name.split(".")[0]
            with zipfile.ZipFile(io.BytesIO(raw)) as z:
                frames[name]=set(
                    Path(path).stem for path in z.namelist() if path.lower().endswith(".csv"))
        files.append({"name":name,"original_figshare_file_id":f["id"],"md5":digest,
                      "original_bytes":len(raw)})
        print("GESUI_SOURCE_FILE_PINNED",json.dumps({"file":name,"original_bytes":len(raw)},sort_keys=True),flush=True)
    stations=[]
    source_stage_counts={}
    for par in ("TN","TP"):
        df=frames[f"{par}.2_meta.csv"].copy()
        years=frames[f"{par}.2_years.csv"].copy()
        df["start"]=pd.to_datetime(df["start"],errors="coerce",utc=True)
        df["end"]=pd.to_datetime(df["end"],errors="coerce",utc=True)
        years["start"]=pd.to_datetime(years["start"],errors="coerce",utc=True)
        years["end"]=pd.to_datetime(years["end"],errors="coerce",utc=True)
        validyears=set(years.loc[years["start"].le("2024-12-31")&years["end"].ge("2011-01-01"),
                                  "WQ.ID"])
        df["WQ.lat"]=pd.to_numeric(df["WQ.lat"],errors="coerce")
        df["WQ.lon"]=pd.to_numeric(df["WQ.lon"],errors="coerce")
        df["sum"]=pd.to_numeric(df["sum"],errors="coerce")
        df["WQ.available"]=pd.to_numeric(df["WQ.available"],errors="coerce")
        full=len(df)
        base=df[
            df["WQ.lat"].between(BBOX[1],BBOX[3])
            &df["WQ.lon"].between(BBOX[0],BBOX[2])
            &df["start"].le("2024-12-31")&df["end"].ge("2011-01-01")
            &df["WQ.ID"].isin(validyears)&df["WQ.available"].eq(1)
            &df["sum"].ge(MIN_OBSERVATIONS_REPORTED)
            &df["WQ.ID"].isin(frames[f"{par}.2_-_Harmonised.zip"])
        ].drop_duplicates("WQ.ID").copy()
        source_stage_counts[par]={"all_original_metadata":full,
            "spatiotemporal_available_and_zip_confirmed":len(base)}
        for s in base.itertuples(index=False):
            stations.append({"parameter":par,"station_id":s[0],
                "lat":float(getattr(s,"_1")) if False else float(s[1]),
                "lon":float(s[2]),"start":str(s[4]) if False else None})
        # Avoid fragile namedtuple fields with dots; explicitly regenerate.
        stations=[x for x in stations if x["parameter"]!=par]
        stations.extend({"parameter":par,"station_id":str(row["WQ.ID"]),
             "lat":float(row["WQ.lat"]),"lon":float(row["WQ.lon"]),
             "start_date":str(row["start"].date()),"end_date":str(row["end"].date()),
             "total_observations_reported":int(row["sum"])}
            for _,row in base.iterrows())
    s=pd.DataFrame(stations)
    if s.empty:raise ValueError("No source-available TN and TP stations")
    s["key"]=s["parameter"]+":"+s["station_id"]
    if s.key.duplicated().any():raise ValueError("Duplicate monitored station key")
    points=gpd.GeoDataFrame(s[["key"]].copy(),
           geometry=gpd.points_from_xy(s.lon,s.lat),crs="EPSG:4326")
    rivers=gpd.read_parquet(a.river_parquet)
    rivers=rivers.to_crs("EPSG:4326")
    rivers=rivers.cx[BBOX[0]:BBOX[2],BBOX[1]:BBOX[3]].copy()
    if "HYRIV_ID" not in rivers.columns:
        rivers=rivers.rename(columns={c:"HYRIV_ID" for c in rivers if c.lower()=="hyriv_id"})
    snapped=gpd.sjoin_nearest(points.to_crs("EPSG:3035"),
        rivers[["HYRIV_ID","geometry"]].to_crs("EPSG:3035"),how="left",
        max_distance=PRIMARY_STATION_SNAP_M,distance_col="snap_distance_m")
    snapped=snapped.sort_values(["key","snap_distance_m","HYRIV_ID"],
         na_position="last").drop_duplicates("key")
    matchmap={}
    for z in snapped.itertuples():
        if pd.notna(z.HYRIV_ID):
            matchmap[z.key]=(int(z.HYRIV_ID),float(z.snap_distance_m))
    ss=s[s.key.isin(matchmap)].copy()
    if not len(ss):raise ValueError("No source monitoring stations within river network")
    ss["snapped_reach_id"]=ss.key.map(lambda x:matchmap[x][0])
    ss["snap_m"]=ss.key.map(lambda x:matchmap[x][1])
    point_tree=cKDTree(xyz(ss.lat,ss.lon))
    rows=[]
    country_count=Counter()
    for outfall in candid:
        nearby=point_tree.query_ball_point(
            xyz([outfall["outfall_lat"]],[outfall["outfall_lon"]])[0],
            chord(PRIMARY_MAX_NETWORK_KM))
        if not nearby:continue
        options={"TN":{"UPSTREAM":[],"DOWNSTREAM":[]},
                 "TP":{"UPSTREAM":[],"DOWNSTREAM":[]}}
        for i in nearby:
            station=ss.iloc[i]
            reach=int(station.snapped_reach_id);r0=int(outfall["river_id"])
            if reach==r0:continue  # same reach cannot assign side of outfall
            upstream_dist=walk_to(reach,r0,river_graph)
            downstream_dist=walk_to(r0,reach,river_graph)
            if upstream_dist is not None and 0<upstream_dist<=PRIMARY_MAX_NETWORK_KM:
                side="UPSTREAM";d=upstream_dist
            elif downstream_dist is not None and 0<downstream_dist<=PRIMARY_MAX_NETWORK_KM:
                side="DOWNSTREAM";d=downstream_dist
            else:continue
            par=str(station.parameter)
            options[par][side].append({
                "station_id":str(station.station_id),"side":side,
                "source_observation_total":int(station.total_observations_reported),
                "lat":float(station.lat),"lon":float(station.lon),
                "station_hyriv_id":reach,
                "source_snap_distance_m":float(station.snap_m),
                "network_distance_km":float(d),
                "source_start":station.start_date,"source_end":station.end_date})
        if not all(options[p][side] for p in ("TN","TP") for side in ("UPSTREAM","DOWNSTREAM")):
            continue
        chosen={p:{} for p in ("TN","TP")}
        for p in ("TN","TP"):
            for side in ("UPSTREAM","DOWNSTREAM"):
                candidates=options[p][side]
                # Predefine nearest river-network distance, then source station id tie-break.
                sel=min(candidates,key=lambda x:(x["network_distance_km"],x["station_id"]))
                chosen[p][side]=sel
        starts=[q["source_start"] for pp in chosen.values() for q in pp.values()]
        ends=[q["source_end"] for pp in chosen.values() for q in pp.values()]
        if max(starts)>min(ends):continue
        rows.append({**outfall,"measured_source_stations":chosen,
                     "four_station_common_time_envelope":[max(starts),min(ends)],
                     "source_only_monitoring_test_eligible":True})
        country_count[outfall["country"]]+=1
    # Source-only cohort selection; never inspect plant species composition.
    rows.sort(key=lambda x:(x["country"],x["effluent_fraction_proxy"],x["wwtp_id"]))
    selected=[]
    # Within country select max three low and three high water fraction cases;
    # leave original 16 outfalls excluded.
    for country in COUNTRIES:
        c=[x for x in rows if x["country"]==country]
        if len(c)<4:continue
        m=max(2,len(c)//4)
        low=c[:m]
        high=c[-m:][::-1]
        for source_type,group in (("low",low),("high",high)):
            for x in group[:3]:
                if x["wwtp_id"] in {z["wwtp_id"] for z in selected}:continue
                selected.append({**x,"monitoring_water_fraction_stratum":source_type})
    report={
      "schema":"gesui_new_independent_monitored_source_cohort_v01",
      "timestamp_utc":datetime.now(timezone.utc).isoformat(),
      "Peters_source_doi":DOI,
      "Peters_original_verified_source_files":files,
      "HydroWASTE_git_blob":SOURCE_HYDROWASTE_BLOB,
      "HydroRIVERS_GeoParquet_SHA256":RIVER_SHA256,
      "HydroRIVERS_independent_official_validation":False,
      "old16_excluded_from_new_cohort":True,
      "original_source_candidates_nonoverlap":len(source_candidates),
      "new_candidate_outfalls_after_old16_exclusion":len(candid),
      "source_monitoring_metadata_counts":source_stage_counts,
      "stations_snapped_under_500m":dict(Counter(ss.parameter)),
      "source_outfalls_with_all_four_available_TNTP_station_positions_and_overlap":len(rows),
      "cohort_geographic_available_by_country":dict(country_count),
      "source_only_stratified_selected_outfalls":len(selected),
      "selected_by_country":dict(Counter(x["country"] for x in selected)),
      "selection_before_any_new_biological_response":True,
      "observation_values_used_for_selection":False,
      "source_measurement_values_loaded":False,
      "caveats":[
        "No observed TN/TP concentration used for selection; original station metadata envelopes may still hide sparse actual years",
        "Source available=1 and site file exists but individual observation dates have not been checked",
        "Natural river branch/other WWTP exposure and distance imbalance not controlled",
        "Mirror is a third-party converted river network, not independently matched official original",
        "Old 16 outfalls must remain separate from new sample and should not be mixed for unbiased generalization",
        "A new monitoring-based sample selection enriches monitored catchments and is not a population-prevalence estimator"
      ]
    }
    (a.out/"SOURCE_READY_COHORT_QC.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    (a.out/"NEW_SOURCE_BLIND_COHORT.json").write_text(json.dumps({
       "provenance":report,"selected":selected},indent=2,sort_keys=True)+"\n")
    print("GESUI_NEW_MEASURED_TNTP_SOURCE_READY_COHORT",
       json.dumps({k:v for k,v in report.items() if k not in ("Peters_original_verified_source_files","caveats")},
                  sort_keys=True),flush=True)
    if not selected:
        print("GESUI_NEW_COHORT_GATE_NO_SOURCE_SUPPORTED_SITES: source coverage failure, do not query plants",flush=True)
if __name__=="__main__":main()
