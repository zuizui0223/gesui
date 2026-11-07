#!/usr/bin/env python3
"""Measured TN/TP *station metadata* coverage for fixed 48 river-plant sites.

Peters 2026 v5 source files verified MD5. Spatial candidate threshold only;
haversine vicinity is NOT an established hydraulic connection. No plant response
read, no water quality concentration estimates made.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import io
import json
import math
from pathlib import Path
import pandas as pd
import requests
from scipy.spatial import cKDTree
import numpy as np

ARTICLE="https://api.figshare.com/v2/articles/30546371"
DOI="10.25400/lincolnuninz.30546371.v5"
NEEDED=("TN.2_meta.csv","TP.2_meta.csv","TN.2_years.csv","TP.2_years.csv")
RADII_KM=(0.5,1,2,5,10)
YEARS_START="2011-01-01"
YEARS_END="2024-12-31"
EARTH=6371.0088

def download(url,cap=5_000_000):
    r=requests.get(url,timeout=(15,90),stream=True,headers={"User-Agent":"gesui-river-nutrient-station-coverage/0.1"})
    r.raise_for_status()
    buf=bytearray()
    for b in r.iter_content(131072):
        buf.extend(b)
        if len(buf)>cap:
            r.close()
            raise ValueError("Nutrient file larger than declared cap")
    r.close()
    return bytes(buf)

def xyz(lat,lon):
    a=np.deg2rad(np.asarray(lat,dtype=float))
    b=np.deg2rad(np.asarray(lon,dtype=float))
    return np.column_stack([np.cos(a)*np.cos(b),np.cos(a)*np.sin(b),np.sin(a)])

def rad_chord(km):
    return 2*math.sin(km/(2*EARTH))

def great_circle(a,b):
    return float(2*EARTH*math.asin(min(1,np.linalg.norm(a-b)/2)))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sites",required=True,type=Path)
    ap.add_argument("--outdir",required=True,type=Path)
    args=ap.parse_args()
    sites=pd.read_csv(args.sites,dtype={"site_id":str},keep_default_na=False)
    if len(sites)!=48 or sites.site_id.nunique()!=48:
        raise ValueError("48 preselected source-first biological sites changed")
    assert set(sites.river_role)=={"u2","u1","d1"}
    assert sites.wwtp_id.nunique()==16
    for c in ("latitude","longitude"):
        sites[c]=pd.to_numeric(sites[c],errors="raise")
    src=json.loads(download(ARTICLE,500_000))
    if src["doi"]!=DOI:raise ValueError("Peters article version changed")
    fidx={x["name"]:x for x in src["files"]}
    if set(NEEDED)-set(fidx):raise ValueError("missing original files")
    site3d=xyz(sites.latitude,sites.longitude)
    output=args.outdir
    output.mkdir(parents=True,exist_ok=True)
    source=[]
    station_matches=[]
    grid={}
    for p in ("TN","TP"):
        original={}
        for suffix in ("meta","years"):
            name=f"{p}.2_{suffix}.csv"
            f=fidx[name]
            raw=download(f["download_url"],5_000_000)
            md5=hashlib.md5(raw).hexdigest()
            if md5!=f["supplied_md5"]:raise ValueError(name+" MD5 mismatched")
            table=pd.read_csv(io.BytesIO(raw),low_memory=False,dtype={"WQ.ID":str})
            original[suffix]=table
            source.append({"file":name,"file_id":f["id"],"md5":md5,"original_n":len(table),
                          "original_download_size_bytes":len(raw)})
        meta=original["meta"].copy()
        years=original["years"].copy()
        lat=pd.to_numeric(meta["WQ.lat"],errors="coerce")
        lon=pd.to_numeric(meta["WQ.lon"],errors="coerce")
        meta["start_day"]=pd.to_datetime(meta["start"],errors="coerce",utc=True)
        meta["end_day"]=pd.to_datetime(meta["end"],errors="coerce",utc=True)
        years["start_day"]=pd.to_datetime(years["start"],errors="coerce",utc=True)
        years["end_day"]=pd.to_datetime(years["end"],errors="coerce",utc=True)
        subset=meta.loc[
            lat.between(35,65)&lon.between(-15,35)
            &(meta["start_day"]<=pd.Timestamp(YEARS_END,tz="UTC"))
            &(meta["end_day"]>=pd.Timestamp(YEARS_START,tz="UTC"))
        ].copy()
        subset["WQ.lat"]=pd.to_numeric(subset["WQ.lat"],errors="coerce")
        subset["WQ.lon"]=pd.to_numeric(subset["WQ.lon"],errors="coerce")
        # Year table is an extra source-level history gate, not proof of annual measurements.
        yearids=set(years.loc[
            (years["start_day"]<=pd.Timestamp(YEARS_END,tz="UTC"))&
            (years["end_day"]>=pd.Timestamp(YEARS_START,tz="UTC")),"WQ.ID"])
        subset=subset[subset["WQ.ID"].isin(yearids)].copy()
        subset=subset.drop_duplicates("WQ.ID")
        keys=subset["WQ.ID"].tolist()
        coord=xyz(subset["WQ.lat"],subset["WQ.lon"])
        tree=cKDTree(coord) if len(coord) else None
        hits={}
        for km in RADII_KM:
            candidates=tree.query_ball_point(site3d,rad_chord(km)) if tree is not None else [[] for _ in range(48)]
            perout={}
            for row,indices in zip(sites.itertuples(index=False),candidates):
                for i in indices:
                    perout.setdefault(int(row.wwtp_id),set()).add(keys[i])
            hits[str(km)]={
                "number_sites_with_station_within_radius":sum(bool(x) for x in candidates),
                "unique_monitoring_stations_in_any_radius":len(set(keys[i] for ids in candidates for i in ids)),
                "source_outfalls_with_station_near_any_site":len(perout),
                "sites_u2_with_station":sum(bool(ids) and sites.iloc[ii].river_role=="u2" for ii,ids in enumerate(candidates)),
                "sites_u1_with_station":sum(bool(ids) and sites.iloc[ii].river_role=="u1" for ii,ids in enumerate(candidates)),
                "sites_d1_with_station":sum(bool(ids) and sites.iloc[ii].river_role=="d1" for ii,ids in enumerate(candidates)),
                "source_outfalls_with_u1_and_d1_both_near_station":sum(
                    bool(candidates[i]) and bool(candidates[j])
                    for i in range(48) for j in [i+1]
                    if i%3==1
                ),
                "source_outfalls_with_all_three_sites_near_station":sum(
                    all(bool(candidates[idx+k]) for k in (0,1,2))
                    for idx in range(0,48,3)
                ),
            }
        grid[p]=hits
        print("GESUI_NUTRIENT_STATION_NEARNESS",
              json.dumps({"parameter":p,"total_metadata":len(meta),
                  "time_overlap_and_era_years_table_filtered":len(subset),
                  "threshold_counts":hits},sort_keys=True),flush=True)
        if tree is not None:
            # Candidate site x station pairs within 10km; exact coords and source IDs
            # stored only in internal artifact for topology validation, not public Git.
            neighbors=tree.query_ball_point(site3d,rad_chord(10))
            for site,row_indices,point in zip(sites.itertuples(index=False),neighbors,site3d):
                for i in row_indices:
                    r=subset.iloc[i]
                    station_matches.append({
                      "nutrient":p,"site_id":site.site_id,"site_role":site.river_role,
                      "wwtp_id":int(site.wwtp_id),"river_reach_id":int(site.river_reach_id),
                      "source_station_id":str(r["WQ.ID"]),
                      "source_station_lat":float(r["WQ.lat"]),
                      "source_station_lon":float(r["WQ.lon"]),
                      "candidate_haversine_distance_km":great_circle(point,coord[i]),
                      "source_first_observation":str(r["start"]),
                      "source_last_observation":str(r["end"]),
                      "source_observations_reported":int(r["sum"]),
                      "source_total_days":int(r["n.days"]),
                      "era_years_table_overlap":True,
                      "spatial_status":"PROXIMITY_ONLY_NOT_RIVER_CONNECTED",
                    })
    df=pd.DataFrame(station_matches)
    if len(df):
        df=df.sort_values(["nutrient","wwtp_id","site_role","candidate_haversine_distance_km","source_station_id"])
        df.to_csv(output/"nutrient_candidate_station_proximities.csv",index=False,float_format="%.8f")
    result={"schema":"gesui_Peters2026_real_TN_TP_metadata_spatial_cover_gate_v01",
            "date":datetime.now(timezone.utc).isoformat(),
            "source_doi":DOI,
            "source_license":"CC BY 4.0",
            "all_source_files_md5_verified":True,
            "files":source,"study_sites":48,"study_outfalls":16,
            "date_gate":"metadata start/end overlap 2011-2024 PLUS .2_years start/end overlap",
            "proximity_km":[0.5,1,2,5,10],
            "station_candidate_pairs_within_10km":len(df),
            "by_nutrient":grid,
            "site_to_station_river_topology_verified":False,
            "measurement_concentrations_loaded":False,
            "warning":"A station near a plant site may lie in a different water body; metadata time envelope does not establish actual measurements in the requested years. This is a source support census only."
    }
    (output/"SOURCE_COVERAGE.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_NUTRIENT_SITE_SUPPORT_COMPLETE",
          json.dumps({k:v for k,v in result.items() if k not in ("by_nutrient","files")},sort_keys=True),flush=True)
if __name__=="__main__":main()
