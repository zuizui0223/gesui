#!/usr/bin/env python3
"""Exact original Austrian GZÜV PKINT station ID → verified EEA2022 benthic
site-key source crosswalk and genuine AT145 reported N/P outfall locality.

GZÜV FWMST.gpkg is a public national station registry (3788 records);
PKINT=ATFW... is uniquely sourced and matches EEA public site codes. File
BEGIN_DAT appears to be a common source export timestamp, NOT sampling year.
Do not read actual EQR, TNTP chemistry values or infer private WISA access.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,math,sys
from pathlib import Path
import numpy as np
import geopandas as gpd
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_austria_wisa_official_source"))
from download_ubagzuev_and_match_AT145 import (
    original_public_download,original_station_points,originals,haversine_distance_vector,
    RADIUS_M,ORIGINAL_DATASET_ID
)
import requests

BQE=ROOT/"results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv"
BQE_BLOB="325f199e985d9862ffee5b0edb205c61a7be3833"
SCHEMA="GESUI_20261116_GZUEV_TRUE_PUBLIC_PKINT_TO_EEA2022_BQE_68_AND_AT145_SOURCE_V1"

def blob(p):
    b=p.read_bytes()
    return hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()

def EEA_sites():
    if blob(BQE)!=BQE_BLOB:
        raise RuntimeError("Original merged EEA2022 monitored benthos source file Gitblob changed")
    rows=[]
    with BQE.open(newline="",encoding="utf-8-sig") as f:
        reader=csv.reader(f)
        cols=next(reader)
        wanted=("site_key","countryCode","lat","lon","HYRIV_ID",
                  "procedureClassificationSystem","observedPropertyDeterminandBiologyEQRCode")
        if not set(wanted).issubset(cols):
            raise RuntimeError("Official EEA benthic original source ID/geom/method columns missing")
        idx={k:cols.index(k) for k in wanted}
        n=0
        for raw in reader:
            n+=1
            if raw[idx["countryCode"]]!="AT":continue
            key=raw[idx["site_key"]]
            parts=key.split("|")
            if len(parts)!=3 or parts[0]!="AT" or parts[2]!="euMonitoringSiteCode":
                raise RuntimeError("Original exact Austrian EEA BQE station identity altered")
            rows.append({
              "original_EEA_station_code":parts[1],"EEA_site_key":key,
              "source_EEA_2022_lat":float(raw[idx["lat"]]),
              "source_EEA_2022_lon":float(raw[idx["lon"]]),
              "source_EEA_receiver_HYRIV_ID":raw[idx["HYRIV_ID"]],
              "source_EEA_2022_method":raw[idx["procedureClassificationSystem"]]})
    if n!=2475 or len(rows)!=68 or len({x["original_EEA_station_code"] for x in rows})!=68:
        raise RuntimeError("Original EEA2022 2475 benthic / 68 Austrian distinct station source population changed")
    return rows

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    session=requests.Session()
    session.headers.update({"User-Agent":"gesui-20261116-public-GZUEV-station-EEA-ID-academic/1.0"})
    gpkg,source=original_public_download(session,a.output)
    gdf,station_PK,qa=original_station_points(gpkg)
    if len(gdf)!=3788 or qa["source_primary_unique_station_identifier_field"]!="PK":
        raise RuntimeError("Genuine 3788-source GZÜV PK station registry not reproduced")
    if "PKINT" not in gdf or gdf["PKINT"].isna().any():
        raise RuntimeError("Actual original GZÜV PKINT identifier missing")
    PKINT=gdf["PKINT"].astype(str).str.strip()
    if PKINT.nunique()!=len(gdf) or not PKINT.str.startswith("ATFW").all():
        raise RuntimeError("Austrian GZÜV PKINT original station identifiers not exactly unique ATFW")
    # Original GeoPackage BEGIN_DAT is constant in source (2025-12-09),
    # cannot use as monitoring start or 2022 field sample confirmation.
    raw_begin=gdf["BEGIN_DAT"].dropna().astype(str)
    if raw_begin.nunique()!=1 or raw_begin.iloc[0]!="2025-12-09":
        raise RuntimeError("Original GZÜV station metadata common date changed; do not use as sample year")
    sites=EEA_sites()
    stations={pk:i for i,pk in enumerate(PKINT)}
    matched=[];notmatched=[]
    for z in sites:
        q=z["original_EEA_station_code"]
        if q not in stations:
            notmatched.append(q);continue
        i=stations[q];geo=gdf.geometry.iloc[i]
        m=float(haversine_distance_vector(
              z["source_EEA_2022_lat"],z["source_EEA_2022_lon"],
              np.asarray([geo.y]),np.asarray([geo.x]))[0])
        matched.append({"EEA2022_site_code":q,
            "GZUEV_source_PKG_ID":str(station_PK[i]),
            "GZUEV_source_PKINT_ID":q,
            "EEA2022_to_GZUEV_station_distance_m":m,
            "true_original_same_source_station_identifier":True,
            "EEA2022_benthic_method":z["source_EEA_2022_method"],
            "GZUEV_BEOB_ART_source_value_not_interpreted":str(gdf["BEOB_ART"].iloc[i]),
            "does_not_prove_original_2022_TNTP_same_station":True})
    sourceout,originalgeo=originals(a.outfalls)
    all_lat=np.array([q.geometry.y for q in gdf.itertuples()],float)
    all_lon=np.array([q.geometry.x for q in gdf.itertuples()],float)
    EEA2022_pkint={x["original_EEA_station_code"] for x in sites}
    sourceall={}
    for dist in RADIUS_M:
        n_station_any=0;n_station_BQE=0;ids_all=set();ids_eqr=set()
        for q in sourceout:
            lat=float(q["original_source_outfall_lat"])
            lon=float(q["original_source_outfall_lon"])
            a=haversine_distance_vector(lat,lon,all_lat,all_lon)
            keep=np.flatnonzero(a<=dist)
            if len(keep):n_station_any+=1
            ids_all.update(PKINT.iloc[keep].tolist())
            observed=[pk for pk in PKINT.iloc[keep] if pk in EEA2022_pkint]
            if observed:n_station_BQE+=1
            ids_eqr.update(observed)
        sourceall[str(dist)]={
           "original_AT2022_annualNP_outfalls_near_any_GZUEV_historical_station":n_station_any,
           "distinct_public_GZUEV_station_IDs":len(ids_all),
           "AT2022_outfalls_near_2022_EEA_MONITORED_benthic_site":n_station_BQE,
           "unique_2022_EEA_benthic_field_station_PKINT_IDs_near_outfalls":len(ids_eqr)}
    result={"schema":SCHEMA,
       "official_BMLUK_original_dataset_id":ORIGINAL_DATASET_ID,
       "original_BMLUK_source_GZUEV_GPKG_SHA256":source["original_GZUEV_GeoPackage_source_SHA256"],
       "source_gzuev_registry_station_rows":len(gdf),
       "genuine_source_GZUEV_station_PK_code_field":"PK",
       "genuine_source_GZUEV_EEA_interoperable_station_identity_field":"PKINT",
       "source_GZUEV_PKINT_unique":True,
       "source_GZUEV_BEGIN_DAT_uniform":raw_begin.iloc[0],
       "source_GZUEV_BEGIN_DAT_is_not_2022_original_field_sampling_date":True,
       "EEA2022_original_benthic_AT_source_total":len(sites),
       "EEA2022_68_benthic_sites_matched_exact_PKINT_count":len(matched),
       "EEA2022_BQE_pkint_not_in_GZUEV_original_registry":notmatched,
       "matching_original_PKINT_station_EEA_GZUEV_original_geodesic_offset_m":{
          "max":max(x["EEA2022_to_GZUEV_station_distance_m"] for x in matched) if matched else None,
          "median":float(np.median([x["EEA2022_to_GZUEV_station_distance_m"] for x in matched])) if matched else None},
       "same_2022_station_source_AT145_outfall_proximity":sourceall,
       "published_source_joined_station_metadata_only":matched,
       "no_2022_WISA_benthos_chemistry_sample_source_verified":True,
       "no_numeric_EQR_TP_or_plant_greenness_outcome_fitted":True,
       "original_geo_file_not_rehosted":True}
    (a.output/"GZUEV_EEA2022_AT68_EXACT_PKINT_AT145_STATION_SOURCE_JOIN.json").write_text(
       json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_GZUEV_PKINT_EEA2022_SOURCE_IDENTITY="+json.dumps({
       "official_station_count":len(gdf),
       "EEA2022_AT_site_count":len(sites),
       "PKINT_exact_overlap":len(matched),
       "unmatched_EEA_site_count":len(notmatched),
       "source_station_distance_offsets_m":result["matching_original_PKINT_station_EEA_GZUEV_original_geodesic_offset_m"],
       "AT145_original_source_proximities":sourceall,
       "same_2022_TNTP_MZB_not_yet_verified":True},ensure_ascii=False,sort_keys=True),flush=True)
    for p in (a.output/"ORIGINAL_GZUEV_STATION_SOURCE_DO_NOT_REHOST.zip",gpkg):
        p.unlink(missing_ok=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--outfalls",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
