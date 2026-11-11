#!/usr/bin/env python3
"""2016 biological site/assessment source registry, before 2016 EQR outcomes.

Third replication year selected after 2018/2022 findings; 2016 monitoring-site geography and assessment methods are frozen BEFORE reading any 2016 biological EQR value.
Uses only official EEA Biology2026 original monitoring IDs, geography, year,
biological element code and national method. Numeric result-normalized EQR
is NOT read; this is a source-only sample-selection/coverage step.
"""
from __future__ import annotations
import argparse,hashlib,json,zipfile
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd
import numpy as np

ROOT="eea_t_waterbase-biology-2026_p_1990-2024_v01_r00/"
EQR=ROOT+"Waterbase_v2025_1_T_WISE2_BiologyEQRData.csv"
SPATIAL=ROOT+"Waterbase_v2025_1_S_WISE2_SpatialObject_DerivedData.csv"
SHA_EQR="f11001a57172b618bee95f58fd0ec51c8e71be7b9c01a70927dde21b7b9c01a70927dde21b7ca6e2a9"
CORRECT_SHA_EQR="f11001a57172b618bee95f58fd0ec51c8e71be7b9c01a70927dde21b7ca6e2a9"
CORRECT_SHA_SPATIAL="f6e2448ca2d99f528c21cdb9496673711b1d349e04576cf9024c6165c9a52008"
BQE={"EEA_13-01-4":"benthic_invertebrates",
     "EEA_124-04-9":"phytobenthos",
     "EEA_14-02-8":"fish"}
YEAR="2016"
BBOX=(-12,44,26,64)
def sha_entry(archive,name):
    h=hashlib.sha256()
    with archive.open(name) as f:
        for part in iter(lambda:f.read(1024*1024),b""):h.update(part)
    return h.hexdigest()
def categories(df,column,max_items=40):
    return {str(k):int(v) for k,v in df[column].fillna("MISSING").value_counts().head(max_items).items()}

def freeze_group_folds(frame,code):
    regional=frame[frame.observedPropertyDeterminandBiologyEQRCode.eq(code)].copy()
    by=regional.groupby("region_3deg_block",dropna=False).size().to_dict()
    assignments={};load=[0,0,0,0,0]
    for block,n in sorted(by.items(),key=lambda x:(-x[1],str(x[0]))):
        fold=min(range(5),key=lambda i:(load[i],i))
        assignments[str(block)]=fold+1
        load[fold]+=int(n)
    return assignments,load

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--archive",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    cols=[
      "countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
      "parameterWaterBodyCategory","observedPropertyDeterminandBiologyEQRCode",
      "observedPropertyDeterminandLabel","procedureClassificationSystem",
      "phenomenonTimeReferenceYear","resultObservationStatus"]
    with zipfile.ZipFile(a.archive) as z:
        if sha_entry(z,EQR)!=CORRECT_SHA_EQR or sha_entry(z,SPATIAL)!=CORRECT_SHA_SPATIAL:
            raise ValueError("EEA Biology 2026 source 2016 inner scientific ZIP members changed")
        e=pd.read_csv(z.open(EQR),usecols=cols,dtype="string",low_memory=False)
        locations=pd.read_csv(z.open(SPATIAL),dtype="string",low_memory=False)
    e=e[
        e.phenomenonTimeReferenceYear.eq(YEAR)&
        e.parameterWaterBodyCategory.eq("RW")&
        e.observedPropertyDeterminandBiologyEQRCode.isin(BQE)
    ].copy()
    original_2016_rows=len(e)
    keys=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme"]
    for f in keys:
        e[f]=e[f].astype("string").str.strip()
        locations[f]=locations[f].astype("string").str.strip()
    e["site_key"]=e[keys].fillna("").agg("|".join,axis=1)
    locations["site_key"]=locations[keys].fillna("").agg("|".join,axis=1)
    e.loc[e.monitoringSiteIdentifier.isna()|e.monitoringSiteIdentifier.eq(""),"site_key"]=pd.NA
    locations.loc[locations.monitoringSiteIdentifier.isna()|locations.monitoringSiteIdentifier.eq(""),"site_key"]=pd.NA
    if e.duplicated(["site_key","observedPropertyDeterminandBiologyEQRCode"]).any():
        raise ValueError("Original EQR includes repeated site/BQE/year records")
    for c in ("lat","lon"):
        locations[c]=pd.to_numeric(locations[c],errors="coerce")
    locations["coordinate_string"]=locations[["lat","lon"]].astype("string").fillna("MISSING").agg("|".join,axis=1)
    agg=locations.dropna(subset=["site_key"]).groupby("site_key").agg(
        original_positions=("coordinate_string","nunique"),
        valid_geo=("lat",lambda x:bool(x.between(-90,90).any())))
    usable_keys=set(agg.index[(agg.original_positions==1)&agg.valid_geo])
    loc=locations[locations.site_key.isin(usable_keys)].drop_duplicates("site_key")
    loc=loc[loc.lat.between(BBOX[1],BBOX[3]) & loc.lon.between(BBOX[0],BBOX[2])]
    subset=e.merge(loc[["site_key","countryCode","lat","lon",
        "waterBodyIdentifier","surfaceWaterBodyTypeCode"]],
        on=["site_key","countryCode"],how="inner",validate="many_to_one")
    subset["region_3deg_block"]=(
      np.floor(subset.lat.astype(float)/3).astype(int).astype(str)+"_"+
      np.floor(subset.lon.astype(float)/3).astype(int).astype(str))
    group_reports={}
    for code,label in BQE.items():
        g=subset[subset.observedPropertyDeterminandBiologyEQRCode.eq(code)]
        mapping,foldloads=freeze_group_folds(subset,code)
        subset.loc[subset.observedPropertyDeterminandBiologyEQRCode.eq(code),
                   "pre_2016_source_only_spatial_fold"]=(
            subset.loc[subset.observedPropertyDeterminandBiologyEQRCode.eq(code),
                       "region_3deg_block"].map(mapping))
        group_reports[label]={
          "pre_2016_EQR_outcome_source_only_station_rows":len(g),
          "original_2016_RW_BQE_rows_before_geometry":int(
            (e.observedPropertyDeterminandBiologyEQRCode==code).sum()),
          "monitoring_site_unique_ids":int(g.site_key.nunique()),
          "reporting_countries":int(g.countryCode.nunique()),
          "source_reporting_by_country":categories(g,"countryCode"),
          "source_original_observation_status_codes":categories(g,"resultObservationStatus"),
          "fixed_3degree_spatial_blocks":int(g.region_3deg_block.nunique()),
          "fixed_source_only_fold_1_through_5":{
              str(i+1):int(n) for i,n in enumerate(foldloads)},
          "source_only_block_to_fold":mapping,
          "no_normalized_EQR_fields_loaded":True
        }
        print("GESUI_EEA2016_EQR_PREOUTCOME_SOURCE_SAMPLE",json.dumps(
            {"BQE":label,**{k:v for k,v in group_reports[label].items()
              if k not in ("source_only_block_to_fold",)}},sort_keys=True),flush=True)
    a.out.mkdir(parents=True,exist_ok=True)
    keep=["site_key","countryCode","monitoringSiteIdentifier",
      "monitoringSiteIdentifierScheme","parameterWaterBodyCategory",
      "observedPropertyDeterminandBiologyEQRCode","procedureClassificationSystem",
      "phenomenonTimeReferenceYear","resultObservationStatus",
      "lat","lon","waterBodyIdentifier","surfaceWaterBodyTypeCode",
      "region_3deg_block","pre_2016_source_only_spatial_fold"]
    subset[keep].to_csv(a.out/"2016_FIXED_BQE_MONITORING_SITE_METADATA_NO_EQR.csv",
                        index=False,float_format="%.10g")
    report={
       "schema":"gesui_EEA2016_river_Biology_2026_source_only_BQE_preoutcome_sample_v01",
       "original_uncompressed_EQR_SHA256":CORRECT_SHA_EQR,
       "original_uncompressed_site_geometry_SHA256":CORRECT_SHA_SPATIAL,
       "retrieved_utc":datetime.now(timezone.utc).isoformat(),
       "independent_year_fixed":YEAR,
       "hydrology_extent_geo_bbox":BBOX,
       "original_2016_RW_all_BQE_rows":original_2016_rows,
       "source_only_geolocated_2016_site_element_rows":len(subset),
       "original_2016_BQE_by_group":group_reports,
       "per_group_source_only_deterministic_5fold_fixed":True,
       "pre_2016_Biology_EQR_values_read":False,
       "comparison_prior_2022_BQE_results_not_used_for_2016_site_filter":True,
       "all_source_rows_before_num_EQR_stage_retained":True,
       "caveats":["No EEA 2016 measured annual WWTP P/N mass assigned yet",
          "The original Biological EQR release uses 2025_1 CSV tables within 2026 ZIP",
          "Source-only sample 2016 and prior 2022 may overlap, so year replications are not independent field sites",
          "The 2016 EQR numeric status is unknown, country/source missingness must be counted before model fitting"]
    }
    (a.out/"2016_BQE_PREOUTCOME_SOURCE_QC.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_2016_TEMPORAL_BIOLOGICAL_SOURCE_ONLY_COMPLETE",
        json.dumps({"all_site_group_rows":len(subset),
          "original_2016_BQE_no_numeric_reuse":True,
          "distinct_original_source_countries":{
            p:r["reporting_countries"] for p,r in group_reports.items()}},sort_keys=True),flush=True)

if __name__=="__main__":main()
