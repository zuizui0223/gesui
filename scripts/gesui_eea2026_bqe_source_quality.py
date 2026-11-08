#!/usr/bin/env python3
"""EEA Waterbase Biology 2026 source support gate, 2022 river monitoring.

Read official source-locked EQR/spatial CSVs. Count eligible original monitoring
stations and BQE reporting methods BEFORE reading nutrient-source values or
fitting EQR outcomes. No causal WWTP effect, no GBIF taxa.
"""
from __future__ import annotations
import argparse,hashlib,io,json,math,zipfile
from datetime import datetime,timezone
from pathlib import Path
from collections import Counter
import numpy as np,pandas as pd

SOURCE_SHA="e621b1f5c43e9f2005d3a6fa40adefff6fd30d7030e76c4c8801d7717194769f"
SOURCENAME="eea_t_waterbase-biology-2026_p_1990-2024_v01_r00"
EQR=SOURCENAME+"/Waterbase_v2025_1_T_WISE2_BiologyEQRData.csv"
SPATIAL=SOURCENAME+"/Waterbase_v2025_1_S_WISE2_SpatialObject_DerivedData.csv"
DEFINITIONS=SOURCENAME+"/Waterbase_v2022_1_dataset_definition_T_WISE2_BiologyEQRData.csv"
PRIMARY_YEAR="2022"

def sha_file(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def composite(frame,kind):
    cols=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme"]
    if not set(cols)<=set(frame.columns):
        raise ValueError("Missing official monitoring site key columns")
    frame=frame.copy()
    for c in cols:frame[c]=frame[c].astype("string").str.strip()
    frame["monitor_key"]=frame[cols].fillna("").agg("|".join,axis=1)
    invalid=frame.monitoringSiteIdentifier.isna()|frame.monitoringSiteIdentifier.eq("")
    frame.loc[invalid,"monitor_key"]=pd.NA
    return frame

def top_counts(df,col,max_items=80):
    if col not in df.columns:return {}
    return {str(k):int(v) for k,v in df[col].fillna("MISSING").value_counts().head(max_items).items()}

def inspect(archive,output):
    if sha_file(archive)!=SOURCE_SHA:raise ValueError("Official EEA Biology 2026 ZIP changed")
    with zipfile.ZipFile(archive) as z:
        index=set(z.namelist())
        if not {EQR,SPATIAL,DEFINITIONS}<=index:
            raise ValueError("Official EQR CSV, spatial CSV or definitions are missing")
        frame=pd.read_csv(z.open(EQR),dtype="string",low_memory=False)
        points=pd.read_csv(z.open(SPATIAL),dtype="string",low_memory=False)
        definitions=pd.read_csv(z.open(DEFINITIONS),dtype="string")
    original_n=len(frame);spatial_n=len(points)
    required=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
        "parameterWaterBodyCategory","observedPropertyDeterminandLabel",
        "observedPropertyDeterminandBiologyEQRCode","procedureClassificationSystem",
        "phenomenonTimeReferenceYear","resultNormalisedEQRValue","resultObservationStatus"]
    if set(required)-set(frame.columns):raise ValueError("Original EEA EQR fields changed")
    frame=composite(frame,"EQR")
    points=composite(points,"spatial")
    for c in ("lat","lon"):
        points[c]=pd.to_numeric(points[c],errors="coerce")
    points["valid_coordinate"]=points.lat.between(-90,90)&points.lon.between(-180,180)
    points["location_sig"]=points[["lat","lon"]].astype("string").agg("|".join,axis=1)
    p=points.dropna(subset=["monitor_key"])
    sizes=p.groupby("monitor_key").agg(n_position=("location_sig","nunique"),
                                     with_coords=("valid_coordinate","any"),
                                     waterbody_count=("waterBodyIdentifier","nunique"))
    reliable=set(sizes.loc[sizes.n_position.eq(1)&sizes.with_coords,"monitor_key"])
    points_unique=p[p.monitor_key.isin(reliable)].drop_duplicates("monitor_key")
    year=frame[frame.phenomenonTimeReferenceYear.eq(PRIMARY_YEAR)].copy()
    year["normalizedEQR_numeric"]=pd.to_numeric(year.resultNormalisedEQRValue,errors="coerce")
    allriver=year[year.parameterWaterBodyCategory.astype(str).str.upper().isin(
        ("RW","RIVER","RIVERS","R"))].copy()
    raw=year
    # A strict river-category assumption is a *source* hypothesis. If water
    # code unexpectedly differs, report all categories and do not silently
    # reclassify nonriver water bodies as rivers.
    loc=allriver.copy()
    loc["spatial_ID_covered"]=loc.monitor_key.isin(reliable)
    loc["normalized_EQR_parseable"]=np.isfinite(loc.normalizedEQR_numeric)
    valid=loc[loc.spatial_ID_covered & loc.normalized_EQR_parseable].copy()
    unique_count=valid.monitor_key.nunique()
    duplicates=int(valid.duplicated(["monitor_key","observedPropertyDeterminandBiologyEQRCode"]).sum())
    detail=[]
    for (code,label),g in valid.groupby(
        ["observedPropertyDeterminandBiologyEQRCode","observedPropertyDeterminandLabel"],
        dropna=False):
        detail.append({
            "biology_code":str(code),"biology_label":str(label),
            "original_2022_rows":len(g),
            "unique_monitoring_sites":g.monitor_key.nunique(),
            "countries_reporting":g.countryCode.nunique(),
            "different_classification_systems":g.procedureClassificationSystem.nunique(),
            "rows_with_normalized_EQR":int(g.normalized_EQR_parseable.sum()),
            "unique_site_biology_code_pairs":len(g.drop_duplicates(["monitor_key","observedPropertyDeterminandBiologyEQRCode"])),
            "status_code_counts":top_counts(g,"resultObservationStatus",15),
        })
    detail=sorted(detail,key=lambda x:(-x["unique_monitoring_sites"],x["biology_code"]))
    summary={
      "schema":"gesui_eea_biology2026_EQR_2022_river_site_coverage_preoutcome_v01",
      "original_source_url":"https://www.eea.europa.eu/en/datahub/datahubitem-view/7bd881f0-59ca-4f0a-a43c-7e0cd0c66d6c",
      "original_archive_sha256":SOURCE_SHA,
      "source_EQR_CSV_name":EQR,"source_spatial_CSV_name":SPATIAL,
      "timestamp_utc":datetime.now(timezone.utc).isoformat(),
      "2022_calendar_year_locked_before_EQR_modeling":PRIMARY_YEAR,
      "EQR_original_all_years_rows":original_n,
      "original_spatial_metadata_rows":spatial_n,
      "distinct_spatial_monitoring_IDs":int(p.monitor_key.nunique()),
      "distinct_one_coordinate_spatial_IDs":len(reliable),
      "spatial_IDs_with_conflicting_coordinate_sets":int(sizes.n_position.gt(1).sum()),
      "2022_EQR_all_categories_rows":len(year),
      "2022_category_codes_in_original":top_counts(year,"parameterWaterBodyCategory"),
      "2022_strict_river_RW_RIVER_R_rows":len(allriver),
      "2022_river_rows_with_unambiguous_coordinates":int(loc.spatial_ID_covered.sum()),
      "2022_river_rows_with_numeric_normalized_EQR":int(loc.normalized_EQR_parseable.sum()),
      "2022_river_rows_with_both_site_geometry_and_EQR":len(valid),
      "2022_river_unique_monitoring_sites_with_both":unique_count,
      "2022_river_site_biology_duplicate_records":duplicates,
      "2022_river_biology_quality_elements":detail,
      "2022_qualified_spatial_countries":top_counts(valid,"countryCode",45),
      "2022_source_quality_status_codes":top_counts(raw,"resultObservationStatus",35),
      "geographic_coding_overlap_has_not_been_validated_against_HydroRIVERS":True,
      "actual_EQR_numeric_values_used_for_EEA_WWTP_impact_model":False,
      "reported_2022_eea_WWTP_loads_not_linked_to_biology":True,
      "limitations":[
        "Monitor_key requires exact country+monitoringSiteIdentifier+monitoringSiteIdentifierScheme",
        "SpatialObject 2025_1 records in 2026 official archive do not guarantee complete year-2022 geolocation",
        "National classification procedures, waterbody types and monitoring design affect normalized EQR comparability",
        "Only actual 2022 class and EQR source support identified; no upstream wastewater association tested",
        "River code use RW/RIVER/R may be adjusted only after checking official metadata, before reading EQR effect values",
        "The presence of a biological monitoring site is not a statistically independent or wastewater-exposed river catchment"
      ]
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_EEA2026_BQE_2022_SOURCE_COUNTS",json.dumps({
      k:summary[k] for k in ("EQR_original_all_years_rows","original_spatial_metadata_rows",
       "distinct_one_coordinate_spatial_IDs","2022_EQR_all_categories_rows",
       "2022_category_codes_in_original","2022_strict_river_RW_RIVER_R_rows",
       "2022_river_rows_with_both_site_geometry_and_EQR",
       "2022_river_unique_monitoring_sites_with_both")},sort_keys=True),flush=True)
    print("GESUI_EEA2026_BQE_2022_LABEL_SUPPORT",json.dumps(detail[:40],ensure_ascii=False),flush=True)
    return summary

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--archive",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    inspect(a.archive,a.out)
