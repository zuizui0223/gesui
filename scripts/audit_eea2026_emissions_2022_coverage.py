#!/usr/bin/env python3
"""Original Waterbase Emissions 2026: 2022 diffuse agricultural/urban N/P data coverage.

Source-only QA: validates EEA ZIP members against the first source inventory;
reports data grain, source categories, pollutant codes, measurement procedures,
spatialUnitIdentifier scheme, and exact waterbody/RBD/subunit ID coverage among
source-only 2022 BQE sites. NO EQR response values or ecological model accessed.
"""
from __future__ import annotations
import argparse,hashlib,json,re,zipfile
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd
import numpy as np

ROOT="eea_t_waterbase-emissions_p_1977-2025_v01_r00/"
TABLE=ROOT+"Waterbase_v2025_1_T_WISE1_Emissions.csv"
SPATIAL=ROOT+"Waterbase_v2025_1_S_WISE6_SpatialObject_DerivedData.csv"
DEFINITION=ROOT+"Waterbase_v2022_1_dataset_definition_T_WISE1_Emissions.csv"
EXPECTED_BIO_SOURCE_SHA="e12df511d64c4c580ef2422dccdbbd11a1309c84d29d467fd4527206e41fa636"

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def members_checked(archive,first_manifest):
    manifest=json.loads(first_manifest.read_text())
    trusted={x["file"]:(int(x["size_uncompressed"]),str(x["CRC32"]))
         for x in manifest["members"]}
    with zipfile.ZipFile(archive) as z:
        live={x.filename:(int(x.file_size),format(x.CRC,"08x"))
              for x in z.infolist() if not x.is_dir()}
        if live!=trusted:
            raise ValueError("Original EEA emissions ZIP member sizes/CRC differ from first inventory")
        for name in (TABLE,SPATIAL,DEFINITION):
            if name not in live:raise ValueError("Expected EEA emission source table absent")
        hashes={}
        for name in (TABLE,SPATIAL,DEFINITION):
            h=hashlib.sha256()
            with z.open(name) as f:
                for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
            hashes[name]=h.hexdigest()
    return hashes,len(live)

def cats(df,k,top=45):
    return {str(x):int(n) for x,n in df[k].fillna("MISSING").value_counts().head(top).items()}

def distinct(values):
    return {str(s).strip() for s in values.dropna().unique()
            if str(s).strip() and str(s).strip().lower() not in ("nan","none")}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--archive",type=Path,required=True)
    ap.add_argument("--first-inventory",type=Path,required=True)
    ap.add_argument("--source-biology-registry",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    source_hashes,n_files=members_checked(a.archive,a.first_inventory)
    if sha(a.source_biology_registry)!=EXPECTED_BIO_SOURCE_SHA:
        raise ValueError("Original 8245 source-only 2022 EQR monitoring sites changed")
    bqe=pd.read_csv(a.source_biology_registry,
                    dtype={"site_key":str,"waterBodyIdentifier":str},low_memory=False)
    if len(bqe)!=8245:raise ValueError("Wrong source-only biological denominator")
    with zipfile.ZipFile(a.archive) as z:
        emissions=pd.read_csv(z.open(TABLE),dtype="string",low_memory=False)
        places=pd.read_csv(z.open(SPATIAL),dtype="string",low_memory=False)
        dictionary=pd.read_csv(z.open(DEFINITION),dtype="string",low_memory=False)
    required=("countryCode","spatialUnitIdentifier","spatialUnitIdentifierScheme",
              "phenomenonTimeReferencePeriod","observedPropertyDeterminandCode",
              "observedPropertyDeterminandLabel","parameterEmissionsSourceCategory",
              "resultEmissionsValue","resultEmissionsUom","procedureEmissionsMethod")
    if set(required)-set(emissions.columns):
        raise ValueError("Official Emissions source schema changed")
    period=emissions.phenomenonTimeReferencePeriod.fillna("")
    is_2022=period.str.contains("2022",regex=False)
    rows=emissions.loc[is_2022].copy()
    rows["numeric_emission"]=pd.to_numeric(rows.resultEmissionsValue,errors="coerce")
    rows["has_numeric_positive_emission"]=rows.numeric_emission.gt(0)
    pollutant_labels=cats(rows,"observedPropertyDeterminandLabel",70)
    pollutant_codes=cats(rows,"observedPropertyDeterminandCode",70)
    source_categories=cats(rows,"parameterEmissionsSourceCategory",80)
    report={
      "schema":"gesui_original_EEA_waterbase_emissions2026_2022_source_and_grain_qa_v01",
      "original_version":"eea_t_waterbase-emissions_p_1977-2025_v01_r00",
      "source_reference":"https://www.eea.europa.eu/en/datahub/datahubitem-view/5bd950b1-52b9-42f3-b0d7-bf51ff3056b7",
      "retrieved_utc":datetime.now(timezone.utc).isoformat(),
      "ZIP_member_original_inventory_matched":True,
      "source_uncompressed_table_sha256":source_hashes,
      "original_EEA_data_archive_bytes":a.archive.stat().st_size,
      "ZIP_member_count":n_files,
      "original_emissions_rows_all_years":len(emissions),
      "source_period_contains_2022_rows":len(rows),
      "2022_emissions_numeric_positive_rows":int(rows.has_numeric_positive_emission.sum()),
      "2022_countries":cats(rows,"countryCode",50),
      "2022_spatial_identifier_schemes":cats(rows,"spatialUnitIdentifierScheme"),
      "2022_distinct_spatialUnitIdentifier":int(rows.spatialUnitIdentifier.nunique()),
      "2022_source_categories":source_categories,
      "2022_determinands_by_label":pollutant_labels,
      "2022_determinands_by_code":pollutant_codes,
      "2022_reported_methods":cats(rows,"procedureEmissionsMethod"),
      "2022_units":cats(rows,"resultEmissionsUom"),
      "2022_observation_status":cats(rows,"resultObservationStatus"),
      "2022_time_period_string_examples":list(sorted(distinct(rows.phenomenonTimeReferencePeriod)))[:35],
      "original_spatial_reference_rows":len(places),
      "original_spatial_columns":list(places.columns),
      "data_dictionary_original_columns":list(dictionary.columns),
      "data_dictionary_published_source_rows":len(dictionary),
      "data_dictionary_emissions_source_category_definition":{
          str(r["Element name"]):str(r["Element definition"])[:650]
          for _,r in dictionary.iterrows() if "EmissionsSourceCategory" in str(r.get("Element name",""))
      },
    }
    for grain in ("waterBodyIdentifier","subUnitIdentifier","rbdIdentifier"):
        if grain not in places.columns:
            report["spatial_join_"+grain]={"status":"MISSING_SOURCE_COLUMN"};continue
        location=places[places[grain].notna()].copy()
        map_=location[["countryCode","waterBodyIdentifier",grain]].copy()
        map_=map_.dropna(subset=["waterBodyIdentifier"])
        map_["waterBodyIdentifier"]=map_.waterBodyIdentifier.str.strip()
        map_[grain]=map_[grain].str.strip()
        map_=map_.drop_duplicates(["countryCode","waterBodyIdentifier",grain])
        source_units=rows[["countryCode","spatialUnitIdentifier"]].dropna().drop_duplicates()
        source_units["spatialUnitIdentifier"]=source_units.spatialUnitIdentifier.str.strip()
        # exact string and exact country; no fuzzy waterbody/hydro id matching.
        candidates=bqe[["countryCode","waterBodyIdentifier","site_key"]].dropna(subset=["waterBodyIdentifier"])
        candidates=candidates.drop_duplicates("site_key").copy()
        x=candidates.merge(map_,how="left",on=["countryCode","waterBodyIdentifier"],validate="many_to_many")
        y=x.merge(source_units,how="left",left_on=["countryCode",grain],
                  right_on=["countryCode","spatialUnitIdentifier"],indicator=True)
        matched=set(y.loc[y._merge.eq("both"),"site_key"])
        report["spatial_join_"+grain]={
           "source_EEA_spatial_reference_nonmissing_identifiers":int(location[grain].nunique()),
           "original_WISE_BQE_monitoring_site_ids":len(candidates),
           "BQE_stations_where_2022_emission_spatialUnit_matches_exact_id":len(matched),
           "by_country_exact_metadata_joined_2022":{
             str(k):int(v) for k,v in
             y[y.site_key.isin(matched)].drop_duplicates("site_key").countryCode.value_counts().items()},
           "spatial_grain_is_annual_reported_source_unit_not_site_level":True
        }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_2022_DIFFUSE_EMISSIONS_SOURCE_COVERAGE",
          json.dumps({k:v for k,v in report.items() if k in (
            "original_emissions_rows_all_years","source_period_contains_2022_rows",
            "2022_emissions_numeric_positive_rows","2022_countries",
            "2022_spatial_identifier_schemes","2022_source_categories",
            "2022_determinands_by_label","2022_reported_methods","2022_units",
            "spatial_join_waterBodyIdentifier","spatial_join_subUnitIdentifier",
            "spatial_join_rbdIdentifier")},sort_keys=True,ensure_ascii=False),flush=True)
if __name__=="__main__":
    main()
