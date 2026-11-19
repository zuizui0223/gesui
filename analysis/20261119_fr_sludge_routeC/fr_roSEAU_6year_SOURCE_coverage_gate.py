#!/usr/bin/env python3
"""Six-year ROSEAU sludge SOURCE gate: coverage and record grain, not biota.

Only validated government source-metadata JSON. Never infer application
nonexistence from missing records, never count assay rows as fields or events.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from datetime import datetime,timezone

YEARS=tuple(str(y) for y in range(2017,2023))
SCHEMA="gesui_FR_ROSEAU_sludge_routeC_2017to22_SOURCE_only_record_grain_20261119_v1"

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def summarize(paths):
    if len(paths)!=6:raise ValueError("Six original frozen source receipts required")
    byyear={}
    shas={}
    for p in paths:
        data=json.loads(Path(p).read_text())
        if data.get("schema")!=SCHEMA or not data.get("no_biotic_EQR_NUTRI_HNS_values_read"):
            raise ValueError("Unexpected source report/outcome exposure")
        year=str(data.get("year"))
        if year not in YEARS or year in byyear:raise ValueError("Missing or duplicated year")
        byyear[year]=data.get("result",{})
        shas[year]=sha(p)
    if set(byyear)!=set(YEARS):raise ValueError("Full source year census missing")
    good={y:byyear[y] for y in YEARS if byyear[y].get("status")=="ORIGINAL_CSV_RETRIEVED_QA"}
    early=any(y in good for y in YEARS[:3])
    late=any(y in good for y in YEARS[3:])
    valid_codes=0
    rows=0
    stats=[]
    for y in YEARS:
        a=byyear[y]
        v=a.get("source_row_grain",{})
        n=v.get("raw_records")
        if n is not None:rows+=int(n)
        stats.append({
           "year":y,"status":a.get("status"),
           "raw_assay_rows":n,"unique_STEU":v.get("official_unique_stations_nonempty"),
           "unique_station_commune_literal_pairs":v.get("official_distinct_station_commune_literal_pairs"),
           "unique_commune_literals":v.get("distinct_source_commune_literal_strings"),
           "unique_confirmed_INSEE_shape_from_unvalidated_strings":v.get("distinct_explicit_INSEE_codes_only_not_validated_against_official_INSEE_geometry"),
           "positive_station_commune_literal_pairs":v.get("distinct_positive_qty_station_commune_pairs"),
           "year_source_CSV_SHA":a.get("raw_csv_sha256"),
           "original_2018_sha_matches_archive":a.get("earlier_2018_sha_matched"),
        })
    # INSEE-code format does not establish a true application parcel or an
    # independent spreading event; never aggregate year sums as locations.
    return {
      "source_year_reports":stats,
      "years_retrieved":sorted(good),
      "years_stopped":sorted(set(YEARS)-set(good)),
      "record_rows_sum_across_years_not_unique_events":rows,
      "independent_communes_ACROSS_years_not_computable_from_counts":True,
      "each_period_at_least_one_public_export_retrieved":early and late,
      "actual_spreading_event_year_NOT_confirmed_by_current_CSV_columns":True,
      "INSEE_commune_code_verification_requires_external_official_mapping":True,
      "france_H7_source_GO":False,
      "reason_no_H7_GO":"Public exporter has station/assay grain and application commune strings, not verified dated field-application records or independent municipality INSEE crosswalk.",
      "year_receipt_sha256":shas,
      "old_2018_water_2017_plant_NUTS2_2x2_and_JRC_native_source_STOP_remain":True
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--input",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    v=a.parse_args()
    files=sorted(v.input.rglob("FR_ROSEAU_BOUES_????_SOURCE_GRAIN_QA.json"))
    gate=summarize(files)
    v.out.mkdir(parents=True,exist_ok=True)
    output={"schema":"gesui_FR_ROSEAU_2017to22_source_grain_coverage_GATE_v1",
            "UTC":datetime.now(timezone.utc).isoformat(),
            "source_only":gate,"no_biological_outcomes_read":True}
    (v.out/"FR_ROSEAU_2017TO22_SOURCE_COVERAGE_AND_SPREADING_ID_GO_STOP.json").write_text(
       json.dumps(output,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_FR_ROSEAU_SIX_YEAR_SOURCE_COVERAGE",json.dumps({
      k:gate[k] for k in ("years_retrieved","years_stopped",
                         "each_period_at_least_one_public_export_retrieved",
                         "france_H7_source_GO","record_rows_sum_across_years_not_unique_events")
    },ensure_ascii=False),flush=True)
if __name__=="__main__":main()
