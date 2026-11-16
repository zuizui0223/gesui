#!/usr/bin/env python3
"""Post-2026-Part2 original source-only chemistry availability audit.

Read *only* previously source-locked official Part2 3525 EEA site-years
artifact, not source science responses. Diagnose whether the negative
same-station positive-report gate is due to chemistry, geography,
sampling quality or missing Part1. Outputs no observed TNTP numeric levels.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path

EXPECTED_ZIP_SHA="315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a"
SOURCE_SCHEMA="GESUI_20261116_EEA_ICM2026_OFFICIAL_CSV_P2_SOURCE_ONLY"

def main(args):
    qa=json.loads((args.input/"EEA2026_EXACT_STATION_YEAR_TNTP_SOURCE_QA.json").read_text())
    if qa.get("schema")!=SOURCE_SCHEMA:
        raise RuntimeError("Official 2026 original Part2 source schema altered")
    if qa.get("official_publisher_2026_original",{}).get("official_2026_original_part2_ZIP_SHA256")!=EXPECTED_ZIP_SHA:
        raise RuntimeError("Official Part2 original source ZIP SHA mismatch")
    file=args.input/"EEA2026_ONLY_MATCHED_3525_STATION_YEAR_TNTP_SOURCE.csv"
    with file.open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=3525 or len({(r["monitoringSiteIdentifier"],r["monitoringSiteIdentifierScheme"],r["year"]) for r in rows})!=3525:
        raise RuntimeError("Original source 3525 exact monitoringSiteIdentifier+scheme+year drift")
    if sum(r["strict_same_station_year_TNTP"]=="True" for r in rows)!=20:
        raise RuntimeError("Original strict 20 site/year total-N and total-P source changed")
    expected=qa.get("actual_annual_CSV_source",{})
    stage=expected.get("original_exact_station_source_gates",{})
    failures=expected.get("source_ineligible_reasons",{})
    row_year=Counter()
    together=Counter()
    reported_both=Counter()
    tn_year=Counter()
    tp_year=Counter()
    group=defaultdict(set)
    for r in rows:
        year=r["year"]
        row_year[year]+=1
        if r["strict_single_TN"]=="True":tn_year[year]+=1
        if r["strict_single_TP"]=="True":tp_year[year]+=1
        if r["strict_same_station_year_TNTP"]=="True":
            together[year]+=1
            group[(r["country"],year)].add((r["h"],r["period"]))
            if r["original_site_reported_upstream_annual_N_or_P_positive"]=="True":
                reported_both[year]+=1
    out={
      "schema":"GESUI_20261116_OFFICIAL_EEA2026_PART2_SOURCE_NEGATIVE_DIAGNOSTIC_V1",
      "official_2026_publisher_fullPart2_ZIP_SHA256":EXPECTED_ZIP_SHA,
      "original_CSV_scanned_rows":expected.get("original_total_source_csv_rows_scanned"),
      "original_source_2018_2022_record_rows":expected.get("source_selected_2018_2022_rows_by_year"),
      "original_3525_biology_station_years_by_era":dict(row_year),
      "source_match_stages_from_full_CSV":stage,
      "source_ineligible_qualified_TNTP_reasons_from_full_CSV":failures,
      "native_tn_year":dict(tn_year),"native_tp_year":dict(tp_year),
      "native_both_year":dict(together),"native_both_report_positive_year":dict(reported_both),
      "countries_and_h7periods_with_any_strict_same_site_TNTP":{
        country+"_"+year:len(h) for (country,year),h in sorted(group.items())},
      "no_2026_DisaggregatedData_part1_parsed":True,
      "no_wastewater_causal_mediation_or_biological_response_fitted":True,
      "no_report_absence_as_pristine_zero":True,
      "source_part2_only_cannot_prove_missing_NP_samples":True}
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"EEA2026_PART2_3525_SITE_YEAR_FAILURE_DIAGNOSTIC.json").write_text(
        json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_EEA2026_PART2_SITE_YEAR_DIAGNOSTIC="+json.dumps(out,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--input",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
