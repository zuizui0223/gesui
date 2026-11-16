#!/usr/bin/env python3
"""Official 2026 full Part1+Part2 source-only TN-only / TP-only opportunity.

Secondary POSTHOC study after Part1+Part2 TNTP BOTH source gate STOP.
Preserve original 3525 genuine EEA EQR station+scheme+year and original
2018/22 upstream WWTP *reporting* state; no EQR/HNS/EIV/W observations.
No lower limits; no invented no-sewage control, no causal mediation.
"""
from __future__ import annotations
import argparse,csv,json,math
from collections import Counter,defaultdict
from pathlib import Path

ORIGINAL_SCHEMA="GESUI_20261116_EEA2026_PART1_PART2_ORIGINAL_3525_NATIVE_SITE_YEAR_TNTP_UNION_GATE_V1"
P1_SHA="ff779521749b40b952cce35a4671240c3b6c47d87f32a8fef9240a4ac503a082"
P2_SHA="315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a"
THRESH={"H7period":30,"H7":20,"H6":10,"countries":3}
LOCAL_MIN=10
YEARS=("2018","2022")

def yes(v):
    if v not in ("True","False","true","false"):
        raise RuntimeError("Original site-year source report state invalid: "+repr(v))
    return v.lower()=="true"

def data(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def main(a):
    source=json.loads((a.union/"EEA2026_PART1_PART2_ORIGINAL_3525_SITE_YEAR_UNION_SOURCE_GATE.json").read_text())
    if (source.get("schema")!=ORIGINAL_SCHEMA or
        source.get("original_2026_publisher_part1_SHA256")!=P1_SHA or
        source.get("original_2026_publisher_part2_SHA256")!=P2_SHA or
        source.get("source_gate_pass") is not False):
        raise RuntimeError("Original prespecified TNTP site/year STOP and publisher source SHAs required")
    rr=data(a.union/"EEA2026_ONLY_3525_NATIVE_SITE_YEAR_SOURCE_UNION.csv")
    if len(rr)!=3525:raise RuntimeError("Original source station-year count changed")
    keys={(x["original_monitoringSiteIdentifier"],x["original_monitoringSiteIdentifierScheme"],x["year"]) for x in rr}
    if len(keys)!=3525:
        raise RuntimeError("Original native station+scheme+year duplicate; no pseudo-observations")
    if len({(x["h"],x["period"]) for x in rr})!=454:
        raise RuntimeError("Original source 454 H7-period geography changed")
    count=Counter()
    result={}
    for kind in ("TN","TP","TNTP"):
        pos=[]
        all_valid=[]
        group_positive=defaultdict(set)
        peryear=Counter()
        for x in rr:
            pTN=x["TN_mg_element_per_L"] not in ("",None)
            pTP=x["TP_mg_element_per_L"] not in ("",None)
            good=(pTN if kind=="TN" else pTP if kind=="TP" else (pTN and pTP))
            if not good:continue
            for nutrient in (("TN","TP") if kind=="TNTP" else (kind,)):
                val=x[nutrient+"_mg_element_per_L"]
                if val in ("",None):
                    raise RuntimeError("Site-year false chemical indicator")
                number=float(val)
                if not math.isfinite(number) or number<0:
                    raise RuntimeError("Invalid published total nutrient source value")
            all_valid.append(x)
            if yes(x["original_reported_upstream_positive"]):
                pos.append(x)
                group_positive[x["h6"]].add(x["h"])
                peryear[x["year"]]+=1
        geography={
            "H7period":len({(x["h"],x["period"]) for x in pos}),
            "H7":len({x["h"] for x in pos}),
            "H6":len({x["h6"] for x in pos}),
            "countries":len({x["country"] for x in pos})}
        local=len([h6 for h6,hs in group_positive.items() if len(hs)>=2])
        sourcepass=all(geography[k]>=THRESH[k] for k in THRESH)
        true_parentpass=local>=LOCAL_MIN
        result[kind]={
           "strict_all_EEA_biology_native_station_years_with_chemical":len(all_valid),
           "strict_with_upstream_positive_annual_source_report_station_years":len(pos),
           "positive_report_H7period_H7_H6_countries":geography,
           "original_report_positive_year_site_counts":dict(peryear),
           "same_true_H6_parent_with_two_distinct_positive_measured_H7":local,
           "original_same_source_replication_minimum":THRESH,
           "original_true_H6_multiH7_source_minimum":LOCAL_MIN,
           "general_source_site_gate_pass":bool(sourcepass),
           "true_parent_positive_vs_positive_gate_pass":bool(true_parentpass),
           "source_only_not_direct_WWTP_effect":True}
    if result["TNTP"]["strict_all_EEA_biology_native_station_years_with_chemical"]!=101 or \
        result["TNTP"]["strict_with_upstream_positive_annual_source_report_station_years"]!=12:
        raise RuntimeError("Cannot independently reproduce original TNTP union result")
    if result["TN"]["strict_all_EEA_biology_native_station_years_with_chemical"]!=186 or \
       result["TP"]["strict_all_EEA_biology_native_station_years_with_chemical"]!=262:
        raise RuntimeError("Original TN/TP totals differ from source-verified union")
    out={"schema":"GESUI_20261116_EEA2026_TOTALN_TOTALP_SEPARATE_SOURCE_REPLICATION_POSTHOC_V1",
         "original_Part1_2026_ZIP_SHA256":P1_SHA,
         "original_Part2_2026_ZIP_SHA256":P2_SHA,
         "original_exact_EEA_biology_station_years":len(rr),
         "study_target_2018_and_2022":True,
         "tests_fixed_before_this_single_nutrient_source_result":["TN","TP","TNTP"],
         "separate_tests_are_posthoc_after_previous_TNTP_STOP":True,
         "source_coverage_by_total_nutrient":result,
         "original_TNTP_MAIN_SOURCE_GATE_STOP_unmodified":True,
         "report_unmatched_does_not_imply_no_wastewater":True,
         "2026_official_ambient_TNTP_is_NOT_outfall_source_identification":True,
         "no_biological_outcomes_or_causal_mediation_fitted":True}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"EEA2026_TOTALN_TOTALP_SEPARATE_SOURCE_OPPORTUNITY.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_EEA2026_TN_TP_SEPARATE_SOURCE_GATE="+json.dumps(out,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--union",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
