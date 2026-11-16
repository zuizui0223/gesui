#!/usr/bin/env python3
"""Source-first EEA 2026 same-station/year TNTP annual from official Parts 1+2.

Only after each original publisher Part is SHA-verified and independently
matched to genuine EEA Biology site 3525 rows. Part1 valid annual arithmetic
mean from >=3 distinct water-sample results. Part2 valid publisher annual
site mean with >=3 publisher samples. Never count as independent biology.
For each station/year/nutrient prefer actual valid Part1 over Part2; *no*
averaging across datasets, no imputation of nonmatches, no wastewater causality.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path

P1_SHA="ff779521749b40b952cce35a4671240c3b6c47d87f32a8fef9240a4ac503a082"
P2_SHA="315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a"
P1_SCHEMA="GESUI_20261116_OFFICIAL_EEA2026_DISAGGREGATED_PART1_ORIGINAL_MATCHED_EEA_STATION_TNTP_V1"
P2_SCHEMA="GESUI_20261116_EEA_ICM2026_OFFICIAL_CSV_P2_SOURCE_ONLY"
OUTPUT_SCHEMA="GESUI_20261116_EEA2026_PART1_PART2_ORIGINAL_3525_NATIVE_SITE_YEAR_TNTP_UNION_GATE_V1"
LIMITS={"H7period":30,"H7":20,"H6":10,"country":3}

def get_rows(p):
    with p.open(newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))
def b(v):
    if str(v).lower() not in ("true","false"):
        raise RuntimeError("Source original site-year has unexpected boolean "+str(v))
    return str(v).lower()=="true"
def valid_num(v):
    if v in (None,"","nan","None"):return None
    z=float(v)
    if not math.isfinite(z) or z<0:raise RuntimeError("Original EEA source chemistry nonfinite/negative")
    return z
def main(a):
    p1=json.loads((a.part1/"EEA2026_PART1_EXACT_STATION_YEAR_SOURCE_QA.json").read_text())
    p2=json.loads((a.part2/"EEA2026_EXACT_STATION_YEAR_TNTP_SOURCE_QA.json").read_text())
    if p1.get("schema")!=P1_SCHEMA or p2.get("schema")!=P2_SCHEMA:
        raise RuntimeError("Original 2026 source Part1/Part2 schema unverified")
    if p1.get("original_2026_public_part1_source",{}).get("original_part1_SHA256")!=P1_SHA:
        raise RuntimeError("2026 source Part1 Deflate64 ZIP SHA mismatch")
    if p2.get("official_publisher_2026_original",{}).get("official_2026_original_part2_ZIP_SHA256")!=P2_SHA:
        raise RuntimeError("2026 source Part2 original ZIP SHA mismatch")
    one=get_rows(a.part1/"EEA2026_ONLY_MATCHED_3525_STATION_YEAR_TNTP_PART1.csv")
    two=get_rows(a.part2/"EEA2026_ONLY_MATCHED_3525_STATION_YEAR_TNTP_SOURCE.csv")
    if len(one)!=3525 or len(two)!=3525:
        raise RuntimeError("Original 3525 genuine EEA biological station-year source changed")
    key1=lambda z:(z["original_site"],z["original_scheme"],z["year"])
    key2=lambda z:(z["monitoringSiteIdentifier"],z["monitoringSiteIdentifierScheme"],z["year"])
    a1={key1(z):z for z in one}
    a2={key2(z):z for z in two}
    if len(a1)!=3525 or len(a2)!=3525 or set(a1)!=set(a2):
        raise RuntimeError("Source Part1 and Part2 do not match exact original 3525 EQR site/year population")
    statuses=Counter();groups=defaultdict(set);overlap={}
    source_rows=[]
    for k in sorted(a1):
        q=a1[k];r=a2[k]
        if (q["h7"],q["h6"],q["country"])!=(r["h"],r["h6"],r["country"]):
            raise RuntimeError("Genuine EEA source H7/H6/country geographic mismatch")
        if b(q["original_upstream_reported_annual_N_or_P_positive"])!=b(r["original_site_reported_upstream_annual_N_or_P_positive"]):
            raise RuntimeError("Original EEA source positive site-year label inconsistency")
        line={"h":q["h7"],"h6":q["h6"],"country":q["country"],"year":q["year"],"period":r["period"],
              "original_monitoringSiteIdentifier":k[0],
              "original_monitoringSiteIdentifierScheme":k[1],
              "original_reported_upstream_positive":b(q["original_upstream_reported_annual_N_or_P_positive"])}
        for nutrient in ("TN","TP"):
            qok=b(q["valid_annual_"+nutrient])
            rok=b(r["strict_single_"+nutrient])
            p1val=valid_num(q["EEA2026_part1_"+nutrient+"_mean_mg"+("N" if nutrient=="TN" else "P")+"_L"]) if qok else None
            p2val=valid_num(r["source_original_"+nutrient+"_mg"+("N" if nutrient=="TN" else "P")+"_L"]) if rok else None
            if qok and p1val is None:raise RuntimeError("Genuine Part1 quality flag claims missing chemical value")
            if rok and p2val is None:raise RuntimeError("Genuine Part2 quality flag claims missing chemical value")
            if qok and rok:
                statuses[nutrient+"_both_tables_valid"]+=1
                diff=abs(p1val-p2val)
                # Diagnostic only: do not filter by observed agreement.
                if diff>0.20*max(1e-8,abs(p1val),abs(p2val)):
                    statuses[nutrient+"_crosspart_relative_discordance_gt20pct"]+=1
            if qok:
                v=p1val;which="Part1_individual_mean"
            elif rok:
                v=p2val;which="Part2_annual_original"
            else:
                v=None;which="NOT_VERIFIED"
            line[nutrient+"_mg_element_per_L"]=v
            line[nutrient+"_official_source_table"]=which
            if v is not None:
                statuses[nutrient+"_union_valid"]+=1
                groups[nutrient+"_H7period"].add((line["h"],line["period"]))
        complete=all(line[t+"_mg_element_per_L"] is not None for t in ("TN","TP"))
        line["strict_same_station_year_both_TNTP"]=complete
        if complete:
            statuses["TNTP_any_union"]+=1
            groups["all_H7period"].add((line["h"],line["period"]))
            if line["original_reported_upstream_positive"]:
                statuses["TNTP_and_report_positive"]+=1
                groups["positive_H7period"].add((line["h"],line["period"]))
                groups["positive_H7"].add(line["h"])
                groups["positive_H6"].add(line["h6"])
                groups["positive_country"].add(line["country"])
        source_rows.append(line)
    support={k:len(groups["positive_"+k]) for k in LIMITS}
    passgate=all(support[k]>=LIMITS[k] for k in LIMITS)
    out={
      "schema":OUTPUT_SCHEMA,
      "original_2026_publisher_part1_SHA256":P1_SHA,
      "original_2026_publisher_part2_SHA256":P2_SHA,
      "genuine_EEA_biology_native_station_years":3525,
      "original_station_year_source_positive":132,
      "valid_source_by_table_and_union_counts":dict(statuses),
      "strict_same_station_year_both_TNTP_H7periods_total":len(groups["all_H7period"]),
      "strict_report_positive_both_source_support":support,
      "original_precommitted_source_gate_limits":LIMITS,
      "source_gate_pass":passgate,
      "source_disposition":"SOURCE_READY_FOR_NONCAUSAL_ASSOCIATION_MODEL_IF_OTHER_CONTROLS_PASS"
           if passgate else "STOP_SAME_SITE_YEAR_REPORT_POSITIVE_MEASURED_TNTP_SUPPORT",
      "Part1_priority_per_source_nutrient":True,
      "cross_part_diff_flag_not_used_for_after_outcome_selection":True,
      "source_part1_Part2_complementary_not_independent_biological_studies":True,
      "report_unmatched_never_treated_as_undischarged":True,
      "measured_ambient_TNTP_does_not_identify_sewage_origin":True,
      "no_new_EQR_NUTRI_HNS_or_mediation_effect_fitted":True}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"EEA2026_PART1_PART2_ORIGINAL_3525_SITE_YEAR_UNION_SOURCE_GATE.json").write_text(
        json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    with (a.output/"EEA2026_ONLY_3525_NATIVE_SITE_YEAR_SOURCE_UNION.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(source_rows[0]))
        w.writeheader();w.writerows(source_rows)
    print("GESUI_EEA2026_PART1_PART2_SOURCE_UNION="+json.dumps(out,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--part1",required=True,type=Path)
    p.add_argument("--part2",required=True,type=Path);p.add_argument("--output",required=True,type=Path);main(p.parse_args())
