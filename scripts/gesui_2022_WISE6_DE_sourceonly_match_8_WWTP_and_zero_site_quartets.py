#!/usr/bin/env python3
"""DE-only balanced 2022 WISE6 TP original monitor BEFORE/AFTER WWTP path matched
with strictly no-registered-source monitor path, selected by source-only
hydro-distance and 2degree region, before any numeric river TP outcome.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

MATCH_N=8
MIN_SAME_BLOCK=4
MIN_DISTINCT_BLOCKS=6
GEO_MISMATCH_COST=0.4
EXPECTED_DE_SOURCE_PAIRS=33
EXPECTED_DE_ZERO_PAIRS=40

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def sites(x):return (x["upstream_station_SHA256"],x["downstream_station_SHA256"])
def classified(x):
    s=x["real_exact_source_between"]
    if s["strictly_between_total_original_Fw_outfalls_n"]>0:return "SOURCE_INTERVENES"
    if s["strictly_between_total_original_Fw_outfalls_n"]==0 and s["unresolved_same_reach_or_first_join_into_downstream_reach_source_n"]==0:return "CLEAN_ZERO_REGISTERED_SOURCE"
    return "UNKNOWN_INTERVENTION"

def score(s,c):
    d=abs(math.log1p(float(s["full_reach_path_to_downstream_reach_start_km"]))-
          math.log1p(float(c["full_reach_path_to_downstream_reach_start_km"])))
    g=s["upstream_country2degree_block"]==c["upstream_country2degree_block"]
    return d+(0 if g else GEO_MISMATCH_COST),d,g

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--national-source-between",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    src=json.loads(a.national_source_between.read_text())
    if not src["original_national_WWTP_between_8_and_clean_zero_8_source_gate_pass"] or not src["replayed_all_original_79_photo_cell_2022_WWTP_source_interposition_records_exactly"]:
        raise ValueError("Original SOURCE-ONLY national stage B not completed and reproducible")
    base=[r for r in src["all_national_2022_same_river_source_only_pairs_with_verified_outlet_between_reaches"] if r["countryCode"]=="DE"]
    sources=[r for r in base if classified(r)=="SOURCE_INTERVENES"]
    controls=[r for r in base if classified(r)=="CLEAN_ZERO_REGISTERED_SOURCE"]
    if len(sources)!=EXPECTED_DE_SOURCE_PAIRS or len(controls)!=EXPECTED_DE_ZERO_PAIRS:
        raise ValueError("Original exact DE registered WWTP-between source/zero strata changed")
    candidates=[]
    for s in sources:
        for c in controls:
            four=sites(s)+sites(c)
            if len(set(four))<4:continue
            value,dist,same=score(s,c)
            candidates.append({
              "value":value,"length_log_difference_abs":dist,"matched_2degree_upstream_block":same,
              "source":s,"registered_source_zero":c,
              "tiebreak":hashlib.sha256("|".join(four).encode()).hexdigest()})
    candidates.sort(key=lambda x:(x["value"],x["tiebreak"]))
    used=set();selected=[]
    for x in candidates:
        four=sites(x["source"])+sites(x["registered_source_zero"])
        if any(k in used for k in four):continue
        selected.append(x)
        used.update(four)
        if len(selected)>=MATCH_N:break
    blockset={x["source"]["upstream_country2degree_block"] for x in selected}|{
        x["registered_source_zero"]["upstream_country2degree_block"] for x in selected}
    same=sum(x["matched_2degree_upstream_block"] for x in selected)
    gate=len(selected)==MATCH_N and same>=MIN_SAME_BLOCK and len(blockset)>=MIN_DISTINCT_BLOCKS and len(used)==MATCH_N*4
    report={
       "schema":"gesui_Germany2022_EEA_original_WWTP_strict_between_WISE6_same_river_TP_SOURCEONLY_matched_quartets_v01",
       "UTC":datetime.now(timezone.utc).isoformat(),
       "original_NATIONAL_2022_EEA_WWTP_source_reach_label_QC_SHA256":sha(a.national_source_between),
       "source_2022_EEA_original_both_NP_one_FW_outlet_provenance_reproduced":True,
       "country_only":"DE",
       "original_WISE6_2022_total_TP_numeric_sample_values_OPENED":False,
       "original_all_DE_candidate_WWTP_intervening_monitor_paths":len(sources),
       "original_all_DE_candidate_clean_zero_registered_intervening_source_paths":len(controls),
       "all_four_station_disjoint_source_vs_control_quartet_OPTIONS_sourceonly_n":len(candidates),
       "matching_cost_definition":"abs(log1p(source_full_reach_path_km)-log1p(control_full_reach_path_km)) + 0.4 if different original 2degree upstream station geographical block",
       "matching_source_only_sort_order":"matching_cost, original four pre-hashed monitoring site IDs SHA256",
       "matching_greedy_strict_no_repeat_monitor_site":True,
       "country2degree_same_upstream_block_matched_quartets_n":same,
       "country2degree_unique_geographical_blocks_n":len(blockset),
       "predeclared_8_quartet_sourceonly_quality_GATE":gate,
       "selected_original_8_source_vs_clean_zero_monitor_quartets":selected,
       "selected_original_2022_WISE6_distinct_hashed_monitor_sites_n":len(used),
       "original_selected_quartet_matching_cost_median":sorted([x["value"] for x in selected])[len(selected)//2] if selected else None,
       "original_source_pairs_with_confirmed_intervening_outfall_are_NOT_causal_experiments":True,
       "no_numeric_river_TP_or_Lycaenidae_butterfly_response_accessed":True,
       "next_if_gate_PASS":"Read only preselected 32 WISE6 2022 river TP monitoring sites; exact 2022 TP, LOQ/date sample quality, and same-country noncausal upstream-downstream delta conditional on 6 valid disjoint quartets and >=3 same block",
       "next_if_gate_FAIL":"STOP no numerical TP effect analysis",
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_DE_2022_SOURCEONLY_WISE6_TP_MATCHED_8_GERMAN_SOURCE_VS_CLEAN_CONTROL_QUARTETS",
       json.dumps({"candidate_pair_of_pairs":len(candidates),"matched_quartets":len(selected),
        "matched_same_2degree_block":same,"regions":len(blockset),
        "all_32_independent_sites":len(used),"sourceonly_gate":gate,
        "source_DE_pair_candidates":len(sources),"control_DE_pair_candidates":len(controls),
        "NO_TP_NUMERIC_OUTCOMES_READ":True},sort_keys=True),flush=True)
if __name__=="__main__":main()
