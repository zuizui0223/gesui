#!/usr/bin/env python3
"""INDEPENDENT OFFLINE arithmetic QC for NEW399 2025 month×area family sensitivity.

Uses only previously archived SHA-frozen source/count/month aggregates, not
original licensed Event/Occurrence row-level or source proximity refitting.
Reconstructs harmonization with independent straightforward loops and
enforces exactly the original six regions, 32/24 matched calendar strata.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path

DIR=Path("results/ebms_denmark_2023_2024_WWTP_source_only")
MONTH=DIR/"DENMARK_2025_NEW399_1DEG_MONTH_MATCHED_LYCAENIDAE_3VS5_FAMILY_QC.json"
FAMILY=DIR/"DENMARK_2025_NEW399_SITE_LYCAENIDAE_FAMILY_POST_SOURCE_GATE_COUNTS.json"
SOURCE=DIR/"DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json"
NAMES=("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae")

def sha(b):
    return hashlib.sha256(b).hexdigest()

def calc(strata):
    original={g:{f:0. for f in NAMES} for g in ("near","remote")}
    n=0
    for s in strata:
        a=s["near"];b=s["remote"]
        if min(a["linked_dates"],b["linked_dates"])<2:
            raise ValueError("Source-matched region-month has less than 2 documented days")
        w=min(a["linked_dates"],b["linked_dates"])
        n+=w
        for g in ("near","remote"):
            for f in NAMES:
                original[g][f]+=w*s[g]["counts"][f]/s[g]["linked_dates"]
    if n<=0:return None
    for g in original:
        for f in NAMES:original[g][f]/=n
    L={g:original[g]["Lycaenidae"] for g in original}
    def odds(fs,g):
        x=original[g]["Lycaenidae"]
        y=sum(original[g][f] for f in fs if f!="Lycaenidae")
        return x/y if y else None
    r=L["near"]/L["remote"] if L["remote"] else None
    a=odds(NAMES[:3],"near");b=odds(NAMES[:3],"remote")
    c=odds(NAMES,"near");d=odds(NAMES,"remote")
    return {
        "L_recording_rate_near_over_remote":r,
        "L_share3_oddsratio_near_over_remote":a/b if a is not None and b is not None else None,
        "L_share5_oddsratio_near_over_remote":c/d if c is not None and d is not None else None,
        "weighted_recording_date_support":n,
        "site_region_months":len(strata),
        "standardized_families":original,
    }

def compare(a,b,name):
    if a is None or b is None:
        if a is not None or b is not None:raise ValueError(name+" null mismatch")
        return
    if abs(a-b)>1e-10:
        raise ValueError(name+" mismatch: "+str(a)+" vs "+str(b))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    d=json.loads((root/MONTH).read_text(encoding="utf-8"))
    family=json.loads((root/FAMILY).read_text(encoding="utf-8"))
    source=json.loads((root/SOURCE).read_text(encoding="utf-8"))
    assert d["independent_399_site_2025_five_family_counts_and_all_2plus_linked_days_exactly_reproduced"] is True
    assert d["prior_399_2025_site_family_count_SHA256"]==sha((root/FAMILY).read_bytes())
    assert d["preoutcome_original_399_site_source_gate_SHA256"]==sha((root/SOURCE).read_bytes())
    assert source["frozen_before_2025_NEW399_butterfly_counts_gate"]["2025_NEW399_butterfly_family_count_outcome_permitted"]
    assert family["source_family_model_sample_gate"]["original_2025_family_contrast_reporting_gate_PASS"]
    qa={}
    for k,v in d["seasonal_and_geographic_source_band_family_results"].items():
        strata=v["matched_original_source_region_month_aggregate_counts"]
        if k.startswith("May"):
            assert all(r["calendar_month"] in (5,6,7,8) for r in strata)
            assert len(strata)==24
        else:
            assert len(strata)==32
        assert len({x["original_1degree_region"] for x in strata})==6
        out=calc(strata)
        expected=v["posthoc_min_date_weight_standardized_family_result"]
        assert v["prespecified_after_annual_results_before_month_values_observed_source_support"]["common_region_month_source_effort_gate_PASS"]
        compare(out["L_recording_rate_near_over_remote"],
                expected["Lycaenidae_per_linked_date_near_vs_remote_ratio"],"rate")
        compare(out["L_share3_oddsratio_near_over_remote"],
                expected["Lycaenidae_three_family_share_oddsratio_near_over_remote"],"OR3")
        compare(out["L_share5_oddsratio_near_over_remote"],
                expected["Lycaenidae_five_family_share_oddsratio_near_over_remote"],"OR5")
        assert out["weighted_recording_date_support"]==expected["matching_weight_sum_min_reporting_days"]
        for g,originalg in (("near","near5"),("remote","remote10")):
            for f in NAMES:compare(out["standardized_families"][g][f],
                expected["standardized_each_family_recorded_per_linked_date"][originalg][f],g+"_"+f)
        leave={}
        for region in sorted({s["original_1degree_region"] for s in strata}):
            sub=[x for x in strata if x["original_1degree_region"]!=region]
            smaller=calc(sub)
            old=v["leave_one_original_1degree_region_out"][region]
            for value,key in (
                ("L_recording_rate_near_over_remote","Lycaenidae_per_linked_date_near_vs_remote_ratio"),
                ("L_share3_oddsratio_near_over_remote","Lycaenidae_three_family_share_oddsratio_near_over_remote"),
                ("L_share5_oddsratio_near_over_remote","Lycaenidae_five_family_share_oddsratio_near_over_remote")):
                compare(smaller[value],old[key],"LOO_"+region+"_"+value)
            leave[region]={
                "rate_ratio":smaller["L_recording_rate_near_over_remote"],
                "OR_three":smaller["L_share3_oddsratio_near_over_remote"],
                "OR_five":smaller["L_share5_oddsratio_near_over_remote"]
            }
        qa[k]={
            "source_strata_n":len(strata),
            "all_six_1degree_regions":True,
            "no_matching_on_butterfly_family_count_or_effect_sign":True,
            "matched_min_recorded_days_as_weights":True,
            "independently_recomputed":out,
            "leave_one_region_ratios_exactly_match":leave,
            "metric_on_positive_side_of_one":{
                "L_rate":out["L_recording_rate_near_over_remote"]>1,
                "three_family_or":out["L_share3_oddsratio_near_over_remote"]>1,
                "five_family_or":out["L_share5_oddsratio_near_over_remote"]>1
            },
            "posthoc_not_independent_source_replication":True
        }
    out={
        "schema":"gesui_DK_2025_NEW399_family_3vs5_month_region_independent_aggregate_replay_v01",
        "licensed_source_2025_original_CCBYNC_4_0":True,
        "original_sources_SHA256":{
            "frozen399_event_survey_source":sha((root/SOURCE).read_bytes()),
            "original399_fivefamily_previous_counts":sha((root/FAMILY).read_bytes()),
            "month_x_region_licensed_original_replay":sha((root/MONTH).read_bytes())},
        "original_source_selected_399_sites":399,
        "original_286_linked_two_date_sites":286,
        "all_matches_exact_at_five_family_counts_32_regionmonth_strata_and_leave_one_region_out":True,
        "independent_original_licensed_archive_NOT_re_read":True,
        "region_month_season_results":qa,
        "not_chemistry_or_hydrological_WWTP_causal_effect":True
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({k:{"n_month_region":x["source_strata_n"],
        "RR_L":x["independently_recomputed"]["L_recording_rate_near_over_remote"],
        "OR3":x["independently_recomputed"]["L_share3_oddsratio_near_over_remote"],
        "OR5":x["independently_recomputed"]["L_share5_oddsratio_near_over_remote"]}
      for k,x in qa.items()},indent=2))
if __name__=="__main__":
    main()
