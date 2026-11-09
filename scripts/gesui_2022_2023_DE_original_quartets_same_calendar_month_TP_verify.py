#!/usr/bin/env python3
"""Independent offline calendar-month-of-year 2022 vs 2023 river TP
four-arm WWTP/control gradient persistence check, same *original* sites.

Source cohort and four quartet/year overlap months were fixed on source-
only monitoring dates before 2023 numeric values were examined.
Pure JSON input; no network, no new chemistry, no butterfly data.
Historical 2022 WWTP P source reporter is not year-2023 treatment data.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json,math,statistics
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

B=9999
SEED=20261109
MIN_OVERLAP_MONTHS=3
MIN_QUARTETS=4
MIN_BLOCKS=3
EXPECTED_SELECTED_SOURCE_ONLY_QUARTETS={2,5,6,7}
OFFSET=.001

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for q in iter(lambda:f.read(1048576),b""):h.update(q)
    return h.hexdigest()

def source_only_preselection(meta22,meta23):
    date22={x["quartet_index"]:x for x in meta22["all_original_8_no_replacement_quartets_and_date_coverage"]}
    if len(date22)!=8:raise ValueError("Original 2022 source quartets altered")
    candidates=[]
    for p in meta23["original8_unreplaced_quartets_2023_fourarm_sample_time_coverage"]:
        idx=p["original_2022_quartet_index"];a=date22[idx]
        m22={m[5:] for m in a["four_arm_same_YYYY_MM_months"]}
        m23={m[5:] for m in p["original_four_arm_shared_2023_calendar_months"]}
        shared=sorted(m22&m23)
        if len(shared)>=MIN_OVERLAP_MONTHS:
            candidates.append({"quartet":idx,"block":p["original_2022_same_river_country2degree_source_block"],"month_of_year":shared,
                "exact_2022_source_frozen_arms":a["four_arm_original_WISE_country_scheme_site_hashes"],
                "exact_2023_source_frozen_arms":p["original_WISE_source_only_fixed_4site_hashes"]})
    if {q["quartet"] for q in candidates}!=EXPECTED_SELECTED_SOURCE_ONLY_QUARTETS:
        raise ValueError("Original 2023 source-only sample changed after source preselection")
    if len(candidates)!=4 or len(set(x["block"] for x in candidates))!=3:
        raise ValueError("Original 4quartet/3block prior-year source comparability gate changed")
    for p in candidates:
        if p["exact_2022_source_frozen_arms"]!=p["exact_2023_source_frozen_arms"]:
            raise ValueError("Original 2022/23 same four WISE site IDs not source-matched")
    return candidates

def measured_delta(period):
    if "matched_gradient_source_minus_control" in period:
        old=float(period["matched_gradient_source_minus_control"])
        # Independent algebra checks from four original source-up/down, ctrl-up/down medians.
        def v(k):return float(period[k])
        computed=math.log(v("source_downstream_TP_mgP_L")+OFFSET)-math.log(v("source_upstream_TP_mgP_L")+OFFSET)
        computed-=math.log(v("control_downstream_TP_mgP_L")+OFFSET)-math.log(v("control_upstream_TP_mgP_L")+OFFSET)
    else:
        old=float(period["source_minus_reference_logTP_gradient"])
        def v(k):return float(period[k])
        computed=math.log(v("source_downstream_mgP_L")+OFFSET)-math.log(v("source_upstream_mgP_L")+OFFSET)
        computed-=math.log(v("registered_zero_downstream_mgP_L")+OFFSET)-math.log(v("registered_zero_upstream_mgP_L")+OFFSET)
    if abs(computed-old)>1e-10:
        raise ValueError("Independent 4-arm TP gradient differs from stored source output")
    return old

def exact_signflip(vals):
    q=np.asarray(vals,float);a=float(np.mean(q))
    perms=np.asarray([float(np.mean(q*np.asarray(signs,float)))
          for signs in itertools.product((-1.,1.),repeat=len(q))])
    return float(np.mean(np.abs(perms)>=abs(a)-1e-12))

def boot(vals):
    v=np.asarray(vals,float)
    rng=np.random.default_rng(SEED)
    draws=np.mean(rng.choice(v,size=(B,len(v)),replace=True),axis=1)
    return [float(x) for x in np.quantile(draws,[.025,.975])]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-2022-source-dates",type=Path,required=True)
    ap.add_argument("--original-2023-source-dates",type=Path,required=True)
    ap.add_argument("--original-2022-numeric-months",type=Path,required=True)
    ap.add_argument("--original-2023-numeric-months",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    d22=json.loads(a.original_2022_source_dates.read_text())
    d23=json.loads(a.original_2023_source_dates.read_text())
    x22=json.loads(a.original_2022_numeric_months.read_text())
    x23=json.loads(a.original_2023_numeric_months.read_text())
    if not(d23["strict_2023_frozen_4arm_same_MONTH_numeric_gate_PASS"] and
           x23["2023_numeric_same_month_min5_quartets_min4_blocks_GATE"] and
           x22["source_4site_numeric_same_month_geographic_gate_pass"]):
        raise ValueError("Original source/date/numeric quality gates must pass before 2022–2023 matched-month study")
    if x23["original_2023_sampling_date_source_QC_SHA256"]!=sha(a.original_2023_source_dates):
        raise ValueError("2023 numeric source not correctly tied to original date-only SHA")
    if x23["original_2022_source_year_matched_TP_file_SHA256"]!=sha(a.original_2022_numeric_months):
        raise ValueError("2023 numeric check not matched to original 2022 source year")
    qs=source_only_preselection(d22,d23)
    prev={x["original_quartet_index"]:x for x in x22["all_original_8_quartet_source_month_numeric_QC"]}
    future={x["source_original_2022_quartet_index"]:x for x in x23["all_original8_2023_numeric_fourarm_quartet_source_comparisons"]}
    records=[];kept=[]
    for item in qs:
        idx=item["quartet"];old=prev[idx];new=future[idx]
        if old["original_source_only_four_site_hashes"]!=item["exact_2022_source_frozen_arms"] or new["original_station_hashes"]!=item["exact_2023_source_frozen_arms"]:
            raise ValueError("Numerical same-year source arms not original frozen source-only 8 quartet stations")
        if not(old["numeric_month_gate_this_quartet"] and new["2023_numeric_fourarm_month_gate"]):
            records.append({"original_2022_quartet_index":idx,"source_2022_2023_numeric_LOQ_gate":False,
                            "source_2022_2023_3same_month_overlap":False})
            continue
        mo22={p["YYYY_MM"][5:]:measured_delta(p) for p in old["month_specific_all4_TP_measurements_and_source_control_deltas"]}
        mo23={p["original_2023_sampling_period"][5:]:measured_delta(p) for p in new["2023_samemonth_numeric_source_vs_zero_logTP_gradients"]}
        common=sorted(set(mo22)&set(mo23)&set(item["month_of_year"]))
        if len(common)<MIN_OVERLAP_MONTHS:
            records.append({"original_2022_quartet_index":idx,"source_2022_2023_numeric_LOQ_gate":True,
                            "source_2022_2023_3same_month_overlap":False,"months_n":len(common)})
            continue
        diffs=[mo23[m]-mo22[m] for m in common]
        row={
          "original_2022_quartet_index":idx,"original_country_2deg_source_block":item["block"],
          "strict_original_same_calendar_months_in_both_years":common,
          "original_2022_mean_gradient_over_identical_month_of_year":float(np.mean([mo22[m] for m in common])),
          "original_2023_mean_gradient_over_identical_month_of_year":float(np.mean([mo23[m] for m in common])),
          "original_2023_minus_2022_mean_gradient_on_same_months":float(np.mean(diffs)),
          "independently_recomputed_monthwise_difference":dict(zip(common,diffs)),
          "no_year2023_WWTP_mass_or_true_dose_measured":True,
        }
        kept.append(row);records.append(row)
    blocks={r["original_country_2deg_source_block"] for r in kept}
    passgate=len(kept)>=MIN_QUARTETS and len(blocks)>=MIN_BLOCKS
    out={
      "schema":"gesui_2022_2023_DE_WISE6_realTP_original_same4arms_same_calendar_month_of_year_source_geoQC_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "original_2022_sample_calendar_file_SHA256":sha(a.original_2022_source_dates),
      "original_2023_sample_calendar_file_SHA256":sha(a.original_2023_source_dates),
      "original_2022_real_TP_numeric_months_file_SHA256":sha(a.original_2022_numeric_months),
      "original_2023_real_TP_numeric_months_file_SHA256":sha(a.original_2023_numeric_months),
      "posthoc_to_known_2022_positive_and_2023_year_temporal_test":True,
      "source_only_four_original_quartets_three_geoblocks_fixed_prior2023_numeric":True,
      "strict_numeric_same_calendar_month_quartet_gate_pass":passgate,
      "original_2022_2023_same_exact_4site_and_calendar_month_quartet_records":records,
      "n_eligible_fourarm_quartets_after_numeric_QC":len(kept),
      "n_original_country2deg_geographic_blocks":len(blocks),
      "no_true_wastewater_source_treatment_intervention_or_2023_reported_P_mass":True,
    }
    if passgate:
        q22=np.asarray([r["original_2022_mean_gradient_over_identical_month_of_year"] for r in kept])
        q23=np.asarray([r["original_2023_mean_gradient_over_identical_month_of_year"] for r in kept])
        diff=q23-q22
        out["strict_original_same_month_year_to_year_comparison"]={
          "original_same_quartets_n":len(kept),
          "original_2022_same_month_source_minus_control_mean":float(np.mean(q22)),
          "original_2023_same_month_source_minus_control_mean":float(np.mean(q23)),
          "original_2023_minus2022_same_month_paired_quartet_mean_change":float(np.mean(diff)),
          "n_year_to_year_gradient_increased":int(np.sum(diff>0)),
          "n_year_to_year_gradient_decreased":int(np.sum(diff<0)),
          "year_change_9999_quartet_resampled_95pct_descriptive_interval":boot(diff),
          "original_exact_quartet_year_change_signflip_two_sided_fraction":exact_signflip(diff),
          "no_causal_effluent_change_inference":True,
        }
    else:
        out["stage"]="STOP_BEFORE_AGGREGATE_2022_2023_COMPARISON_AFTER_SOURCE_GATED_MONTH_DATA"
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_GERMANY_2022_2023_SAME_CALENDAR_MONTH_REAL_TP_QUARTET_YEARCHECK",
      json.dumps({"source_frozen_QC_4quartets_3blocks":passgate,
         "eligible_numeric_quartets":len(kept),"blocks":len(blocks),
         "year_change":out.get("strict_original_same_month_year_to_year_comparison"),"ecological_causality":False},sort_keys=True),flush=True)
if __name__=="__main__":main()
