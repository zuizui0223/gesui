#!/usr/bin/env python3
"""Offline independent arithmetic replay of original 2022 same-month WWTP TP quartets.

Reads only existing source-locked JSON: original source-only sampling dates
and derived site/month median TP results. NO source API queries, person data,
butterfly records, or additional sampling. Verifies all four arm differences,
time window, equal-quartet mean, 9999 descriptive bootstrap and sign-flips.
"""
from __future__ import annotations
import argparse,hashlib,json,math,itertools
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

OFFSET=0.001
SEED=20261109
NBOOT=9999
TOL=1e-10

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def close(a,b):
    if abs(a-b)>TOL:raise ValueError(f"Original monthly contrast mismatch {a} != {b}")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--date-source",type=Path,required=True)
    p.add_argument("--monthly-derived",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    d=json.loads(a.date_source.read_text())
    r=json.loads(a.monthly_derived.read_text())
    if d["source_exact_fourstation_calendar_date_gate_PASS"] or not d["source_same_calendar_month_fourstation_gate_PASS"]:
        raise ValueError("Original pre-numeric same-day/month source gate changed")
    if r["original_source_first_sample_date_gate_SHA256"]!=sha(a.date_source):
        raise ValueError("Original source-only monthly dates don't match numeric results")
    if r["source_4site_numeric_same_day_geographic_gate_pass"] or not r["source_4site_numeric_same_month_geographic_gate_pass"]:
        raise ValueError("Original 2022 monthly numeric source-support gate invalid")
    if r["source_same_day_valid_match_results"] is not None or r["source_same_month_valid_match_results"] is None:
        raise ValueError("Same-day effect inferred despite failed original source gate")
    a_quartets={int(x["quartet_index"]):x for x in d["all_original_8_no_replacement_quartets_and_date_coverage"]}
    b_quartets=r["all_original_8_quartet_source_month_numeric_QC"]
    if len(a_quartets)!=8 or len(b_quartets)!=8:raise ValueError("Original 8 quartet source population changed")
    vals=[];groups=set();meta=[]
    for q in b_quartets:
        idx=int(q["original_quartet_index"])
        if idx not in a_quartets:raise ValueError("Unknown original source frozen quartet")
        original=a_quartets[idx]
        if original["source_country2degree_upstream_block"]!=q["source_country2degree_upstream_block"]:
            raise ValueError("Original quartet country 2degree group changed")
        pvals=q["month_specific_all4_TP_measurements_and_source_control_deltas"]
        if q["numeric_month_gate_this_quartet"]:
            if len(pvals)<3:raise ValueError("Accepted quartet not >=3 common four-arm monthly observations")
            originalmonths=set(original["four_arm_same_YYYY_MM_months"])
            deltas=[]
            for x in pvals:
                month=x["YYYY_MM"]
                if month not in originalmonths:raise ValueError("Numeric monthly TP from month not prequalified prior to TP outcome")
                v=[float(x[k]) for k in (
                    "source_upstream_TP_mgP_L","source_downstream_TP_mgP_L",
                    "control_upstream_TP_mgP_L","control_downstream_TP_mgP_L")]
                if any(t<0 or not math.isfinite(t) for t in v):raise ValueError("Invalid original real river TP from QA")
                s=math.log(v[1]+OFFSET)-math.log(v[0]+OFFSET)
                z=math.log(v[3]+OFFSET)-math.log(v[2]+OFFSET)
                contrast=s-z
                close(float(x["source_gradient_logTP"]),s)
                close(float(x["control_gradient_logTP"]),z)
                close(float(x["matched_gradient_source_minus_control"]),contrast)
                deltas.append(contrast)
            qmean=float(np.mean(deltas))
            close(float(q["mean_matched_source_minus_control_monthly_logTP_gradient"]),qmean)
            vals.append(qmean)
            groups.add(q["source_country2degree_upstream_block"])
            meta.append({"original_quartet_index":idx,"n_4station_sampled_months":len(pvals),
                         "source_minus_control_monthly_logTP_mean_reproduced":qmean})
        else:
            if q["mean_matched_source_minus_control_monthly_logTP_gradient"] is not None:
                raise ValueError("Unqualified original quartet retained TP comparison")
    if len(vals)!=r["source_same_month_valid_match_results"]["original_independent_matched_quartets_n"]:
        raise ValueError("Original number of monthly TP quartets differs")
    if len(vals)<5 or len(groups)<4:raise ValueError("Original month support no longer passes gate")
    arr=np.asarray(vals,float)
    mean=float(np.mean(arr))
    rr=r["source_same_month_valid_match_results"]
    close(mean,rr["quartet_equal_weight_mean_matched_source_minus_sourcezero_logTP_gradient"])
    close(float(np.median(arr)),rr["quartet_equal_weight_median_matched_gradient"])
    if (int(np.sum(arr>0)),int(np.sum(arr<0)))!=(rr["positive_contrast_quartets"],rr["negative_contrast_quartets"]):
        raise ValueError("Source quartet direction counts don't reproduce")
    rng=np.random.default_rng(SEED)
    sim=np.mean(rng.choice(arr,size=(NBOOT,len(arr)),replace=True),axis=1)
    ci=[float(x) for x in np.quantile(sim,[0.025,0.975])]
    for x,y in zip(ci,rr["descriptive_9999_quartet_resample95pct_interval"]):close(x,y)
    flips=[float(np.mean(arr*np.array(s,float))) for s in itertools.product((-1,1),repeat=len(arr))]
    tail=sum(abs(x)>=abs(mean)-1e-12 for x in flips)/len(flips)
    close(tail,rr["source_only_equal_quartet_signflip_2sided_reference_fraction"])
    loo=[float(np.mean(np.delete(arr,i))) for i in range(len(arr))]
    for x,y in zip((min(loo),max(loo)),rr["original_quartet_leave_one_out_mean_range"]):close(x,y)
    status_codes={}
    for v in r["all_original_source_32monitoring_site_raw_LOQ_and_status_qc"].values():
        for k,n in v["original_resultObservationStatus_source_codes"].items():
            status_codes[k]=status_codes.get(k,0)+n
    result={
      "schema":"gesui_2022_DE_8source_vs_registeredzero_2022_WISE6_monthlyTP_offline_reproduction_v01",
      "utc":datetime.now(timezone.utc).isoformat(),
      "exact_source_only_dates_sha256":sha(a.date_source),
      "official_measured_TP_derived_original_sha256":sha(a.monthly_derived),
      "original_2022_all8_quartets_untouched":True,
      "original_source_quartets_permitted_monthly":len(vals),
      "original_country2degree_blocks_permitted":len(groups),
      "original_2022_monthly_summary_all_recomputed_from_four_station_TP":meta,
      "original_equal_quartet_mean":mean,
      "original_monthly_95pct_descriptive_bootstrap":ci,
      "original_2sided_signflip_reference_fraction":tail,
      "leave_one_out_mean_range":[min(loo),max(loo)],
      "all_original_32_public_water_data_resultObservationStatus_source_counts":status_codes,
      "no_same_day_numeric_inference":True,
      "original_WISE6_numeric_TP_is_not_butterfly_exposure_or_causal_WWTP_effect":True,
      "no_new_observations_or_2022_TP_measurements_downloaded":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_2022_WISE6_4ARM_MONTHLY_TP_OFFLINE_INDEPENDENT_REPLAY_SUCCESS",
      json.dumps({"quartets":len(vals),"blocks":len(groups),"mean":mean,
         "95pct":ci,"signflip":tail,"leaveoneout":[min(loo),max(loo)],
         "source_original_resultObservationStatus_codes":status_codes},sort_keys=True),flush=True)

if __name__=="__main__":main()
