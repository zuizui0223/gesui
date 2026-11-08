#!/usr/bin/env python3
"""POSTHOC regional block bootstrap of already archived 2021/2022 Lycaenidae
photo reporting 2022 WWTP-P source model predictive errors, no new outcome.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

N=1999
SEED=20261108

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for v in iter(lambda:f.read(1048576),b""):h.update(v)
    return h.hexdigest()
def bootstrap(model,rng):
    preds=model["per_country2deg_heldout_cell_predictions"]
    grouped=defaultdict(list)
    for v in preds:grouped[v["group"]].append(v)
    groups=sorted(grouped)
    baseline=float(model["baseline_MSE"])
    source=float(model["source_P_added_MSE"])
    observed=100*(baseline-source)/baseline
    draws=[]
    for _ in range(N):
        sampled=rng.choice(groups,size=len(groups),replace=True)
        rows=[x for key in sampled for x in grouped[key]]
        b=float(np.mean([x["baseline_sqerr"] for x in rows]))
        s=float(np.mean([x["P_source_sqerr"] for x in rows]))
        if b>0:draws.append(100*(b-s)/b)
    if len(draws)!=N:raise ValueError("Unexpected zero MSE in source-heldout error bootstrap")
    return {
      "actual_MSE_improvement_percent_replayed":observed,
      "original_predicted_cells":len(preds),
      "unique_original_country2degree_heldout_groups":len(groups),
      "original_baseline_MSE":baseline,"original_P_source_MSE":source,
      "regional_block_bootstrap_N":len(draws),
      "seed":SEED,
      "median_MSE_improvement_percent":float(np.median(draws)),
      "descriptive_95pct_MSE_improvement_interval":[float(v) for v in np.quantile(draws,[.025,.975])],
      "fraction_regional_resamples_gain_le0":float(np.mean(np.asarray(draws)<=0)),
      "NOT_confidence_interval_for_wastewater_causality_or_independent_field_experiment":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--prior-result",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    j=json.loads(a.prior_result.read_text())
    if j["status"]!="POSTHOC_HISTORICAL_PHOTO_COMPOSITION_AND_SOURCE_PREDICTIVE_INCREMENT_COMPLETED":
        raise ValueError("Required prior historical photo + source CV not completed")
    if j["same_original_2021_and_2022_complete_family_photo_supported_cells"]!=42:
        raise ValueError("No new source photo sample selection allowed")
    rng=np.random.default_rng(SEED)
    reports={}
    for name,m in j["fixed_2degree_block_heldout_source_P_models"].items():
        reports[name]=bootstrap(m,rng)
        old=m["EEA2022_reported_P_source_incremental_MSE_improvement_percent"]
        if abs(old-reports[name]["actual_MSE_improvement_percent_replayed"])>1e-9:
            raise ValueError("Original historical source photo MSE does not reproduce")
        print("GESUI_2021_PRIORYEAR_BUTTERFLY_2022_P_SOURCE_POSTHOC_REGION_BOOTSTRAP",
          json.dumps({"model":name,
           "source_P_gain":old,
           "bootstrap_2deg95_range":reports[name]["descriptive_95pct_MSE_improvement_interval"],
           "fraction_gain_le_zero":reports[name]["fraction_regional_resamples_gain_le0"],
           "original_model_n":reports[name]["original_predicted_cells"]},sort_keys=True),flush=True)
    out={
      "schema":"gesui_2022_2021_butterfly_historical_source_P_posthoc_2deg_block_prediction_error_resample_v01",
      "utc":datetime.now(timezone.utc).isoformat(),
      "POSTHOC_after_2021_history_P_source_3_70pct_gain_seen":True,
      "original_2021_history_2022_outcome_source_SHA256":sha(a.prior_result),
      "original_region_error_model_results":reports,
      "original_2022_source_P_historical_baseline_MSE_improvement_is_not_causal":True,
      "no_photo_or_2021_2022_API_data_requested":True,
      "no_retroactive_outcome_year_sample_selection":True,
      "important_statistical_limitations":[
         "2deg blocks are spatial group units but training fitted models overlap",
         "Bootstrap of heldout residuals does not capture full source data collection uncertainty",
         "Broad ecological land-use and observer effort confounds are not removed",
         "2021 WWTP source operation/history not known and 2022 annual reported mass is not environmental dose",
         "Photo submissions are not butterflies, behavior or larval performance"
      ]
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_POSTHOC_PRIORYEAR_HISTORY_REGIONAL_ERROR_UNCERTAINTY_COMPLETE",flush=True)
if __name__=="__main__":main()
