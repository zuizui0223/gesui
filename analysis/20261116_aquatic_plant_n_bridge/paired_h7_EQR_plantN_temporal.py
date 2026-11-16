#!/usr/bin/env python3
"""Paired true-H7 temporal EQR-PlantN cochange; observed, NOT causal.

The same native watershed is compared across source windows 2017-19 and
2020-22; original EEA EQR evaluations represent calendar 2018/2022.
Primary cohort requires unchanged majority EQR assessment method.
Country and genuine H6 remain fixed, original W is static and cannot be
used to claim a temporal wastewater cause. Source opportunity counts
and eligible cohorts are frozen before opening any fitted coefficients.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT / "analysis/20261116_aquatic_plant_n_bridge"))
from same_H7_EQR_plantN import source_cohort, EQR, P_N, HNS
sys.path.insert(0,str(ROOT / "analysis/20261115_eu454_plant_nutri"))
from h6_model_if_gate import source_cluster_regression,cv_h6

EARLY="2017-2019"
LATE="2020-2022"
SEED=20261116
MIN_PAIRS=100
MIN_PARENTS=40
CORE=["d_nEEA","d_nPlant","d_nButterfly","d_nAnimal",
      "d_plant_coverage","d_scored_plant_species"]

def pair_source(source):
    df=source_cohort(source)
    rows=[]
    for h,g in df.groupby("h",sort=True):
        if len(g)==1:
            continue
        if len(g)!=2 or set(g.period)!={EARLY,LATE}:
            raise RuntimeError("True original site pair nonunique/unexpected periods "+str(h))
        x=g.loc[g.period==EARLY].iloc[0]
        y=g.loc[g.period==LATE].iloc[0]
        if x.h6!=y.h6 or x.country_holdout!=y.country_holdout:
            raise RuntimeError("Same H7 changed true parent/country over two source periods")
        if not np.isclose(float(x.W),float(y.W),rtol=0,atol=1e-12):
            raise RuntimeError("Static source W unexpectedly changed")
        def delta(col):
            return float(y[col])-float(x[col])
        a={"h":str(h),"h6":str(x.h6),"country_method":str(x.country_holdout),
           "period":"delta_late_minus_early",
           "same_majority_EQR_method":bool(str(x.method)==str(y.method)),
           "d_EQR":delta(EQR),"d_plant_N":delta(P_N),
           "d_plant_N_record":delta("Plant_N_record"),
           "d_HNS_incidence":delta(HNS),
           "d_HNS_record":delta("primary_record_mean"),
           "d_nEEA":delta("nEEA"),"d_nPlant":delta("nPlant"),
           "d_nButterfly":delta("nButterfly"),"d_nAnimal":delta("nAnimal"),
           "d_plant_coverage":delta("NUTRI_coverage_all_Plantae_speciesKeys"),
           "d_scored_plant_species":delta("log1p_eiv_scored_taxa"),
           "d_hns_coverage":delta("primary_record_coverage_raw"),
           "d_hns_completeness":delta("primary_host_completeness")}
        rows.append(a)
    d=pd.DataFrame(rows)
    if len(d)!=219 or d.h6.nunique()!=86:
        raise RuntimeError("Previously source-frozen 219 paired H7/86 H6 drift: "+str((len(d),d.h6.nunique())))
    if d.h.duplicated().any() or d[["d_EQR","d_plant_N","d_HNS_incidence"]].isna().any().any():
        raise RuntimeError("Original actual EQR NUTRI/HNS missing for paired H7")
    # Fixed pre-outcome main: same original majority EQR assessment method.
    same=d.loc[d.same_majority_EQR_method].reset_index(drop=True)
    if len(same)<MIN_PAIRS or same.h6.nunique()<MIN_PARENTS:
        return d,same,False
    return d,same,True


def eval_model(df,y,x,controls):
    sub=df.copy().reset_index(drop=True)
    # The generic fitted routine only reads these two fields to describe
    # underlying design; neither is added as the temporal covariate.
    sub["W"]=0.
    sub["Plant_N_incidence"]=sub["d_plant_N"]
    for k in [y,x]+controls:
        if not np.isfinite(sub[k]).all() or sub[k].std(ddof=0)<1e-12:
            raise RuntimeError("Fixed original time-difference source has invalid/constant model column: "+k)
    regression=source_cluster_regression(sub,y,controls+[x],x)
    prediction=cv_h6(sub,y,controls,x)
    return {"adjusted_beta_per_1SD_delta":regression["adjusted_beta_per_1SD_focal"],
            "H6_descriptive95":regression["H6_cluster_descriptive95"],
            "H6_out_gain_percent":prediction["H6_equal_added_MSE_gain_pct"],
            "H6_out_gain95_fixed_error":prediction["fixed_out_of_fold_error_H6_resampling95_pct"],
            "n_true_H6":int(df.h6.nunique()),"n_true_H7_pairs":len(df),
            "regression":regression,"holdout":prediction,
            "not_intervention_or_causal_change":True}


def h6_raw_corr(df,a,b):
    # Geographic units are H6 means; fixed-value group bootstrap.
    grouped=df.groupby("h6",sort=True)[[a,b]].mean()
    x=grouped[a].to_numpy(float);y=grouped[b].to_numpy(float)
    r=float(np.corrcoef(x,y)[0,1])
    rng=np.random.default_rng(SEED)
    sample=rng.integers(0,len(x),size=(2000,len(x)))
    arr=np.array([np.corrcoef(x[z],y[z])[0,1] for z in sample])
    arr=arr[np.isfinite(arr)]
    if len(arr)<1000:raise RuntimeError("H6 group uncertainty simulation unstable")
    return {"H6_groupmean_pearson_r":r,
            "H6_fixed_group_values_bootstrap95":np.quantile(arr,[.025,.975]).tolist(),
            "not_ecological_mediation":True}


def main(a):
    df,same,supported=pair_source(a.source)
    out={"schema":"GESUI_20261116_PAIRED_TRUE_H7_EQR_PLANT_NUTRI_TEMPORAL_BRIDGE_V1",
         "status":"PREOUTCOME_GATE_PASS" if supported else "STOP_UNCHANGED_EQR_METHOD_SOURCE_SUPPORT",
         "source_unchanged_451_232_90":True,
         "candidate_two_period_pairs":len(df),
         "candidate_true_H6":int(df.h6.nunique()),
         "primary_exact_majority_EQR_method_unchanged_pairs":len(same),
         "primary_exact_majority_EQR_method_unchanged_H6":int(same.h6.nunique()),
         "original_periods":[EARLY,LATE],
         "frozen_min_pairs":MIN_PAIRS,"frozen_min_H6":MIN_PARENTS,
         "reporting_source_EEA_2018_2022_not_same_physical_butterfly_sites":True,
         "original_static_W_not_time_varying":True,
         "no_temporal_discharge_assignment":True,
         "fixed_controls_CORE":CORE,
         "source_only_before_fits":True,
         "posthoc_exploratory":True}
    if supported:
        # Primary original same-method stations geographically aggregated to H7.
        prim=eval_model(same,"d_EQR","d_plant_N",CORE)
        # Secondary: butterfly resource-habitat affinity difference, beyond
        # plant N difference and original recording-quality shifts.
        hns=eval_model(same,"d_HNS_incidence","d_EQR",
                       CORE+["d_plant_N","d_hns_coverage","d_hns_completeness"])
        # Reversal must not be interpreted as causal direction.
        rev=eval_model(same,"d_plant_N","d_EQR",CORE)
        out["primary_same_EQR_method_EQR_on_plant_N_delta"]=prim
        out["secondary_same_method_HNS_on_EQR_after_plantN_delta"]=hns
        out["secondary_same_method_plantN_on_EQR_delta"]=rev
        out["descriptive_H6_raw_mean_delta_EQR_plantN"]=h6_raw_corr(same,"d_EQR","d_plant_N")
        brief={"status":out["status"],"n_H7":len(same),"n_H6":int(same.h6.nunique()),
               "EQR_dPlantN":{k:prim[k] for k in ("adjusted_beta_per_1SD_delta","H6_descriptive95","H6_out_gain_percent","H6_out_gain95_fixed_error")},
               "HNS_dEQR_after_dPlantN":{k:hns[k] for k in ("adjusted_beta_per_1SD_delta","H6_descriptive95","H6_out_gain_percent","H6_out_gain95_fixed_error")},
               "H6_raw_delta_corr":out["descriptive_H6_raw_mean_delta_EQR_plantN"]}
    else:
        brief={"status":out["status"],"n_H7":len(same),"n_H6":int(same.h6.nunique()),"outcome_not_fit":True}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"GESUI_TRUE_H7_PAIRED_EQR_PLANTN_TEMPORAL_BRIDGE.json").write_text(
        json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_TRUE_H7_EQR_PLANTN_TEMPORAL_BRIDGE="+json.dumps(brief,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
