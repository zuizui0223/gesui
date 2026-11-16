#!/usr/bin/env python3
"""Fixed posthoc five pair-specific EIV specificity contrasts for butterfly HNS.

Earlier stricter six-indicator same-species gate STOP is retained. Each of
five source-approved (NUTRI, rival) pairs uses exactly the same source
plant species and same H7×period observations for the two trait scores.
No multi-control single-panel claim, wastewater dose, mediation or
species abundance interpretation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261115_eu454_plant_nutri"))
from h6_model_if_gate import (ORIGINAL_PANEL_SHA,EIV_SHA,FACETS_SHA,
    H6_weights, design_full, source_cluster_regression, cv_h6)
from adversarial_host_exclusion_country import country_predictive

RIVALS=["LIGHT","TEMPE","MOIST","REACT","SALIN"]
BASE=["crop","urban","riverflow","naturalout","area",
      "nEEA","nAnimal","nPlant","nButterfly",
      "primary_record_coverage_raw","primary_host_completeness",
      "W","pair_coverage","log1p_pair_scored"]
PRIMARY="primary_incidence_mean"
ALPHA=10.0

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def vif_N(d):
    X,names=design_full(d,BASE+["pair_other","pair_N"])
    j=names.index("pair_N")
    z=X[:,j]
    X0=np.delete(X,j,axis=1)
    w=H6_weights(d)
    root=np.sqrt(w)
    coeff=np.linalg.lstsq(root[:,None]*X0,root*z,rcond=None)[0]
    res=z-X0@coeff
    zm=np.average(z,weights=w)
    tss=float(np.sum(w*(z-zm)**2))
    rss=float(np.sum(w*res**2))
    if rss<1e-10:return float("inf")
    return tss/rss

def bh(pvalues):
    sorted_keys=sorted(pvalues,key=lambda k:pvalues[k])
    m=len(sorted_keys)
    output={}
    previous=1.
    for j in reversed(range(m)):
        k=sorted_keys[j]
        previous=min(previous,float(pvalues[k])*m/(j+1))
        output[k]=previous
    return output

def main(args):
    gate=json.loads((args.source/"EU454_EIV_NUTRI_FIVE_PAIR_SOURCE_GATE.json").read_text())
    if (gate.get("schema")!="GESUI_EIV_5_NUTRI_PAIRWISE_SAME_SCORED_SPECIES_SUPPORT_V1"
        or not gate.get("original_six_alltrait_source_gate_STOP")
        or gate.get("original_EIV_source_SHA")!=EIV_SHA
        or gate.get("original_2808_facets_SHA")!=FACETS_SHA
        or gate.get("original_454_panel_SHA")!=ORIGINAL_PANEL_SHA
        or gate.get("all_five_controls_prespecified")!=RIVALS
        or not gate.get("each_pair_on_exact_same_plants_for_both_N_and_other")):
        raise RuntimeError("Source-first genuine five original EIV pair opportunity gate missing")
    original=ROOT/"results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"
    if sha(original)!=ORIGINAL_PANEL_SHA:raise RuntimeError("Original 454 H7 science panel SHA drift")
    df=pd.read_csv(original,dtype={"h":str,"h6":str})
    if (len(df),df.h.nunique(),df.h6.nunique())!=(454,232,90):
        raise RuntimeError("True geography identity altered")
    source=pd.read_csv(args.source/"EU454_EIV_NUTRI_FIVE_PAIRS_SOURCE.csv",dtype={"h":str,"h6":str})
    if len(source)!=454 or source.duplicated(["h","period"]).any():
        raise RuntimeError("Source fifth EIV gate row ambiguity")
    output={
      "schema":"GESUI_EU454_NUTRI_PAIRED_5_EIV_SPECIFICITY_POSTHOC_V1",
      "scope":"posthoc five separate pairwise same-scored-plant tests; six-trait joint-model STOP",
      "source_original_eiv_sha":EIV_SHA,
      "source_original_h7_sha":ORIGINAL_PANEL_SHA,
      "source_original_2808_facets_sha":FACETS_SHA,
      "original_4011_exact_NUTRI_species_no_imputation":True,
      "all_five_rivals_prespecified":RIVALS,
      "source_gate":gate["pairwise_source"],
      "individual_models":{},
      "multivariate_all_five_control_model_was_not_run":True,
      "waterborne_effluent_N_not_measured":True,
      "interaction_or_causal_mediation_not_identified":True}
    summaries={}
    main_p={}
    control_p={}
    for rival in RIVALS:
        sg=gate["pairwise_source"][rival]
        if not sg.get("support_pass"):
            output["individual_models"][rival]={"status":"STOP_PAIR_SOURCE_SUPPORT",
                                                "source":sg}
            continue
        use=source.loc[source["source_pair_eligible_"+rival].astype(str)=="True",
            ["h","period","n_pair_"+rival,"cov_pair_"+rival,
             "pair_N_"+rival,"pair_"+rival]].copy()
        use=use.rename(columns={
            "n_pair_"+rival:"n_pair_scored",
            "cov_pair_"+rival:"pair_coverage",
            "pair_N_"+rival:"pair_N",
            "pair_"+rival:"pair_other"})
        d=df.merge(use,on=["h","period"],validate="one_to_one",how="inner")
        if len(d)!=sg["paired_support_H7period"] or d.h6.nunique()!=sg["paired_support_true_H6"]:
            raise RuntimeError("Original pair-specific source/H6 cohort changed "+rival)
        if d.country_holdout.nunique()!=sg["paired_support_country_groups"]:
            raise RuntimeError("Original country cohort changed "+rival)
        if (d.n_pair_scored.min()<10 or d.pair_coverage.min()<.4 or
            not np.isfinite(d[["pair_N","pair_other"]]).all().all()):
            raise RuntimeError("Invalid paired exact scored source quality "+rival)
        d["log1p_pair_scored"]=np.log1p(d.n_pair_scored.astype(float))
        d["country_method"]=d.country.fillna("__UNKNOWN__").astype(str)+"::"+d.method.fillna("__UNKNOWN__").astype(str)
        # The frozen generic estimator reports Plant_N_incidence source SD only
        # as metadata. Here that exact value is the pair-shared-species NUTRI.
        d["Plant_N_incidence"]=d.pair_N.astype(float)
        # Four balanced within-same-source comparisons, never compare their
        # *coefficients between rival-specific different cohorts*.
        N_single=source_cluster_regression(d,PRIMARY,BASE+["pair_N"],"pair_N")
        O_single=source_cluster_regression(d,PRIMARY,BASE+["pair_other"],"pair_other")
        N_after_other=source_cluster_regression(d,PRIMARY,BASE+["pair_other","pair_N"],"pair_N")
        other_after_N=source_cluster_regression(d,PRIMARY,BASE+["pair_N","pair_other"],"pair_other")
        Ncv=cv_h6(d,PRIMARY,BASE,"pair_N")
        Ocv=cv_h6(d,PRIMARY,BASE,"pair_other")
        Nextra=cv_h6(d,PRIMARY,BASE+["pair_other"],"pair_N")
        Oextra=cv_h6(d,PRIMARY,BASE+["pair_N"],"pair_other")
        outside=country_predictive(d,PRIMARY,BASE+["pair_other"],"pair_N")
        corr=float(np.corrcoef(d.pair_N,d.pair_other)[0,1])
        vif=vif_N(d)
        output["individual_models"][rival]={
            "status":"SOURCE_GATE_PASS_POSTHOC_PAIRWISE_OBSERVATION",
            "n_H7period":len(d),"n_H6":int(d.h6.nunique()),
            "n_country_groups":int(d.country_holdout.nunique()),
            "fixed_common_scored_species_keys_in_source_union":sg["scored_species_keys_original_11511_union"],
            "NUTRI_vs_other_pair_CWM_correlation":corr,
            "NUTRI_regressor_weighted_VIF_after_other_and_source_controls":vif,
            "VIF_gt_10_warning":bool(vif>10),
            "NUTRI_single":{"WLS_H6_beta":N_single,"H6_pred_extra":Ncv},
            "other_single":{"WLS_H6_beta":O_single,"H6_pred_extra":Ocv},
            "NUTRI_added_after_other":{"WLS_H6_beta":N_after_other,
                                       "H6_pred_extra":Nextra,
                                       "H6_country_out_extra":outside},
            "other_added_after_NUTRI":{"WLS_H6_beta":other_after_N,
                                        "H6_pred_extra":Oextra},
            "these_beta_MSE_are_not_independent_studies":True}
        main_p[rival]=N_after_other["nominal_p_descriptive"]
        control_p[rival]=other_after_N["nominal_p_descriptive"]
        summaries[rival]={
            "n_H7period":len(d),"n_H6":int(d.h6.nunique()),
            "NUTRI_vs_other_CWM_r":corr,"NUTRI_VIF":vif,
            "NUTRI_after_other_beta":N_after_other["adjusted_beta_per_1SD_focal"],
            "NUTRI_after_other_H6_cluster95":N_after_other["H6_cluster_descriptive95"],
            "NUTRI_after_other_H6_MSE_gain":Nextra["H6_equal_added_MSE_gain_pct"],
            "NUTRI_after_other_H6_gain95":Nextra["fixed_out_of_fold_error_H6_resampling95_pct"],
            "NUTRI_after_other_country_gain":outside.get("country_equal_added_MSE_gain_pct"),
            "NUTRI_after_other_country_gain95":outside.get("country_equal_fixed_errors_bootstrap95"),
            "other_after_NUTRI_beta":other_after_N["adjusted_beta_per_1SD_focal"],
            "other_after_NUTRI_H6_cluster95":other_after_N["H6_cluster_descriptive95"],
            "other_after_NUTRI_H6_MSE_gain":Oextra["H6_equal_added_MSE_gain_pct"],
            "other_after_NUTRI_H6_gain95":Oextra["fixed_out_of_fold_error_H6_resampling95_pct"]}
    q_N=bh(main_p) if main_p else {}
    q_other=bh(control_p) if control_p else {}
    for r in summaries:
        summaries[r]["NUTRI_after_other_BH5_q"]=q_N[r]
        summaries[r]["other_after_NUTRI_BH5_q"]=q_other[r]
        output["individual_models"][r]["NUTRI_added_after_other"]["BH5_q_for_beta_N"]=q_N[r]
        output["individual_models"][r]["other_added_after_NUTRI"]["BH5_q_for_beta_other"]=q_other[r]
    output["five_N_after_other_BH_q"]=q_N
    output["five_other_after_N_BH_q"]=q_other
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"EU454_EIV_NUTRI_FIVE_PAIRWISE_DIAGNOSTICS.json").write_text(
        json.dumps(output,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_EIV_FIVE_PAIR_NONCAUSAL_RESULT="+json.dumps(
        {"scope":output["scope"],"comparison":summaries},sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
