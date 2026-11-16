#!/usr/bin/env python3
"""Post-result adversarial robustness on original 183 H7 / 40 actual H6 EIV/HNS.

Full 13-controls H6 FE source slope is reproduced exactly before comparing:
- old 30 multi-H7 parents only, same source controls and 9999 H6-within null
- 40-H6 full model reduced/no controls, H6 vs H7 equal weights
- 10 new parent H6 external out-of-parent contrast predictions FROM old30 H6,
  holding the non-H6 coefficients fixed and cancelling unseen parent intercepts
  by new-minus-old H7 differences (10 geographic pairs).

Not independent preregistration; no EEA water/river TN/TP, local host feeding,
butterfly individual abundance or causal waste N mediation is observed.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
from pathlib import Path
from collections import Counter
import math
import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261115_h6_source_expansion"))
from test_183_40_within_H6 import (
    clean_source,features,fit_once,weights,blocked_fl_null,
    PREDICTOR,Y_PRIMARY,Y_SECONDARY,CONTROLS,
    PANEL_SOURCE_SHA,SOURCE_SCHEMA,SEED,N_PERM
)

THREE=["W","crop","urban"]
MODEL_SETS={"all_13":CONTROLS,"environment_three":THREE,"H6_only":[]}
SOURCE_FILE="GENUINE_TERRESTRIAL_ONLY_183_H7_40_H6_PAIRED_H7_MEANS.csv"
SOURCE_QA="GENUINE_TERRESTRIAL_183_H7_40_H6_SOURCE_RECONSTRUCTION_QA.json"
BETA_CHECK=0.11807615189599471
N_EXTERNAL=10
N_BOOT=2000

def fixed_full_centers(d):
    centers={}
    for c in CONTROLS+[PREDICTOR]:
        v=d[c].to_numpy(float)
        s=float(v.std(ddof=0))
        if s<1e-12:raise RuntimeError("Original full source predictor constant "+c)
        centers[c]=(float(v.mean()),s)
    return centers

def model(d,target,controls,center,weight_method="H6_equal",cluster=True):
    # Force within-genuine H6 fixed intercepts; do not put old country×method
    # EEA dummies into the new single-era terrestrial cohort.
    if len(d)<20 or d.h6.nunique()<10 or d.h.duplicated().any():
        raise RuntimeError("Invalid H7/H6 source design for FE model")
    if set(controls)-set(CONTROLS):
        raise RuntimeError("Non-precommitted covariate entered model")
    X=[np.ones(len(d))];names=["intercept"]
    for c in list(controls)+[PREDICTOR]:
        mu,sd=center[c]
        X.append((d[c].to_numpy(float)-mu)/sd)
        names.append(c)
    h=d.h6.astype(str)
    for g in sorted(h.unique())[1:]:
        X.append((h==g).to_numpy(float))
        names.append("H6_FE:"+g)
    X=np.column_stack(X)
    rank=int(np.linalg.matrix_rank(X))
    if rank!=X.shape[1]:
        return {"status":"STOP_COLLINEAR_DESIGN","rank":rank,"ncols":X.shape[1]}
    w=weights(d) if weight_method=="H6_equal" else np.ones(len(d))
    if cluster:
        result=sm.WLS(d[target].to_numpy(float),X,weights=w).fit(
            cov_type="cluster",cov_kwds={"groups":h.to_numpy(),"use_correction":True},use_t=False)
    else:
        result=sm.WLS(d[target].to_numpy(float),X,weights=w).fit()
    j=names.index(PREDICTOR)
    beta=float(result.params[j]);se=float(result.bse[j])
    return {"status":"FIT",
            "beta_full_source_plant_N_1sd":beta,
            "H6_cluster_descriptive95":[beta-1.96*se,beta+1.96*se] if cluster else None,
            "H6_cluster_nominal_p":float(result.pvalues[j]) if cluster else None,
            "design_condition":float(np.linalg.cond(X*np.sqrt(w)[:,None])),
            "n_original_H7":len(d),"n_true_H6":int(d.h6.nunique()),
            "source_controls":list(controls),"weight":weight_method,
            "names":names,"X":X,"w":w}

def fit_metrics(d,target,controls,center,weight_method):
    result=model(d,target,controls,center,weight_method)
    if result["status"]!="FIT":return result
    for key in ("names","X","w"):result.pop(key)
    return result

def old30_permutation(old,center):
    fit=model(old,Y_PRIMARY,CONTROLS,center,"H6_equal")
    if fit["status"]!="FIT":return {"model":fit}
    X,names,w=fit.pop("X"),fit.pop("names"),fit.pop("w")
    null=blocked_fl_null(old,old[Y_PRIMARY].to_numpy(float),X,names,w,SEED,N_PERM)
    if abs(null["source_residualized_plant_N_WLS_beta"]-fit["beta_full_source_plant_N_1sd"])>1e-8:
        raise RuntimeError("Independent old30 FWL slope does not match WLS")
    return {"cluster_fit":fit,"H6_within_residual_label_null":null}

def withinH6_vif(d,center):
    full=model(d,Y_PRIMARY,CONTROLS,center)
    if full["status"]!="FIT":return {"status":"STOP_COLLINEAR"}
    X,names,w=full["X"],full["names"],full["w"]
    j=names.index(PREDICTOR)
    x=X[:,j]
    root=np.sqrt(w)
    groups=d.h6.astype(str)
    fixed=np.column_stack([np.ones(len(d))]+[(groups==g).to_numpy(float)
        for g in sorted(groups.unique())[1:]])
    xf=x-fixed@np.linalg.lstsq(fixed*root[:,None],x*root,rcond=None)[0]
    X0=np.delete(X,j,axis=1)
    xr=x-X0@np.linalg.lstsq(X0*root[:,None],x*root,rcond=None)[0]
    dfe=float(np.sum(w*xf**2))
    drem=float(np.sum(w*xr**2))
    vif=dfe/drem if drem>1e-12 else float("inf")
    return {"post_H6_fixed_effect_source_VIF_for_plantN_after_13controls":vif,
            "original_design_condition":full["design_condition"],
            "source_within_H6_x_variance":dfe,
            "source_after_H6_plus13_x_variance":drem,
            "warning_VIF_gt10":bool(vif>10),
            "warning_condition_gt100":bool(full["design_condition"]>100)}

def regression_train(old,target,controls,add_plant):
    # Train-only transformations. FE coefficients never require extrapolated
    # categorical labels on new H6: test is pair DIFFERENCE within H6.
    features=list(controls)+([PREDICTOR] if add_plant else [])
    scaling={}
    A=[np.ones(len(old))];names=["intercept"]
    for name in features:
        v=old[name].to_numpy(float)
        mu=float(v.mean());sd=float(v.std(ddof=0))
        if sd<1e-12:raise RuntimeError("Training controls unexpectedly constant "+name)
        scaling[name]=(mu,sd)
        A.append((v-mu)/sd);names.append(name)
    groups=old.h6.astype(str)
    for h in sorted(groups.unique())[1:]:
        A.append((groups==h).to_numpy(float));names.append("H6_FE:"+h)
    X=np.column_stack(A);w=weights(old)
    if np.linalg.matrix_rank(X)!=X.shape[1]:
        raise RuntimeError("Old30 163 H7 source model not full rank in transfer")
    coef=np.linalg.lstsq(X*np.sqrt(w)[:,None],old[target].to_numpy(float)*np.sqrt(w),rcond=None)[0]
    return {n:float(coef[i]) for i,n in enumerate(names) if n in features},scaling

def numeric_corr(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    if np.std(a)<1e-12 or np.std(b)<1e-12:return None
    return float(np.corrcoef(a,b)[0,1])

def external_pairs(old,new,target):
    if len(new)!=20 or new.h6.nunique()!=10:raise RuntimeError("Expected ten unseen paired H6 with exactly two H7")
    by=new.groupby("h6")
    if any(len(g)!=2 for _,g in by):
        raise RuntimeError("A new parent H6 does not have one genuine new H7 and one old H7")
    per=[]
    for h6,g in by:
        oldrow=g.loc[g.original_source_cohort=="original_old_H7"]
        newrow=g.loc[g.original_source_cohort!="original_old_H7"]
        if len(oldrow)!=1 or len(newrow)!=1:
            raise RuntimeError("New H6 does not contain one frozen old/new original H7 pair")
        o=oldrow.iloc[0];n=newrow.iloc[0]
        per.append({"h6":str(h6),"old_H7":str(o.h),"new_H7":str(n.h),
            "delta_plant_N":float(n[PREDICTOR]-o[PREDICTOR]),
            "delta_y":float(n[target]-o[target]),
            "delta_source_controls":{c:float(n[c]-o[c]) for c in CONTROLS}})
    if len(per)!=10 or len({p["new_H7"] for p in per})!=10:
        raise RuntimeError("External pairs not ten unique native H7")
    c0,s0=regression_train(old,target,CONTROLS,False)
    c1,s1=regression_train(old,target,CONTROLS,True)
    pred0=[];pred1=[];yy=[];x=[]
    for p in per:
        y=p["delta_y"];yy.append(y);x.append(p["delta_plant_N"])
        a=sum(p["delta_source_controls"][k]*c0[k]/s0[k][1] for k in CONTROLS)
        b=sum(p["delta_source_controls"][k]*c1[k]/s1[k][1] for k in CONTROLS)
        b+=p["delta_plant_N"]*c1[PREDICTOR]/s1[PREDICTOR][1]
        pred0.append(float(a));pred1.append(float(b))
        p["pred_delta_baseline"]=float(a);p["pred_delta_with_plantN"]=float(b)
        p["squared_error_baseline"]=float((y-a)**2)
        p["squared_error_with_plantN"]=float((y-b)**2)
    y=np.array(yy);a=np.array(pred0);b=np.array(pred1);x=np.array(x)
    e0=(y-a)**2;e1=(y-b)**2
    m0=float(e0.mean());m1=float(e1.mean())
    if m0<1e-10:raise RuntimeError("External 10pair baseline MSE trivial")
    improvement=100*(m0-m1)/m0
    rng=np.random.default_rng(SEED)
    indices=rng.integers(0,len(per),size=(N_BOOT,len(per)))
    e0draw=e0[indices].mean(axis=1)
    e1draw=e1[indices].mean(axis=1)
    bt=100*(e0draw-e1draw)/e0draw
    sign_potential=int(np.sum((x*(y-a))>0))
    notzero=int(np.sum(np.abs(x*(y-a))>1e-10))
    return {"n_source_independent_parent_H6":len(per),"n_test_distinct_H7":2*len(per),
            "source_original_external_true_H6_train_overlap":0,
            "n_train_H6":int(old.h6.nunique()),"n_train_H7":len(old),
            "external_pair_differences_outcome_mean":float(y.mean()),
            "external_pair_differences_plantN_mean":float(x.mean()),
            "out_of_parent_pair_MSE_baseline":m0,
            "out_of_parent_pair_MSE_with_plantN":m1,
            "out_of_parent_pair_added_predictive_MSE_gain_percent":float(improvement),
            "fixed_10_pair_resampling_descriptive95":np.quantile(bt,[.025,.975]).tolist(),
            "deltaPlantN_deltaHNS_raw_pearson":numeric_corr(x,y),
            "deltaPlantN_baseline_control_residual_pearson":numeric_corr(x,y-a),
            "pair_direction_concordance_after_controls":sign_potential,
            "n_nonzero_pair_differences_after_controls":notzero,
            "pairs":per,
            "baseline_train_only_13_controls":True,
            "source_full_model_train_only_13_controls_plus_EIV":True,
            "genuine_H6_intercept_cancels_in_new_minus_old_pair_differences":True,
            "all_10_pairs_frozen_before_this_posthoc_outcome_check":True,
            "fixed_error_interval_not_model_refit":True,
            "not_an_independent_experiment_or_true_wastewater_effect":True}

def main(args):
    panel=args.source/SOURCE_FILE
    qa=json.loads((args.source/SOURCE_QA).read_text())
    if qa.get("schema")!=SOURCE_SCHEMA:
        raise RuntimeError("Source reconstruction did not carry original 183 H7 SHA")
    d=clean_source(panel,qa)
    d=d.sort_values(["h6","h"],kind="stable").reset_index(drop=True)
    n_old=(d.original_source_cohort=="original_old_H7")
    if n_old.sum()!=173 or (~n_old).sum()!=10:
        raise RuntimeError("Source-old versus added H7 identities altered")
    newH6=set(d.loc[~n_old,"h6"])
    if len(newH6)!=10:raise RuntimeError("Ten new H7 do not identify ten distinct H6")
    old=d.loc[~d.h6.isin(newH6)].copy().reset_index(drop=True)
    external=d.loc[d.h6.isin(newH6)].copy().reset_index(drop=True)
    if (len(old),old.h6.nunique(),len(external),external.h6.nunique())!=(163,30,20,10):
        raise RuntimeError("Original 30 old multi-H7 H6 and 10 old/new parent groups drift")
    if not set(old.h6).isdisjoint(set(external.h6)):
        raise RuntimeError("Old30/new10 H6 train-test geography leaks")
    center=fixed_full_centers(d)
    original=fit_metrics(d,Y_PRIMARY,CONTROLS,center,"H6_equal")
    if original["status"]!="FIT" or abs(original["beta_full_source_plant_N_1sd"]-BETA_CHECK)>1e-8:
        raise RuntimeError("Previously authenticated full source 183-H7/40-H6 β changed")
    cases={}
    for name,cols in MODEL_SETS.items():
        cases[name]=fit_metrics(d,Y_PRIMARY,cols,center,"H6_equal")
    cases["all_13_H7_equal"]=fit_metrics(d,Y_PRIMARY,CONTROLS,center,"H7_equal")
    report={
        "schema":"GESUI_183H7_40H6_ADVERSARIAL_OLD30_VS_NEW10_TRUEH6_TRANSFER_V1",
        "source_panel_SHA256":PANEL_SOURCE_SHA,
        "source_original_QA_schema":SOURCE_SCHEMA,
        "n_all_source_distinct_H7":len(d),
        "n_all_true_H6":int(d.h6.nunique()),
        "n_old30_H7":len(old),"n_old30_true_H6":int(old.h6.nunique()),
        "n_new10_pairs_H7":len(external),"n_new10_true_H6":int(external.h6.nunique()),
        "source_original_H6_FE_beta_reproduced":original["beta_full_source_plant_N_1sd"],
        "full_40_H6_model_sensitivity":cases,
        "full_40_H6_partial_N_collinearity_diagnostic":withinH6_vif(d,center),
        "old30_only_same_13_controls":old30_permutation(old,center),
        "external_new10_H6_true_group_holdout":{
            Y_PRIMARY:external_pairs(old,external,Y_PRIMARY),
            Y_SECONDARY:external_pairs(old,external,Y_SECONDARY)},
        "within_H6_original_9999_perm_result_retained_not_retuned":True,
        "not_new_independent_confirmatory_study":True,
        "true_EEA_2018_and_2022_same_location_not_assumed_for_newH7":True,
        "not_plant_herbivore_or_waterborne_N_causal_mediation":True,
        "original_40H6_main_preoutcome_support_not_rewritten":True}
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"183H7_40H6_ADVERSARIAL_OLD30_NEW10_TRANSFER_RESULTS.json").write_text(
         json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    brief={"source_genuine_H7":len(d),"genuine_H6":int(d.h6.nunique()),
         "source_original_40H6_beta":original["beta_full_source_plant_N_1sd"],
         "full_models":{k:{"status":v["status"],"beta":v.get("beta_full_source_plant_N_1sd"),
                            "beta95":v.get("H6_cluster_descriptive95")} for k,v in cases.items()},
         "collinearity":report["full_40_H6_partial_N_collinearity_diagnostic"],
         "old30":{"beta":report["old30_only_same_13_controls"]["cluster_fit"].get("beta_full_source_plant_N_1sd"),
                  "95":report["old30_only_same_13_controls"]["cluster_fit"].get("H6_cluster_descriptive95"),
                  "blocked_perm":report["old30_only_same_13_controls"]["H6_within_residual_label_null"].get("two_sided_H6_blocked_Freedman_Lane_permutation_p")},
         "new10_external":{k:{"gain":v["out_of_parent_pair_added_predictive_MSE_gain_percent"],
                             "gain95":v["fixed_10_pair_resampling_descriptive95"],
                             "n_source_pairs":v["n_source_independent_parent_H6"],
                             "sign_concordance":v["pair_direction_concordance_after_controls"],
                             "raw_r":v["deltaPlantN_deltaHNS_raw_pearson"],
                             "residual_r":v["deltaPlantN_baseline_control_residual_pearson"]}
                            for k,v in report["external_new10_H6_true_group_holdout"].items()}}
    print("GESUI_ADVERSARIAL_OLD30_NEW10_WITHIN_H6="+json.dumps(brief,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    main(ap.parse_args())
