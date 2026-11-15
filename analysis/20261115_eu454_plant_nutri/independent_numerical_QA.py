#!/usr/bin/env python3
"""Numerically independent FWL weighted-source and augmented-LS Ridge audit.

Reads original true 454 source data and source-verified 451 EIV scores. Compares
all four predeclared coefficients and out-H6 prediction MSEs to original
statsmodels/normal-equations model artifact; neither fitting approach is reused.
No new ecological effects selected, no H6 leakage, no cause/mediator inference.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
SOURCE_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
P=["crop","urban","riverflow","naturalout","area","nEEA","nPlant",
   "NUTRI_coverage_all_Plantae_speciesKeys","log1p_eiv_scored_taxa"]
B=["crop","urban","riverflow","naturalout","area","nEEA","nAnimal","nPlant",
   "nButterfly","primary_record_coverage_raw","primary_host_completeness",
   "NUTRI_coverage_all_Plantae_speciesKeys","log1p_eiv_scored_taxa","W"]
CAT=("country_method","period")
TOL=1e-7

def numerical(df, cols, ref, fit):
    h=df.h6.astype(str).to_numpy()
    counts=df.h6.astype(str).value_counts()
    w=np.array([len(df)/(len(counts)*counts[z]) for z in h],dtype=float)
    def zz(column):
        a=df[column].to_numpy(float);std=np.std(a)
        if std<=1e-12:raise RuntimeError("Constant independent QA source "+column)
        return (a-np.mean(a))/std
    keep=["intercept"]+[x for x in cols if x!=ref]
    array=[np.ones(len(df))]+[zz(k) for k in cols if k!=ref]
    for cat in CAT:
        values=df[cat].astype(str)
        for name in sorted(values.unique().tolist())[1:]:
            keep.append(cat+"="+name);array.append((values==name).to_numpy(float))
    X=np.column_stack(array)
    r=zz(ref);y=df[fit].to_numpy(float);sw=np.sqrt(w)
    res_r=r-X@np.linalg.lstsq(X*sw[:,None],r*sw,rcond=None)[0]
    res_y=y-X@np.linalg.lstsq(X*sw[:,None],y*sw,rcond=None)[0]
    beta=np.dot(w*res_r,res_y)/np.dot(w*res_r,res_r)
    return float(beta)

def make_pair(tr,te,cols):
    xtr=[np.ones(len(tr))];xte=[np.ones(len(te))]
    for col in cols:
        u=tr[col].to_numpy(float);v=te[col].to_numpy(float)
        mu=float(np.mean(u));sd=float(np.std(u))
        if sd<1e-12:sd=1
        xtr.append((u-mu)/sd);xte.append((v-mu)/sd)
    for cat in CAT:
        observed=tr[cat].astype(str)
        for label in sorted(observed.unique().tolist())[1:]:
            xtr.append((observed==label).to_numpy(float))
            xte.append((te[cat].astype(str)==label).to_numpy(float))
    return np.column_stack(xtr),np.column_stack(xte)

def augmented_lstsq_predict(tr,te,ycol,columns):
    X,T=make_pair(tr,te,columns)
    g=tr.h6.astype(str).to_numpy()
    sizes=tr.h6.astype(str).value_counts()
    weights=np.array([len(tr)/(len(sizes)*sizes[z]) for z in g],float)
    sw=np.sqrt(weights)
    penalty=np.eye(X.shape[1])*np.sqrt(10.0)
    penalty[0,0]=0.
    # augmented penalized least squares, not the main model's normal equations
    aug=np.vstack((X*sw[:,None],penalty))
    response=np.r_[tr[ycol].to_numpy(float)*sw,np.zeros(X.shape[1])]
    coef=np.linalg.lstsq(aug,response,rcond=None)[0]
    return T@coef

def crossval(df,ycol,base,plus):
    h6=df.h6.astype(str).to_numpy()
    indices=np.array([int(hashlib.sha256(str(g).encode()).hexdigest()[:8],16)%5 for g in h6])
    p0=np.full(len(df),np.nan);p1=np.full(len(df),np.nan)
    for fold in range(5):
        tr=np.flatnonzero(indices!=fold);te=np.flatnonzero(indices==fold)
        if set(h6[tr])&set(h6[te]):raise RuntimeError("True H6 leaked!")
        a=df.iloc[tr].reset_index(drop=True)
        b=df.iloc[te].reset_index(drop=True)
        p0[te]=augmented_lstsq_predict(a,b,ycol,base)
        p1[te]=augmented_lstsq_predict(a,b,ycol,base+[plus])
    if not np.isfinite(p0).all() or not np.isfinite(p1).all():
        raise RuntimeError("Incomplete original H6 predictions")
    y=df[ycol].to_numpy(float)
    by=[]
    for g in sorted(set(h6)):
        ix=h6==g
        by.append([np.mean((y[ix]-p0[ix])**2),np.mean((y[ix]-p1[ix])**2)])
    a,b=np.array(by).T
    return float(100*(a.mean()-b.mean())/a.mean())

def main(args):
    original=ROOT/"results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"
    if hashlib.sha256(original.read_bytes()).hexdigest()!=SOURCE_SHA:
        raise RuntimeError("Original 454 source panel checksum mismatch")
    gate=json.loads((args.source/"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if not gate.get("hard_source_support_pass") or gate.get("schema")!="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA":
        raise RuntimeError("No full EIV source gate")
    main=json.loads((args.main/"EU454_PLANT_HNS_MODELS.json").read_text())
    if main.get("n_analysis_H7period")!=451 or main.get("n_true_H6")!=90:
        raise RuntimeError("Expected original complete H7 model source")
    p=pd.read_csv(original,dtype={"h":str,"h6":str})
    e=pd.read_csv(args.source/"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv",dtype={"h":str})
    v=pd.read_csv(args.source/"EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv",dtype={"h":str})
    e=e.loc[e.full_source_trait_eligible.astype(str)=="True"]
    if len(e)!=451 or len(v)!=451:raise RuntimeError("Original source eligibility changed")
    df=p.merge(e[["h","period","NUTRI_coverage_all_Plantae_speciesKeys",
                     "exact_NUTRI_scored_vascular_accepted_species"]],
               on=["h","period"],validate="one_to_one",how="inner")
    df=df.merge(v[["h","period","Plant_N_incidence","Plant_N_record"]],
                on=["h","period"],validate="one_to_one",how="inner")
    df["log1p_eiv_scored_taxa"]=np.log1p(df.exact_NUTRI_scored_vascular_accepted_species.to_numpy(float))
    df["country_method"]=(df.country.fillna("__UNKNOWN__").astype(str)+"::"+
                          df.method.fillna("__UNKNOWN__").astype(str))
    if len(df)!=451 or df.h6.nunique()!=90:raise RuntimeError("Source cohort mismatch")
    tests=[
      ("W_to_plant_N","Plant_N_incidence",P,"W"),
      ("W_to_plant_N","Plant_N_record",P,"W"),
      ("plant_N_added_to_butterfly_HNS","primary_incidence_mean",B,"Plant_N_incidence"),
      ("plant_N_added_to_butterfly_HNS","primary_record_mean",B,"Plant_N_incidence")]
    results=[]
    for section,outcome,controls,focus in tests:
        if section=="W_to_plant_N":
            m=main[section][outcome]; expected=m["association"]["adjusted_beta_per_1SD_focal"]
            expected_cv=m["H6_out_W_increment"]["H6_equal_added_MSE_gain_pct"]
        else:
            m=main[section][outcome]; expected=m["association"]["adjusted_beta_per_1SD_focal"]
            expected_cv=m["H6_out_PlantN_increment"]["H6_equal_added_MSE_gain_pct"]
        beta=numerical(df,controls+[focus],focus,outcome)
        gain=crossval(df,outcome,controls,focus)
        errb=abs(beta-expected);errcv=abs(gain-expected_cv)
        if errb>TOL or errcv>TOL:
            raise RuntimeError(f"Independent weighted-FWL/Ridge audit failed: {section}/{outcome}: {errb}, {errcv}")
        results.append({"model":section,"outcome":outcome,"beta_independent_FWL":beta,
                        "H6_add_MSE_gain_independent_augmentedLS":gain,
                        "beta_abs_error":errb,"H6_added_gain_abs_error":errcv,
                        "no_H6_train_test_overlap":True})
    out={"schema":"GESUI_EU454_FOURMODELS_FWL_AUGMENTEDRIDGE_INDEPENDENT_QA_V1",
         "status":"PASSED_ALL_FOUR_ORIGINAL_SOURCE_NUMERIC_AUDITS",
         "n_H7period":451,"n_true_H6":90,
         "max_abs_beta_error":max(x["beta_abs_error"] for x in results),
         "max_abs_holdout_percent_gain_error":max(x["H6_added_gain_abs_error"] for x in results),
         "independent_reconstruction_checks":results,
         "no_new_causal_or_ecology_inference":True}
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"EU454_PLANT_HNS_NUMERIC_INDEPENDENT_QA.json").write_text(
        json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_EU454_FOUR_MODEL_INDEPENDENT_QA="+json.dumps({
        k:out[k] for k in ["status","n_H7period","n_true_H6","max_abs_beta_error","max_abs_holdout_percent_gain_error"]
        },ensure_ascii=False),flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",required=True,type=Path)
    parser.add_argument("--main",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    main(parser.parse_args())

