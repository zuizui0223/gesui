#!/usr/bin/env python3
"""Posthoc, uniformly specified 3-compartment EU H7 period models.

Single original native 451 H7-periods / 90 independent H6 with complete
aquatic BQE-EQR, plant EIV NUTRI and butterfly candidate-host HNS.

Source-positive WWTP reporting is NOT a measured receiving-water dose.
W is modeled static wastewater geography. This standardized screen is
a descriptive, noncausal SAME-GRAIN comparison, not a third-party new
independent site sampling, true interaction or mediation model.

H6 disjoint fivefold Ridge(alpha=10), H6-equal WLS and cluster normal CI,
H6 fixed out-of-fold prediction-error resampling 2000 draws (not refit CI).
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

UNIT_SOURCE="results/20261116_evidence_plane/EU_H7_PERIOD_454_WATER_PLANT_BUTTERFLY_WWTP_COMMON_PLANE.csv"
QA_SOURCE="results/20261116_evidence_plane/EU_H7_PERIOD_454_COMMON_PLANE_SOURCE_QA.json"
N=451
N_H6=90
N_FOLDS=5
BOOTSTRAP=2000
RIDGE=10.
SEED=20261116
BASE=[
    "crop","urban","riverflow","nEEA","nPlant","nButterfly","nAnimal",
    "primary_record_coverage_raw","primary_host_completeness",
    "NUTRI_coverage_all_Plantae_speciesKeys","log1p_source_scored_plant_species"
]
CATEGORICAL=["country_method","period"]
W="modelled_static_wastewater_geography_W"
EQR="aquatic_EEA_original_EQR"
PLANT="plant_EIV_NUTRI_equal_species"
HNS="butterfly_hostN_HNS_equal_species"
MODELS=[
  ("W_to_aquatic_EQR",EQR,BASE,W),
  ("W_to_plant_NUTRI",PLANT,BASE,W),
  ("W_to_butterfly_HNS",HNS,BASE,W),
  ("plant_NUTRI_added_to_EQR_after_W",EQR,BASE+[W],PLANT),
  ("plant_NUTRI_added_to_butterfly_after_W",HNS,BASE+[W],PLANT),
  ("aquatic_EQR_added_to_butterfly_after_plantN_and_W",HNS,BASE+[W,PLANT],EQR)
]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def check_input():
    panel=Path(UNIT_SOURCE)
    receipt=json.loads(Path(QA_SOURCE).read_text())
    if (receipt.get("schema")!="GESUI_20261116_TRUE_EU_H7xPERIOD_WATER_PLANT_BUTTERFLY_PLANE_V1"
        or sha(panel)!=receipt.get("generated_common_plane_csv_SHA256")):
        raise RuntimeError("Original source-integrity-gated 454 H7-period plane drift")
    df=pd.read_csv(panel,dtype={"h":str,"h6":str,"period":str,"country":str,"method":str})
    if (len(df),df.h.nunique(),df.h6.nunique())!=(454,232,90):
        raise RuntimeError("Original 454×232×90 source geography drift")
    if df.duplicated(["h","period"]).any():raise RuntimeError("H7-period duplicates")
    cols=[EQR,PLANT,HNS]
    use=df.loc[df[cols].notna().all(axis=1)].copy().reset_index(drop=True)
    if len(use)!=N or use.h6.nunique()!=N_H6:
        raise RuntimeError("Original 451 complete three-layer H7-periods/90 H6 drift")
    if not use.plant_EIV_original_source_eligible.astype(bool).all():
        raise RuntimeError("EIV missingness artificially included")
    if int(df.EEA_upstream_positive_report_matched_flag.sum())!=56:
        raise RuntimeError("Original report positive source count changed")
    use["log1p_source_scored_plant_species"]=np.log1p(
        use["exact_NUTRI_scored_vascular_accepted_species"].to_numpy(float))
    use["country_method"]=use.country.fillna("__UNKNOWN__").astype(str)+"::"+use.method.fillna("__UNKNOWN__").astype(str)
    required=BASE+[W,EQR,PLANT,HNS]
    if not np.isfinite(use[required].to_numpy(float)).all():
        raise RuntimeError("Nonfinite required original site controls/outcomes")
    for name in required:
        if use[name].std(ddof=0)<1e-11:
            raise RuntimeError("Prespecified common source covariate lacks variation: "+name)
    return use,receipt

def weights(data):
    group=data.h6.astype(str)
    counts=group.map(group.value_counts()).to_numpy(float)
    return len(data)/(counts*group.nunique())

def fixed_design(train,test,numeric):
    a=[np.ones(len(train))];b=[np.ones(len(test))];names=["intercept"]
    for name in numeric:
        tr=train[name].to_numpy(float);te=test[name].to_numpy(float)
        mu=tr.mean();sd=tr.std(ddof=0)
        if sd<1e-11:raise RuntimeError("Training-only numeric covariate constant "+name)
        a.append((tr-mu)/sd);b.append((te-mu)/sd);names.append(name)
    for name in CATEGORICAL:
        tr=train[name].astype(str);te=test[name].astype(str)
        for lev in sorted(tr.unique())[1:]:
            a.append((tr==lev).to_numpy(float));b.append((te==lev).to_numpy(float))
            names.append(name+":"+lev)
    return np.column_stack(a),np.column_stack(b),names

def full_fit(data,target,base,focus):
    X,_,names=fixed_design(data,data,base+[focus])
    if np.linalg.matrix_rank(X)!=X.shape[1]:
        raise RuntimeError("Full-rank source model invalid: "+target+" / "+focus)
    y=data[target].to_numpy(float)
    ysd=y.std(ddof=0)
    ystandard=(y-y.mean())/ysd
    w=weights(data)
    model=sm.WLS(ystandard,X,weights=w).fit(
        cov_type="cluster",cov_kwds={"groups":data.h6.to_numpy(),"use_correction":True},
        use_t=False)
    idx=names.index(focus)
    beta=float(model.params[idx]);se=float(model.bse[idx])
    return {"beta_outcome_SD_per_predictor_SD":beta,
            "H6_cluster_descriptive95":[beta-1.96*se,beta+1.96*se],
            "H6_cluster_nominal_p_exploratory":float(model.pvalues[idx]),
            "outcome_source_SD":float(ysd),
            "predictor_source_SD":float(data[focus].std(ddof=0)),
            "n_controls_and_dummy_columns":len(names)-2,
            "WLS_H6_equal":True,
            "joint_design_condition_number":float(np.linalg.cond(X*np.sqrt(w)[:,None]))}

def fitted_ridge(train,test,target,numeric):
    A,B,names=fixed_design(train,test,numeric)
    w=weights(train)
    gram=A.T@(w[:,None]*A)
    penalty=np.eye(A.shape[1])*RIDGE
    penalty[0,0]=0
    if np.linalg.cond(gram+penalty)>1e12:
        raise RuntimeError("Source ridge matrix nearly singular")
    coef=np.linalg.solve(gram+penalty,A.T@(w*train[target].to_numpy(float)))
    return B@coef

def fold_id(h):
    return int(hashlib.sha256(str(h).encode()).hexdigest()[:8],16)%N_FOLDS

def h6_cv(data,target,base,focus):
    groups=data.h6.astype(str).to_numpy()
    fold=np.array([fold_id(x) for x in groups])
    for group in sorted(set(groups)):
        if len(set(fold[groups==group]))!=1:raise RuntimeError("H6 leakage")
    y=data[target].to_numpy(float)
    pred0=np.full(len(data),np.nan);pred1=np.full(len(data),np.nan)
    fold_support=[]
    for j in range(N_FOLDS):
        te=np.flatnonzero(fold==j);tr=np.flatnonzero(fold!=j)
        if len(te)<5 or len(tr)<50:raise RuntimeError("Unbalanced original H6 CV split")
        if set(groups[te])&set(groups[tr]):raise RuntimeError("H6 parent train test overlap")
        fit=data.iloc[tr].reset_index(drop=True)
        test=data.iloc[te].reset_index(drop=True)
        pred0[te]=fitted_ridge(fit,test,target,base)
        pred1[te]=fitted_ridge(fit,test,target,base+[focus])
        fold_support.append({"fold":j,"test_original_H7period":len(te),
                             "test_distinct_H6":len(set(groups[te]))})
    if not np.isfinite(pred0).all() or not np.isfinite(pred1).all():
        raise RuntimeError("Missing H6-out predictions")
    err0=(y-pred0)**2;err1=(y-pred1)**2
    h6=np.array(sorted(set(groups)))
    per0=np.array([err0[groups==h].mean() for h in h6])
    per1=np.array([err1[groups==h].mean() for h in h6])
    b0=per0.mean();b1=per1.mean()
    if b0<=1e-12:raise RuntimeError("Baseline H6 squared error too small")
    gain=float(100*(b0-b1)/b0)
    rng=np.random.default_rng(SEED)
    ix=rng.integers(0,len(h6),size=(BOOTSTRAP,len(h6)))
    e0=per0[ix].mean(axis=1);e1=per1[ix].mean(axis=1)
    interval=np.quantile(100*(e0-e1)/e0,[.025,.975])
    return {"H6_equal_baseline_MSE":float(b0),
            "H6_equal_added_predictor_MSE":float(b1),
            "H6_out_Ridge10_MSE_gain_pct":gain,
            "fixed_H6_error_bootstrap_descriptive95_pct":interval.tolist(),
            "H6_nonoverlap_fivefold_source_counts":fold_support,
            "n_real_H6":len(h6),
            "no_model_refit_bootstrap":True,
            "train_only_predictor_scaling":True}

def main(args):
    df,receipt=check_input()
    results=[]
    for label,target,base,focus in MODELS:
        fit=full_fit(df,target,base,focus)
        cv=h6_cv(df,target,base,focus)
        results.append({"id":label,"outcome":target,"focal_source":focus,
            "baseline_predictors":base,
            "observational_H6_cluster_standardized":fit,"heldout_H6_CV":cv,
            "not_causal_or_independent_confirmation":True})
    out={"schema":"GESUI_20261116_SAME451_UNIFORM_THREE_COMPARTMENT_H6_SCREEN_V1",
         "source_input_SHA256":sha(Path(UNIT_SOURCE)),
         "source_QA_SHA256":sha(Path(QA_SOURCE)),
         "n_exact_same_H7period":len(df),
         "n_original_H7":int(df.h.nunique()),
         "n_true_H6":int(df.h6.nunique()),
         "original_WWTP_annual_positive_report_periods_in_full454":56,
         "all_models_exact_same_original_rows":True,
         "unified_source_controls":BASE,
         "categorical_controls":CATEGORICAL,
         "models":results,
         "posthoc_source_aligned_description_not_new_independent_field_data":True,
         "W_is_static_modelled_wastewater_geography_NOT_river_N":True,
         "source_EEA_actual_reported_NP_is_NOT_exposure_in_these_standardized_models":True,
         "US_NRSA_is_NOT_included_as_if_H7":True,
         "H6_holdout_not_year_independent":True,
         "not_biotic_interaction_or_causal_mediation":True}
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"EU_H7_451_UNIFORM_W_AND_PAIRWISE_H6_COMPARISON.json").write_text(
        json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    rows=[]
    for item in results:
        w=item["observational_H6_cluster_standardized"]
        h=item["heldout_H6_CV"]
        rows.append({"id":item["id"],"outcome":item["outcome"],"focal_source":item["focal_source"],
            "n_H7period":len(df),"true_H6":int(df.h6.nunique()),
            "beta_outcomeSD_per_predictorSD":w["beta_outcome_SD_per_predictor_SD"],
            "H6_cluster_95_low":w["H6_cluster_descriptive95"][0],
            "H6_cluster_95_high":w["H6_cluster_descriptive95"][1],
            "H6_out_extra_predictor_MSE_gain_pct":h["H6_out_Ridge10_MSE_gain_pct"],
            "H6_fixed_errors_gain95_low":h["fixed_H6_error_bootstrap_descriptive95_pct"][0],
            "H6_fixed_errors_gain95_high":h["fixed_H6_error_bootstrap_descriptive95_pct"][1],
            "same_unit_not_sewage_causal":True})
    pd.DataFrame(rows).to_csv(args.output/"EU_H7_451_UNIFORM_W_AND_PAIRWISE_H6_COMPARISON.csv",index=False)
    print("GESUI_ONE_PLANE_SAME451_JOINT_STANDARDIZED_H6="+json.dumps(
        {"source_H7period":len(df),"source_true_H6":int(df.h6.nunique()),
         "outcomes":[{"name":r["id"],
          "beta":r["observational_H6_cluster_standardized"]["beta_outcome_SD_per_predictor_SD"],
          "beta95":r["observational_H6_cluster_standardized"]["H6_cluster_descriptive95"],
          "CV_gain":r["heldout_H6_CV"]["H6_out_Ridge10_MSE_gain_pct"],
          "CV_95":r["heldout_H6_CV"]["fixed_H6_error_bootstrap_descriptive95_pct"]}
         for r in results]},ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
