#!/usr/bin/env python3
"""Posthoc fixed-effect within-true-H6 terrestrial plant N vs butterfly HNS.

Original 183 independent source-locked H7 (40 genuine HydroBASINS H6),
each aggregated over exactly two GBIF years/period bins, genuine original
EIV NUTRI and butterfly host HNS. This is NOT an EEA water model, actual
wastewater-N dose, ecological biomass, species occupancy, or causal mediation.

Primary: H6-equal H6-fixed-intercept WLS, cluster robust 40 genuine H6,
9,999 within-H6 null-model residual-label Freedman–Lane draws.
Secondary: occurrence record-weighted HNS, 40 leave-one-parent-out fits.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

SOURCE_SCHEMA="GESUI_ORIGINAL_NEW_183H7_40H6_TERRESTRIAL_SOURCE_COHORT_AND_TRAIT_REPLAY_V1"
PANEL_SOURCE_SHA="c0840b1c5c6e85342fea41d57dc40fa5ec295fabf6bcd6908a4072557db9a736"
GEO_SOURCE_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
PREDICTOR="plant_n_incidence"
Y_PRIMARY="hns_incidence"
Y_SECONDARY="hns_record"
CONTROLS=["W","crop","urban","naturalout","area","nPlant","nButterfly",
          "nAnimal","plant_coverage","hns_coverage","hns_host_completeness",
          "plant_scored_species","hns_scored_species"]
SEED=20261115
N_PERM=9999

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def clean_source(p,qa):
    source=pd.read_csv(p,dtype={"h":str,"h6":str})
    if (len(source)!=183 or source.h.nunique()!=183 or
        source.h6.nunique()!=40 or source.h.duplicated().any()):
        raise RuntimeError("Not exactly frozen 183 H7/40 H6 genuine geography")
    if sha(p)!=PANEL_SOURCE_SHA or qa.get("data_source_panel_SHA256")!=PANEL_SOURCE_SHA:
        raise RuntimeError("Original reconstituted H7 plant/butterfly source score panel SHA has changed")
    if qa.get("old_native_panel_sha256")!=GEO_SOURCE_SHA or qa.get("original_old_score_mismatch_n")!=0:
        raise RuntimeError("Previously observed HNS / NUTRI source result mismatches")
    if qa.get("n_new_H7")!=10 or qa.get("n_old_H7")!=173:
        raise RuntimeError("New vs prior source cohort composition changed")
    if not set(["h","h6",PREDICTOR,Y_PRIMARY,Y_SECONDARY]+CONTROLS).issubset(source.columns):
        raise RuntimeError("Frozen terrestrial predictor and control column missing")
    vals=source[[PREDICTOR,Y_PRIMARY,Y_SECONDARY]+CONTROLS].to_numpy(dtype=float)
    if not np.isfinite(vals).all():
        raise RuntimeError("Genuine selected source numeric data contain missing/bad fields")
    counts=source.h6.value_counts()
    if counts.min()<2 or len(counts)!=40:
        raise RuntimeError("At least two distinct H7 within each true H6 is required")
    if (source.W<0).any() or (source.plant_coverage<.4).any():
        raise RuntimeError("Source W and original 40% plant coverage criterion failed")
    return source

def features(d,center=None):
    """Use fixed full-sample standardization even in leave-H6-out sensitivity."""
    if center is None:
        center={}
        for k in CONTROLS+[PREDICTOR]:
            x=d[k].to_numpy(float);mu=x.mean();sd=x.std(ddof=0)
            if sd<1e-12:raise RuntimeError("Fixed source covariate constant "+k)
            center[k]=(float(mu),float(sd))
    X=[np.ones(len(d))]
    names=["intercept"]
    for key in CONTROLS+[PREDICTOR]:
        mu,sd=center[key]
        X.append((d[key].to_numpy(float)-mu)/sd);names.append(key)
    h6=d.h6.astype(str)
    groups=sorted(h6.unique().tolist())
    for h in groups[1:]:
        X.append((h6==h).to_numpy(float));names.append("native_H6_FE:"+h)
    return np.column_stack(X),names,center

def weights(d):
    h=d.h6.astype(str)
    num=h.map(h.value_counts()).to_numpy(float)
    return len(d)/(h.nunique()*num)

def fit_once(d,target,center=None,cluster=True):
    X,names,center=features(d,center=center)
    yw=d[target].to_numpy(float)
    w=weights(d)
    if np.linalg.matrix_rank(X)!=X.shape[1]:
        raise RuntimeError(f"H6 within-parent true source model not full rank; {np.linalg.matrix_rank(X)}/{X.shape[1]} columns")
    if cluster:
        fit=sm.WLS(yw,X,weights=w).fit(
            cov_type="cluster",
            cov_kwds={"groups":d.h6.astype(str).to_numpy(),"use_correction":True},
            use_t=False)
    else:
        fit=sm.WLS(yw,X,weights=w).fit()
    j=names.index(PREDICTOR)
    beta=float(fit.params[j]);se=float(fit.bse[j])
    return fit,X,names,center,w,{
       "beta_NUTRI_per_full_sample_1SD":beta,
       "standard_error":se,
       "descriptive_H6_cluster_95":[beta-1.96*se,beta+1.96*se] if cluster else None,
       "nominal_two_sided_p":float(fit.pvalues[j]) if cluster else None,
       "true_H6_cluster_n":int(d.h6.nunique()),
       "different_H7_n":len(d),
       "source_model_design_columns":len(names),
       "source_model_design_rank":int(np.linalg.matrix_rank(X)),
       "source_design_condition_number":float(np.linalg.cond(X*np.sqrt(w)[:,None])),
       "not_causal":True}

def blocked_fl_null(d,y,Xfull,names,w,seed,nperm):
    """Freedman–Lane residual shuffle within each genuine H6 only.

    Original H6 group weights are identical within each parent, preserving
    weighting under permutation. Approximate null, not randomized exposure.
    """
    j=names.index(PREDICTOR)
    X0=np.delete(Xfull,j,axis=1)
    x=Xfull[:,j]
    root=np.sqrt(w)
    coefx=np.linalg.lstsq(X0*root[:,None],x*root,rcond=None)[0]
    rx=x-X0@coefx
    den=np.sum(w*rx**2)
    if den<=1e-10:raise RuntimeError("Plant N score has no identifiable within-H6 residual variation")
    beta=float(np.dot(w*rx,y)/den)
    coef0=np.linalg.lstsq(X0*root[:,None],y*root,rcond=None)[0]
    resid=y-X0@coef0
    rng=np.random.default_rng(seed)
    groups=[np.flatnonzero(d.h6.astype(str).to_numpy()==g) for g in sorted(d.h6.unique())]
    shuffle=np.empty_like(resid)
    obs=[]
    for t in range(nperm):
        for gi in groups:
            shuffle[gi]=rng.permutation(resid[gi])
        obs.append(float(np.dot(w*rx,shuffle)/den))
    null=np.asarray(obs,float)
    ext=int(np.sum(np.abs(null)>=abs(beta)-1e-12))
    return {
       "source_residualized_plant_N_WLS_beta":beta,
       "two_sided_H6_blocked_Freedman_Lane_permutation_p":float((ext+1)/(nperm+1)),
       "number_of_null_draws":nperm,
       "extreme_draws_abs_beta":ext,
       "null_2p5_to_97p5":np.quantile(null,[.025,.975]).tolist(),
       "null_median":float(np.median(null)),
       "null_SD":float(np.std(null)),
       "original_H6_groups_permute_residuals_only_within":len(groups),
       "null_approximate_residual_exchangeability_not_experiment":True,
       "no_geographic_or_temporal_H6_crossleakage":True}

def describe_one(d,target):
    fit,X,names,center,w,stats=fit_once(d,target)
    y=d[target].to_numpy(float)
    j=names.index(PREDICTOR)
    null=np.delete(X,j,axis=1)
    baseline=sm.WLS(y,null,weights=w).fit()
    full_loss=float(np.sum(w*(y-fit.fittedvalues)**2))
    base_loss=float(np.sum(w*(y-baseline.fittedvalues)**2))
    if not base_loss>0:raise RuntimeError("Null model cannot provide nonzero source residual loss")
    stats["in_sample_H6FE_partial_weighted_SSE_improvement_pct"]=100*(base_loss-full_loss)/base_loss
    stats["same_H6_FE_is_not_out_of_H6_predictive_transfer"]=True
    blocked=blocked_fl_null(d,y,X,names,w,SEED,N_PERM)
    if abs(blocked["source_residualized_plant_N_WLS_beta"]-stats["beta_NUTRI_per_full_sample_1SD"])>1e-8:
        raise RuntimeError("Independent within-H6 WLS FWL slope mismatch")
    influence=[]
    for group in sorted(d.h6.unique()):
        sub=d.loc[d.h6!=group].copy().reset_index(drop=True)
        m,_,_,_,_,s=fit_once(sub,target,center=center,cluster=False)
        influence.append({"removed_original_H6":str(group),
            "source_retained_H6":int(sub.h6.nunique()),
            "beta_per_fixed_original_N_1SD":s["beta_NUTRI_per_full_sample_1SD"]})
    signs=[x["beta_per_fixed_original_N_1SD"] for x in influence]
    return {
      "source_cohort_183H7_40H6":True,
      "response":target,
      "focal_source_index":"mean of original two-period species-equal plant EIV NUTRI",
      "model":"genuine H6 parent fixed intercept + 13 source-original terrestrial controls",
      "H6_equal_original_WLS":stats,
      "source_within_true_H6_residual_label_null":blocked,
      "leave_one_true_H6_out_40_fullsource_coefficient_range":[float(min(signs)),float(max(signs))],
      "number_LOO_positive":sum(x>0 for x in signs),
      "number_LOO_negative":sum(x<0 for x in signs),
      "full_40_H6_source_influence":influence,
      "not_an_out_of_new_H6_generalization_test":True,
      "not_wastewater_causal_mediation":True}

def main(a):
    p=a.source/"GENUINE_TERRESTRIAL_ONLY_183_H7_40_H6_PAIRED_H7_MEANS.csv"
    qa=json.loads((a.source/"GENUINE_TERRESTRIAL_183_H7_40_H6_SOURCE_RECONSTRUCTION_QA.json").read_text())
    if qa.get("schema")!=SOURCE_SCHEMA:
        raise RuntimeError("Original source restoration must pass before analysis")
    d=clean_source(p,qa)
    d=d.sort_values(["h6","h"],kind="stable").reset_index(drop=True)
    co=d.h6.value_counts()
    output={
      "schema":"GESUI_SOURCE_LOCKED_183H7_40H6_WITHIN_H6_TERRESTRIAL_NUTRI_HNS_FE_V1",
      "original_183_H7_source_sha256":PANEL_SOURCE_SHA,
      "source_origin_NUTRI":"Tichy 2023 European exact accepted vascular EIV",
      "source_origin_HNS":"fixed 186 accepted GBIF butterfly species larval-host habitat N preference",
      "model_conditions_fixed_before_source_effect_open":True,
      "original_183_H7":len(d),
      "original_true_parent_H6":int(d.h6.nunique()),
      "min_H7_per_original_H6":int(co.min()),
      "max_H7_per_original_H6":int(co.max()),
      "source_old_H7":int((d.original_source_cohort=="original_old_H7").sum()),
      "source_new_H7":int((d.original_source_cohort!="original_old_H7").sum()),
      "source_controls":CONTROLS,
      "two_periods_averaged_unweighted_per_H7_before_40H6_analysis":True,
      "W_is_static_modelled_outlet_not_measured_sewage_N_at_flora":True,
      "additional_EEA_waterstation_measures_not_source_compatible":True,
      "primary":describe_one(d,Y_PRIMARY),
      "secondary_record_weighted":describe_one(d,Y_SECONDARY),
      "existing_within_H7_temporal_non_synchrony_unchanged":True,
      "not_true_plant_butterfly_feeding_or_wastewater_heterotrophic_energy_flux":True}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"GENUINE_183_H7_40_H6_WITHIN_PARENT_EIV_HNS_RESULTS.json").write_text(
        json.dumps(output,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    brief={"n_H7":len(d),"n_H6":d.h6.nunique(),"controls":len(CONTROLS)}
    for key in ("primary","secondary_record_weighted"):
        v=output[key]
        brief[key]={
           "beta":v["H6_equal_original_WLS"]["beta_NUTRI_per_full_sample_1SD"],
           "beta95":v["H6_equal_original_WLS"]["descriptive_H6_cluster_95"],
           "beta_nominal_p":v["H6_equal_original_WLS"]["nominal_two_sided_p"],
           "H6_blocked_residual_perm_p":v["source_within_true_H6_residual_label_null"]["two_sided_H6_blocked_Freedman_Lane_permutation_p"],
           "source_null_extreme":v["source_within_true_H6_residual_label_null"]["extreme_draws_abs_beta"],
           "source_null95":v["source_within_true_H6_residual_label_null"]["null_2p5_to_97p5"],
           "within_H6_in_sample_SSE_gain_pct":v["H6_equal_original_WLS"]["in_sample_H6FE_partial_weighted_SSE_improvement_pct"],
           "leave_H6_out_beta_range":v["leave_one_true_H6_out_40_fullsource_coefficient_range"],
           "leave_H6_out_positive":v["number_LOO_positive"],
           "source_design_condition_number":v["H6_equal_original_WLS"]["source_design_condition_number"]}
    print("GESUI_183_H7_40_H6_WITHIN_PARENT_SOURCE_TEST="+json.dumps(brief,sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
