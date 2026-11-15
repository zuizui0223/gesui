#!/usr/bin/env python3
"""Post-result adversarial within-H6 source association diagnostics.

- modelled static wastewater-geography W -> species-equal Plant NUTRI
- W -> butterfly equal-incidence HNS after Plant N
- Plant N -> HNS after old/new GBIF source-cohort indicator
- old-only 163 H7 in 30 genuine multiH7 H6 descriptive sensitivity
- H6 source W variation, N VIF and ten paired old/new H7 descriptions

No sewage-N delivered dose, insect counts, local feeding, causal mediation,
or independent geographic out-of-H6 predictive transfer is estimated.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import defaultdict
from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm

SRC_SHA="c0840b1c5c6e85342fea41d57dc40fa5ec295fabf6bcd6908a4072557db9a736"
CONTROLS=["W","crop","urban","naturalout","area","nPlant","nButterfly","nAnimal",
          "plant_coverage","hns_coverage","hns_host_completeness",
          "plant_scored_species","hns_scored_species"]
PLANT="plant_n_incidence"
BUTTERFLY="hns_incidence"
SEED=20261115
N_PERM=9999

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def weight(d):
    h=d.h6.astype(str)
    sz=h.map(h.value_counts()).to_numpy(float)
    return len(d)/(h.nunique()*sz)
def construct(d,control,focal,scales,source=False):
    columns=[np.ones(len(d))]
    names=["intercept"]
    for k in list(control)+[focal]:
        mu,sd=scales[k]
        columns.append((d[k].to_numpy(float)-mu)/sd)
        names.append(k)
    if source:
        columns.append((d.original_source_cohort!="original_old_H7").astype(float).to_numpy())
        names.append("source_new_H7_fixed_indicator")
    g=d.h6.astype(str)
    for h in sorted(g.unique().tolist())[1:]:
        columns.append((g==h).to_numpy(float))
        names.append("true_H6_FE="+h)
    X=np.column_stack(columns)
    if np.linalg.matrix_rank(X)!=X.shape[1]:
        raise RuntimeError("Source-locked H6 fixed-effect design not full rank "+str((np.linalg.matrix_rank(X),X.shape[1])))
    return X,names

def model(d,target,focal,controls,scales,source=False,blocked=True):
    X,names=construct(d,controls,focal,scales,source)
    j=names.index(focal)
    y=d[target].to_numpy(float)
    w=weight(d)
    res=sm.WLS(y,X,weights=w).fit(
        cov_type="cluster",cov_kwds={"groups":d.h6.astype(str).to_numpy(),"use_correction":True},
        use_t=False)
    b,se=float(res.params[j]),float(res.bse[j])
    output={"outcome":target,"focal":focal,
       "source_cohort_indicator_controlled":source,
       "beta_per_original_full_183_H7_1SD":b,
       "source_H6_cluster_descriptive95":[b-1.96*se,b+1.96*se],
       "nominal_p_descriptive":float(res.pvalues[j]),
       "n_H7":len(d),"n_true_H6":int(d.h6.nunique()),
       "n_design_columns":X.shape[1],
       "source_weighted_condition_number":float(np.linalg.cond(X*np.sqrt(w[:,None]))),
       "not_causal":True}
    X0=np.delete(X,j,axis=1);x=X[:,j]
    root=np.sqrt(w)
    betax=np.linalg.lstsq(X0*root[:,None],x*root,rcond=None)[0]
    rx=x-X0@betax
    den=float(np.dot(w*rx,rx))
    tss=float(np.dot(w*(x-np.average(x,weights=w)),x-np.average(x,weights=w)))
    output["source_residualized_focal_VIF"]=float(tss/den) if den>1e-12 else None
    if den<=1e-12:raise RuntimeError("Focal variable no residual within-H6 variation")
    beta_independent=float(np.dot(w*rx,y)/den)
    if abs(beta_independent-b)>1e-8:raise RuntimeError("Independent original FWL beta not reproduced")
    if blocked:
        b0=np.linalg.lstsq(X0*root[:,None],y*root,rcond=None)[0]
        residual=y-X0@b0
        h=d.h6.astype(str).to_numpy()
        members=[np.flatnonzero(h==g) for g in sorted(set(h))]
        rng=np.random.default_rng(SEED)
        shuffled=np.empty_like(residual)
        null=np.empty(N_PERM)
        for t in range(N_PERM):
            for group in members:
                shuffled[group]=rng.permutation(residual[group])
            null[t]=np.dot(w*rx,shuffled)/den
        ext=int(np.sum(np.abs(null)>=abs(b)-1e-12))
        output["source_H6_blocked_residual_permutation"]={
          "n_draws":N_PERM,"extreme_at_least_abs_observed":ext,
          "two_sided_p":(ext+1)/(N_PERM+1),
          "null_q025_q975":np.quantile(null,[.025,.975]).tolist(),
          "null_not_an_assigned_experiment":True}
    return output

def main(a):
    p=a.source/"GENUINE_TERRESTRIAL_ONLY_183_H7_40_H6_PAIRED_H7_MEANS.csv"
    qa=json.loads((a.source/"GENUINE_TERRESTRIAL_183_H7_40_H6_SOURCE_RECONSTRUCTION_QA.json").read_text())
    if (sha(p)!=SRC_SHA or qa.get("data_source_panel_SHA256")!=SRC_SHA
       or qa.get("n_old_H7")!=173 or qa.get("n_new_H7")!=10):
        raise RuntimeError("Original 183-H7 genuine original EIV/HNS source population not intact")
    d=pd.read_csv(p,dtype={"h":str,"h6":str})
    if len(d)!=183 or d.h.nunique()!=183 or d.h6.nunique()!=40:
        raise RuntimeError("Original strict 183 H7 40 H6 spatial frame changed")
    all_keys=CONTROLS+[PLANT,BUTTERFLY]
    if not np.isfinite(d[all_keys].to_numpy(float)).all():
        raise RuntimeError("Original all-source response fields contain nonfinite values")
    scales={}
    for k in all_keys:
        v=d[k].to_numpy(float)
        mu,sd=float(np.mean(v)),float(np.std(v))
        if sd<=1e-12:raise RuntimeError("Original common source numeric constant "+k)
        scales[k]=(mu,sd)
    group=d.groupby("h6")
    n_wvar=int(sum(bool((g.W.max()-g.W.min())>1e-10) for _,g in group))
    n_wpositive=int(sum(bool(g.W.max()>0) for _,g in group))
    n_wallzero=int(sum(bool(g.W.max()==0) for _,g in group))
    # Primary prior coefficient reproduction in the same source prior frozen
    # model, before adding posthoc tests.
    prior=model(d,BUTTERFLY,PLANT,CONTROLS,scales,source=False,blocked=False)
    if abs(prior["beta_per_original_full_183_H7_1SD"]-.11807615189599471)>1e-8:
        raise RuntimeError("Original separate H6-positive benchmark not reproduced exactly")
    results={
      "schema":"GESUI_183_H7_40_H6_POSTHOC_W_SOURCE_AND_COHORT_BIAS_ADVERSARIAL_V1",
      "analysis_stage":"posthoc after original within-H6 positive",
      "source_original_183H7_panel_SHA256":SRC_SHA,
      "original_H7":183,"original_H6":40,
      "within_H6_with_some_W_gradient_H6":n_wvar,
      "within_H6_at_least_one_nonzero_W_H6":n_wpositive,
      "within_H6_all_W_zero_H6":n_wallzero,
      "original_prior_positive_beta_source_recovered":prior["beta_per_original_full_183_H7_1SD"],
      "original_prior_plantN_residualized_weighted_VIF":prior["source_residualized_focal_VIF"],
      "original_prior_design_condition_number":prior["source_weighted_condition_number"],
      "modeled_W_to_plant_N_within_true_H6":model(
          d,PLANT,"W",[z for z in CONTROLS if z!="W"],scales),
      "modeled_W_to_butterfly_HNS_after_plantN_within_H6":model(
          d,BUTTERFLY,"W",[z for z in CONTROLS if z!="W"]+[PLANT],scales),
      "plant_N_to_butterfly_HNS_with_old_new_source_indicator":model(
          d,BUTTERFLY,PLANT,CONTROLS,scales,source=True),
      "original_source_empty_W_group_cannot_identify_W_effect":True,
      "no_new_researcher_choice_of_H7_or_H6":True,
      "not_delivered_sewage_nutrient_exposure":True,
      "not_plant_mediated_causal_effect":True}
    old=d[d.original_source_cohort=="original_old_H7"].copy()
    old_group=old.h6.value_counts()
    keep=set(old_group.index[old_group>=2])
    old=old[old.h6.isin(keep)].reset_index(drop=True)
    if len(old)!=163 or old.h6.nunique()!=30:
        raise RuntimeError("Former independent old-only multi H7 H6 group count altered")
    results["old_only_below_primary_40H6_support"]={
       "interpretation":"posthoc descriptive only, frozen forty H6 main gate not passed",
       "test":model(old,BUTTERFLY,PLANT,CONTROLS,scales,blocked=False)}
    ten=[]
    for parent,g in d.groupby("h6"):
        oldg=g[g.original_source_cohort=="original_old_H7"]
        newg=g[g.original_source_cohort!="original_old_H7"]
        if len(newg)>0:
            if len(newg)!=1 or len(oldg)!=1 or len(g)!=2:
                raise RuntimeError("New source true H6 was not original singleton plus one selected H7")
            delta_n=float(newg.plant_n_incidence.iloc[0]-oldg.plant_n_incidence.iloc[0])
            delta_b=float(newg.hns_incidence.iloc[0]-oldg.hns_incidence.iloc[0])
            ten.append({"native_H6":parent,"H7_old":str(oldg.h.iloc[0]),"H7_new":str(newg.h.iloc[0]),
              "new_minus_old_Plant_N":delta_n,
              "new_minus_old_Butterfly_HNS":delta_b,
              "concordant_direction":bool(delta_n*delta_b>0)})
    if len(ten)!=10:raise RuntimeError("Original ten true H6 source new-old pairs missing")
    results["source_new_10_H6_original_paired_direction_only"]={
       "n_true_H6":len(ten),"concordant_direction_n":int(sum(q["concordant_direction"] for q in ten)),
       "new_vs_old_paired_mean_DeltaPlantN":float(np.mean([q["new_minus_old_Plant_N"] for q in ten])),
       "new_vs_old_paired_mean_DeltaButterflyHNS":float(np.mean([q["new_minus_old_Butterfly_HNS"] for q in ten])),
       "source_ten_pair_not_a_full_adjusted_or_experimental_test":True,
       "original_each_pair":ten}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"POSTHOC_WITHIN_H6_W_MODEL_AND_SOURCE_COHORT_QA.json").write_text(
        json.dumps(results,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    short={
      "original_183_H7_40_H6":True,
      "W_varying_H6":n_wvar,"W_some_positive_H6":n_wpositive,"W_all_zero_H6":n_wallzero,
      "original_NUTRI_residualized_VIF":results["original_prior_plantN_residualized_weighted_VIF"],
      "ten_new_H6_paired_concordant_count":results["source_new_10_H6_original_paired_direction_only"]["concordant_direction_n"]}
    for name in ("modeled_W_to_plant_N_within_true_H6",
                 "modeled_W_to_butterfly_HNS_after_plantN_within_H6",
                 "plant_N_to_butterfly_HNS_with_old_new_source_indicator"):
        z=results[name]
        short[name]={"beta":z["beta_per_original_full_183_H7_1SD"],
          "H6_cluster95":z["source_H6_cluster_descriptive95"],
          "H6_nominal_p":z["nominal_p_descriptive"],
          "H6_blocked_perm_p":z["source_H6_blocked_residual_permutation"]["two_sided_p"],
          "focal_weighted_VIF":z["source_residualized_focal_VIF"],
          "condition_number":z["source_weighted_condition_number"]}
    origonly=results["old_only_below_primary_40H6_support"]["test"]
    short["original_only_below40_H6"]={"n_H7":old.shape[0],"n_H6":int(old.h6.nunique()),
        "beta":origonly["beta_per_original_full_183_H7_1SD"],
        "cluster95":origonly["source_H6_cluster_descriptive95"]}
    print("GESUI_WITHIN_H6_W_AND_SOURCE_BIAS_POSTHOC="+json.dumps(short,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
