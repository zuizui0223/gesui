#!/usr/bin/env python3
"""Adversarial H6 fixed-effect EIV NUTRI vs HNS eight covariate block grid.

All eight block combinations on the SAME genuine original 183 H7, 40 H6.
Explore suppression, high partial VIF and source-driven conditional estimates.
This is posthoc after positive full adjusted β and not a prospective test.
"""
from __future__ import annotations
import argparse,itertools,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm

SOURCE_SHA="c0840b1c5c6e85342fea41d57dc40fa5ec295fabf6bcd6908a4072557db9a736"
FOCAL="plant_n_incidence"
OUTCOME="hns_incidence"
SEED=20261115
BLOCKS={
  "environment":["W","crop","urban","naturalout","area"],
  "recording":["nPlant","nButterfly","nAnimal"],
  "trait_support":["plant_coverage","hns_coverage","hns_host_completeness",
                   "plant_scored_species","hns_scored_species"]}
BLOCK_ORDER=["environment","recording","trait_support"]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def one(d,name):
    a=list(d.rglob(name))
    if len(a)!=1:raise RuntimeError("Original 183-H7 source file not unique "+name)
    return a[0]
def weight(d):
    size=d.h6.map(d.h6.value_counts()).to_numpy(float)
    return len(d)/(d.h6.nunique()*size)

def model(d,controls,scale):
    y=d[OUTCOME].to_numpy(float)
    feat=[np.ones(len(d))]
    names=["intercept"]
    for key in controls+[FOCAL]:
        mu,sd=scale[key]
        feat.append((d[key].to_numpy(float)-mu)/sd)
        names.append(key)
    g=d.h6.astype(str)
    for h in sorted(g.unique().tolist())[1:]:
        feat.append((g==h).to_numpy(float));names.append("H6="+h)
    X=np.column_stack(feat)
    rk=int(np.linalg.matrix_rank(X))
    if rk!=X.shape[1]:raise RuntimeError("Full H6 model block rank deficiency "+str((controls,rk,X.shape[1])))
    w=weight(d)
    fitted=sm.WLS(y,X,weights=w).fit(
        cov_type="cluster",cov_kwds={"groups":g.to_numpy(),"use_correction":True},use_t=False)
    j=names.index(FOCAL)
    beta=float(fitted.params[j]);se=float(fitted.bse[j])
    X0=np.delete(X,j,axis=1);x=X[:,j]
    root=np.sqrt(w)
    mx=np.linalg.lstsq(X0*root[:,None],x*root,rcond=None)[0]
    rx=x-X0@mx
    den=float(np.dot(w*rx,rx))
    mean=np.average(x,weights=w)
    tss=float(np.dot(w*(x-mean),x-mean))
    if den<=1e-12:raise RuntimeError("No identifiable source within-H6 N after covariates")
    fwl=float(np.dot(w*rx,y)/den)
    if abs(fwl-beta)>1e-8:raise RuntimeError("Source WLS FWL independent coefficient mismatch")
    return {
        "source_environment_controls":controls,
        "beta_per_original_full183_1SD_PlantN":beta,
        "H6_cluster_normal_descriptive95":[beta-1.96*se,beta+1.96*se],
        "nominal_source_p":float(fitted.pvalues[j]),
        "partial_PlantN_weighted_VIF":float(tss/den),
        "weighted_design_condition_number":float(np.linalg.cond(X*root[:,None])),
        "number_design_columns":int(X.shape[1]),
        "number_true_H6":int(g.nunique()),
        "not_causal":True}

def withinraw(d):
    # Within each H6 mean-centering and then H6-equal weighting, without
    # geography or effort controls, purely descriptive.
    w=weight(d)
    x=d[FOCAL].to_numpy(float)
    y=d[OUTCOME].to_numpy(float)
    groups=d.h6.astype(str).to_numpy()
    xm=np.empty_like(x);ym=np.empty_like(y)
    for h in sorted(set(groups)):
        ids=groups==h
        xm[ids]=x[ids]-x[ids].mean()
        ym[ids]=y[ids]-y[ids].mean()
    ux=np.dot(w*xm,xm);uy=np.dot(w*ym,ym)
    xy=np.dot(w*xm,ym)
    return {"H6_equal_within_parent_demeaned_pearson_r":float(xy/np.sqrt(ux*uy)),
            "H6_equal_unadjusted_beta_per_raw_unit_PlantN":float(xy/ux),
            "frozen_n_H7":len(d),"frozen_n_H6":len(set(groups)),
            "not_source_nutrient_dose_or_mediation":True}

def main(a):
    p=one(a.source,"GENUINE_TERRESTRIAL_ONLY_183_H7_40_H6_PAIRED_H7_MEANS.csv")
    qa=json.loads(one(a.source,"GENUINE_TERRESTRIAL_183_H7_40_H6_SOURCE_RECONSTRUCTION_QA.json").read_text())
    if sha(p)!=SOURCE_SHA or qa.get("data_source_panel_SHA256")!=SOURCE_SHA:
        raise RuntimeError("Immutable genuine source-population SHA mismatch")
    d=pd.read_csv(p,dtype={"h":str,"h6":str})
    d=d.sort_values(["h6","h"]).reset_index(drop=True)
    if len(d)!=183 or d.h.nunique()!=183 or d.h6.nunique()!=40 or d.h6.value_counts().min()<2:
        raise RuntimeError("True 183 H7 /40 H6 paired regional source frame mismatch")
    required=[z for keys in BLOCKS.values() for z in keys]+[FOCAL,OUTCOME]
    if not np.isfinite(d[required].to_numpy(float)).all():
        raise RuntimeError("Original H6 plant/HNS source or environmental data not finite")
    scales={}
    for key in required:
        v=d[key].to_numpy(float)
        sd=v.std(ddof=0)
        if sd<=1e-12:raise RuntimeError("Original source block numeric variable constant "+key)
        scales[key]=(float(v.mean()),float(sd))
    models={}
    for n in range(4):
        for subset in itertools.combinations(BLOCK_ORDER,n):
            key="+".join(subset) if subset else "H6_FE_ONLY"
            cov=[c for group in subset for c in BLOCKS[group]]
            models[key]=model(d,cov,scales)
    if len(models)!=8:
        raise RuntimeError("All eight posthoc block combinations not evaluated")
    full=models["environment+recording+trait_support"]
    if abs(full["beta_per_original_full183_1SD_PlantN"]-.11807615189599471)>1e-8:
        raise RuntimeError("Original frozen fully adjusted main coefficient failed source replay")
    out={"schema":"GESUI_183_H7_40_H6_POSTHOC_ALL8_COVARIATE_BLOCKS_V1",
         "source_183_H7_panel_SHA256":SOURCE_SHA,
         "source_original_H7":183,"source_original_true_H6":40,
         "focal_original_NUTRI_source":"species-equal EIV habitat nutrient affinity",
         "response_original_HNS_source":"species-equal butterfly potential larval host N affinity",
         "source_block_definitions":BLOCKS,
         "all_eight_subsets_reported_no_selection":True,
         "posthoc_after_original_positive_and_highVIF_not_preregistered":True,
         "full_regression_independently_reproduced":True,
         "within_H6_unadjusted_correlation":withinraw(d),
         "every_source_model":models,
         "not_wastewater_mediation_or_true_Ecological_survey_abundance":True}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"ALL_8_H6_FIXED_EIV_HNS_CONTROL_BLOCK_RESULTS.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    summary={k:{"beta":v["beta_per_original_full183_1SD_PlantN"],
                "H6_cluster95":v["H6_cluster_normal_descriptive95"],
                "source_VIF_NUTRI":v["partial_PlantN_weighted_VIF"],
                "condition_number":v["weighted_design_condition_number"]}
             for k,v in models.items()}
    print("GESUI_ALL8_H6_CONTROL_BLOCK_ADVERSARIAL="+json.dumps({
        "raw_withinH6":out["within_H6_unadjusted_correlation"],
        "all8":summary},ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
