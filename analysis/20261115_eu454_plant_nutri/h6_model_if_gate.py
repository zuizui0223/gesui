#!/usr/bin/env python3
"""Gate-controlled SAME-H7 plant EIV NUTRI vs W and butterfly HNS models.

No outcome opened if original 11511 taxonomy/source eligibility gate is false.
Posthoc source association, NOT wastewater treatment, larval use, plant biomass
or a true mediation analysis.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

ORIGINAL_PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
FACETS_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
SEED=20261115
N_BOOT=2000
ALPHA=10.0
P_GEOG=["crop","urban","riverflow","naturalout","area","nEEA","nPlant",
        "NUTRI_coverage_all_Plantae_speciesKeys","log1p_eiv_scored_taxa"]
B_GEOG=["crop","urban","riverflow","naturalout","area","nEEA","nAnimal","nPlant",
        "nButterfly","primary_record_coverage_raw","primary_host_completeness",
        "NUTRI_coverage_all_Plantae_speciesKeys","log1p_eiv_scored_taxa","W"]
CATS=["country_method","period"]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def H6_weights(df):
    group=df.h6.astype(str)
    counts=group.map(group.value_counts()).to_numpy(float)
    w=len(df)/(group.nunique()*counts)
    if (w<=0).any() or not np.isfinite(w).all():
        raise RuntimeError("H6 weights invalid")
    return w

def design_full(df,numeric):
    # For WLS inference, fixed in-cohort standardization and source reference.
    columns=[np.ones(len(df))];labels=["intercept"]
    for k in numeric:
        x=df[k].to_numpy(float)
        if not np.isfinite(x).all():
            raise RuntimeError("Nonfinite original source covariate: "+k)
        mu=x.mean();sd=x.std(ddof=0)
        if sd<1e-12:raise RuntimeError("Constant fixed covariate: "+k)
        columns.append((x-mu)/sd);labels.append(k)
    for k in CATS:
        s=df[k].astype(str)
        for name in sorted(s.unique())[1:]:
            columns.append((s==name).to_numpy(float))
            labels.append(f"{k}={name}")
    X=np.column_stack(columns)
    if np.linalg.matrix_rank(X)<X.shape[1]:
        raise RuntimeError(f"Full-rank original source WLS design invalid: {len(df)} rows / {X.shape[1]} columns")
    return X,labels

def source_cluster_regression(df,target,numeric,focus):
    X,labels=design_full(df,numeric)
    y=df[target].to_numpy(float)
    if not np.isfinite(y).all():raise RuntimeError("Response is not fully supported "+target)
    w=H6_weights(df)
    groups=df.h6.astype(str).to_numpy()
    m=sm.WLS(y,X,weights=w).fit(
        cov_type="cluster",cov_kwds={"groups":groups,"use_correction":True},use_t=False)
    j=labels.index(focus)
    b=float(m.params[j]);se=float(m.bse[j])
    return {"outcome":target,"focal_predictor":focus,"adjusted_beta_per_1SD_focal":b,
            "H6_cluster_standard_error":se,"H6_cluster_descriptive95":[b-1.96*se,b+1.96*se],
            "nominal_p_descriptive":float(m.pvalues[j]),"n_source_H7period":len(df),
            "n_true_H6":int(df.h6.nunique()),"n_design_parameters":X.shape[1],
            "effective_W_SD":float(df.W.std(ddof=0)),
            "normalized_eiv_N_SD":float(df.Plant_N_incidence.std(ddof=0)),
            "model_not_causal_or_mediation":True}

def matrix_pair(train,test,numeric):
    a=[np.ones(len(train))];b=[np.ones(len(test))]
    for k in numeric:
        v=train[k].to_numpy(float);u=test[k].to_numpy(float)
        if not np.isfinite(v).all() or not np.isfinite(u).all():
            raise RuntimeError("Nonfinite CV numeric: "+k)
        mu=v.mean();sd=v.std(ddof=0)
        if sd<1e-12:sd=1.
        a.append((v-mu)/sd);b.append((u-mu)/sd)
    for c in CATS:
        names=sorted(train[c].astype(str).unique().tolist())
        for cat in names[1:]:
            a.append((train[c].astype(str)==cat).to_numpy(float))
            b.append((test[c].astype(str)==cat).to_numpy(float))
    return np.column_stack(a),np.column_stack(b)

def predict_ridge(train,test,target,numeric):
    A,B=matrix_pair(train,test,numeric)
    y=train[target].to_numpy(float)
    w=H6_weights(train)
    gram=A.T@(w[:,None]*A)
    penalty=np.eye(A.shape[1])*ALPHA
    penalty[0,0]=0
    coef=np.linalg.solve(gram+penalty,A.T@(w*y))
    return B@coef

def fold_index(group):
    return int(hashlib.sha256(str(group).encode()).hexdigest()[:8],16)%5

def cv_h6(df,target,base,extra):
    h6=df.h6.astype(str).to_numpy()
    folds=np.array([fold_index(h) for h in h6],dtype=int)
    pred0=np.full(len(df),np.nan);pred1=np.full(len(df),np.nan)
    for fold in range(5):
        test=np.flatnonzero(folds==fold);train=np.flatnonzero(folds!=fold)
        if len(test)<1 or len(train)<30:
            raise RuntimeError("Original deterministic H6 folds empty")
        if set(h6[test])&set(h6[train]):
            raise RuntimeError("H6 parent split leakage")
        d0=df.iloc[train].reset_index(drop=True)
        d1=df.iloc[test].reset_index(drop=True)
        pred0[test]=predict_ridge(d0,d1,target,base)
        pred1[test]=predict_ridge(d0,d1,target,base+[extra])
    if not (np.isfinite(pred0).all() and np.isfinite(pred1).all()):
        raise RuntimeError("Incomplete original H6 out-of-fold prediction")
    y=df[target].to_numpy(float)
    r0=(y-pred0)**2;r1=(y-pred1)**2
    by0=[];by1=[]
    for h in sorted(set(h6)):
        ok=h6==h
        by0.append(float(np.mean(r0[ok])))
        by1.append(float(np.mean(r1[ok])))
    a=np.array(by0);b=np.array(by1)
    if a.mean()<1e-12:
        raise RuntimeError("Trivial zero original group baseline MSE")
    improvement=100*(a.mean()-b.mean())/a.mean()
    rng=np.random.default_rng(SEED)
    draws=rng.integers(0,len(a),size=(N_BOOT,len(a)))
    a_bs=a[draws].mean(axis=1);b_bs=b[draws].mean(axis=1)
    boots=100*(a_bs-b_bs)/a_bs
    return {"n_true_H6":len(a),"n_original_H7period":len(df),
            "original_holdout_H6_baseline_MSE":float(a.mean()),
            "original_holdout_H6_added_MSE":float(b.mean()),
            "H6_equal_added_MSE_gain_pct":float(improvement),
            "fixed_out_of_fold_error_H6_resampling95_pct":np.quantile(boots,[.025,.975]).tolist(),
            "fixed_source_folds":5,"deterministic_H6_hash":True,
            "train_only_scaling_and_country_method_dummies":True,
            "ridge_alpha_fixed":ALPHA,
            "bootstrap_resampling_H6_groups":N_BOOT,
            "uncertainty_is_fixed_errors_not_model_refit":True,
            "not_causal":True}

def main(args):
    gate=json.loads((args.source_dir/"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if gate.get("schema")!="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V1" or gate.get("source_facets_SHA256")!=FACETS_SHA or gate.get("source_EIV",{}).get("original_source_SHA256")!=EIV_SHA:
        raise RuntimeError("Full-original source taxonomy/EIV gate incorrect")
    args.output.mkdir(parents=True,exist_ok=True)
    if not gate.get("hard_source_support_pass"):
        report={"schema":"GESUI_EU454_PLANT_HNS_SAMEBASIN_MODEL_RESULT_V1",
                "status":"STOP_ORIGINAL_SOURCE_EIV_SUPPORT",
                "reason":"Pre-outcome fixed trait coverage/independent H6/country/source API gate did not pass",
                "original_gate_eligible_H7period":gate.get("eligible_H7_period"),
                "original_gate_eligible_true_H6":gate.get("eligible_true_H6"),
                "original_gate_eligible_country":gate.get("eligible_majority_country_groups"),
                "NO_BIOLOGICAL_OUTCOME_FITTED":True}
        (args.output/"EU454_PLANT_HNS_MODELS.json").write_text(json.dumps(report,indent=2)+"\n")
        print("GESUI_H7_PLANT_HNS_MODELS="+json.dumps(report,sort_keys=True),flush=True)
        return
    source=args.source_dir/"EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv"
    if not source.exists():
        raise RuntimeError("Full source gate passed but original trait means absent")
    original=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")
    if sha(original)!=ORIGINAL_PANEL_SHA:raise RuntimeError("Original 454 scientific panel SHA mismatch")
    df=pd.read_csv(original,dtype={"h":str,"h6":str})
    scores=pd.read_csv(source,dtype={"h":str})
    elig=pd.read_csv(args.source_dir/"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv",dtype={"h":str,"h6":str})
    if (len(df),df.h.nunique(),df.h6.nunique())!=(454,232,90):
        raise RuntimeError("Analysis H7/H6 source cohort drift")
    if len(elig)!=454 or elig.duplicated(["h","period"]).any() or scores.duplicated(["h","period"]).any():
        raise RuntimeError("Source NUTRI unique geography drift")
    eligible=elig[elig.full_source_trait_eligible.astype(str)=="True"].copy()
    if len(eligible)!=gate["eligible_H7_period"]:
        raise RuntimeError("Source quality threshold rows not exactly matched")
    merged=df.merge(eligible[["h","period","NUTRI_coverage_all_Plantae_speciesKeys",
                            "exact_NUTRI_scored_vascular_accepted_species"]],
                   on=["h","period"],validate="one_to_one",how="inner")
    merged=merged.merge(scores[["h","period","Plant_N_incidence","Plant_N_record","n_scored_species"]],
                         on=["h","period"],validate="one_to_one",how="inner")
    if len(merged)!=len(eligible) or merged.h6.nunique()!=gate["eligible_true_H6"]:
        raise RuntimeError("Species trait table changed original eligible grain/H6")
    if merged.country_holdout.nunique()!=gate["eligible_majority_country_groups"]:
        raise RuntimeError("Trait-eligibility prior country group count changed")
    if merged.n_scored_species.min()<10 or merged.NUTRI_coverage_all_Plantae_speciesKeys.min()<.4:
        raise RuntimeError("Pre-outcome trait support gate bypassed")
    merged["log1p_eiv_scored_taxa"]=np.log1p(merged.exact_NUTRI_scored_vascular_accepted_species.to_numpy(float))
    merged["country_method"]=(merged.country.fillna("__UNKNOWN__").astype(str)+"::"+
                              merged.method.fillna("__UNKNOWN__").astype(str))
    if merged[["Plant_N_incidence","Plant_N_record","primary_record_mean",
               "primary_incidence_mean","W"]].isna().any().any():
        raise RuntimeError("Missing NUTRI/HNS response despite source quality gate")
    # Two observational W->Plant N models; both use exactly frozen source covariates.
    A={}
    for plant in ["Plant_N_incidence","Plant_N_record"]:
        A[plant]={"association":source_cluster_regression(merged,plant,P_GEOG+["W"],"W"),
                  "H6_out_W_increment":cv_h6(merged,plant,P_GEOG,"W")}
    # Additional *same-basin* noncausal compositional model: predictive addition
    # of observed plant N index to an otherwise W/geography/effort-adjusted
    # butterfly HNS. This is not an identified mediator or dietary flow.
    B={}
    for butterfly in ["primary_record_mean","primary_incidence_mean"]:
        B[butterfly]={"association":source_cluster_regression(
            merged,butterfly,B_GEOG+["Plant_N_incidence"],"Plant_N_incidence"),
            "H6_out_PlantN_increment":cv_h6(merged,butterfly,B_GEOG,"Plant_N_incidence")}
    out={"schema":"GESUI_EU454_PLANT_HNS_SAMEBASIN_MODEL_RESULT_V1",
         "status":"SOURCE_QUALITY_GATE_PASS_POSTHOC_OBSERVATIONAL_ASSOCIATIONS",
         "original_panel_sha256":ORIGINAL_PANEL_SHA,
         "true_original_source_plant_facets_SHA256":FACETS_SHA,
         "original_EIV_source_SHA256":EIV_SHA,
         "source_gate_hard_pass":True,
         "source_gate_H7period":gate["eligible_H7_period"],
         "source_gate_H6":gate["eligible_true_H6"],
         "source_gate_country_holdout_groups":gate["eligible_majority_country_groups"],
         "n_analysis_H7period":len(merged),
         "n_true_H6":int(merged.h6.nunique()),
         "n_country_holdout_groups":int(merged.country_holdout.nunique()),
         "source_model_adjustment_P_GEOG":P_GEOG,
         "source_model_adjustment_B_GEOG":B_GEOG,
         "plant_NUTRI_eiv_range_source_only":[float(merged.Plant_N_incidence.min()),float(merged.Plant_N_incidence.max())],
         "W_to_plant_N":A,"plant_N_added_to_butterfly_HNS":B,
         "original_scored_records_are_GBIF_NOT_plant_cover":True,
         "environment_W_is_modelled_static_geography_not_site_effluent_N":True,
         "PlantN_to_HNS_is_not_causal_mediation":True,
         "H6_resampling_intervals_are_descriptive":True,
         "source_trait_selection_is_nonrandom_missingness_risk":True}
    (args.output/"EU454_PLANT_HNS_MODELS.json").write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    desc={k:out[k] for k in ["status","n_analysis_H7period","n_true_H6","n_country_holdout_groups"]}
    desc["W_to_plant_N"]={key:{"beta":v["association"]["adjusted_beta_per_1SD_focal"],
              "H6_CI":v["association"]["H6_cluster_descriptive95"],
              "H6_MSE_gain":v["H6_out_W_increment"]["H6_equal_added_MSE_gain_pct"],
              "gain95":v["H6_out_W_increment"]["fixed_out_of_fold_error_H6_resampling95_pct"]}
              for key,v in A.items()}
    desc["plant_N_added_to_butterfly_HNS"]={key:{"beta":v["association"]["adjusted_beta_per_1SD_focal"],
              "H6_CI":v["association"]["H6_cluster_descriptive95"],
              "H6_MSE_gain":v["H6_out_PlantN_increment"]["H6_equal_added_MSE_gain_pct"],
              "gain95":v["H6_out_PlantN_increment"]["fixed_out_of_fold_error_H6_resampling95_pct"]}
              for key,v in B.items()}
    print("GESUI_H7_PLANT_HNS_MODELS="+json.dumps(desc,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--source-dir",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    main(parser.parse_args())

