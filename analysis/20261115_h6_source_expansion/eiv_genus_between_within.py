#!/usr/bin/env python3
"""Strict algebraic within/between-GENUS NUTRI composition decomposition.

Post-result, noncausal: source-authentic 183 H7 x exactly two periods and
40 genuine HydroBASINS H6. All 4011 accepted binomial plant scores retained.
A(H7,species) constructed and checked in plant_eiv_global_genus_label_null.
Original parent FE and thirteen source covariates NOT modified.
"""
from __future__ import annotations
import argparse
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

from plant_eiv_global_genus_label_null import source_frame
from test_183_40_within_H6 import (
    features,weights,CONTROLS,PREDICTOR,Y_PRIMARY,Y_SECONDARY
)

SOURCE_EIV_KEYS=4011
EXPECTED_BETA=0.11807615189599471

def one(folder,name):
    paths=list(folder.rglob(name))
    if len(paths)!=1:raise RuntimeError("Original source exact file is missing/duplicated "+name)
    return paths[0]

def group_eiv_genera(cross,values):
    t=pd.read_csv(one(cross,"EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv"),
                  dtype={"speciesKey":str})
    if len(t)!=11511 or t.speciesKey.duplicated().any():
        raise RuntimeError("Original 11511 speciesKey original source lost")
    t=t.loc[t.EIV_NUTRI_exact.notna()].sort_values("speciesKey",kind="stable")
    if len(t)!=SOURCE_EIV_KEYS or np.max(np.abs(t.EIV_NUTRI_exact.to_numpy(float)-values))>1e-12:
        raise RuntimeError("Original EIV score order/values do not exactly reproduce plant incidence matrix")
    genus=np.array([str(n).split()[0].casefold() for n in t.gbif_canonical_name],dtype=str)
    if any(not re.fullmatch(r"[a-z][a-z-]+",g) for g in genus):
        raise RuntimeError("Not all 4011 original accepted GBIF source taxa have verified genus")
    ref=np.zeros(len(genus))
    allg={}
    for g in sorted(set(genus)):
        idx=np.flatnonzero(genus==g)
        ref[idx]=values[idx].mean()
        allg[g]={"n_scored_species":len(idx),"NUTRI_numeric_mean":float(ref[idx[0]]),
                 "NUTRI_numeric_SD":float(np.std(values[idx]))}
    return genus,ref,allg

def residuals(d,v,Q,root):
    weighted=root*np.asarray(v,float)
    return weighted-Q@(Q.T@weighted)

def numeric_summarize(u,v):
    if u@u<1e-15 or v@v<1e-15:return None
    return float((u@v)/np.sqrt(float((u@u)*(v@v))))

def cluster_model(d,y,baseX,bres,wres,wb):
    # Include both genuine genus components as distinct predictors while
    # keeping exact source full-13 + genuine 40 H6 fixed design.
    scale_b=float(np.std(bres,ddof=0));scale_w=float(np.std(wres,ddof=0))
    if scale_b<1e-8 or scale_w<1e-8:raise RuntimeError("Genus component has no score variance")
    X=np.column_stack([baseX,(bres-np.mean(bres))/scale_b,
                         (wres-np.mean(wres))/scale_w])
    if np.linalg.matrix_rank(X)!=X.shape[1]:
        return {"status":"STOP_GENUS_COMPONENTS_COLLINEAR",
                "rank":int(np.linalg.matrix_rank(X)),"columns":X.shape[1]}
    res=sm.WLS(y,X,weights=wb).fit(
        cov_type="cluster",cov_kwds={"groups":d.h6.to_numpy(),"use_correction":True},
        use_t=False)
    ci={}
    for j,name in ((-2,"between_genus"),(-1,"within_genus")):
        beta=float(res.params[j]);se=float(res.bse[j])
        ci[name]={"beta_per_component_1SD":beta,
                  "H6_cluster_description95":[beta-1.96*se,beta+1.96*se],
                  "nominal_p":float(res.pvalues[j])}
    return {"status":"FIT",
            "source_components":ci,
            "weighted_source_design_condition":float(np.linalg.cond(X*np.sqrt(wb)[:,None])),
            "same_original_H7_and_H6_as_total":True,
            "component_beta_scaling_is_DIFFERENT_not_additive":True}

def model_partial_stats(d,sourceA,sourceScore,genusMean,originalQa):
    srcTotal=sourceA@sourceScore
    srcB=sourceA@genusMean
    srcW=sourceA@(sourceScore-genusMean)
    if max(np.abs(srcTotal-(srcB+srcW)))>1e-12:
        raise RuntimeError("Original NUTRI genus decomposition identity fails")
    if max(np.abs(srcTotal-d[PREDICTOR].to_numpy(float)))>1e-10:
        raise RuntimeError("Original 183 H7 plant N source differs from decomposition")
    Xfull,names,_=features(d)
    j=names.index(PREDICTOR)
    X0=np.delete(Xfull,j,axis=1)
    w=weights(d)
    root=np.sqrt(w)
    Q,_=np.linalg.qr(root[:,None]*X0,mode="reduced")
    ut=residuals(d,srcTotal,Q,root)
    ub=residuals(d,srcB,Q,root)
    uw=residuals(d,srcW,Q,root)
    if np.max(np.abs(ut-(ub+uw)))>1e-10:
        raise RuntimeError("Within genus and between genus weighted residual FWL identity fails")
    den=float(ut@ut)
    if den<1e-10:raise RuntimeError("No original source within-H6 NUTRI variation")
    orig_sd=float(np.std(srcTotal,ddof=0))
    youtputs={}
    for y_name in (Y_PRIMARY,Y_SECONDARY):
        yr=residuals(d,d[y_name].to_numpy(float),Q,root)
        denom_y=float(yr@yr)
        if denom_y<1e-10:raise RuntimeError("No original butterfly source variance in H6")
        sourceBeta=orig_sd*float(ut@yr)/den
        betweenContribution=orig_sd*float(ub@yr)/den
        withinContribution=orig_sd*float(uw@yr)/den
        if abs(sourceBeta-betweenContribution-withinContribution)>1e-9:
            raise RuntimeError("Original beta fixed-denominator algebraic genus decompo did not reproduce")
        if y_name==Y_PRIMARY and abs(sourceBeta-EXPECTED_BETA)>1e-8:
            raise RuntimeError("Original reported within-H6 primary EIV source beta mismatch")
        corrs={
            "total":numeric_summarize(yr,ut),
            "between_genus":numeric_summarize(yr,ub),
            "within_genus":numeric_summarize(yr,uw)
        }
        # In-sample H6 FE + same original controls partial gain.
        SSE0=denom_y
        Bcoef=float((ub@yr)/(ub@ub))
        SSEB=float((yr-ub*Bcoef)@(yr-ub*Bcoef))
        Z=np.column_stack([ub,uw])
        if np.linalg.matrix_rank(Z)<2:
            raise RuntimeError("Genus components not separable after same source controls")
        coefs=np.linalg.lstsq(Z,yr,rcond=None)[0]
        resid_full=yr-Z@coefs
        SSEBW=float(resid_full@resid_full)
        if SSEB<SSEBW-1e-10 or SSE0<SSEB-1e-10:
            raise RuntimeError("Nested genotype source sums of squares inconsistent")
        after_between_y=yr-ub*Bcoef
        after_between_w=uw-ub*float((ub@uw)/(ub@ub))
        adjcorr=numeric_summarize(after_between_y,after_between_w)
        youtputs[y_name]={
           "original_full_source_plantN_beta_scaled_1SD":sourceBeta,
           "algebraic_between_genus_covariance_component_original_denominator":betweenContribution,
           "algebraic_within_genus_covariance_component_original_denominator":withinContribution,
           "beta_components_sum_equals_full":True,
           "full_source_original_FWL_partial_r":corrs["total"],
           "component_partial_r_after_original_H6_FE_and_13controls":corrs,
           "between_component_H6_in_sample_partial_SSE_gain_pct":100*(SSE0-SSEB)/SSE0,
           "within_added_after_between_H6_in_sample_SSE_gain_pct":100*(SSEB-SSEBW)/SSEB,
           "both_genus_components_H6_in_sample_partial_SSE_gain_pct":100*(SSE0-SSEBW)/SSE0,
           "within_genus_given_between_component_partial_r":adjcorr,
           "both_component_model_H6_cluster":cluster_model(d,d[y_name].to_numpy(float),X0,srcB,srcW,w),
           "no_experiment_or_out_of_region_prediction_fitted":True}
    return {"source_genus_between_vs_within_CWM_correlation":float(np.corrcoef(srcB,srcW)[0,1]),
            "source_component_ranges":{"between":[float(srcB.min()),float(srcB.max())],
               "within":[float(srcW.min()),float(srcW.max())]},
            "source_full_original_NUTRI_CWM_variance":float(np.var(srcTotal)),
            "source_between_genus_NUTRI_CWM_variance":float(np.var(srcB)),
            "source_within_genus_NUTRI_CWM_variance":float(np.var(srcW)),
            "model_targets":youtputs,
            "n_fixed_original_H7":len(d),
            "n_fixed_original_H6":int(d.h6.nunique())}

def main(a):
    frame,A,values,_,qa=source_frame(a)
    genus,genusmean,groups=group_eiv_genera(a.eiv,values)
    out=model_partial_stats(frame,A,values,genusmean,qa)
    out.update({
        "schema":"GESUI_183H7_40H6_PLANT_EIV_NUTRI_BETWEEN_WITHIN_GENUS_DIAGNOSTIC_V1",
        "original_source_receipts":qa,
        "n_scored_plant_species":len(values),
        "n_original_accepted_genus":len(groups),
        "n_scored_singleton_genera":sum(v["n_scored_species"]==1 for v in groups.values()),
        "n_scored_multispecies_genera":sum(v["n_scored_species"]>=2 for v in groups.values()),
        "n_original_scored_species_in_genus_size_ge2":int(sum(v["n_scored_species"] for v in groups.values() if v["n_scored_species"]>=2)),
        "genus_nutri_numeric_SD_nonzero_groups":sum(v["NUTRI_numeric_SD"]>1e-10 for v in groups.values()),
        "source_genus_reference_equal_scored_species_average":True,
        "original_species_by_H7_incidence_and_two_periods_preserved":True,
        "posthoc_full_W_FE_controls_and_all_site_records_fixed":True,
        "original_9999_group_residual_and_4999_label_perm_results_not_rewritten":True,
        "no_original_local_wastewater_N_or_larval_consumption_observed":True,
        "not_a_genus_phylogeny_or_causal_mediation":True})
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"183H7_40H6_PLANT_EIV_BETWEEN_WITHIN_GENUS_RESULTS.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    focus={}
    for key,val in out["model_targets"].items():
        focus[key]={k:val[k] for k in ("original_full_source_plantN_beta_scaled_1SD",
          "algebraic_between_genus_covariance_component_original_denominator",
          "algebraic_within_genus_covariance_component_original_denominator",
          "component_partial_r_after_original_H6_FE_and_13controls",
          "between_component_H6_in_sample_partial_SSE_gain_pct",
          "within_added_after_between_H6_in_sample_SSE_gain_pct",
          "both_genus_components_H6_in_sample_partial_SSE_gain_pct",
          "within_genus_given_between_component_partial_r",
          "both_component_model_H6_cluster")}
    summary={"n_H7":len(frame),"n_H6":int(frame.h6.nunique()),
      "original_n_scored_species":len(values),"genus_n":len(groups),
      "source_genus_corr":out["source_genus_between_vs_within_CWM_correlation"],
      "targets":focus}
    print("GESUI_EIV_GENUS_COMPONENTS_183_H7_40_H6="+json.dumps(
        summary,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source183",required=True,type=Path)
    p.add_argument("--oldfacets",required=True,type=Path)
    p.add_argument("--newfacets",required=True,type=Path)
    p.add_argument("--eiv",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
