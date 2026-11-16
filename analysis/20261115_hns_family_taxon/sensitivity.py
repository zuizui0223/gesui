#!/usr/bin/env python3
"""Post-result source-locked sensitivities for family HNS association.

Tests taxon leverage, record effort, eligibility cohorts and trait bounds.
Not preregistered, not causal, no new data or field observations.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
import statsmodels.api as sm
from analyse import (load_reconstruction, build_family, cluster_fit, country_prediction,
                     FAMS, COMMON, MIN_SCORED_RECORDS, MIN_SCORED_SPECIES, MIN_SCORE_COVERAGE,
                     design, cv_add_value)


def eligible(df, family, min_n, min_species, min_coverage):
    return ((df["count_records_"+family].to_numpy()>=min_n)&
            (df["count_species_"+family].to_numpy()>=min_species)&
            (df["source_score_coverage_"+family].to_numpy()>=min_coverage))


def fit_extra(df, y, extras):
    """Fitted solely as robustness; not controlled-direct causal regression."""
    X, names=design(df,include_family=False)
    n0=X.shape[1]
    for k in extras:
        v=df[k].to_numpy(float)
        if not np.isfinite(v).all(): raise ValueError("Missing extra source field "+k)
        sd=float(np.std(v))
        if sd < 1e-12: raise ValueError("Constant effort covariate "+k)
        X=np.column_stack([X,(v-v.mean())/sd])
        names.append(k)
    if np.linalg.matrix_rank(X)!=X.shape[1]:
        raise ValueError("Effort sensitivity collinear design: "+str((n0,len(names))))
    h6=df.h6.astype(str); ct=h6.value_counts()
    w=len(df)/(h6.nunique()*h6.map(ct).to_numpy(dtype=float))
    res=sm.WLS(y,X,weights=w).fit(cov_type="cluster",
         cov_kwds={"groups":h6.to_numpy(),"use_correction":True},use_t=False)
    j=names.index("W");beta,se=float(res.params[j]),float(res.bse[j])
    return {"beta_W_per_1sd_of_supported_cohort":beta,
            "H6_cluster_descriptive95":[beta-1.96*se,beta+1.96*se],
            "n_H6":h6.nunique(),"n_rows":len(df),
            "additional_effort_covariates":list(extras),
            "source_W_SD_in_supported_cohort":float(df.W.std(ddof=0)),
            "not_an_identified_direct_effect":True}


def trait_bounds(df, B, t, known, score, family, fams, mask):
    """Best/worst W coefficient for all family-known species independently allowed
    to vary between archived per-species documented-host lower/upper N bounds.
    H6 source fixed, controls fixed; not a probability or confidence interval.
    """
    sub=df.loc[mask].reset_index(drop=True)
    bix=B[np.flatnonzero(mask)]
    fidx=np.flatnonzero((fams==family)&known)
    sums=bix[:,fidx].sum(1)
    if np.any(sums<=0):raise ValueError("Empty supported scored-family denominator")
    X,names=design(sub,include_family=False)
    h6=sub.h6.astype(str); w=len(sub)/(h6.nunique()*h6.map(h6.value_counts()).to_numpy(float))
    j=names.index("W");x0=np.delete(X,j,axis=1);wz=X[:,j]
    rw=wz-x0@np.linalg.lstsq(x0*np.sqrt(w[:,None]),wz*np.sqrt(w),rcond=None)[0]
    coeff=(w*rw)@(bix[:,fidx]/sums[:,None])/np.sum(w*rw**2)
    lo=t.documented_host_missing_lower.to_numpy(float)[fidx]
    hi=t.documented_host_missing_upper.to_numpy(float)[fidx]
    if not np.isfinite(lo).all() or not np.isfinite(hi).all() or np.any(lo>hi):
        return {"status":"bound_source_missing_or_invalid"}
    baseline=float(coeff@score[fidx])
    base=cluster_fit(sub,sub["record_mean_"+family].to_numpy(float))["beta_W_per_1sd_this_cohort"]
    if abs(baseline-base)>1e-7:raise ValueError("Trait linear combination disagrees with main association")
    worst=float(np.sum(coeff*np.where(coeff>=0,lo,hi)))
    best=float(np.sum(coeff*np.where(coeff>=0,hi,lo)))
    return {"baseline_beta_recovered":baseline,
            "worst_allowed_fixed_design_beta":worst,
            "best_allowed_fixed_design_beta":best,
            "sign_can_change_with_documented_host_N_uncertainty":worst<=0<=best,
            "not_sampling_confidence_or_probability":True,
            "source_species_variants_do_not_cover_unlisted_hosts_or_117_unknown_keys":True}


def run(args):
    df,taxa,B,score,known,families,raw_fam,q,receipt=load_reconstruction(
        args.panel,args.scores,args.inputdir)
    df,means=build_family(df,B,score,known,families,raw_fam,q)
    result={"schema":"GESUI_FIVEFAMILY_HNS_POSTHOC_EFFORT_TAXON_TRAIT_SENSITIVITY_V1",
      "source_receipts":receipt,
      "posthoc_exploratory_not_corrected_for_multiple_comparisons":True,
      "original_support_thresholds":{"min_records":MIN_SCORED_RECORDS,
        "min_species":MIN_SCORED_SPECIES,"min_coverage":MIN_SCORE_COVERAGE},
      "families":{}}
    allloo=[]
    for family in FAMS:
        mask=eligible(df,family,MIN_SCORED_RECORDS,MIN_SCORED_SPECIES,MIN_SCORE_COVERAGE)
        sub=df.loc[mask].copy().reset_index(drop=True)
        y=sub["record_mean_"+family].to_numpy(float)
        x=B[np.flatnonzero(mask)]
        idx=np.flatnonzero((families==family)&known)
        c=x[:,idx].astype(float)
        ys=c@score[idx];total=c.sum(1)
        np.testing.assert_allclose(ys/total,y,atol=1e-10,rtol=0)
        orig=cluster_fit(sub,y)
        loorows=[]
        for k,i in enumerate(idx):
            n=c[:,k];remaining=total-n
            if np.any(remaining<=0):
                raise ValueError("Original two-species support not preserved after one species exclusion")
            yl=(ys-n*score[i])/remaining
            model=cluster_fit(sub,yl)
            loorows.append({"family":family,"butterflyKey":str(taxa.butterflyKey.iloc[i]),
                            "accepted_name":str(taxa.accepted_name.iloc[i]),
                            "score_fixed":float(score[i]),
                            "source_records_in_supported_periods":int(n.sum()),
                            "beta_without_species":model["beta_W_per_1sd_this_cohort"],
                            "change_beta_vs_full":model["beta_W_per_1sd_this_cohort"] -
                                                 orig["beta_W_per_1sd_this_cohort"],
                            "beta_without_species_ci_low":model["cluster95"][0],
                            "beta_without_species_ci_high":model["cluster95"][1]})
        loorows.sort(key=lambda x:abs(x["change_beta_vs_full"]),reverse=True)
        allloo.extend(loorows)
        sub["family_log1p_records"]=np.log1p(sub["count_records_"+family])
        sub["family_log1p_scored_species"]=np.log1p(sub["count_species_"+family])
        sub["family_scored_record_coverage"]=sub["source_score_coverage_"+family]
        additions=["family_log1p_records","family_log1p_scored_species",
                   "family_scored_record_coverage"]
        adjusted=fit_extra(sub,y,additions)
        cv=cv_add_value(sub,"record_mean_"+family,COMMON+additions,["W"])
        country=country_prediction(sub,"record_mean_"+family,numeric=COMMON+additions)
        b=trait_bounds(df,B,taxa,known,score,family,families,mask)
        family_result={
           "baseline_supported_rows":len(sub),"baseline_H6":sub.h6.nunique(),
           "baseline_W_beta":orig["beta_W_per_1sd_this_cohort"],
           "baseline_beta95":orig["cluster95"],
           "n_scored_species_in_family":len(idx),
           "all_leave_one_species_out_beta_min":float(min(x["beta_without_species"] for x in loorows)),
           "all_leave_one_species_out_beta_max":float(max(x["beta_without_species"] for x in loorows)),
           "n_leave_one_out_beta_positive":sum(x["beta_without_species"]>0 for x in loorows),
           "n_leave_one_out_ci_above_zero":sum(x["beta_without_species_ci_low"]>0 for x in loorows),
           "top_five_single_taxon_leverage_posthoc":loorows[:5],
           "family_effort_adjusted_model":{
             "W_association":adjusted,
             "H6_out_W_added_value":cv,
             "country_leave_one_out_H6_purged":country},
           "documented_host_uncertainty_fixed_design":b}
        if family=="Lycaenidae":
            eligibility=[]
            for name,n_min,sp_min,cov_min in [
               ("inclusive",1,1,.50),("original",10,2,.70),("strict",20,3,.80)]:
                z=eligible(df,family,n_min,sp_min,cov_min)
                sub2=df.loc[z].reset_index(drop=True)
                if sub2.h6.nunique()<25 or len(sub2)<100:
                    eligibility.append({"cohort":name,"n":len(sub2),"not_enough_support":True})
                    continue
                a=cluster_fit(sub2,sub2["record_mean_"+family].to_numpy(float))
                g=cv_add_value(sub2,"record_mean_"+family,COMMON,["W"])
                eligibility.append({"cohort":name,"n":len(sub2),
                  "n_H6":sub2.h6.nunique(),
                  "beta":a["beta_W_per_1sd_this_cohort"],
                  "beta95":a["cluster95"],
                  "H6_gain":g["H6_equal_MSE_gain_percent"],
                  "H6_gain95":g["H6_heldout_fixed_errors_bootstrap95_percent"]})
            family_result["posthoc_support_threshold_sensitivity"]=eligibility
            top3=[x["butterflyKey"] for x in loorows[:3]]
            drop=[i for i in idx if str(taxa.butterflyKey.iloc[i]) in top3]
            sel_mask=~np.isin(idx,drop)
            remaining=c[:,sel_mask].sum(1)
            if np.any(remaining<=0):
                family_result["remove_top3_influence"]="unsupported_nonzero_requirement"
            else:
                y2=(c[:,sel_mask]@score[idx[sel_mask]])/remaining
                fit=cluster_fit(sub,y2)
                sub3=sub.copy()
                sub3["HNS_without_3_posthoc"]=y2
                h6gain=cv_add_value(sub3,"HNS_without_3_posthoc",COMMON,["W"])
                family_result["remove_top3_influence_posthoc"]={
                  "species_keys":top3,"fixed_original_support":True,
                  "beta":fit["beta_W_per_1sd_this_cohort"],
                  "beta95":fit["cluster95"],
                  "H6_gain":h6gain["H6_equal_MSE_gain_percent"],
                  "H6_gain95":h6gain["H6_heldout_fixed_errors_bootstrap95_percent"]}
        result["families"][family]=family_result
        print("FAMILY_STRESS="+json.dumps({
            "family":family,"original_beta":family_result["baseline_W_beta"],
            "effort_adjusted_beta":adjusted["beta_W_per_1sd_of_supported_cohort"],
            "effort_added_H6_gain":cv["H6_equal_MSE_gain_percent"],
            "LOO_beta_min":family_result["all_leave_one_species_out_beta_min"],
            "LOO_beta_max":family_result["all_leave_one_species_out_beta_max"],
            "trait_bounds":b},allow_nan=False,separators=(',',':')),flush=True)
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"FAMILY_TAXON_EFFORT_SENSITIVITY.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    allloo.sort(key=lambda x:(x["family"],-abs(x["change_beta_vs_full"])))
    pd.DataFrame(allloo).to_csv(args.output/"ALL_186_TAXON_LEAVE_ONE_OUT.csv",index=False)
    print("FULL_FAMILY_TAXON_SENSITIVITY_JSON="+json.dumps(
        result,ensure_ascii=False,separators=(",",":"),allow_nan=False),flush=True)
    print("FULL_FAMILY_TAXON_LOO_CSV="+(args.output/"ALL_186_TAXON_LEAVE_ONE_OUT.csv").read_text().replace("\n","\\n"),flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--panel",type=Path,default=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"))
    ap.add_argument("--scores",type=Path,default=Path("results/20261115_sync/GESUI_303_HOST_NUTRIENT_SCORES.csv"))
    ap.add_argument("--inputdir",type=Path,default=Path("analysis/20261115_hns_family_taxon/download"))
    ap.add_argument("--output",type=Path,default=Path("results/20261115_hns_family_taxon"))
    run(ap.parse_args())
