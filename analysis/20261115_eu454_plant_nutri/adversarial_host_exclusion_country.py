#!/usr/bin/env python3
"""Adversarial SAME-H7 transfer & catalogued-larval-host exclusion diagnostics.

Use unchanged official v2 plant N traits and original GBIF 454 H7 flora,
and independent PRE-EXISTING 990-key potential host catalogue from PR86.
Remove known candidate hosts from the plant trait index; do not relabel other
plants as true nonhosts. Group-country holdout purges shared genuine H6.
This is explicitly POSTHOC source-consistent stress testing; not a mediation.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parent))
from h6_model_if_gate import (P_GEOG,B_GEOG,SEED,ALPHA,N_BOOT,
    ORIGINAL_PANEL_SHA,FACETS_SHA,EIV_SHA,
    source_cluster_regression,cv_h6,predict_ridge)

ROOT=Path(__file__).resolve().parents[2]
HOST_PATH=ROOT/"results/20261115_sync/network/host_plant_candidate_contributions.csv"
HOST_GIT_BLOB="6ee4dfa7b57137cc61966a5346eca63137dc9f05"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def gitsha(p):
    b=p.read_bytes()
    return hashlib.sha1(f"blob {len(b)}\0".encode()+b).hexdigest()
def read_csv(p,dtype=None):return pd.read_csv(p,dtype=dtype)
def source_records(original_dir):
    arr=list(original_dir.rglob("GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json"))
    if len(arr)!=1 or sha(arr[0])!=FACETS_SHA:
        raise RuntimeError("Original authentic 2808 source-facet JSON checksum mismatch")
    rr=json.loads(arr[0].read_text())
    if len(rr)!=2808 or any(x.get("source_status")!="ok" for x in rr):
        raise RuntimeError("Original 2808-facet source status changed")
    plants={(str(x["HYBAS_ID_L07"]),x["period"]):x
            for x in rr if x["group"]=="Plantae"}
    if len(plants)!=468:raise RuntimeError("Original 468 source Plantae facets absent")
    return plants

def prepare(original_dir,score_dir):
    gate=json.loads((score_dir/"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if gate.get("schema")!="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA" or \
        not gate.get("hard_source_support_pass") or \
        not gate.get("original_EIV_NUTRI_scored_keys_unchanged") or \
        gate.get("source_hybrid_score_assigned")!=0 or \
        gate.get("official_confirmed_valid_unscorable_nonbinomial_hybrids")!=175:
        raise RuntimeError("No original-source hybrid/refinement gate before outcome")
    if gate["source_EIV"]["original_source_SHA256"]!=EIV_SHA:raise RuntimeError("EIV source mismatch")
    if gitsha(HOST_PATH)!=HOST_GIT_BLOB:raise RuntimeError("Original independent 990 source host source drift")
    hosts=read_csv(HOST_PATH,dtype={"plantKey":str})
    if len(hosts)!=990 or hosts.plantKey.duplicated().any():
        raise RuntimeError("Original 990 potential host plants source changed")
    hostkeys=set(hosts.plantKey)
    scores=read_csv(score_dir/"EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv",
                    dtype={"speciesKey":str})
    if len(scores)!=11511 or scores.speciesKey.duplicated().any():
        raise RuntimeError("Original plant species crosswalk source changed")
    numeric={row.speciesKey:float(row.EIV_NUTRI_exact) for row in scores.itertuples()
             if pd.notna(row.EIV_NUTRI_exact)}
    if len(numeric)!=4011 or any(not 1<=v<=12 for v in numeric.values()):
        raise RuntimeError("Original EIV 4011 unchanged exact numeric scores invalid")
    original=ROOT/"results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"
    if sha(original)!=ORIGINAL_PANEL_SHA:
        raise RuntimeError("Original HNS 454 panel checksum mismatched")
    panel=read_csv(original,dtype={"h":str,"h6":str})
    if len(panel)!=454 or panel.h.nunique()!=232 or panel.h6.nunique()!=90:
        raise RuntimeError("Original true H7/H6 panel altered")
    elig=read_csv(score_dir/"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv",
                  dtype={"h":str,"h6":str})
    trait=read_csv(score_dir/"EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv",
                   dtype={"h":str})
    if len(elig)!=454 or len(trait)!=451:
        raise RuntimeError("Original 454 and original 451 gate-supported source trait rows inconsistent")
    eligible=elig.loc[elig.full_source_trait_eligible.astype(str)=="True"]
    if len(eligible)!=451:raise RuntimeError("Not original source eligible 451 H7-periods")
    df=panel.merge(eligible[["h","period","NUTRI_coverage_all_Plantae_speciesKeys",
                             "exact_NUTRI_scored_vascular_accepted_species"]],
                   on=["h","period"],validate="one_to_one",how="inner")
    df=df.merge(trait[["h","period","Plant_N_incidence","Plant_N_record"]],
                on=["h","period"],validate="one_to_one",how="inner")
    if len(df)!=451 or df.h6.nunique()!=90:raise RuntimeError("Original 451 source grain mismatch")
    df["log1p_eiv_scored_taxa"]=np.log1p(df.exact_NUTRI_scored_vascular_accepted_species.astype(float))
    df["country_method"]=(df.country.fillna("__UNKNOWN__").astype(str)+"::"+
                          df.method.fillna("__UNKNOWN__").astype(str))
    plants=source_records(original_dir)
    nonhosts=[]
    for row in df.itertuples(index=False):
        k=(row.h,row.period)
        facet=plants.get(k)
        if facet is None:raise RuntimeError("Absent original H7 plant facet "+repr(k))
        allcounts={str(key):int(v) for key,v in facet["source_speciesKey_counts"].items() if int(v)>0}
        other_keys={key:v for key,v in allcounts.items() if key not in hostkeys}
        scored=[(numeric[key],v) for key,v in other_keys.items() if key in numeric]
        nonhosts.append({"h":row.h,"period":row.period,
            "n_all_original_source_plants":len(allcounts),
            "n_catalogue_host_plant_records":sum(v for key,v in allcounts.items() if key in hostkeys),
            "n_original_no_catalogue_species":len(other_keys),
            "n_non_catalogue_eiv_scored_species":len(scored),
            "fraction_non_catalogue_scored":len(scored)/len(other_keys) if other_keys else 0.,
            "non_catalogue_incidence_N":float(np.mean([x[0] for x in scored])) if scored else np.nan,
            "non_catalogue_record_N":float(np.average([x[0] for x in scored],
                 weights=[x[1] for x in scored])) if scored else np.nan})
    df=df.merge(pd.DataFrame(nonhosts),on=["h","period"],validate="one_to_one",how="inner")
    unionkeys=set(numeric)
    return df,hostkeys,unionkeys,numeric,gate

def country_predictive(df,target,baseline,extra):
    countries=df.country_holdout.astype(str).to_numpy()
    h6=df.h6.astype(str).to_numpy()
    y=df[target].to_numpy(float)
    rows=[]
    for c in sorted(set(countries)):
        test=np.flatnonzero(countries==c)
        if len(test)<10 or len(set(h6[test]))<3:continue
        train=np.flatnonzero((countries!=c)&~np.isin(h6,h6[test]))
        if len(train)<50 or len(set(h6[train]))<10:continue
        if set(h6[train])&set(h6[test]):raise RuntimeError("country H6 leak")
        tr=df.iloc[train].reset_index(drop=True);te=df.iloc[test].reset_index(drop=True)
        base=predict_ridge(tr,te,target,baseline)
        add=predict_ridge(tr,te,target,baseline+[extra])
        loss0=(y[test]-base)**2;loss1=(y[test]-add)**2
        rows.append({"country_holdout":c,"n_H7period":len(test),
                     "n_true_H6":len(set(h6[test])),
                     "baseline_MSE":float(loss0.mean()),
                     "added_MSE":float(loss1.mean()),
                     "country_MSE_improvement_pct":float(100*(loss0.mean()-loss1.mean())/loss0.mean())})
    if not rows:return {"status":"NO_ELIGIBLE_COUNTRY_REPLICATION"}
    a=np.array([x["baseline_MSE"] for x in rows]);b=np.array([x["added_MSE"] for x in rows])
    rng=np.random.default_rng(SEED)
    ix=rng.integers(0,len(rows),size=(N_BOOT,len(rows)))
    sim=100*(a[ix].mean(axis=1)-b[ix].mean(axis=1))/a[ix].mean(axis=1)
    return {"n_eligible_countries":len(rows),
            "country_equal_added_MSE_gain_pct":float(100*(a.mean()-b.mean())/a.mean()),
            "country_equal_fixed_errors_bootstrap95":np.quantile(sim,[.025,.975]).tolist(),
            "country_holdout_source_H6_purged":True,
            "source_country_results":rows,
            "country_resampling_descriptive_not_study_meta_analysis":True}

def main(a):
    df,hosts,unionscores,numeric,gate=prepare(a.facets,a.source)
    support=((df.n_non_catalogue_eiv_scored_species>=10) &
             (df.fraction_non_catalogue_scored>=.40))
    sub=df.loc[support].copy().reset_index(drop=True)
    external={
        "W_to_original_plant_N":country_predictive(df,"Plant_N_incidence",P_GEOG,"W"),
        "original_plant_N_to_incidence_butterfly_HNS":country_predictive(
            df,"primary_incidence_mean",B_GEOG,"Plant_N_incidence"),
        "original_plant_N_to_record_butterfly_HNS":country_predictive(
            df,"primary_record_mean",B_GEOG,"Plant_N_incidence")}
    out={
        "schema":"GESUI_EU454_H6_COUNTRY_TRANSFER_AND_NO_CATALOGUED_HOST_PLANTS_V1",
        "posthoc_adversarial_not_causal_or_preregistered":True,
        "original_H7period":451,"original_true_H6":90,
        "original_EIV_numeric_species":len(unionscores),
        "independent_previous_catalogued_host_plant_keys":len(hosts),
        "original_scored_plant_species_in_host_catalogue":len(unionscores&hosts),
        "original_scored_plant_species_not_in_host_catalogue":len(unionscores-hosts),
        "country_holdout":external,
        "noncatalogue_score_gates":{"min_scored_non_catalogue_plants":10,
           "min_scored_non_catalogue_fraction":.4},
        "non_catalogue_supported_H7period":len(sub),
        "non_catalogue_supported_true_H6":int(sub.h6.nunique()),
        "non_catalogue_supported_country_groups":int(sub.country_holdout.nunique()),
        "strict_non_catalogue_score_cohort_support_pass":bool(
           len(sub)>=100 and sub.h6.nunique()>=40 and sub.country_holdout.nunique()>=5),
        "known_990_plants_are_potential_hosts_not_verified_local_host_use":True,
        "remaining_plants_are_unknown_host_status_NOT_true_nonhosts":True,
        "host_exclusion_does_not_prove_ecological_food_web_mediation":True,
        "original_4011_trait_scores_unmodified":True}
    if out["strict_non_catalogue_score_cohort_support_pass"]:
        out["same_eligible_cohort_models"]={}
        for target in ("primary_incidence_mean","primary_record_mean"):
            models={}
            for feature in ("Plant_N_incidence","non_catalogue_incidence_N"):
                models[feature]={
                    "W_adjusted_H6_cluster_association":source_cluster_regression(
                        sub,target,B_GEOG+[feature],feature),
                    "H6_out_addition":cv_h6(sub,target,B_GEOG,feature),
                    "entire_country_out_transfer":country_predictive(sub,target,B_GEOG,feature)}
            out["same_eligible_cohort_models"][target]=models
        out["W_to_non_catalogue_plant_N"]={
            "association":source_cluster_regression(sub,"non_catalogue_incidence_N",P_GEOG+["W"],"W"),
            "H6_out_W_increment":cv_h6(sub,"non_catalogue_incidence_N",P_GEOG,"W")}
        out["full_versus_non_catalogue_N_correlation"]=float(
            np.corrcoef(sub.Plant_N_incidence,sub.non_catalogue_incidence_N)[0,1])
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"EU454_HNS_HOST_EXCLUSION_COUNTRY_TRANSFER.json").write_text(json.dumps(
        out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    summary={"n_original_H7period":451,
             "n_valid_non_catalogue_H7period":len(sub),
             "n_valid_non_catalogue_H6":int(sub.h6.nunique()),
             "scored_total":len(unionscores),"scored_excluded_host_catalogue":len(unionscores&hosts),
             "country_transfer_original":{k:{
                "country_gain":v.get("country_equal_added_MSE_gain_pct"),
                "country_gain95":v.get("country_equal_fixed_errors_bootstrap95"),
                "n_country":v.get("n_eligible_countries")}
                for k,v in external.items()},
             "posthoc_causal_mediation":False}
    if "same_eligible_cohort_models" in out:
        summary["noncatalogue_vs_original"]={
            k:{t:{"beta":v["W_adjusted_H6_cluster_association"]["adjusted_beta_per_1SD_focal"],
                  "beta95":v["W_adjusted_H6_cluster_association"]["H6_cluster_descriptive95"],
                  "h6_gain":v["H6_out_addition"]["H6_equal_added_MSE_gain_pct"],
                  "gain95":v["H6_out_addition"]["fixed_out_of_fold_error_H6_resampling95_pct"],
                  "country_gain":v["entire_country_out_transfer"].get("country_equal_added_MSE_gain_pct"),
                  "country_gain95":v["entire_country_out_transfer"].get("country_equal_fixed_errors_bootstrap95")}
                  for t,v in val.items()}
                 for k,val in out["same_eligible_cohort_models"].items()}
        summary["noncatalogue_W_to_plant_N_beta"]=out["W_to_non_catalogue_plant_N"]["association"]["adjusted_beta_per_1SD_focal"]
        summary["full_noncat_N_r"]=out["full_versus_non_catalogue_N_correlation"]
    print("GESUI_H7_NO_CATALOGUE_HOST_AND_COUNTRY_TRANSFER="+json.dumps(
        summary,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--facets",required=True,type=Path)
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())

