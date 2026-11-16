#!/usr/bin/env python3
"""Actual 454 source GBIF butterflies with TRUE same-individual HOST taxa EIV N+F.

Rebuild original source B (454x303) from publisher-verified 2808 GBIF facets.
Join 303-butterfly strictly exact original 2014 Ferrer-Paris host plant names
paired on identical host plant species with original European EIV N+MOIST,
n162 eligible accepted butterflies from original source SHA.

Compare plant CWM MOIST and aquatic EQR source association on exact same H7
period and authentic H6 parent holdout. No actual local butterfly larval
feeding, same-day stream P, sewage effluent exposure or water→land cause.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
from rebuild_H7_hostN_F_and_water_land_provisional import (
    genuine_original_data,source_cluster_regression,cv_h6,country_predictive,
    SCORES,PANEL,SCORES as ORIGINAL_ROSTER,FACET_SHA,FACET_NAME,
    FAMS,SCORES as SCORES_SOURCE
)

SCHEMA="GESUI_20261117_H7_BUTTERFLY_SAME_ACTUAL_FERRERPARIS_LARVALHOST_EIV_NF_MATCHED_PLANT_AND_WATER_EQR_V1"
STRICT_HOST_SCHEMA="GESUI_20261117_FERRERPARIS2014_EIV2023_SAME_LARVAL_HOST_N_MOIST_EXACT_SOURCE_GATE_V1"
HOST_DERIVED_CSV_SHA="a59eac46718b426bf8c057ab1d7197c4860bf4f8e9e4ba5aba7f5895e6e83268"
ORIGINAL_NUMERIC_N_HOSTS=162
ORIGINAL_FAMILIES={"Hesperiidae":22,"Lycaenidae":29,"Nymphalidae":81,"Papilionidae":10,"Pieridae":20}
MIN_SITE_SCORING=2
MIN_FACET_SCORING_RECORDS=10
MIN_MODEL_H7_PERIOD=200
MIN_MODEL_PARENT_H6=40
MIN_COUNTRIES=5

BASE=[
 "crop","urban","riverflow","naturalout","area","nEEA",
 "nAnimal","nPlant","nButterfly","W",
 "primary_record_coverage_raw","primary_host_completeness",
 "source_matched_host_butterfly_record_coverage",
 "source_matched_host_log_scored_species",
 "source_matched_2014_host_double_scored_coverage",
 "plant_pair_EIV_coverage","plant_pair_log_scored_species",
 "Bfly_same_host_N_incidence"]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def find(folder,name):
    paths=list(folder.rglob(name))
    if len(paths)!=1:raise RuntimeError("Expected exactly one historical original source file "+name)
    return paths[0]
def load_original_2014_host_score(source_dir):
    gatepath=find(source_dir,"GESUI_ORIGINAL_HOST_N_F_SAME_PLANT_SOURCE_GATE.json")
    gate=json.loads(gatepath.read_text())
    if (gate["schema"]!=STRICT_HOST_SCHEMA or
        not gate["source_min40_accepted_bfly_and4_families"] or
        gate["source_derived_303_original_taxon_trait_CSV_SHA256"]!=HOST_DERIVED_CSV_SHA):
        raise RuntimeError("Original true SAME HOST EIV N/F certified source is not available")
    path=find(source_dir,"GESUI_ORIGINAL_303_FERRERPARIS_EIV_SAME_HOST_N_F_SOURCE_ONLY.csv")
    if sha(path)!=HOST_DERIVED_CSV_SHA:
        raise RuntimeError("Authenticated original derived exact 2014 source N and F paired scores changed")
    scores=pd.read_csv(path,dtype={"butterflyKey":str})
    if len(scores)!=303 or scores.butterflyKey.duplicated().any():
        raise RuntimeError("Original 303 butterfly accepted-name source per HOSTS not unique")
    subset=scores.loc[scores.eligible_strict_same_HOST_N_F_and_accepted_butterfly.astype(str)=="True"]
    if len(subset)!=ORIGINAL_NUMERIC_N_HOSTS:
        raise RuntimeError("Original strictly matched same HOST source 162 butterfly roster changed")
    if collections.Counter(subset.original_butterfly_family)!=ORIGINAL_FAMILIES:
        raise RuntimeError("Original 5 family strictly scored butterfly quality changed")
    if ((subset.source_2014_host_taxa_with_BOTH_EIV_N_and_MOIST<1).any() or
        (subset.source_2014_host_taxa_paired_trait_coverage<=0).any()):
        raise RuntimeError("Strict host pair denominator EIV species availability invalid")
    return scores,gate

def reconstructed_original_gbif_community(root,panel,oldroster,score):
    facet=find(root,FACET_NAME)
    if sha(facet)!=FACET_SHA:raise RuntimeError("Original 2808 GBIF species source SHA failure")
    source=json.loads(facet.read_text())
    if len(source)!=2808:
        raise RuntimeError("True original source 2808 family/plant facets missing")
    lookup={(str(x["HYBAS_ID_L07"]),str(x["period"]),x["group"]):x for x in source}
    if len(lookup)!=2808:
        raise RuntimeError("Duplicated source H7×period family facets")
    keys=oldroster.butterflyKey.astype(str).tolist()
    ix={k:i for i,k in enumerate(keys)}
    mask=np.array([str(k) in score for k in keys],dtype=bool)
    if int(mask.sum())!=162:
        raise RuntimeError("Original 303 source roster does not fully match 162 scored original butterflyKeys")
    N=np.array([float(score[k]["N"]) for k in keys if k in score])
    F=np.array([float(score[k]["F"]) for k in keys if k in score])
    HC=np.array([float(score[k]["HC"]) for k in keys if k in score])
    B=np.zeros((len(panel),len(keys)),dtype=np.int64)
    total=[]
    for j,row in enumerate(panel.itertuples(index=False)):
        count=0
        for fam in FAMS:
            z=lookup[(row.h,row.period,fam)]
            count+=int(z["source_query_occurrence_count"])
            for k,v in z["source_speciesKey_counts"].items():
                if k not in ix or oldroster.loc[ix[k],"family"]!=fam:
                    raise RuntimeError("Original key/family mismatch in stored historical 4.318m GBIF source")
                if B[j,ix[k]]:raise RuntimeError("Original one butterfly species counted more than once in basin")
                B[j,ix[k]]=int(v)
        total.append(count)
    if (int(B.sum()),sum(total))!=(4318269,4400147):
        raise RuntimeError("Original butterfly GBIF 454 source reconstructed denominator not authentic")
    C=B[:,mask];P=(C>0);n=P.sum(axis=1);nr=C.sum(axis=1)
    with np.errstate(divide="ignore",invalid="ignore"):
        panel["Bfly_same_host_N_incidence"]=np.divide(P@N,n,
            out=np.full(len(panel),np.nan),where=n>0)
        panel["Bfly_same_host_F_incidence"]=np.divide(P@F,n,
            out=np.full(len(panel),np.nan),where=n>0)
        panel["Bfly_same_host_F_record"]=np.divide(C@F,nr,
            out=np.full(len(panel),np.nan),where=nr>0)
        panel["source_matched_host_butterfly_record_coverage"]=nr/np.asarray(total,float)
        panel["source_matched_2014_host_double_scored_coverage"]=np.divide(
            P@HC,n,out=np.full(len(panel),np.nan),where=n>0)
    panel["source_matched_host_butterfly_scored_species"]=n
    panel["source_matched_host_butterfly_scored_records"]=nr
    panel["source_matched_host_log_scored_species"]=np.log1p(n)
    return panel,{"original_B_matrix_454_speciesKey_record_total":int(B.sum()),
        "source_bfly_family_raw_record_total":sum(total),
        "source_strict_same_hostN_F_plant_accepted_butterfly_species":len(score),
        "source_strict_2014_host_all_original_true":True}

def regress(df,y,base,focus,country=False):
    out={"source_H6_WLS":source_cluster_regression(df,y,base+[focus],focus),
       "source_H6_disjoint_Ridge_added_prediction":cv_h6(df,y,base,focus)}
    if country:out["source_country_H6_purged_prediction"]=country_predictive(df,y,base,focus)
    return out

def paired_time(df):
    cols=["Bfly_same_host_F_incidence","Bfly_same_host_N_incidence",
          "Plant_F_incidence","Plant_N_incidence",
          "EQR_original_site_equal_basin_period_mean"]
    t=df.pivot(index="h",columns="period",values=cols).dropna()
    if len(t)<40:return {"status":"STOP_TIME_PAIR_SOURCE_SUPPORT","n":len(t)}
    times=["2017-2019","2020-2022"]
    shift={col:t[(col,times[1])].to_numpy(float)-t[(col,times[0])].to_numpy(float) for col in cols}
    return {"source_two_periods_H7_same_match_n":len(t),
      "delta_fly_hostF_vs_plantF_Spearman":float(spearmanr(shift["Bfly_same_host_F_incidence"],shift["Plant_F_incidence"]).statistic),
      "delta_fly_hostF_vs_plantN_Spearman":float(spearmanr(shift["Bfly_same_host_F_incidence"],shift["Plant_N_incidence"]).statistic),
      "delta_fly_hostF_vs_EQR_Spearman":float(spearmanr(shift["Bfly_same_host_F_incidence"],shift["EQR_original_site_equal_basin_period_mean"]).statistic),
      "source_3year_periods_not_same_days_or_2022_actual_larval_use":True}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    # Previous independently certified original H7×plant N/F source and all303
    # old HNS reconstruction (the old 180 score study is not used as output).
    merged,oldvalid,origqa=genuine_original_data(a)
    original=pd.read_csv(ORIGINAL_ROSTER,dtype={"butterflyKey":str},encoding="utf-8-sig")
    original=original.reset_index(drop=True)
    scores,scoregate=load_original_2014_host_score(a.samehosts)
    score={
      str(z.butterflyKey):{"N":float(z.source_same_exact_host_NUTRI_mean),
                           "F":float(z.source_same_exact_host_MOIST_mean),
                           "HC":float(z.source_2014_host_taxa_paired_trait_coverage)}
      for z in scores.itertuples(index=False)
      if str(z.eligible_strict_same_HOST_N_F_and_accepted_butterfly).lower()=="true"}
    if len(score)!=162:raise RuntimeError("Original exact same-host EIV 2014 source score set not 162")
    # Avoid silent misalignment: original panel fixed source H7 row ordering
    # and original paired plant scoring restricted to n447 from pre-response gates.
    pairsource=merged[["h","period","h6","Plant_N_incidence","Plant_F_incidence",
                     "plant_pair_EIV_coverage","plant_pair_log_scored_species",
                     "source_bf_record_coverage"]]
    df=merged.drop(columns=["source_bf_record_coverage"],errors="ignore").copy()
    df,traitsQA=reconstructed_original_gbif_community(a.facets,df,original,score)
    # BF source read old site count doesn't determine strict 162 source sites.
    quality=(df.source_matched_host_butterfly_scored_species>=MIN_SITE_SCORING)&(
             df.source_matched_host_butterfly_scored_records>=MIN_FACET_SCORING_RECORDS)&(
             np.isfinite(df[["Bfly_same_host_N_incidence","Bfly_same_host_F_incidence",
                              "Bfly_same_host_F_record",
                              "Plant_N_incidence","Plant_F_incidence",
                              "EQR_original_site_equal_basin_period_mean"]]).all(axis=1))
    selected=df.loc[quality].copy().reset_index(drop=True)
    nH6=int(selected.h6.nunique())
    nC=int(selected.country_holdout.nunique())
    passed=len(selected)>=MIN_MODEL_H7_PERIOD and nH6>=MIN_MODEL_PARENT_H6 and nC>=MIN_COUNTRIES
    result={"schema":SCHEMA,
      "status":"PASS_SAME_PUBLISHED_HOST_PLANT_NF_SOURCE_BUT_CAUSAL_WATER_LAND_UNIDENTIFIED" if passed else
         "STOP_SAME_ORIGINAL_HOST_NF_COVERAGE_FOR_WATERLAND_MODELS",
      "original_source_facet_H7_and_previous_provisionalQA":origqa,
      "source_true_same_plant_host_N_F_FerrerParis_gate":scoregate,
      "original_strict_162_butterfly_species_source":traitsQA,
      "original_454_siteperiods_with_paired_plant_source":len(df),
      "source_siteperiods_with_geographic_species_support":len(selected),
      "source_true_H6_of_supported_sites":nH6,
      "source_countries_of_supported_sites":nC,
      "source_scoring_min_two_original_butterfly_species_and10_GBIF_records":True,
      "model_fitted_only_after_source_gate_pass":passed,
      "observed_source_larval_host_plants_from_world_literature_not_actual_local_use":True,
      "observed_butterfly_GBIF_species_keys_not_actual_individual_abundance":True,
      "observed_2022_EQR_basin_average_not_same_local_butterfly_site":True,
      "modeled_static_wastewater_W_not_measured_2022_site_dose":True,
      "no_true_water_to_land_causal_mediation_or_sewage_dependency":True}
    if passed:
        baseline=[
          "crop","urban","riverflow","naturalout","area","nEEA","nAnimal",
          "nPlant","nButterfly","W","primary_record_coverage_raw",
          "primary_host_completeness",
          "source_matched_host_butterfly_record_coverage",
          "source_matched_host_log_scored_species",
          "source_matched_2014_host_double_scored_coverage",
          "plant_pair_EIV_coverage","plant_pair_log_scored_species",
          "Bfly_same_host_N_incidence"]
        plantN=baseline+["Plant_N_incidence"]
        plantF=baseline+["Plant_F_incidence"]
        BOTH=baseline+["Plant_N_incidence","Plant_F_incidence"]
        res={
         "original_paired_plant_MOIST_added_after_N_and_bfly_samehostN":regress(
             selected,"Bfly_same_host_F_incidence",plantN,"Plant_F_incidence",country=True),
         "original_paired_plant_NUTRI_added_after_MOIST_and_bfly_samehostN":regress(
             selected,"Bfly_same_host_F_incidence",plantF,"Plant_N_incidence"),
         "actual_aquatic_EQR_after_BOTH_plant_traits_and_bfly_samehostN":regress(
             selected,"Bfly_same_host_F_incidence",BOTH,
             "EQR_original_site_equal_basin_period_mean",country=True),
         "wastewater_related_W_after_BOTH_plant_traits_and_bfly_samehostN":regress(
             selected,"Bfly_same_host_F_incidence",[x for x in BOTH if x!="W"],"W"),
         "GBIF_record_weighted_MOIST_after_samehost_N_and_plant_N":regress(
             selected,"Bfly_same_host_F_record",plantN,"Plant_F_incidence"),
         "two_period_source_correlations":paired_time(selected)
        }
        result["source_same_host_family_score_water_land_2022_observational_results"]=res
        result["source_strict_double_scored_host_vs_old_precomputed_butterfly_score_limitation"]=(
          "different exact accepted species sample and potentially different original host source; not independent field validation")
    (a.output/"GESUI_STRICT_SAME_HOST_NF_H7_PLANTF_EQR_CROSSDOMAIN_RESULTS.json").write_text(
      json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    summary={
       "status":result["status"],
       "strict_butterfly_species_162":len(score),
       "original_454_joint_plant_N_F_source":len(df),
       "same_host_source_H7_period_support":len(selected),
       "same_host_source_H6_support":nH6,
       "source_n_countries":nC,
       "not_a_realized_aquatic_terrestrial_trophic_effect":True}
    for key,val in result.get("source_same_host_family_score_water_land_2022_observational_results",{}).items():
        if "source_H6_WLS" in val:
            summary[key]={
             "conditional_beta_per_source_SD":val["source_H6_WLS"]["adjusted_beta_per_1SD_focal"],
             "conditional_beta95":val["source_H6_WLS"]["H6_cluster_descriptive95"],
             "H6_oof_MSE_gain_pct":val["source_H6_disjoint_Ridge_added_prediction"]["H6_equal_added_MSE_gain_pct"],
             "fixed_H6_error_gain95":val["source_H6_disjoint_Ridge_added_prediction"]["fixed_out_of_fold_error_H6_resampling95_pct"],
             "country_H6_purged_gain":val.get("source_country_H6_purged_prediction",{}).get("country_equal_added_MSE_gain_pct")}
        else:summary[key]=val
    print("GESUI_TRUE_SAME_HOST_NF_H7_WATER_LAND_SOURCE_RESULT="+json.dumps(
      summary,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--samehosts",required=True,type=Path)
    p.add_argument("--facets",required=True,type=Path)
    p.add_argument("--plants",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
