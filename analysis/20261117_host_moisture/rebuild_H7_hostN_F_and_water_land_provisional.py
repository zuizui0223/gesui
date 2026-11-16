#!/usr/bin/env python3
"""Exploratory same-H7/period butterfly host moisture X original water EQR
and paired-plant EIV MOIST; ALL scores from official source-locked GBIF/EIV.

Requires original 454×303 five-family butterfly speciesKey occurrence counts
with original 2808 successful facets. Scores N,MOIST both finite on SAME
observed butterfly species subset n180. Their underlying individual host
plants may differ between N/MOIST score means; hence this is only
a provisional cross-system score association, NOT confirmed functional
host-trait co-use or local wastewater→butterfly mechanism.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,math,sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261115_eu454_plant_nutri"))
from h6_model_if_gate import (cv_h6,source_cluster_regression)
from adversarial_host_exclusion_country import country_predictive

FACET_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
ROSTER_SHA="00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb"
PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
SCHEMA="GESUI_20261117_EU454_H7_BFLY_HOST_F_vs_PLANT_EIV_F_vs_REAL_AQUATIC_EQR_SAMESPECIES_PROVISIONAL_V1"
FACET_NAME="GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json"
SCORES=ROOT/"results/20261115_sync/GESUI_303_HOST_NUTRIENT_SCORES.csv"
PANEL=ROOT/"results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"
ORIG_PLANTS="EU454_EIV_NUTRI_FIVE_PAIRS_SOURCE.csv"
PAIR_GATE="EU454_EIV_NUTRI_FIVE_PAIR_SOURCE_GATE.json"
FAMS=["Hesperiidae","Lycaenidae","Nymphalidae","Papilionidae","Pieridae"]
# Primary comparison on common taxon denominator; NO trait search.
MIN_SOURCE_BFLY_SPECIES=2
MIN_SOURCE_BFLY_RECORDS=10
MIN_MODEL_H7PERIOD=200
MIN_MODEL_PARENT_H6=40
MIN_COUNTRY=5
BASE=[
 "crop","urban","riverflow","naturalout","area","nEEA","nAnimal","nPlant","nButterfly",
 "W","primary_record_coverage_raw","primary_host_completeness",
 "source_bf_record_coverage",
 "source_bf_log_scored_species",
 "plant_pair_EIV_coverage","plant_pair_log_scored_species",
 "Bfly_hostN_incidence"
]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def single_source(root,name):
    files=list(root.rglob(name))
    if len(files)!=1:
        raise RuntimeError("Genuine original source file missing or duplicated: "+name+" "+str(len(files)))
    return files[0]

def genuine_original_data(a):
    if sha(SCORES)!=ROSTER_SHA or sha(PANEL)!=PANEL_SHA:
        raise RuntimeError("Original immutable 303 taxon and H7 science panel source changed")
    taxa=pd.read_csv(SCORES,dtype={"butterflyKey":str},encoding="utf-8-sig")
    panel=pd.read_csv(PANEL,dtype={"h":str,"h6":str})
    if len(taxa)!=303 or taxa.butterflyKey.duplicated().any():
        raise RuntimeError("Original 303 butterfly taxon IDs changed")
    if len(panel)!=454 or panel.h.nunique()!=232 or panel.h6.nunique()!=90 or panel.country.nunique()!=15:
        raise RuntimeError("Original genuine 454 panel geography changed")
    f=single_source(a.facets,FACET_NAME)
    if sha(f)!=FACET_SHA:
        raise RuntimeError("Original GBIF source 2808 facet SHA256 mismatch")
    data=json.loads(f.read_text())
    if len(data)!=2808 or any(x["source_status"]!="ok" for x in data):
        raise RuntimeError("Original 2808 GBIF family/plant facets source compromised")
    look={(str(x["HYBAS_ID_L07"]),str(x["period"]),x["group"]):x for x in data}
    if len(look)!=2808:raise RuntimeError("Original authentic butterfly/plant H7×period facets duplicate")
    c=np.zeros((len(panel),len(taxa)),np.int64)
    raw=[]
    ixs={key:i for i,key in enumerate(taxa.butterflyKey)}
    for rid,rec in enumerate(panel.itertuples(index=False)):
        source_total=0
        seen=set()
        for fam in FAMS:
            x=look[(str(rec.h),str(rec.period),fam)]
            source_total+=int(x["source_query_occurrence_count"])
            for name,val in x["source_speciesKey_counts"].items():
                if name in seen or name not in ixs:raise RuntimeError("Species counted twice or key not in original 303")
                seen.add(name)
                col=ixs[name]
                if taxa.loc[col,"family"]!=fam:raise RuntimeError("Source family-mapped species differs")
                c[rid,col]=int(val)
        raw.append(source_total)
    if (int(c.sum()),sum(raw))!=(4318269,4400147):
        raise RuntimeError("Actual original 4.318m source butterfly key records not reconstructed")
    # Independent reproduction of the actual ORIGINAL 186-numeric HNS source,
    # before changing the butterfly scored set to 180 N and F same taxa.
    clean=(taxa.taxon_status=="ACCEPTED")&(
         ~taxa.actual_used_alias_conflict.astype(str).str.lower().isin(["true","yes","1"]))
    nscore=taxa.host_nutrient_mean.to_numpy(float)
    fscore=taxa.host_moisture_mean.to_numpy(float)
    origN=clean.to_numpy()&np.isfinite(nscore)
    if int(origN.sum())!=186:
        raise RuntimeError("Original primary nutrient HNS source 186 numeric accepted taxa missing")
    cp=c[:,origN]
    meanR=(cp@nscore[origN])/cp.sum(axis=1)
    meanI=((cp>0)@nscore[origN])/(cp>0).sum(axis=1)
    err=max(
      float(np.max(np.abs(meanR-panel.primary_record_mean.to_numpy(float)))),
      float(np.max(np.abs(meanI-panel.primary_incidence_mean.to_numpy(float)))))
    if err>1e-10:raise RuntimeError("Source matrix nutrient HNS reconstruction failed original primary scientific result "+str(err))
    # Same butterflies, but the underlying source host-plant list may differ by trait.
    valid=clean.to_numpy()&np.isfinite(nscore)&np.isfinite(fscore)
    if int(valid.sum())!=180:
        raise RuntimeError("Original 180 same-BUTTERFLY-species N and MOIST numeric roster support drift")
    fam=collections.Counter(taxa.loc[valid,"family"])
    if fam!={"Hesperiidae":23,"Lycaenidae":44,"Nymphalidae":85,"Papilionidae":10,"Pieridae":18}:
        raise RuntimeError("Original 180 scored butterfly family distribution drift")
    count=c[:,valid]
    inc=(count>0)
    den=inc.sum(axis=1);records=count.sum(axis=1)
    with np.errstate(divide="ignore",invalid="ignore"):
        panel["Bfly_hostN_incidence"]=np.divide(inc@nscore[valid],den,
                                                 out=np.full(len(panel),np.nan),where=den>0)
        panel["Bfly_hostF_incidence"]=np.divide(inc@fscore[valid],den,
                                                 out=np.full(len(panel),np.nan),where=den>0)
        panel["Bfly_hostF_record"]=np.divide(count@fscore[valid],records,
                                              out=np.full(len(panel),np.nan),where=records>0)
        panel["source_bf_record_coverage"]=records/np.asarray(raw,dtype=float)
    panel["source_bf_scored_species"]=den
    panel["source_bf_scored_records"]=records
    panel["source_bf_log_scored_species"]=np.log1p(den)
    # Host completeness was calculated for N; don't fabricate F coverage.
    # Same normalized HNS population but possibly different original EIV host sets.
    sg=json.loads(single_source(a.plants,PAIR_GATE).read_text())
    if (sg.get("schema")!="GESUI_EIV_5_NUTRI_PAIRWISE_SAME_SCORED_SPECIES_SUPPORT_V1" or
        sg.get("original_454_panel_SHA")!=PANEL_SHA or
        sg.get("original_2808_facets_SHA")!=FACET_SHA or
        sg.get("original_EIV_source_SHA")!=EIV_SHA or
        not sg.get("each_pair_on_exact_same_plants_for_both_N_and_other")):
        raise RuntimeError("Original plant EIV same-plant source gate not source-authenticated")
    mo=sg["pairwise_source"]["MOIST"]
    if (mo["paired_support_H7period"],mo["paired_support_true_H6"],
        mo["scored_species_keys_original_11511_union"])!=(447,89,3722):
        raise RuntimeError("Original verified same-plant N/F EIV source opportunity changed")
    pair=pd.read_csv(single_source(a.plants,ORIG_PLANTS),dtype={"h":str,"h6":str})
    if len(pair)!=454 or pair.duplicated(["h","period"]).any():
        raise RuntimeError("Original fixed same-H7×period paired plant trait source changed")
    if pair["source_pair_eligible_MOIST"].astype(str).eq("True").sum()!=447:
        raise RuntimeError("Original 447 eligible same-plant N/F site×period source cohort changed")
    p=pair.loc[pair.source_pair_eligible_MOIST.astype(str)=="True",
        ["h","period","pair_N_MOIST","pair_MOIST",
         "n_pair_MOIST","cov_pair_MOIST"]].rename(columns={
            "pair_N_MOIST":"Plant_N_incidence",
            "pair_MOIST":"Plant_F_incidence",
            "n_pair_MOIST":"n_pair_scored",
            "cov_pair_MOIST":"plant_pair_EIV_coverage"})
    merged=panel.merge(p,on=["h","period"],validate="one_to_one",how="inner")
    if len(merged)!=447 or merged.h6.nunique()!=89:
        raise RuntimeError("Source fixed paired plant EIV n447/n89 changed by butterfly merge")
    merged["plant_pair_log_scored_species"]=np.log1p(merged.n_pair_scored.astype(float))
    merged["country_method"]=merged.country.fillna("__UNKNOWN__").astype(str)+"::"+merged.method.fillna("__UNKNOWN__").astype(str)
    passed=(merged.source_bf_scored_species>=MIN_SOURCE_BFLY_SPECIES)&(
        merged.source_bf_scored_records>=MIN_SOURCE_BFLY_RECORDS)&(
        np.isfinite(merged[["Bfly_hostN_incidence","Bfly_hostF_incidence",
                    "Bfly_hostF_record","Plant_F_incidence","Plant_N_incidence",
                    "EQR_original_site_equal_basin_period_mean"]]).all(axis=1))
    support=merged.loc[passed].copy().reset_index(drop=True)
    detail={
     "original_facet_source_SHA256":FACET_SHA,"original_panel_source_SHA256":PANEL_SHA,
     "original_303_trait_roster_SHA256":ROSTER_SHA,
     "original_EIV_plant_source_SHA256":EIV_SHA,
     "original_454x303_butterfly_counts_original_sum":int(c.sum()),
     "original_454_family_raw_counts_sum":int(sum(raw)),
     "old_186score_HNS_both_original_response_reconstruction_max_absolute_error":err,
     "same_scored_butterfly_species_N_MOIST_count":int(valid.sum()),
     "source_five_family_N_F_same_butterfly_species":dict(fam),
     "plant_N_MOIST_paired_source_species_union":3722,
     "original_447_plant_pair_eligible_site_period":len(merged),
     "source_plant_EIV_eligible_H6":int(merged.h6.nunique()),
     "butterfly_score_quality_gate":{"min_same_bfly_species":MIN_SOURCE_BFLY_SPECIES,
                                     "min_scored_bfly_records":MIN_SOURCE_BFLY_RECORDS},
     "source_butterfly_joint_score_supported_rows":len(support),
     "source_butterfly_joint_score_supported_H6":int(support.h6.nunique()),
     "source_butterfly_joint_score_supported_country":int(support.country_holdout.nunique()),
     "source_butterfly_host_individual_plant_taxa_not_verified_same_between_N_and_MOIST":True,
     "source_butterfly_site_incidence_not_individual_density":True,
     "source_EQR_basin_period_aggregate_not_actual_same_butterfly_site":True}
    return merged,support,detail

def analyze_one(df,outcome,baseline,feature,country=False):
    kwargs={"WLS_H6_cluster":source_cluster_regression(df,outcome,baseline+[feature],feature),
        "H6_parent_purged_incremental_out_of_parent_prediction":cv_h6(df,outcome,baseline,feature)}
    if country:
        kwargs["country_holdout"]=country_predictive(df,outcome,baseline,feature)
    return kwargs

def paired_temporal(support):
    # Compare TWO periods within same H7 ONLY if both sources are measured,
    # and do not imply an exact 2022 source or direct river→land flux.
    vals=["Bfly_hostF_incidence","Plant_F_incidence",
          "Bfly_hostN_incidence","Plant_N_incidence",
          "EQR_original_site_equal_basin_period_mean"]
    paired=support.pivot(index="h",columns="period",values=vals)
    cols=sorted(support.period.unique().tolist())
    if cols!=["2017-2019","2020-2022"]:
        raise RuntimeError("Original two 3-year observational periods changed")
    paired=paired.dropna()
    times={"n_distinct_H7_two_periods_with_both_source_support":len(paired)}
    if len(paired)<30:
        times["status"]="STOP_TEMPORAL_H7_PAIR_SUPPORT"
        return times
    shifts={v:paired[(v,cols[1])].to_numpy(float)-paired[(v,cols[0])].to_numpy(float) for v in vals}
    for k in ("Plant_F_incidence","Plant_N_incidence","EQR_original_site_equal_basin_period_mean"):
        sp=spearmanr(shifts["Bfly_hostF_incidence"],shifts[k])
        times["source_delta_hostF_vs_delta_"+k+"_Spearman"]=float(sp.statistic)
    times["source_two_periods_are_not_2022_station_same_day_observations"]=True
    return times

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    merged,use,source=genuine_original_data(a)
    source_pass=(len(use)>=MIN_MODEL_H7PERIOD and
             use.h6.nunique()>=MIN_MODEL_PARENT_H6 and
             use.country_holdout.nunique()>=MIN_COUNTRY)
    result={
       "schema":SCHEMA,"status":"PASS_POSTHOC_PROVISIONAL_SAME_BUTTERFLY_SPECIES_TRAIT_CO_RESPONSE" if source_pass else "STOP_SAME_SCORED_BUTTERFLY_TRAIT_GEO_COVERAGE",
       "origin_original_source_integrity":source,
       "scientific_interpretation":"numeric butterfly species host habitat N and F derived possibly from different individual host plants; aquatic EQR is another observational compartment in a basin, not a causal input",
       "observational_after_previous_wastewater_findings_and_posthoc":True,
       "source_model_support_pass":bool(source_pass),
       "no_butterfly_larval_actual_host_eating_or_abundance":True,
       "no_true_measured_WWTP_delivered_N_P_to_host_plants":True,
       "no_direct_water_to_land_trophic_mediation_or_ecosystem_services":True}
    if source_pass:
        # Design is the SAME subset of source H7×period and same 180
        # butterfly species for both response scores.
        source_base=[
           "crop","urban","riverflow","naturalout","area","nEEA",
           "nAnimal","nPlant","nButterfly","W",
           "primary_record_coverage_raw","primary_host_completeness",
           "source_bf_record_coverage","source_bf_log_scored_species",
           "plant_pair_EIV_coverage","plant_pair_log_scored_species",
           "Bfly_hostN_incidence"]
        f_both=["Plant_N_incidence","Plant_F_incidence"]
        source_base_plantN=source_base+["Plant_N_incidence"]
        source_base_plantF=source_base+["Plant_F_incidence"]
        source_base_both=source_base+f_both
        if len(source_base_both)!=len(set(source_base_both)):
            raise RuntimeError("Original plant/host covariates repeated")
        results={
         "plant_f_added_after_hostN_and_paired_plantN":{
            **analyze_one(use,"Bfly_hostF_incidence",source_base_plantN,"Plant_F_incidence",country=True),
            "interpretation":"potential same H7 terrestrial host affinity matching, both plant traits use same scored plant taxa"},
         "plant_n_added_after_hostN_and_paired_plantF":analyze_one(
             use,"Bfly_hostF_incidence",source_base_plantF,"Plant_N_incidence"),
         "same_H7_actual_aquatic_EQR_added_after_plant_N_F_and_hostN":{
             **analyze_one(use,"Bfly_hostF_incidence",source_base_both,
              "EQR_original_site_equal_basin_period_mean",country=True),
             "interpretation":"aquatic/terrestrial CONDITIONAL spatial co-response, never EQR-caused butterfly host use"},
         "static_W_geographic_indicator_added_after_plant_and_hostN":{
            **analyze_one(use,"Bfly_hostF_incidence",
                [x for x in source_base_both if x!="W"],"W"),
            "interpretation":"modelled wastewater-related geography, not observed 2022 nutrient dose"},
         "observed_plantF_added_to_record_weighted_hostF_after_plantN_and_hostN":analyze_one(
            use,"Bfly_hostF_record",source_base_plantN,"Plant_F_incidence"),
         "H7_observed_two_period_delta_trait_correspondence":paired_temporal(use),
        }
        result["observational_provisional_models"]=results
    (a.output/"EU454_HOST_MOISTURE_PLANT_MOISTURE_AQUATIC_EQR_PROVISIONAL.json").write_text(
      json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    synopsis={"status":result["status"],"original_454":len(merged),"supported_rows":len(use),
      "supported_trueH6":int(use.h6.nunique()),"source_butterfly_joint180_species":180,
      "original_preexisting_HNS_numeric_reconstruction_error":source[
         "old_186score_HNS_both_original_response_reconstruction_max_absolute_error"],
      "not_same_individual_host_PLANTS_for_N_and_F":True}
    for label,value in result.get("observational_provisional_models",{}).items():
        if "WLS_H6_cluster" in value:
            synopsis[label]={
              "beta":value["WLS_H6_cluster"]["adjusted_beta_per_1SD_focal"],
              "cluster95":value["WLS_H6_cluster"]["H6_cluster_descriptive95"],
              "H6_MSE_gain":value["H6_parent_purged_incremental_out_of_parent_prediction"]["H6_equal_added_MSE_gain_pct"],
              "fixed_H6_errors95":value["H6_parent_purged_incremental_out_of_parent_prediction"]["fixed_out_of_fold_error_H6_resampling95_pct"],
              "country_gain":value.get("country_holdout",{}).get("country_equal_added_MSE_gain_pct")}
        else:synopsis[label]=value
    print("GESUI_H7_PROVISIONAL_HOSTMOISTURE_PLANTMOIST_EQR_ORIGINAL_SOURCE_RESULT="+json.dumps(
         synopsis,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--facets",required=True,type=Path)
    p.add_argument("--plants",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
