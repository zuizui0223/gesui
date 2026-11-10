#!/usr/bin/env python3
"""Original global WWp × Plantae, Animalia and Fungi GBIF RECORDING comparator.

New two-kingdom source counts were acquired only after exact 1055/752 source
manifest and original positive Plantae result were already known. Distinguish
cross-kingdom indexed recording from direct sewage-driven ecological effects.
No source-site addition or removal using the sign of any biological endpoint.
"""
from __future__ import annotations

import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.formula.api as smf

D=Path("results/global_all_Plantae_WWp2015_20261110")
FROZEN=D/"ORIGINAL_ALL1055_42COUNTRY_WWp2015_BEFORE_NEW808_GBIF_SOURCE.json"
PLANT=D/"ORIGINAL_ALL1055_42COUNTRY_GBIF_PLANTAE_2011_24_50RECORD_RAREFACTION_SOURCE.json"
ORIGINAL=D/"GLOBAL_ALL1055_WWp2015_ALL_PLANTAE_50RAREFIED_42COUNTRY_MACRO_AND_COUNTRY_HOLDOUT.json"
COUNTS=D/"GLOBAL1055_WWp2015_GBIF_ANIMALIA_FUNGI_COUNT_SOURCE_BEFORE_EFFECT.json"

BASE="log_2015_WWTP_count+worldcover2021_built+worldcover2021_crop+abs_lat+sin_lon+cos_lon"
WWP="log_2015_modelled_WWp"
REC="log1p_plant_records"
RARE="log1p_50rarefied_Plantae"
RAW="log1p_raw_Plantae_speciesKey"
ANI="log1p_Animalia_indexed_records"
FUN="log1p_Fungi_indexed_records"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def prepare(root):
    files={key:root/path for key,path in {
        "WPI_1055_WWp_source":FROZEN,
        "original_GBIF_1055_Plantae":PLANT,
        "earlier_full_752_plant_model":ORIGINAL,
        "new_Animalia_Fungi_count_source":COUNTS
    }.items()}
    records={key:json.loads(p.read_text()) for key,p in files.items()}
    source=records["WPI_1055_WWp_source"]
    plant=records["original_GBIF_1055_Plantae"]
    previous=records["earlier_full_752_plant_model"]
    new=records["new_Animalia_Fungi_count_source"]
    if source["schema"]!="gesui_original_1055_WPI_ports42country_808NEW_2011_24_ALLPlantae_PREOUTCOME_frame_v01":
        raise ValueError("Wrong original WWp WPI source site frame")
    if plant["original_all1055_WWp_WPI_country_source_frame_SHA256"]!=sha(files["WPI_1055_WWp_source"]):
        raise ValueError("Original all-Plantae observations no longer fit WPI source")
    if new["original_source1055_WPI2019_WWp2015_manifest_SHA256"]!=sha(files["WPI_1055_WWp_source"]):
        raise ValueError("Cross-kingdom counts not from exactly original 1055 WWp WPI port sites")
    if new["previously_frozen_original_1055_allPlantae_source_SHA256"]!=sha(files["original_GBIF_1055_Plantae"]):
        raise ValueError("Cross-kingdom count site source not from ORIGINAL plant source")
    if previous["original_247old_plus_808new_allPlantae_GBIF_1055_result_SHA256"]!=sha(files["original_GBIF_1055_Plantae"]):
        raise ValueError("Earlier full 752 WWp plant models not on exact same GBIF source")
    original_model=previous["all1055_2015WWp_and_2011_24_all_Plantae_countryFE_models"]
    if previous["original_all_Plantae_full1055_model_actual_sample"]["n_original_port_records"]!=752:
        raise ValueError("Original 752 plant study case count changed")
    gate=new["precommitted_SOURCE_QC_before_effect"]
    if new["original_port_source_N"]!=1055 or new["original_source_countries_N"]!=42:
        raise ValueError("Original 1055 source quality count changed")
    if not gate["cross_kingdom_recording_control_statistical_models_PERMITTED"]:
        return files,source,plant,previous,new,None,"PRECOMMITTED_CROSSKINGDOM_SOURCE_COUNT_AVAILABILITY_FAILED"
    orig={x["original_port_id"]:x for x in source["source_all_1055_public_country_port_site_metadata_no_new_GBIF_response"]}
    floral={x["original_port_id"]:x for x in plant["original_1055_port_all_Plantae_source_GBIF_facet_summaries"]}
    others={x["original_port_id"]:x for x in new["all_original_source_port_animal_and_fungi_indexed_counts"]}
    if len(orig)!=len(floral)!=len(others)!=1055 or set(orig)!=set(floral) or set(orig)!=set(others):
        raise ValueError("Original source full1055 IDs not identical across three taxa")
    rows=[]
    complete=Counter()
    for pid in sorted(orig):
        o=orig[pid];p=floral[pid];n=others[pid]
        if o["original_source_country"]!=p["source_country"] or o["original_source_country"]!=n["original_source_country"]:
            raise ValueError("Original all-WPI source country differs from GBIF counts")
        if not p["source_quality_min50_faceted_records_min10_species_coverage70pct"]:
            continue
        complete["plant752"]+=1
        if not all(n[t]["GBIF_status"]=="COMPLETE" for t in ("Animalia","Fungi")):
            complete["animal_or_fungi_API_missing_not_zero"]+=1
            continue
        animals=n["Animalia"]["original_indexed_record_count"]
        fungi=n["Fungi"]["original_indexed_record_count"]
        if type(animals) is not int or type(fungi) is not int or animals<0 or fungi<0:
            raise ValueError("Original two-kingdom indexed record count invalid")
        if n["source_original_all_Plantae_quality_752_flag"] is not True:
            raise ValueError("Source GBIF comparator altered original plant quality site eligibility")
        lat=float(o["original_port_lat"])
        lon=float(o["original_port_lon"])
        plant_records=int(p["indexed_GBIF_all_Plantae_2011_2024_records"])
        raw_species=int(p["faceted_original_distinct_Plantae_speciesKeys"])
        rare=float(p["original_exact_hypergeometric_50_plant_record_rarefied_speciesKey"])
        if not(plant_records>=50 and raw_species>=10 and 1<=rare<=50):
            raise ValueError("Previously source accepted all-Plantae individualCount/rarefaction QC no longer valid")
        rows.append({
           "port_id":pid,"country":str(o["original_source_country"]),
           WWP:float(o["source_WWp_log1p"]),
           "log_2015_WWTP_count":math.log1p(float(o["source_WWTP_count_25km"])),
           "worldcover2021_built":float(o["source_ESA2021_built_fraction"]),
           "worldcover2021_crop":float(o["source_ESA2021_crop_fraction"]),
           "abs_lat":abs(lat),"sin_lon":math.sin(math.radians(lon)),"cos_lon":math.cos(math.radians(lon)),
           REC:math.log1p(plant_records),
           RARE:math.log1p(rare),
           RAW:math.log1p(raw_species),
           ANI:math.log1p(animals),
           FUN:math.log1p(fungi),
           "log1p_animal_and_fungi_records":math.log1p(animals+fungi),
           "log1p_plants_minus_log1p_animals":math.log1p(plant_records)-math.log1p(animals),
           "log1p_plants_minus_log1p_fungi":math.log1p(plant_records)-math.log1p(fungi),
           "original_recorded_Plantae":plant_records,
           "original_recorded_Animalia":animals,
           "original_recorded_Fungi":fungi,
           "original_50record_rarefied_Plantae_species":rare
        })
    if complete["plant752"]!=752:
        raise ValueError("Original full all-Plantae 752 source sites no longer fixed")
    df=pd.DataFrame(rows).sort_values(["country","port_id"])
    if len(df)<700 or df.country.nunique()<30:
        raise ValueError("Previously passed source gate but model site support too sparse")
    by=df.groupby("country").agg(n=("port_id","size"),k=(WWP,"nunique"))
    usable=by.index[(by.n>=2)&(by.k>=2)]
    df=df[df.country.isin(usable)].copy()
    if df.country.nunique()<25 or len(df)<690:
        return files,source,plant,previous,new,None,"POST_SUCCESS_COMPARATOR_COUNTRY_SUPPORT_FAIL"
    return files,source,plant,previous,new,df,None

def fit_one(df,response,extra=""):
    formula=response+"~"+BASE+"+"+WWP+extra+"+C(country)"
    m=smf.ols(formula,data=df).fit(cov_type="cluster",cov_kwds={"groups":df.country})
    return {
        "response_definition":response,
        "formula":formula,
        "source_port_count":len(df),
        "source_country_count":int(df.country.nunique()),
        "beta_2015_modelled_WWp":float(m.params[WWP]),
        "country_cluster95":list(map(float,m.conf_int().loc[WWP])),
        "nominal_p":float(m.pvalues[WWP]),
        "R2":float(m.rsquared),
        "posthoc_not_causal_sewage_ecology":True
    }

def countries(df,response):
    table=[]
    for country,g in df.groupby("country",sort=True):
        if len(g)>=8 and g[WWP].nunique()>=3 and g[response].nunique()>=3:
            s=float(spearmanr(g[WWP],g[response]).statistic)
            table.append({"country":str(country),"sites":len(g),"Spearman_WWp_recording":s})
    return {
        "source_eligible_countries_with_at_least8_ports":len(table),
        "country_source_spearman_positive_count":sum(r["Spearman_WWp_recording"]>0 for r in table),
        "country_source_spearman_negative_count":sum(r["Spearman_WWp_recording"]<0 for r in table),
        "individual_country_observational_ranks":table,
        "not_independent_randomized_country_treatment_effect":True
    }

def run(root):
    files,source,plant,previous,new,df,stop=prepare(root)
    out={
        "schema":"gesui_WWp2015_original_global1055_Plantae_Animalia_Fungi_recording_comparator_countryFE_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "source_files_SHA256":{k:sha(v) for k,v in files.items()},
        "preoutcome_source_gate_original":new["precommitted_SOURCE_QC_before_effect"],
        "original_752_GBIF_allPlantae_quality_sites_frozen_before_new_Animalia_Fungi_counts":True,
        "crosskingdom_source_model_allowed":stop is None,
        "stop_reason":stop,
        "global_original_WWp2015_production_NOT_released_nitrogen_phosphorus":True,
        "Animalia_Fungi_record_counts_are_NOT_biologically_unaffected_WWp_negative_controls":True,
        "original_2011_24_GBIF_index_varies_with_taxonomic_indexing_and_recording":True,
        "plant_rarefied50_richness_is_not_censused_plant_community_diversity":True,
        "all_macro_models_POSTHOC_after_plant_results":True
    }
    if stop:return out
    out["source_case_support"]={
        "n_sites":len(df),"n_original_countries":int(df.country.nunique()),
        "original_full_plant752_subset_n":new["previous_plant_quality752_ports_with_both_kingdom_queries_successful"],
        "observed_zero_Animalia_record_sites":int((df.original_recorded_Animalia==0).sum()),
        "observed_zero_Fungi_record_sites":int((df.original_recorded_Fungi==0).sum()),
        "median_Plant_GBIF_records":float(df.original_recorded_Plantae.median()),
        "median_Animalia_GBIF_records":float(df.original_recorded_Animalia.median()),
        "median_Fungi_GBIF_records":float(df.original_recorded_Fungi.median()),
        "country_source_site_sizes":dict(Counter(df.country))
    }
    names={
        "Plantae_record_count":(REC,""),
        "Animalia_record_count":(ANI,""),
        "Fungi_record_count":(FUN,""),
        "Combined_Animalia_Fungi_record_count":("log1p_animal_and_fungi_records",""),
        "Plantae_vs_Animalia_log_record_ratio":("log1p_plants_minus_log1p_animals",""),
        "Plantae_vs_Fungi_log_record_ratio":("log1p_plants_minus_log1p_fungi",""),
        "Plantae_raw_speciesKey":(RAW,""),
        "Plantae_rarefied_50record":(RARE,""),
        "Plantae_raw_speciesKey_plus_Animalia_Fungi_records":(RAW,"+"+ANI+"+"+FUN),
        "Plantae_rarefied_50record_plus_Animalia_Fungi_records":(RARE,"+"+ANI+"+"+FUN),
        "Plantae_GBIF_record_count_plus_Animalia_Fungi_records":(REC,"+"+ANI+"+"+FUN),
    }
    results={}
    for key,(y,extra) in names.items():
        results[key]=fit_one(df,y,extra)
    out["same_source_2015_WWp_countryFE_crosskingdom_recording_models"]=results
    out["country_original_WWp_correlations_with_raw_Animalia_Fungi_record_counts"]={
        "Plantae":countries(df,REC),
        "Animalia":countries(df,ANI),
        "Fungi":countries(df,FUN)
    }
    if len(df)==752 and df.country.nunique()==38:
        earlier=previous["all1055_2015WWp_and_2011_24_all_Plantae_countryFE_models"]
        for name,oldname in (("Plantae_record_count","all_Plantae_GBIF_record_count"),
                              ("Plantae_raw_speciesKey","raw_Plantae_speciesKey"),
                              ("Plantae_rarefied_50record","rarefied50_Plantae")):
            got=results[name]
            expected=earlier[oldname]
            if abs(got["beta_2015_modelled_WWp"]-expected["beta_on_log1p_modeled_WWp2015"])>1e-9 or any(abs(a-b)>1e-8 for a,b in zip(got["country_cluster95"],expected["country_cluster95"])):
                raise ValueError("Prior full1055 752-port all-Plantae model not reproduced before cross-kingdom effect")
        out["previous_752_flora_source_model_exactly_reproduced_before_animal_fungi_effect"]=True
    else:
        out["previous_752_flora_source_model_exactly_reproduced_before_animal_fungi_effect"]=None
        out["previous_floristic_model_not_same_sample_due_to_API_completeness"]={"new_n":len(df),"old_n":752}
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=run(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
       "gate":out["preoutcome_source_gate_original"],
       "allowed":out["crosskingdom_source_model_allowed"],
       "sample":out.get("source_case_support"),
       "models":out.get("same_source_2015_WWp_countryFE_crosskingdom_recording_models"),
       "country_signs":{k:{z:v[z] for z in ("source_eligible_countries_with_at_least8_ports",
             "country_source_spearman_positive_count","country_source_spearman_negative_count")}
           for k,v in out.get("country_original_WWp_correlations_with_raw_Animalia_Fungi_record_counts",{}).items()}
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
