#!/usr/bin/env python3
"""Analyze NEW GBIF ALL Plantae species richness after source-locked source gate.

This is independent of old butterfly-host candidate record-type numerator.
Globally selected 247 source sites across original 42 WWp countries remain
unchanged. Country clustered associations are observational posthoc estimates.
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

DIR=Path("results/global_all_Plantae_WWp2015_20261110")
SOURCE=DIR/"ORIGINAL_42COUNTRY_WWp_STRATIFIED_PORTS_PRE_GBIF_ALLFLORA_SOURCE.json"
PLANT=DIR/"ORIGINAL_247PORT_42COUNTRY_GBIF_ALL_PLANTAE_2011_24_RAREFACTION_SOURCE.json"
MAIN_RESPONSE="log1p_50rarefied_allPlantae"
PREDICTOR="log1p_WWp2015"
BASE="log1p_WWTP25+original_built_frac+original_crop_frac+abs_lat+sin_lon+cos_lon"
N_BOOT=2000
SEED=20261110

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def sample(root):
    old=json.loads((root/SOURCE).read_text())
    fresh=json.loads((root/PLANT).read_text())
    if not old["source_selection_is_independent_of_all_new_GBIF_Plantae_records"]:
        raise ValueError("Original 42country source sample was not outcome blind")
    if fresh["source_only_study_manifest_sha256"]!=sha(root/SOURCE):
        raise ValueError("New all-Plantae result not on the frozen source set")
    if old["new_all_Plantae_GBIF_pilot_site_n"]!=247 or old["countries_in_source_pilot"]!=42:
        raise ValueError("Old source country-port count changed")
    if fresh["original_GBIF_result_completion_summary"]["fixed_ports"]!=247:
        raise ValueError("New ALL Plantae source count not all 247 port requests")
    indexed={x["original_port_id"]:x for x in fresh["original_site_GBIF_facet_summary"]}
    psource={x["original_port_id"]:x for x in old["all_pilot_port_georeferenced_original_sites"]}
    if len(indexed)!=247 or len(psource)!=247 or set(indexed)!=set(psource):
        raise ValueError("GBIF all-Plantae sampled a changed or duplicated source port list")
    if not fresh["original_source_coverage_gate"]["ALL_PLANTAE_WWp_statistical_comparison_ALLOWED"]:
        return old,fresh,None,"ORIGINAL_ALL_PLANTAE_GBIF_RECORD_SUPPORT_GATE_FAILED"
    data=[]
    for sid in sorted(psource):
        p=psource[sid]
        z=indexed[sid]
        if p["original_source_country"]!=z["source_country"]:
            raise ValueError("Port source country and GBIF source country mismatch")
        if not z["source_quality_min50_faceted_records_min10_species_coverage70pct"]:
            continue
        if z["GBIF_status"]!="COMPLETE" or z["faceted_original_speciesKey_records"]<50:
            raise ValueError("Original GBIF opportunity gate marked an invalid sample as complete")
        if z["faceted_original_distinct_Plantae_speciesKeys"]<10 or z["facet_record_coverage"]<.7:
            raise ValueError("Accepted original all-Plantae species richness violates its gate")
        val=z["original_exact_hypergeometric_50_plant_record_rarefied_speciesKey"]
        if val is None or val<1 or val>50:
            raise ValueError("Rarefied 50 plant-record species richness outside bounds")
        lat=p["original_port_lat"];lon=p["original_port_lon"]
        data.append({
            "port_id":int(sid),
            "country":p["original_source_country"],
            "lat":float(lat),
            "lon":float(lon),
            PREDICTOR:float(p["source_WWp_log1p"]),
            "log1p_WWTP25":math.log1p(p["source_WWTP_count_25km"]),
            "original_built_frac":float(p["source_ESA2021_built_fraction"]),
            "original_crop_frac":float(p["source_ESA2021_crop_fraction"]),
            "abs_lat":abs(float(lat)),
            "sin_lon":math.sin(math.radians(float(lon))),
            "cos_lon":math.cos(math.radians(float(lon))),
            MAIN_RESPONSE:math.log1p(val),
            "rarefied_50_Plantae_speciesKey":float(val),
            "log1p_raw_Plantae_speciesKey":math.log1p(z["faceted_original_distinct_Plantae_speciesKeys"]),
            "log1p_Plantae_indexed_records":math.log1p(z["indexed_GBIF_all_Plantae_2011_2024_records"]),
            "GBIF_Plantae_facet_coverage":float(z["facet_record_coverage"]),
            "orig_GBIF_all_Plantae_occurrences":int(z["indexed_GBIF_all_Plantae_2011_2024_records"]),
            "orig_GBIF_species_key_count":int(z["faceted_original_distinct_Plantae_speciesKeys"]),
        })
    df=pd.DataFrame(data).sort_values(["country","port_id"])
    if len(df)<120 or df.country.nunique()<25:
        raise ValueError("Original precommit rarefied plant source pass but fewer than 120 valid ports /25 countries")
    valid=Counter(df.country)
    var=df.groupby("country")[PREDICTOR].nunique()
    good=set(c for c,n in valid.items() if n>=2 and var[c]>=2)
    if len(good)<20:
        return old,fresh,None,"WITHIN_COUNTRY_WWp_VARIATION_AFTER_GBIF_REPORTING_FAIL"
    d=df[df.country.isin(good)].copy()
    return old,fresh,d,None

def fit_models(df):
    controls=BASE+"+C(country)"
    result={}
    for outcome,predictors in (
        (MAIN_RESPONSE,[PREDICTOR]),
        ("log1p_raw_Plantae_speciesKey",[PREDICTOR]),
        ("log1p_Plantae_indexed_records",[PREDICTOR]),
        (MAIN_RESPONSE,[PREDICTOR,"log1p_Plantae_indexed_records"])):
        label=outcome+("_plus_reporting_records_control" if len(predictors)>1 else "_main")
        formula=outcome+"~"+controls+"+"+"+".join(predictors)
        model=smf.ols(formula,data=df).fit(cov_type="cluster",cov_kwds={"groups":df.country})
        x=model.params[PREDICTOR]
        ci=list(map(float,model.conf_int().loc[PREDICTOR]))
        result[label]={
            "formula":formula,
            "original_new_ALL_Plantae_ports_n":len(df),
            "original_new_ALL_Plantae_countries_n":int(df.country.nunique()),
            "beta_log_WWp":float(x),
            "country_cluster95":ci,
            "nominal_p":float(model.pvalues[PREDICTOR]),
            "R2":float(model.rsquared),
            "positive_coefficient":bool(x>0),
            "is_precisely_all_Plantae_not_original_butterfly_candidate_hosts":True,
            "country_FE_not_unseen_country_forecast_model":True,
        }
    return result

def country_rank(df):
    rows=[]
    for country,g in df.groupby("country",sort=True):
        if len(g)>=4 and g[PREDICTOR].nunique()>=3 and g[MAIN_RESPONSE].nunique()>=3:
            rho=spearmanr(g[PREDICTOR],g[MAIN_RESPONSE]).statistic
            rows.append({"country":country,"ports":len(g),"rho":float(rho)})
    return {"eligible_countries_n":len(rows),
            "positive_rho_n":sum(v["rho"]>0 for v in rows),
            "negative_rho_n":sum(v["rho"]<0 for v in rows),
            "exact_zero_rho_n":sum(v["rho"]==0 for v in rows),
            "all_reported_supported_country_ranks":rows,
            "country_n_under6_per_country_small_and_weak_power":True}

def country_LOCO(df):
    result=[]
    for left in sorted(df.country.unique()):
        train=df[df.country!=left]
        if train.country.nunique()<18:raise ValueError("Country source support no longer macro")
        model=smf.ols(MAIN_RESPONSE+"~"+BASE+"+"+PREDICTOR+"+C(country)",
                      data=train).fit()
        result.append({"left_out_country":str(left),
                       "n_train_ports":len(train),
                       "n_train_countries":int(train.country.nunique()),
                       "reestimated_new_all_Plantae_rarefied_WWp_coefficient":float(model.params[PREDICTOR])})
    return {"source_countries_left_out":len(result),
            "positive_count":sum(x["reestimated_new_all_Plantae_rarefied_WWp_coefficient"]>0 for x in result),
            "negative_count":sum(x["reestimated_new_all_Plantae_rarefied_WWp_coefficient"]<0 for x in result),
            "min_coefficient":min(x["reestimated_new_all_Plantae_rarefied_WWp_coefficient"] for x in result),
            "max_coefficient":max(x["reestimated_new_all_Plantae_rarefied_WWp_coefficient"] for x in result),
            "all_country_influence_rows":result,
            "not_prediction_to_new_country":True}

def report(root):
    old,fresh,df,why=sample(root)
    output={
        "schema":"gesui_global_2011_2024_all_Plantae_original_247port_42country_rarefied50_WWp_macro_pilot_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_source_port_WWp_selection_SHA256":sha(root/SOURCE),
        "independent_GBIF_allPlantae_2011_24_247port_facet_summary_SHA256":sha(root/PLANT),
        "original_source_247_ports_42_countries_selected_before_new_GBIF":True,
        "GBIF_all_Plantae_taxonKey6_NOT_only_butterfly_candidate_hosts":True,
        "predeclared_all_Plantae_record_support_gate":fresh["original_source_coverage_gate"],
        "source_original_247port_42country_GBIF_allPlantae_response_QC":fresh["original_GBIF_result_completion_summary"],
        "model_possible":why is None,
        "model_stop_reason":why,
        "2015_WWp_generated_model_not_effluent_or_NP_mass":True,
        "recorded_2011_24_Plantae_species_keys_not_true_biodiversity_plant_abundance":True,
        "all_new_regressions_POSTHOC_after_original_candidate_host_positive":True,
        "selected_countries_from_old_host_eligibility_not_globally_random":True,
    }
    if why is not None:
        output["macro_country_models"]=None
        return output
    assert df is not None
    m=fit_models(df)
    output["source_per_country_WWp_variation_model_sample_n"]={
        "countries":int(df.country.nunique()),"port_sites":len(df),
        "source_model_input_country_sizes":dict(Counter(df.country)),
        "min_record_coverage":float(df.GBIF_Plantae_facet_coverage.min()),
        "median_record_coverage":float(df.GBIF_Plantae_facet_coverage.median()),
        "median_rarefied50_speciesKey":float(df.rarefied_50_Plantae_speciesKey.median()),
        "median_original_speciesKey":float(df.orig_GBIF_species_key_count.median()),
        "median_record_count":float(df.orig_GBIF_all_Plantae_occurrences.median()),
        "countries_with_distinct_WWp_sample_requirement":True,
    }
    output["macro_country_models"]=m
    output["new_all_Plantae_country_2015WWp_rank_direction"]=country_rank(df)
    output["leave_one_original_country_out_influence"]=country_LOCO(df)
    output["sampling_effort_and_species_identity_limits"]=[
        "GBIF taxonKey=6 includes all currently indexed Plantae, not exactly only vascular floristic taxa.",
        "50-record hypergeometric rarefaction standardizes records but not observer site selection or GBIF data completeness.",
        "GBIF 2011-24 occurrence records and 2015 modeled WWp production can covary with country, urbanization, plant indexing and uneven coastal reporting.",
        "Even positive rarefied richness would not show nutrient enrichment, terrestrial disposal, establishment or real plant-population biomass.",
        "Existing 1055 original candidate-host 2015 WWp positive beta is a DIFFERENT response, not a replication of this independent all-Plantae test."
    ]
    return output

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    obj=report(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(obj,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"gate":obj["model_possible"],
                      "source_QC":obj["source_original_247port_42country_GBIF_allPlantae_response_QC"],
                      "model_cohort":obj.get("source_per_country_WWp_variation_model_sample_n"),
                      "models":obj["macro_country_models"],
                      "country_signs":obj.get("new_all_Plantae_country_2015WWp_rank_direction"),
                     },indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
