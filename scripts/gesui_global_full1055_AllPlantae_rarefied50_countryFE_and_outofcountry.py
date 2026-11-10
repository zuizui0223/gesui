#!/usr/bin/env python3
"""FULL 1055 original global 42country WWp × independent all-Plantae response.

Only original source-fixed WPI2019 ports, original 2015 Jones WW generated
volume, original 2021 built/crop, new source-frozen 2011–24 GBIF all Plantae
speciesKey & exact 50-record rarefaction. No historic butterfly-host species
are substituted for all Plantae. Country FE, true unseen-country 2° folds
and sampling effort controls are observational non-causal posthoc measures.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import defaultdict,Counter
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.formula.api as smf

D=Path("results/global_all_Plantae_WWp2015_20261110")
FROZEN=D/"ORIGINAL_ALL1055_42COUNTRY_WWp2015_BEFORE_NEW808_GBIF_SOURCE.json"
GBIF=D/"ORIGINAL_ALL1055_42COUNTRY_GBIF_PLANTAE_2011_24_50RECORD_RAREFACTION_SOURCE.json"
PRED="log_2015_modelled_WWp"
COVARS="log_2015_WWTP_count+worldcover2021_built+worldcover2021_crop+abs_lat+sin_lon+cos_lon"
Y="log_rarefied50_all_Plantae_species_keys"
YEAR="2011,2024"
SEED=20261110
RESAMPLES=1500

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def load(root):
    source=json.loads((root/FROZEN).read_text())
    data=json.loads((root/GBIF).read_text())
    if source["schema"]!="gesui_original_1055_WPI_ports42country_808NEW_2011_24_ALLPlantae_PREOUTCOME_frame_v01":
        raise ValueError("Wrong full-source cohort")
    if data["original_all1055_WWp_WPI_country_source_frame_SHA256"]!=sha(root/FROZEN):
        raise ValueError("GBIF all-Plantae source receipt changed from original all1055 WPI sites")
    if data["original_full_42country_all1055_site_count"]!=1055:
        raise ValueError("All original 1055 WPI ports not included")
    if not data["precommitted_full_1055_AllPlantae_response_support_GATE"]["ALL1055_AllPlantae_WWp2015_macro_regression_ALLOWED"]:
        return source,data,None,"PRECOMMITTED_ALL1055_GBIF_PLANTAE_SAMPLE_QUALITY_GATE_FAIL"
    ports={r["original_port_id"]:r for r in source["source_all_1055_public_country_port_site_metadata_no_new_GBIF_response"]}
    result={r["original_port_id"]:r for r in data["original_1055_port_all_Plantae_source_GBIF_facet_summaries"]}
    if len(ports)!=1055 or len(result)!=1055 or set(ports)!=set(result):
        raise ValueError("Original 1055 source-WWTP and all-Plantae site joins inconsistent")
    rows=[]
    for portid in sorted(ports):
        p=ports[portid];g=result[portid]
        if p["original_source_country"]!=g["source_country"]:
            raise ValueError("Original 2015 model and new all Plantae port country mismatch")
        if not g["source_quality_min50_faceted_records_min10_species_coverage70pct"]:continue
        if g["GBIF_status"]!="COMPLETE":
            raise ValueError("Invalid all-Plantae GBIF sample marked original source-quality supported")
        s50=float(g["original_exact_hypergeometric_50_plant_record_rarefied_speciesKey"])
        if not (1<=s50<=50):
            raise ValueError("Original 50record rarefaction result outside physical range")
        lat=float(p["original_port_lat"]);lon=float(p["original_port_lon"])
        rows.append({
           "port_id":int(portid),"country":str(p["original_source_country"]),
           "source_pilot247":bool(p["previously_extracted_247_AllPlantae_pilot_source_port"]),
           "latitude":lat,"longitude":lon,
           PRED:float(p["source_WWp_log1p"]),
           "log_2015_WWTP_count":math.log1p(float(p["source_WWTP_count_25km"])),
           "worldcover2021_built":float(p["source_ESA2021_built_fraction"]),
           "worldcover2021_crop":float(p["source_ESA2021_crop_fraction"]),
           "abs_lat":abs(lat),"sin_lon":math.sin(math.radians(lon)),"cos_lon":math.cos(math.radians(lon)),
           Y:math.log1p(s50),
           "rarefied50_species":s50,
           "log1p_plant_species_keys":math.log1p(int(g["faceted_original_distinct_Plantae_speciesKeys"])),
           "log1p_total_plant_records":math.log1p(int(g["indexed_GBIF_all_Plantae_2011_2024_records"])),
           "GBIF_plant_species_key_coverage":float(g["facet_record_coverage"]),
           "original_plant_speciesKey_n":int(g["faceted_original_distinct_Plantae_speciesKeys"]),
           "original_plant_indexed_records_n":int(g["indexed_GBIF_all_Plantae_2011_2024_records"]),
        })
    df=pd.DataFrame(rows).sort_values(["country","port_id"])
    size=df.groupby("country").size()
    variation=df.groupby("country")[PRED].nunique()
    valid=size.index[(size>=2)&(variation>=2)]
    df=df[df.country.isin(valid)].copy()
    if len(df)<500 or df.country.nunique()<30:
        return source,data,None,"FULL_ALL1055_SOURCE_COUNTRY_VARIATION_MODEL_SAMPLE_INSUFFICIENT"
    return source,data,df,None

def regression(df):
    result={}
    for name,response,effort in [
        ("rarefied50_Plantae",Y,False),
        ("raw_Plantae_speciesKey","log1p_plant_species_keys",False),
        ("all_Plantae_GBIF_record_count","log1p_total_plant_records",False),
        ("rarefied50_plus_original_plant_record_effort",Y,True)]:
        formula=response+"~"+COVARS+"+"+PRED
        if effort:formula+="+log1p_total_plant_records"
        formula+="+C(country)"
        fit=smf.ols(formula,data=df).fit(cov_type="cluster",cov_kwds={"groups":df.country})
        result[name]={
           "formula":formula,
           "source_eligible_original_1055_port_count":int(len(df)),
           "source_eligible_original_countries":int(df.country.nunique()),
           "beta_on_log1p_modeled_WWp2015":float(fit.params[PRED]),
           "country_cluster95":list(map(float,fit.conf_int().loc[PRED])),
           "nominal_p":float(fit.pvalues[PRED]),
           "R2":float(fit.rsquared),
           "positive_slope":bool(fit.params[PRED]>0),
           "not_measured_sewage_discharge_or_plant_growth":True
        }
    return result

def rank_and_country_influence(df):
    bycountry=[]
    for country,d in df.groupby("country",sort=True):
        if len(d)>=6 and d[PRED].nunique()>=3 and d[Y].nunique()>=3:
            r=spearmanr(d[PRED],d[Y]).statistic
            bycountry.append({"country":str(country),"ports":int(len(d)),"WWp_v_rareS50_rho":float(r)})
    leave=[]
    for country in sorted(df.country.unique()):
        sub=df[df.country!=country]
        formula=Y+"~"+COVARS+"+"+PRED+"+C(country)"
        fit=smf.ols(formula,data=sub).fit()
        leave.append({"left_out_country":str(country),"remaining_port_count":int(len(sub)),
                      "beta_modelled_WWp2015":float(fit.params[PRED])})
    return{
      "n_with_6plus_ports_varied_WWp":len(bycountry),
      "positive_country_rho_n":sum(x["WWp_v_rareS50_rho"]>0 for x in bycountry),
      "negative_country_rho_n":sum(x["WWp_v_rareS50_rho"]<0 for x in bycountry),
      "within_country_rank_rows":bycountry,
      "original_countries_left_out_one_at_time":len(leave),
      "country_deleted_positive_countryFE_WWp_n":sum(x["beta_modelled_WWp2015"]>0 for x in leave),
      "country_deleted_WWp_coefficient_min":min(x["beta_modelled_WWp2015"] for x in leave),
      "country_deleted_WWp_coefficient_max":max(x["beta_modelled_WWp2015"] for x in leave),
      "country_LOCO_influence_not_unseen_country_prediction":True,
      "all_country_deleted_beta":leave
    }

def score(rows):
    # Each port has squared errors for source-identical nested predictors.
    if not rows:return None
    by=defaultdict(list)
    for r in rows:by[r["country"]].append(r)
    base=np.mean([r["baseline_error_sq"] for r in rows])
    aug=np.mean([r["plus2015WWp_error_sq"] for r in rows])
    countries=[]
    for c,rs in sorted(by.items()):
        x=np.mean([r["baseline_error_sq"] for r in rs])
        y=np.mean([r["plus2015WWp_error_sq"] for r in rs])
        countries.append({"country":c,"ports":len(rs),"baseline_MSE":float(x),
                          "plus_WWp_MSE":float(y),
                          "gain_percent":float(100*(x-y)/x) if x else None})
    ebase=np.mean([x["baseline_MSE"] for x in countries])
    eaug=np.mean([x["plus_WWp_MSE"] for x in countries])
    return{
       "n_ports":len(rows),"n_countries":len(countries),
       "port_pooled_MSE_baseline":float(base),
       "port_pooled_MSE_plusWWp":float(aug),
       "port_pooled_MSE_relative_gain_percent":float(100*(base-aug)/base),
       "country_equal_MSE_baseline":float(ebase),
       "country_equal_MSE_plusWWp":float(eaug),
       "country_equal_MSE_relative_gain_percent":float(100*(ebase-eaug)/ebase),
       "country_improved":sum(x["plus_WWp_MSE"]<x["baseline_MSE"] for x in countries),
       "country_worsened":sum(x["plus_WWp_MSE"]>x["baseline_MSE"] for x in countries),
       "per_country":countries,
    }

def geographical_2degree_within_country(df):
    d=df.copy()
    d["country_2degree"]=(d.country.astype(str)+"|"+
        np.floor(d.latitude/2).astype(int).astype(str)+"|"+
        np.floor(d.longitude/2).astype(int).astype(str))
    folds={}
    for country,g in d.groupby("country",sort=True):
        for i,block in enumerate(sorted(g.country_2degree.unique())):folds[block]=i%5
    d["fold"]=d.country_2degree.map(folds)
    raw=[]
    invalid=Counter()
    f0=Y+"~"+COVARS+"+C(country)"
    f1=Y+"~"+COVARS+"+"+PRED+"+C(country)"
    for i in range(5):
        train=d[d.fold!=i].copy()
        test=d[d.fold==i].copy()
        valid=set(train.groupby("country").size().loc[lambda s:s>=5].index)
        invalid["country_with_under_5_train_sites"]+=int((~test.country.isin(valid)).sum())
        test=test[test.country.isin(valid)].copy()
        if test.empty:continue
        model0=smf.ols(f0,data=train).fit()
        model1=smf.ols(f1,data=train).fit()
        e0=(test[Y].to_numpy()-model0.predict(test).to_numpy())**2
        e1=(test[Y].to_numpy()-model1.predict(test).to_numpy())**2
        for t,a,b in zip(test.itertuples(),e0,e1):
            raw.append({"port_id":int(t.port_id),"country":str(t.country),
                        "source_2degree_block":t.country_2degree,"fold":i,
                        "baseline_error_sq":float(a),"plus2015WWp_error_sq":float(b)})
    stats=score(raw)
    if stats is None or stats["n_ports"]<300 or stats["n_countries"]<20:
        return {"status":"INSUFFICIENT_COUNTRY_2DEGREE_SOURCE_TRAIN_TEST_SUPPORT",
                "original_2degree_blocks":len(folds),"tested_source_ports":len(raw),"model_results_not_admissible":True}
    return {
      "status":"ESTIMATED_SPATIALLY_OUTOFBLOCK_WITHIN_ALREADY_REPRESENTED_COUNTRIES",
      "original_2degree_WPI_country_blocks":len(folds),
      "source_excluded_test_cases":dict(invalid),
      "original_same_country_2degree_5fold_source_independent_mse":stats,
      "source_test_port_errors_original":raw,
      "country_FE_not_transfer_to_unseen_country":True}

def new_country_transfer(df):
    controls=Y+"~"+COVARS
    aug=controls+"+"+PRED
    out=[]
    for country in sorted(df.country.unique()):
        train=df[df.country!=country]
        test=df[df.country==country]
        m0=smf.ols(controls,data=train).fit()
        m1=smf.ols(aug,data=train).fit()
        a=(test[Y].to_numpy()-m0.predict(test).to_numpy())**2
        b=(test[Y].to_numpy()-m1.predict(test).to_numpy())**2
        for r,e0,e1 in zip(test.itertuples(),a,b):
            out.append({"country":str(country),"port_id":int(r.port_id),
                        "baseline_error_sq":float(e0),
                        "plus2015WWp_error_sq":float(e1)})
    if len(out)!=len(df):raise ValueError("Unseen-country source scoring did not include all original model sites")
    summ=score(out)
    rng=np.random.default_rng(SEED)
    vals=summ["per_country"]
    draws=[]
    for _ in range(RESAMPLES):
        chosen=rng.integers(0,len(vals),size=len(vals))
        v=[vals[i] for i in chosen]
        mse0=np.mean([r["baseline_MSE"] for r in v])
        mse1=np.mean([r["plus_WWp_MSE"] for r in v])
        draws.append(100*(mse0-mse1)/mse0)
    summ["country_equal_descriptive_bootstrap95_MSE_gain_percent"]=list(map(float,np.quantile(draws,[.025,.975])))
    summ["country_resample_probability_of_positive_modelled_WWp_forecast_gain"]=sum(x>0 for x in draws)/len(draws)
    return {"status":"STRICT_UNSEEN_COUNTRY_ALL_PLANTAE_NO_COUNTRY_FE",
            "full_source_42country_true_OOS_models":int(df.country.nunique()),
            "source_test_port_errors_original":out,
            "unseen_country_modelled_WWp_pilot_forecast_gain":summ,
            "new_country_no_countryFE_correct":True}

def subset_reproduction(df):
    res={}
    for name,flag in (("original247",True),("new808",False)):
        part=df[df.source_pilot247==flag].copy()
        eligible_country=part.groupby("country").agg(n=("port_id","size"),k=(PRED,"nunique"))
        good=set(eligible_country.index[(eligible_country.n>=2)&(eligible_country.k>=2)])
        d=part[part.country.isin(good)]
        record={"source_candidate_sites_qualified":len(part),
                "countries_with_plantae_records":int(part.country.nunique()),
                "country_2site_WWp_variation_model_sites":len(d),
                "country_2site_WWp_variation_model_countries":int(d.country.nunique())}
        if len(d)>=100 and d.country.nunique()>=20:
            m=smf.ols(Y+"~"+COVARS+"+"+PRED+"+C(country)",data=d).fit(cov_type="cluster",cov_kwds={"groups":d.country})
            record["rarefied50_modelled_WWp_beta"]=float(m.params[PRED])
            record["rarefied50_WWp_country_cluster95"]=list(map(float,m.conf_int().loc[PRED]))
        else:
            record["rarefied50_modelled_WWp_beta"]=None
            record["reason"]="INSUFFICIENT_SOURCE_SUPPORT_FOR_SAME_CONTROL_COHORT"
        res[name]=record
    return res

def run(root):
    source,gbif,d,reason=load(root)
    out={
        "schema":"gesui_global_all1055_42country_2015WWp_vs_2011_24_GBIF_ALL_Plantae_record_rarefied50_macro_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_1055_WPI_WWp2015_full_source_SHA256":sha(root/FROZEN),
        "original_247old_plus_808new_allPlantae_GBIF_1055_result_SHA256":sha(root/GBIF),
        "original_source_new_GBIF_country_and_species_quality_gate":gbif["precommitted_full_1055_AllPlantae_response_support_GATE"],
        "actual_original_country_2011_24_GBIF_speciesKey_QC":gbif["original_GBIF_allPlantae_query_and_coverage"],
        "ALL_Plantae_original_macro_model_eligible":reason is None,
        "original_macro_model_STOP_reason":reason,
        "all_Plantae_not_old_butterfly_candidate_host_plant_fraction":True,
        "all_2015_modelled_WWp_generated_not_actual_plant_nutrient_discharge":True,
        "all_2026_GBIF_observations_are_recording_effort_and_habitat_selected":True,
        "new_808_comparison_hypothesis_POSTHOC_after_247_pilot_results_known":True,
    }
    if reason is not None:return out
    out["original_all_Plantae_full1055_model_actual_sample"]={
      "n_original_port_records":len(d),
      "n_original_countries":int(d.country.nunique()),
      "old247_qualified_sites":int(d.source_pilot247.sum()),
      "new808_qualified_sites":int((~d.source_pilot247).sum()),
      "median_plant_speciesKey_richness":float(d.original_plant_speciesKey_n.median()),
      "median_50_record_rarefied_Plantae_species":float(d.rarefied50_species.median()),
      "median_GBIF_all_plant_record_count":float(d.original_plant_indexed_records_n.median()),
      "median_original_speciesKey_facet_coverage":float(d.GBIF_plant_species_key_coverage.median()),
      "country_original_port_counts":dict(Counter(d.country)),
    }
    out["all1055_2015WWp_and_2011_24_all_Plantae_countryFE_models"]=regression(d)
    out["country_rank_sign_and_country_LOCO_influence"]=rank_and_country_influence(d)
    out["original_pilot247_vs_new808_same_model_sensitivity"]=subset_reproduction(d)
    out["same_country_2degree_spatial_block_5fold_predictions"]=geographical_2degree_within_country(d)
    out["STRICT_unseen_country_without_country_fixed_effects"]=new_country_transfer(d)
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    result=run(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
      "source":result["actual_original_country_2011_24_GBIF_speciesKey_QC"],
      "sample":result.get("original_all_Plantae_full1055_model_actual_sample"),
      "models":result.get("all1055_2015WWp_and_2011_24_all_Plantae_countryFE_models"),
      "old_vs_new":result.get("original_pilot247_vs_new808_same_model_sensitivity"),
      "geo":result.get("same_country_2degree_spatial_block_5fold_predictions",{}).get("original_same_country_2degree_5fold_source_independent_mse"),
      "unseen":result.get("STRICT_unseen_country_without_country_fixed_effects",{}).get("unseen_country_modelled_WWp_pilot_forecast_gain"),
    },indent=2))
if __name__=="__main__":
    main()
