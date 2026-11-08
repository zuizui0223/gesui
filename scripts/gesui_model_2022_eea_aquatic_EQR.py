#!/usr/bin/env python3
"""Source-frozen 2022 EEA UWWTD measured N/P loads × aquatic biological EQR.

The source-only 8,245-record river station registry and five heldout 3-degree
regional folds were fixed before extracting any numerical biological response.
EEA official original Waterbase Biology 2026 EQR member is SHA256-verified.
One monitoring site/biological element/year is one observation; never treat
different organism classes or repeat samples as independent identical responses.
No causal effluent->biota claims.
"""
from __future__ import annotations
import argparse, hashlib, io, json, math, zipfile
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler,OneHotEncoder

SOURCE_REGISTRY_SHA="e12df511d64c4c580ef2422dccdbbd11a1309c84d29d467fd4527206e41fa636"
SOURCE_GIS_RUN=42432030617
BIO_EQR_SHA="f11001a57172b618bee95f58fd0ec51c8e71be7b9c01a70927dde21b7ca6e2a9"
EQR_MEMBER="eea_t_waterbase-biology-2026_p_1990-2024_v01_r00/Waterbase_v2025_1_T_WISE2_BiologyEQRData.csv"
YEAR="2022"
CODES={
  "EEA_13-01-4":"benthic_invertebrates",
  "EEA_124-04-9":"phytobenthos",
  "EEA_14-02-8":"fish"
}
MIN_UNIQUE_SITES=100
MIN_EXPOSED_SITES=20
MIN_COUNTRIES=5
MIN_EXPOSED_COUNTRIES=3
MIN_BLOCKS=10
RIDGE_ALPHA=10.0
FOLDS=5
SOURCE_N={
  "P":("upstream_2022_EEA_P_measured_t_per_yr",
       "downstream_2022_EEA_P_t_per_yr_placebo"),
  "N":("upstream_2022_EEA_N_measured_t_per_yr",
       "downstream_2022_EEA_N_t_per_yr_placebo")
}

def file_sha(path):
    h=hashlib.sha256()
    with path.open("rb") as fd:
        for x in iter(lambda:fd.read(1024*1024),b""):
            h.update(x)
    return h.hexdigest()

def get_outcome(arc):
    digest=hashlib.sha256()
    with zipfile.ZipFile(arc) as z:
        with z.open(EQR_MEMBER) as f:
            for b in iter(lambda:f.read(1024*1024),b""):
                digest.update(b)
        if digest.hexdigest()!=BIO_EQR_SHA:
            raise ValueError("Official 2026 Biology EQR original uncompressed SHA changed")
        # Need the entire 2022 official response table with source identifiers
        cols=[
          "countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
          "parameterWaterBodyCategory","observedPropertyDeterminandBiologyEQRCode",
          "procedureClassificationSystem","phenomenonTimeReferenceYear",
          "resultNormalisedEQRValue","resultEQRValue","resultObservationStatus",
          "resultEcologicalStatusClassValue","resultNumberOfSamples"]
        df=pd.read_csv(z.open(EQR_MEMBER),usecols=cols,dtype="string",low_memory=False)
    df=df[df.phenomenonTimeReferenceYear.eq(YEAR)&
          df.parameterWaterBodyCategory.eq("RW")&
          df.observedPropertyDeterminandBiologyEQRCode.isin(CODES)].copy()
    keys=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme"]
    for c in keys:df[c]=df[c].str.strip()
    df["site_key"]=df[keys].fillna("").agg("|".join,axis=1)
    if df.duplicated(["site_key","observedPropertyDeterminandBiologyEQRCode"]).any():
        raise ValueError("Repeated 2022 BQE site/outcome key; not safe to choose duplicate response")
    df["normalised_EQR"]=pd.to_numeric(df.resultNormalisedEQRValue,errors="coerce")
    df["raw_EQR_numeric"]=pd.to_numeric(df.resultEQRValue,errors="coerce")
    return df,digest.hexdigest()

def numeric_and_categories(df):
    frame=pd.DataFrame(index=df.index)
    for var in ("lat","lon"):
        frame[var]=pd.to_numeric(df[var],errors="coerce")
    for orig,out in (("DIS_AV_CMS","log_mean_river_Q"),("UPLAND_SKM","log_upstream_area")):
        x=pd.to_numeric(df[orig],errors="coerce")
        if (x<0).any():raise ValueError("Negative hydrological source variable")
        frame[out]=np.log1p(x)
    frame["country"]=df.countryCode.astype("string").fillna("UNKNOWN").astype(str)
    method=df.procedureClassificationSystem.astype("string").fillna("UNKNOWN").astype(str)
    frame["country_method"]=frame.country+"|"+method
    return frame

def single_model(df,source,name,upstream=False,downstream=False):
    x=numeric_and_categories(df)
    numerical=["lat","lon","log_mean_river_Q","log_upstream_area"]
    categorical=["country","country_method"]
    if upstream:
        x["log_reported_upstream_mass"]=np.log1p(
            pd.to_numeric(df[SOURCE_N[source][0]],errors="raise").clip(lower=0))
        numerical.append("log_reported_upstream_mass")
    if downstream:
        x["log_reported_downstream_mass"]=np.log1p(
            pd.to_numeric(df[SOURCE_N[source][1]],errors="raise").clip(lower=0))
        numerical.append("log_reported_downstream_mass")
    if not np.isfinite(x[numerical].to_numpy(float)).all():
        raise ValueError("Nonfinite source numeric predictors; no model fit")
    y=df.normalised_EQR.to_numpy(float)
    folds=pd.to_numeric(df.fixed_spatial_cv_fold,errors="raise").astype(int).to_numpy()
    if set(folds)!={1,2,3,4,5}:raise ValueError("Source pre-assigned spatial fold 1-5 missing")
    pred=np.full(len(df),np.nan)
    fold_metrics=[]
    for i in range(1,FOLDS+1):
        test=np.flatnonzero(folds==i)
        train=np.flatnonzero(folds!=i)
        if len(test)==0 or len(train)==0:raise ValueError("Empty source-only test/train fold")
        transformer=ColumnTransformer([
            ("numeric",StandardScaler(),numerical),
            ("category",OneHotEncoder(handle_unknown="ignore"),categorical)
        ],remainder="drop")
        m=Pipeline([("features",transformer),("ridge",Ridge(alpha=RIDGE_ALPHA))])
        m.fit(x.iloc[train],y[train])
        p=m.predict(x.iloc[test])
        pred[test]=p
        fold_metrics.append({"fold":i,
           "test_site_n":len(test),
           "distinct_3degree_blocks":int(df.iloc[test].region_3deg_block.nunique()),
           "heldout_MSE":float(np.mean((y[test]-p)**2)),
           "heldout_MAE":float(np.mean(np.abs(y[test]-p)))})
    if not np.isfinite(pred).all():raise ValueError("Unpredicted heldout stations")
    err=y-pred
    return {
        "model":name,"n":len(df),
        "overall_heldout_MSE":float(np.mean(err**2)),
        "overall_heldout_MAE":float(np.mean(np.abs(err))),
        "original_source_fold_assignments_used":True,
        "fold_metrics":fold_metrics
    },err

def cv(df,source):
    numeric_up=pd.to_numeric(df[SOURCE_N[source][0]],errors="raise")
    numeric_down=pd.to_numeric(df[SOURCE_N[source][1]],errors="raise")
    if (numeric_up<0).any() or (numeric_down<0).any():
        raise ValueError("Reported measured source loads must not be negative")
    spec={
      "M0_country_method_hydrogeo":(False,False),
      "M1_plus_upstream_measured_load":(True,False),
      "M2_plus_downstream_negative_control":(False,True),
      "M3_both_upstream_and_downstream":(True,True)
    }
    metrics={}
    residuals={}
    for name,(up,down) in spec.items():
        out,err=single_model(df,source,name,up,down)
        metrics[name]=out
        residuals[name]=err
    baseline=metrics["M0_country_method_hydrogeo"]["overall_heldout_MSE"]
    if baseline==0:raise ValueError("Zero baseline MSE, cannot compute change")
    summary={}
    for name,entry in metrics.items():
        if name== "M0_country_method_hydrogeo":continue
        mse=entry["overall_heldout_MSE"]
        groups=df.region_3deg_block.astype(str).to_numpy()
        block_delta={}
        for block in sorted(set(groups)):
            mask=groups==block
            block_delta[block]=float(np.mean(
                residuals["M0_country_method_hydrogeo"][mask]**2-residuals[name][mask]**2))
        summary[name]={
          "MSE_improvement_percent":float(100*(baseline-mse)/baseline),
          "MSE_absolute_reduction":float(baseline-mse),
          "blocks_positive_improvement":sum(x>0 for x in block_delta.values()),
          "blocks_total":len(block_delta)
        }
    return {"actual_measurement_parameter":source,
            "source_upstream_nonzero_site_n":int((numeric_up>0).sum()),
            "source_downstream_nonzero_site_n":int((numeric_down>0).sum()),
            "CV_MSE_by_model":metrics,
            "added_source_predictor_heldout_performance_change":summary}

def study_one(source_registry,original,code,label):
    df=source_registry[source_registry.observedPropertyDeterminandBiologyEQRCode.eq(code)].copy()
    raw=original[original.observedPropertyDeterminandBiologyEQRCode.eq(code)].copy()
    source_keys=set(df["site_key"])
    original_keys=set(raw["site_key"])
    if not source_keys<=original_keys:
        raise ValueError("Some source selected 2022 EQR sites missing from official outcome table")
    join=df.merge(raw[["site_key","observedPropertyDeterminandBiologyEQRCode",
                 "normalised_EQR","raw_EQR_numeric","resultEcologicalStatusClassValue",
                 "resultNumberOfSamples","resultObservationStatus"]],
            on=["site_key","observedPropertyDeterminandBiologyEQRCode"],
            how="inner",validate="one_to_one",suffixes=("_source","_observed"))
    count_before=len(join)
    measured=join.normalised_EQR
    numerical_valid=np.isfinite(measured) & measured.ge(0)
    q=join[numerical_valid].copy().reset_index(drop=True)
    number_bad=int((~numerical_valid).sum())
    country_n=int(q.countryCode.nunique())
    positive=q[q.upstream_2022_unique_EEA_both_load_WWTP_count.gt(0)]
    positive_countries=int(positive.countryCode.nunique())
    blocks=int(q.region_3deg_block.nunique())
    counts=dict(
      source_sites_before_numeric_quality=count_before,
      numeric_original_EQR_sites=len(q),
      missing_or_invalid_original_EQR_sites=number_bad,
      original_normalized_EQR_above_1=int(q.normalised_EQR.gt(1).sum()),
      original_normalized_EQR_exact_zero=int(q.normalised_EQR.eq(0).sum()),
      original_numeric_EQR_min=float(q.normalised_EQR.min()) if len(q) else None,
      original_numeric_EQR_median=float(q.normalised_EQR.median()) if len(q) else None,
      original_numeric_EQR_max=float(q.normalised_EQR.max()) if len(q) else None,
      country_count=country_n,
      EEA_2022_upstream_source_positive_sites=len(positive),
      EEA_2022_upstream_source_positive_countries=positive_countries,
      spatial_3degree_blocks=blocks,
      EEA_2022_source_zero_registered_sites=int(
        q.upstream_2022_unique_EEA_both_load_WWTP_count.eq(0).sum()),
      distinct_original_waterbody_IDs=int(q.waterBodyIdentifier.nunique()),
      original_classification_methods=int(q.procedureClassificationSystem.nunique()),
      source_EQR_observation_status_counts={
        str(k):int(v) for k,v in
        q.resultObservationStatus_observed.fillna("MISSING").value_counts().items()}
        if "resultObservationStatus_observed" in q else {},
      source_reported_country_counts={str(k):int(v) for k,v in q.countryCode.value_counts().items()}
    )
    supported=(len(q)>=MIN_UNIQUE_SITES and len(positive)>=MIN_EXPOSED_SITES
              and country_n>=MIN_COUNTRIES and positive_countries>=MIN_EXPOSED_COUNTRIES
              and blocks>=MIN_BLOCKS)
    counts["predeclared_geographic_and_source_support_pass"]=supported
    results={
      "source_code":code,"biological_group":label,
      "support":counts,
      "biology_meanings":"Original EEA normalized EQR is biological status index, NOT species abundance",
      "underlying_country_and_classification_procedure_controlled":True,
      "actual_source_water_chemical_NP_concentration_not_treatment_effect":True,
      "source_2022_reported_measured_WWTP_loads_tonnes_per_year":True,
      "outcome_year":2022
    }
    if supported:
        r={}
        for nutrient in ("P","N"):
            up=pd.to_numeric(q[SOURCE_N[nutrient][0]],errors="raise")
            rho,p=spearmanr(np.log1p(up),q.normalised_EQR.to_numpy(float))
            r[nutrient]={
              "source_EQR_spearman_rho_descriptive":float(rho),
              "source_EQR_spearman_asymptotic_p_exploratory":float(p),
              "fixed_geographic_block_predictive_comparisons":cv(q,nutrient)}
        results["numerical_prediction"]=r
        results["status"]="ORIGINAL_2022_SOURCE_AND_BIOLOGY_VALID__OBSERVATIONAL_ONLY"
    else:
        results["status"]="NO_POOLED_ECOLOGICAL_TEST__PREDECLARED_SOURCE_OR_GEOGRAPHIC_SUPPORT_FAIL"
    print("GESUI_ACTUAL_EEA2022_BIOLOGICAL_QUALITY_RESULT",
          json.dumps({k:v for k,v in results.items() if k!="numerical_prediction"},sort_keys=True),flush=True)
    if supported:
        for nutrient,r in results["numerical_prediction"].items():
            small=r["fixed_geographic_block_predictive_comparisons"]
            print("GESUI_EEA2022_BQE_HELDOUT_ECOLOGICAL_PREDICTION",
                json.dumps({"BQE":label,"nutrient":nutrient,
                  "sample_sites":len(q),
                  "upstream_source_positive_sites":small["source_upstream_nonzero_site_n"],
                  "observed_BQE_raw_rho":r["source_EQR_spearman_rho_descriptive"],
                  "CV_baseline_MSE":small["CV_MSE_by_model"]["M0_country_method_hydrogeo"]["overall_heldout_MSE"],
                  "CV_upstream_MSE":small["CV_MSE_by_model"]["M1_plus_upstream_measured_load"]["overall_heldout_MSE"],
                  "CV_downstream_placebo_MSE":small["CV_MSE_by_model"]["M2_plus_downstream_negative_control"]["overall_heldout_MSE"],
                  "MSE_improvement_pct":small["added_source_predictor_heldout_performance_change"]["M1_plus_upstream_measured_load"]["MSE_improvement_percent"],
                  "downstream_placebo_improvement_pct":small["added_source_predictor_heldout_performance_change"]["M2_plus_downstream_negative_control"]["MSE_improvement_percent"]},sort_keys=True),flush=True)
    return results, q

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-registry",required=True,type=Path)
    p.add_argument("--original-biology-zip",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    if file_sha(args.source_registry)!=SOURCE_REGISTRY_SHA:
        raise ValueError("Original pre-EQR outcome source-only site registry SHA changed")
    df=pd.read_csv(args.source_registry,dtype={
       "site_key":str,"monitoringSiteIdentifier":str,
       "monitoringSiteIdentifierScheme":str,"countryCode":str,
       "observedPropertyDeterminandBiologyEQRCode":str,
       "procedureClassificationSystem":str,
       "spatial_3deg_block":str},low_memory=False)
    if len(df)!=8245 or df.duplicated(["site_key","observedPropertyDeterminandBiologyEQRCode"]).any():
        raise ValueError("Source-only station or biological quality element selection changed")
    if set(df.observedPropertyDeterminandBiologyEQRCode)!=set(CODES):
        raise ValueError("Source-only fixed BQE classes changed")
    source_registry=SOURCE_GIS_RUN
    original,eqr_hash=get_outcome(args.original_biology_zip)
    results={}
    source_original_rows={}
    args.out.mkdir(parents=True,exist_ok=True)
    for code,label in CODES.items():
        r,q=study_one(df,original,code,label)
        results[label]=r
        if r["support"]["predeclared_geographic_and_source_support_pass"]:
            # Station-level file retains original EQR values, source loads, and
            # original deterministic source folds. Raw EEA all-years records
            # are never copied into research artifact.
            keep=["site_key","countryCode","waterBodyIdentifier",
                  "observedPropertyDeterminandBiologyEQRCode",
                  "procedureClassificationSystem","lat","lon",
                  "HYRIV_ID","DIS_AV_CMS","UPLAND_SKM","region_3deg_block",
                  "fixed_spatial_cv_fold","normalised_EQR","resultObservationStatus_observed",
                  "upstream_2022_EEA_P_measured_t_per_yr",
                  "upstream_2022_EEA_N_measured_t_per_yr",
                  "downstream_2022_EEA_P_t_per_yr_placebo",
                  "downstream_2022_EEA_N_t_per_yr_placebo"]
            keep=[x for x in keep if x in q]
            q[keep].to_csv(args.out/(label+"_EEA2022_EQR_source_matched_monitoring.csv"),
                            index=False,float_format="%.12g")
    summary={
      "schema":"gesui_2022_official_EEA_Biology_EQR_reported_NP_wastewater_model_v01",
      "source_only_original_run":SOURCE_GIS_RUN,
      "source_only_csv_SHA256":SOURCE_REGISTRY_SHA,
      "official_2026_Biology_original_EQR_member_SHA256":eqr_hash,
      "study_reference_year":2022,
      "chemical_source":"EEA UWWTD reported measured annual nutrient discharge load from uniquely matched direct FW outfall",
      "biological_source":"EEA WISE2 normalized EQR original 2022 monitoring-site record",
      "source_GBIF_plant_species_results_used_in_BQE_selection":False,
      "source_2022_full_EQR_all_values_not_redistributed":True,
      "reported_annual_source_P_and_N_are_not_measured_EQR_causal_treatments":True,
      "predeclared_site_gate":{
        "min_sites":MIN_UNIQUE_SITES,"min_upstream_source_sites":MIN_EXPOSED_SITES,
        "min_countries":MIN_COUNTRIES,"min_upstream_source_countries":MIN_EXPOSED_COUNTRIES,
        "min_3degree_blocks":MIN_BLOCKS},
      "models_and_scale":{
        "ridge_alpha":RIDGE_ALPHA,
        "fixed_source_only_5_geographic_folds":True,
        "baseline_numeric":["latitude","longitude","river_mean_discharge","upstream_drainage_area"],
        "baseline_categorical":["countryCode","countryCode×classificationProcedure"],
        "negative_control":"WWTP source loads physically downstream of biological monitoring site"},
      "results":results,
      "limitations":[
        "Original EEA normalized EQR is national method-calibrated biological quality, not species counts",
        "Source loads are 2022 annually measured WWTP discharge mass, not independently attributed EQR treatment",
        "HydroRIVERS third-party EU topology not independently validated against official original",
        "Source-zero has no eligible registered measured WWTP load within 50km, not no wastewater",
        "Spatial residual urbanization, agricultural runoff, channel modification and EQR monitoring selection may confound",
        "Site-level data/river categories are only EEA official RW and unique source-matched site IDs"
      ]}
    (args.out/"SOURCE_LOCKED_2022_AQUATIC_EQR_RESULTS.json").write_text(
      json.dumps(summary,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_EEA2022_BIOLOGICAL_ALL_GROUP_SUMMARY",
          json.dumps({k:{"n":v["support"]["numeric_original_EQR_sites"],
                         "exposed":v["support"]["EEA_2022_upstream_source_positive_sites"],
                         "countries":v["support"]["country_count"],
                         "status":v["status"]}
                for k,v in results.items()},sort_keys=True),flush=True)

if __name__=="__main__":
    main()
