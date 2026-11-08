#!/usr/bin/env python3
"""Posthoc independent mechanism falsification: reported WWTP nutrient amount vs facility presence/count.

No new model selection for the original paper claim. Original source-only
2022 reported measured N/P WWTP counts, annual loads and biological station
EQR outcomes are immutable. All sites and geographic fold assignments
unchanged; compare (a) WWTP-present binary, (b) upstream WWTP count, (c) actual
reported N/P mass, (d) mass+presence/count in the ORIGINAL common-station CV.

This diagnostic does NOT separate N vs P (highly collinear) and does NOT
prove causal wastewater nutrient effects, even if annual source mass predicts
better than binary registration.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd,sklearn
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge

SRC_HASH="e12df511d64c4c580ef2422dccdbbd11a1309c84d29d467fd4527206e41fa636"
ORIG_HASH={
  "benthic_invertebrates":"5e1f3c31b389d0bd769984ab055d3e2c54b63ccaa4f5eedfb2d82bcdc83c28e6",
  "phytobenthos":"7fa5dfff9023337ac0a62a87fe406eb71783979858d9d924ab97019c6e4e5fc3"}
FILES={"benthic_invertebrates":"BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv",
       "phytobenthos":"PHYTOBENTHOS_2022_SOURCE_EQR.csv"}
SOURCE_AMOUNT={"P":"upstream_2022_EEA_P_measured_t_per_yr",
               "N":"upstream_2022_EEA_N_measured_t_per_yr"}
COUNT="upstream_2022_unique_EEA_both_load_WWTP_count"
SPEC={
  "baseline":[],
  "reported_eligible_WWTP_presence":["source_presence"],
  "registered_upstream_WWTP_count":["source_count"],
  "reported_measured_source_load":["source_mass"],
  "source_load_beyond_presence":["source_presence","source_mass"],
  "source_load_beyond_WWTP_count":["source_count","source_mass"]
}
EXPECTED_N={"benthic_invertebrates":2475,"phytobenthos":2784}
EXPECTED_MSE={"benthic_invertebrates":{"P":(0.05833804346777809,0.05769459552059828),
                                      "N":(0.05833804346777809,0.05752327324679533)},
              "phytobenthos":{"P":(0.035934357335444235,0.035875931550570526),
                              "N":(0.035934357335444235,0.03587108895518606)}}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def source_join(outcome,src):
    x=outcome.merge(src[[
        "site_key","observedPropertyDeterminandBiologyEQRCode",COUNT
       ]],on=["site_key","observedPropertyDeterminandBiologyEQRCode"],
       validate="one_to_one",how="inner")
    if len(x)!=len(outcome):raise ValueError("Source-only facility counts missing from original biology panel")
    assert not x.site_key.duplicated().any()
    assert x[COUNT].ge(0).all()
    if not np.array_equal((x[COUNT]>0).to_numpy(),(x[SOURCE_AMOUNT["P"]]>0).to_numpy()):
        raise ValueError("Source count presence disagrees with source reported P mass")
    return x

def matrix(df,nutrient,optional):
    x=pd.DataFrame(index=df.index)
    x["lat"]=df.lat.to_numpy(float)
    x["lon"]=df.lon.to_numpy(float)
    x["logQ"]=np.log1p(df.DIS_AV_CMS.to_numpy(float))
    x["logArea"]=np.log1p(df.UPLAND_SKM.to_numpy(float))
    nums=["lat","lon","logQ","logArea"]
    if "source_presence" in optional:
        x["source_presence"]=(df[COUNT]>0).astype(float)
        nums.append("source_presence")
    if "source_count" in optional:
        x["source_count"]=np.log1p(df[COUNT].astype(float))
        nums.append("source_count")
    if "source_mass" in optional:
        x["source_mass"]=np.log1p(df[SOURCE_AMOUNT[nutrient]].astype(float))
        nums.append("source_mass")
    x["country"]=df.countryCode.astype("string").fillna("UNKNOWN").astype(str)
    x["country_method"]=x.country+"|"+df.procedureClassificationSystem.astype("string").fillna("UNKNOWN").astype(str)
    if not np.isfinite(x[nums].to_numpy(float)).all():
        raise ValueError("Nonfinite model source predictors")
    return x,nums

def heldout(df,nutrient,optional):
    x,nums=matrix(df,nutrient,optional)
    y=df.normalised_EQR.to_numpy(float)
    fold=df.fixed_spatial_cv_fold.to_numpy(dtype=int)
    assert set(fold)=={1,2,3,4,5}
    pred=np.full(len(df),np.nan)
    for k in range(1,6):
        tr=np.flatnonzero(fold!=k)
        te=np.flatnonzero(fold==k)
        pipe=Pipeline([
          ("preprocess",ColumnTransformer([
            ("numeric",StandardScaler(),nums),
            ("categorical",OneHotEncoder(handle_unknown="ignore"),["country","country_method"])
          ])),
          ("model",Ridge(alpha=10.0))])
        pipe.fit(x.iloc[tr],y[tr]);pred[te]=pipe.predict(x.iloc[te])
    if not np.isfinite(pred).all():raise ValueError("Broken model heldout predictions")
    return float(np.mean((y-pred)**2))

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--results-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    srcpath=a.results_dir/"SOURCE_ONLY_BQE_SITE_REGISTRY.csv"
    if sha(srcpath)!=SRC_HASH:raise ValueError("Frozen source registry SHA differs")
    source=pd.read_csv(srcpath,low_memory=False)
    result={}
    for group,filename in FILES.items():
        fp=a.results_dir/filename
        if sha(fp)!=ORIG_HASH[group]:raise ValueError("Original EQR response sample has been modified")
        old=pd.read_csv(fp,low_memory=False)
        if len(old)!=EXPECTED_N[group]:raise ValueError("Source BQE sample count drift")
        df=source_join(old,source)
        models={}
        for nutrient in ("P","N"):
            results={name:heldout(df,nutrient,extras) for name,extras in SPEC.items()}
            expected_baseline,expected_annual_mass=EXPECTED_MSE[group][nutrient]
            if (abs(results["baseline"]-expected_baseline)>4e-6 or
                abs(results["reported_measured_source_load"]-expected_annual_mass)>4e-6):
                raise ValueError("Preregistered model does not replay in posthoc source count check")
            base=results["baseline"]
            pairs={}
            for name,mse in results.items():
                pairs[name]={"heldout_MSE":mse,
                     "percent_improvement_from_baseline":float(100*(base-mse)/base)}
            pairs["incremental_measured_load_beyond_registered_presence"]={
              "absolute_MSE_reduction":results["reported_eligible_WWTP_presence"]-
                                         results["source_load_beyond_presence"],
              "MSE_reduction_pct_relative_presence_only":
                100*(results["reported_eligible_WWTP_presence"]-
                       results["source_load_beyond_presence"])/
                       results["reported_eligible_WWTP_presence"]}
            models[nutrient]={
              "sites":len(df),"original_exposed_sites":int((df[COUNT]>0).sum()),
              "original_BQE_fold_count":5,
              "original_WWTP_facility_presence_vs_count_vs_source_load":pairs}
            print("GESUI_EEA_BQE_REPORTED_MASS_VS_WWTP_PRESENCE",
              json.dumps({"group":group,"nutrient":nutrient,
                          "exposed":models[nutrient]["original_exposed_sites"],
                          "baseline_gain_pct":{
                            k:v["percent_improvement_from_baseline"]
                            for k,v in pairs.items() if "percent_improvement_from_baseline" in v
                          }},sort_keys=True),flush=True)
        result[group]=models
    a.out.parent.mkdir(parents=True,exist_ok=True)
    out={
       "schema":"gesui_EEA2022_BQE_posthoc_mass_presence_count_comparative_validation_v01",
       "status":"PASS_SAME_ORIGINAL_SITES_FOLDS_POSTHOC_LOAD_AMPLITUDE_COMPARISON",
       "original_source_only_run":42432030617,
       "original_biological_outcome_run":42432384196,
       "source_station_registry_SHA256":SRC_HASH,
       "source_EEA_2022_N_P_both_measured_positive_set_only":True,
       "source_WWTP_count_is_registered_both_positive_annual_load_facilities":True,
       "added_model_variants_are_posthoc_and_not_preregistered":True,
       "original_country_method_and_hydrology_predictors_unchanged":True,
       "fixed_3deg_geographic_holdout_unchanged":True,
       "no_plant_or_GBIF_response_reanalysed":True,
       "no_new_source_network_download":True,
       "not_a_causal_exposure_or_P_vs_N_mechanism_test":True,
       "limits":[
         "Presence count and annual reported source mass correlate and are derived from same limited EEA both-positive measured facility population",
         "Incremental prediction gain beyond binary source presence cannot rule out urbanization/agricultural runoff and monitoring selection",
         "Posthoc exploratory new model comparisons must not be promoted to confirmatory preregistered significance",
         "Same source registry excludes other unreported or multiple-outlet wastewater facilities",
         "National EQR scales and source country/method levels differ, P and N mass cannot be disentangled",
       ],
       "results":result}
    a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_EEA2022_BQE_MASS_VS_PRESENCE_AUDIT",
          json.dumps({"status":out["status"],"groups":list(result)},sort_keys=True),flush=True)

if __name__=="__main__":
    main()
