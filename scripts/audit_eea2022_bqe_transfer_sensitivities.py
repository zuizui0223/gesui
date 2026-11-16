#!/usr/bin/env python3
"""Posthoc EEA-2022 aquatic biological source prediction stress tests.

This is NOT the preregistered model. It repeats the ORIGINAL fixed-fold model
to inspect out-of-fold error by upstream-source presence, applies non-linear
geographic covariates without biological effect-selection, and runs an
explicit leave-one-entire-country-out transfer stress test. All original
sites/outcomes and reported N/P source values are unchanged and SHA-locked.
No new source queries or biological mechanisms inferred.
"""
from __future__ import annotations
import argparse,hashlib,json,sys,math
from pathlib import Path
import numpy as np,pandas as pd,sklearn
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline

SOURCE_DIRECTORY=Path("results/eea2022_biology_bqe_20261108")
CODE_GROUPS={
 "benthic_invertebrates":"BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv",
 "phytobenthos":"PHYTOBENTHOS_2022_SOURCE_EQR.csv"}
SRC={"N":"upstream_2022_EEA_N_measured_t_per_yr",
     "P":"upstream_2022_EEA_P_measured_t_per_yr"}
MODELS={"benthic_invertebrates":{"P":{"baseline":0.05833804346777809,
                                     "source":0.05769459552059828},
                                "N":{"baseline":0.05833804346777809,
                                     "source":0.05752327324679533}},
        "phytobenthos":{"P":{"baseline":0.035934357335444235,
                             "source":0.035875931550570526},
                         "N":{"baseline":0.035934357335444235,
                             "source":0.03587108895518606}}}
TOL=4e-6
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def frame(df,nutrient,include_load,nonlinear_geo):
    out=pd.DataFrame(index=df.index)
    out["lat"]=df.lat.to_numpy(float)
    out["lon"]=df.lon.to_numpy(float)
    out["logQ"]=np.log1p(df.DIS_AV_CMS.to_numpy(float))
    out["logArea"]=np.log1p(df.UPLAND_SKM.to_numpy(float))
    nums=["lat","lon","logQ","logArea"]
    if nonlinear_geo:
        out["lat_center_sq"]=(out.lat-50)**2
        out["lon_center_sq"]=(out.lon-10)**2
        out["latlon"]=(out.lat-50)*(out.lon-10)
        nums+=["lat_center_sq","lon_center_sq","latlon"]
    if include_load:
        out["log_upstream_reported_"+nutrient]=np.log1p(df[SRC[nutrient]].to_numpy(float))
        nums.append("log_upstream_reported_"+nutrient)
    out["country"]=df.countryCode.astype("string").fillna("UNKNOWN").astype(str)
    out["country_method"]=out.country+"|"+df.procedureClassificationSystem.astype("string").fillna("UNKNOWN").astype(str)
    if not np.isfinite(out[nums].to_numpy(float)).all():
        raise ValueError("Nonfinite source-only predictors")
    return out,nums

def fit_predict(df,nutrient,include_load,scheme,nonlinear_geo=False):
    x,nums=frame(df,nutrient,include_load,nonlinear_geo)
    y=df.normalised_EQR.to_numpy(float)
    if scheme=="fixed_geo_blocks":
        groups=df.fixed_spatial_cv_fold.astype(int).to_numpy()
        split=[(np.flatnonzero(groups!=f),np.flatnonzero(groups==f))
               for f in range(1,6)]
        if set(groups)!={1,2,3,4,5}:raise ValueError("Frozen fold labels changed")
    elif scheme=="leave_one_reporting_country_out":
        countries=df.countryCode.astype(str).to_numpy()
        split=[(np.flatnonzero(countries!=country),
                np.flatnonzero(countries==country)) for country in sorted(set(countries))]
    else:raise ValueError("Unrecognized sensitivity scheme")
    pred=np.full(len(df),np.nan)
    for train,test in split:
        if not len(train) or not len(test):raise ValueError("Empty heldout sample")
        model=Pipeline([
          ("source_covariate_processing",ColumnTransformer([
            ("numeric",StandardScaler(),nums),
            ("categorical",OneHotEncoder(handle_unknown="ignore"),["country","country_method"])
          ])),
          ("fixed_ridge",Ridge(alpha=10.0))])
        model.fit(x.iloc[train],y[train])
        pred[test]=model.predict(x.iloc[test])
    if not np.isfinite(pred).all():raise ValueError("Missing heldout model predictions")
    return (y-pred)**2

def compare(df,nutrient,scheme,nonlinear_geo,check_original=False):
    baseline=fit_predict(df,nutrient,False,scheme,nonlinear_geo)
    plus=fit_predict(df,nutrient,True,scheme,nonlinear_geo)
    bmean=float(np.mean(baseline))
    cmean=float(np.mean(plus))
    if check_original:
        saved=MODELS["benthic_invertebrates" if len(df)==2475 else "phytobenthos"][nutrient]
        if abs(bmean-saved["baseline"])>TOL or abs(cmean-saved["source"])>TOL:
            raise ValueError("Original source-frozen sample/method not reproduced")
    exposed=df[SRC[nutrient]].to_numpy(float)>0
    def evaluate(mask):
        if not np.any(mask):return None
        mse_b=float(np.mean(baseline[mask]))
        mse_x=float(np.mean(plus[mask]))
        return {"n":int(np.sum(mask)),"baseline_MSE":mse_b,
                "source_MSE":mse_x,
                "improvement_MSE_pct":100*(mse_b-mse_x)/mse_b}
    bycountry={}
    for country in sorted(df.countryCode.astype(str).unique()):
        mask=df.countryCode.astype(str).eq(country).to_numpy()
        bycountry[country]={
           "all_station_n":int(mask.sum()),
           "upstream_source_positive_n":int(np.sum(mask&exposed)),
           "MSE_improvement_pct":evaluate(mask)["improvement_MSE_pct"]}
    return {
      "evaluation_scheme":scheme,
      "posthoc_nonlinear_geo":nonlinear_geo,
      "source":nutrient,
      "all_stations":evaluate(np.ones(len(df),dtype=bool)),
      "registered_upstream_load_positive_stations":evaluate(exposed),
      "source_zero_registered_stations":evaluate(~exposed),
      "country_specific_heldout_contributions":bycountry,
      "countries_with_positive_MSE_gain":sum(x["MSE_improvement_pct"]>0 for x in bycountry.values()),
      "country_count":len(bycountry)
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--results-dir",type=Path,default=SOURCE_DIRECTORY)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    receipt=json.loads((a.results_dir/"SOURCE_RECEIPT.json").read_text())
    original={x["file"]:x["original_file_SHA256"] for x in receipt["all_original_files"]}
    if receipt["original_measured_BQE_run"]!=42432384196 or receipt["original_source_run"]!=42432030617:
        raise ValueError("Original source run not as expected")
    originalsummary=json.loads((a.results_dir/"ORIGINAL_EQR_BIOLOGICAL_RESULTS.json").read_text())
    audited=json.loads((a.results_dir/"INDEPENDENT_AUDIT.json").read_text())
    assert audited["status"]=="PASS_INDEPENDENT_OFFLINE_RECONSTRUCTION_2022_BQE"
    result={}
    for group,filename in CODE_GROUPS.items():
        if sha(a.results_dir/filename)!=original[filename]:
            raise ValueError("Original EEA station-level outcomes modified: "+filename)
        d=pd.read_csv(a.results_dir/filename,low_memory=False)
        expected_n=originalsummary["results"][group]["support"]["numeric_original_EQR_sites"]
        assert len(d)==expected_n
        assert not d.site_key.duplicated().any()
        group_results={}
        for nutrient in ("P","N"):
            # Original pre-registered 3-degree region CV is reproduced exactly.
            original_cv=compare(d,nutrient,"fixed_geo_blocks",False,True)
            nonlinear=compare(d,nutrient,"fixed_geo_blocks",True)
            # Posthoc country-out-of-sample robustness: does not assume an unseen
            # country's EQR classification scale is calibrated to training countries.
            country=compare(d,nutrient,"leave_one_reporting_country_out",False)
            group_results[nutrient]={
              "original_fixed_3degree_fold_exposure_split_diagnostic":original_cv,
              "posthoc_fixed_fold_quadratic_latlon_sensitivity":nonlinear,
              "posthoc_leave_one_country_out_transfer_stress_test":country}
            print("GESUI_BQE_POSTHOC_COUNTRY_TRANSPORT",json.dumps({
              "group":group,"nutrient":nutrient,
              "original_gain_pct":original_cv["all_stations"]["improvement_MSE_pct"],
              "exposed_site_gain_pct":original_cv["registered_upstream_load_positive_stations"]["improvement_MSE_pct"],
              "source_zero_site_gain_pct":original_cv["source_zero_registered_stations"]["improvement_MSE_pct"],
              "quadratic_geo_gain_pct":nonlinear["all_stations"]["improvement_MSE_pct"],
              "country_holdout_gain_pct":country["all_stations"]["improvement_MSE_pct"]},sort_keys=True),flush=True)
        result[group]=group_results
    report={
      "schema":"gesui_2022_BQE_EEA_source_posthoc_country_and_exposure_specific_stress_tests_v01",
      "status":"PASS_ORIGINAL_FIXED_MODELS_AND_EXPLORATORY_TRANSFER_STRESS_TESTS",
      "original_source_only_run":42432030617,
      "original_observed_BQE_run":42432384196,
      "independent_source_original_audit_run":42436814278,
      "Ridge_fixed_alpha":10,
      "no_new_sites_or_source_variables_added":True,
      "original_geographic_folds_unchanged":True,
      "posthoc_diagnostics_not_preregistered":True,
      "models_refitted_only_to_test_preexisting_predictions_generalization":True,
      "scikit_learn_version":sklearn.__version__,
      "country_leaveout_new_country_EQR_method_levels_unseen_in_training":True,
      "same_source_only_2022_EEA_N_and_P_facilities_collinear":True,
      "no_causal_effluent_pollution_effect_claimed":True,
      "limits":[
        "All exposure-positive subsets evaluated only after observing primary overall MSE: exploratory, not a prospective test",
        "Country leaveout holds out the country-specific EQR national method calibration, making transfer nonexchangeable",
        "Nonlinear latitude/longitude is a geographic proxy, not independent agricultural or urban land cover",
        "Site monitoring and national classification quality vary; source P and N are highly correlated",
        "No 2022 measured river flow seasonality, adjacent other P/N agricultural loads or actual receiving reach transport"
      ],
      "results":result}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_BQE_TRANSFER_SENSITIVITY_COMPLETE",json.dumps({
      "status":report["status"],"groups":list(result),
      "original_source_run_verified":True,"external_queries":0},sort_keys=True),flush=True)

if __name__=="__main__":main()
