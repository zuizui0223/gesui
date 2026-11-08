#!/usr/bin/env python3
"""Posthoc strict two-way temporal AND monitoring-site-disjoint transport audit.

Fit fixed original country/method/hydrology Ridge(alpha=10) on ALL authentic
2018 EQR sites, predict ONLY 2022 sites lacking the same official site_key
in 2018; reverse 2022 -> 2018-exclusive. No ecological site from one year is
ever evaluated after appearing in the training year. The 2018 and 2022
EEA WWTP source loads are their *own actual measured report-year* values.

THIS TEST IS POSTHOC (the 2018/2022 EQR signals were already observed),
NOT a prespecified/confirmatory replication, and cannot prove WWTP causation.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge

FILES={
 "benthic_invertebrates":{
    "2018":"benthic_invertebrates_2018_WWTP_source_real_EQR.csv",
    "2018sha":"2a1a72d43656aca588a832dfae13860054c01a651e76bf25f984c5a374709b02",
    "2022":"BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv",
    "2022sha":"5e1f3c31b389d0bd769984ab055d3e2c54b63ccaa4f5eedfb2d82bcdc83c28e6",
    "2018n":2634,"2022n":2475},
 "phytobenthos":{
    "2018":"phytobenthos_2018_WWTP_source_real_EQR.csv",
    "2018sha":"8302daf2afa3f3556debc5ca769ee667bd06816d0de2e63afd4cb6b538587ad7",
    "2022":"PHYTOBENTHOS_2022_SOURCE_EQR.csv",
    "2022sha":"7fa5dfff9023337ac0a62a87fe406eb71783979858d9d924ab97019c6e4e5fc3",
    "2018n":2670,"2022n":2784}
}
EXPOSURES={
 "2018":{
  "P":["2018_measured_upstream_P_t_per_year",
       "2018_measured_downstream_P_t_per_year_placebo"],
  "N":["2018_measured_upstream_N_t_per_year",
       "2018_measured_downstream_N_t_per_year_placebo"]},
 "2022":{
  "P":["upstream_2022_EEA_P_measured_t_per_yr",
       "downstream_2022_EEA_P_t_per_yr_placebo"],
  "N":["upstream_2022_EEA_N_measured_t_per_yr",
       "downstream_2022_EEA_N_t_per_yr_placebo"]}
}
MODELS={
 "M0_baseline":(False,False),
 "M1_upstream_reported_mass":(True,False),
 "M2_downstream_negative_control":(False,True),
 "M3_upstream_and_downstream":(True,True)}
BASELINE_NUMS=["latitude","longitude","log_Q","log_upstream_area"]
ALPHA=10.0
YEAR_RESPONSE={"2018":"numeric_EQR_2018","2022":"normalised_EQR"}
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def matrix(df,year,chem,add_up,add_down):
    data=pd.DataFrame(index=df.index)
    data["latitude"]=df.lat.to_numpy(float)
    data["longitude"]=df.lon.to_numpy(float)
    data["log_Q"]=np.log1p(df.DIS_AV_CMS.to_numpy(float))
    data["log_upstream_area"]=np.log1p(df.UPLAND_SKM.to_numpy(float))
    nums=BASELINE_NUMS.copy()
    if add_up:
        data["log_upstream_source"]=np.log1p(df[EXPOSURES[year][chem][0]].to_numpy(float))
        nums.append("log_upstream_source")
    if add_down:
        data["log_downstream_source"]=np.log1p(df[EXPOSURES[year][chem][1]].to_numpy(float))
        nums.append("log_downstream_source")
    data["country"]=df.countryCode.fillna("UNKNOWN").astype(str)
    data["country_method"]=data.country+"|"+df.procedureClassificationSystem.fillna("UNKNOWN").astype(str)
    if not np.isfinite(data[nums].to_numpy(float)).all():
        raise ValueError("Invalid numeric source features in time-transport design")
    return data,nums

def transport(train,test,train_year,test_year,chem):
    if set(train.site_key)&set(test.site_key):
        raise ValueError("Original site reused between training year and test year")
    ytrain=train[YEAR_RESPONSE[train_year]].to_numpy(float)
    ytest=test[YEAR_RESPONSE[test_year]].to_numpy(float)
    if not np.isfinite(ytrain).all() or not np.isfinite(ytest).all():
        raise ValueError("Unusable EQR biological responses")
    source_test=test[EXPOSURES[test_year][chem][0]].to_numpy(float)>0
    source_train=train[EXPOSURES[train_year][chem][0]].to_numpy(float)>0
    outcomes={}
    errors={}
    train_countries=set(train.countryCode.astype(str))
    train_methods=set(train.countryCode.astype(str)+"|"+train.procedureClassificationSystem.fillna("UNKNOWN").astype(str))
    test_countries=test.countryCode.astype(str)
    test_methods=test_countries+"|"+test.procedureClassificationSystem.fillna("UNKNOWN").astype(str)
    for name,(up,down) in MODELS.items():
        xtrain,nums=matrix(train,train_year,chem,up,down)
        xtest,nums_test=matrix(test,test_year,chem,up,down)
        if nums!=nums_test:raise ValueError("Cross-year model schema changed")
        model=Pipeline([
          ("source_normalizer",ColumnTransformer([
            ("numeric",StandardScaler(),nums),
            ("categorical",OneHotEncoder(handle_unknown="ignore"),["country","country_method"])
          ],remainder="drop")),
          ("ridge_fixed_10",Ridge(alpha=ALPHA))])
        model.fit(xtrain,ytrain)
        predicted=model.predict(xtest)
        square=(ytest-predicted)**2
        mse=float(np.mean(square))
        errors[name]=square
        outcomes[name]={
           "test_MSE":mse,
           "test_MAE":float(np.mean(np.abs(ytest-predicted)))}
    baseline=outcomes["M0_baseline"]["test_MSE"]
    improvement={}
    for name in ("M1_upstream_reported_mass",
                 "M2_downstream_negative_control","M3_upstream_and_downstream"):
        delta=(baseline-outcomes[name]["test_MSE"])/baseline*100
        improvement[name]={"MSE_improvement_vs_M0_percent":float(delta),
                           "predeclared_direction_upstream_outperforms_M0_and_placebo":
                             bool(name=="M1_upstream_reported_mass" and
                                  outcomes[name]["test_MSE"]<baseline and
                                  outcomes[name]["test_MSE"]<outcomes["M2_downstream_negative_control"]["test_MSE"])}
    def split_errors(mask):
        mse_b=float(errors["M0_baseline"][mask].mean()) if mask.any() else None
        mse_1=float(errors["M1_upstream_reported_mass"][mask].mean()) if mask.any() else None
        return {
          "sites":int(mask.sum()),
          "baseline_MSE":mse_b,"upstream_MSE":mse_1,
          "MSE_improvement_pct":100*(mse_b-mse_1)/mse_b if mse_b else None}
    return {
      "training_year":train_year,"evaluation_year":test_year,
      "training_site_n":len(train),"heldout_other_year_disjoint_site_n":len(test),
      "source_positive_training_n":int(source_train.sum()),
      "source_positive_test_n":int(source_test.sum()),
      "test_country_count":test.countryCode.nunique(),
      "test_countries_not_in_training":sorted(set(test_countries)-train_countries),
      "test_sites_unknown_country_n":int((~test_countries.isin(train_countries)).sum()),
      "test_sites_unknown_country_method_n":int((~test_methods.isin(train_methods)).sum()),
      "all_test_EQR_scores_site_disjoint_from_training":True,
      "MSE_and_MAE_by_model":outcomes,
      "MSE_changes":improvement,
      "source_positive_only_posthoc_QC":split_errors(source_test),
      "source_zero_only_posthoc_QC":split_errors(~source_test),
      "cross_year_reported_NP_mass_is_year_matched_not_constant":True,
      "method_calibration_and_country_coding_are_not_identical_in_years":True,
    }
def main():
    aparser=argparse.ArgumentParser()
    aparser.add_argument("--year2018",required=True,type=Path)
    aparser.add_argument("--year2022",required=True,type=Path)
    aparser.add_argument("--out",required=True,type=Path)
    args=aparser.parse_args()
    result={}
    for group,meta in FILES.items():
        data={}
        for yr in ("2018","2022"):
            original=args.year2018 if yr=="2018" else args.year2022
            file=original/meta[yr]
            if sha(file)!=meta[yr+"sha"]:
                raise ValueError("Original year EQR/site file altered: "+group+yr)
            data[yr]=pd.read_csv(file,dtype={"site_key":str},low_memory=False)
            if len(data[yr])!=meta[yr+"n"] or data[yr].site_key.duplicated().any():
                raise ValueError("Changed original source only biological sites")
        original_overlap=len(set(data["2018"].site_key)&set(data["2022"].site_key))
        disjoint22=data["2022"][~data["2022"].site_key.isin(set(data["2018"].site_key))].copy()
        disjoint18=data["2018"][~data["2018"].site_key.isin(set(data["2022"].site_key))].copy()
        result[group]={
          "site_year_original_overlap":original_overlap,
          "site_disjoint_test2022_n":len(disjoint22),
          "site_disjoint_test2018_n":len(disjoint18),
          "2018_training_to_2022_NEW_SITES":{
             chem:transport(data["2018"],disjoint22,"2018","2022",chem) for chem in ("P","N")},
          "2022_training_to_2018_NEW_SITES":{
             chem:transport(data["2022"],disjoint18,"2022","2018",chem) for chem in ("P","N")}}
        for design in ("2018_training_to_2022_NEW_SITES","2022_training_to_2018_NEW_SITES"):
            for chem in ("P","N"):
                q=result[group][design][chem]
                print("GESUI_2018_2022_SITE_DISJOINT_MODEL_YEAR_TRANSPORT",
                  json.dumps({"group":group,"transport":design,"chem":chem,
                    "train_n":q["training_site_n"],"test_n":q["heldout_other_year_disjoint_site_n"],
                    "test_exposed_n":q["source_positive_test_n"],
                    "test_new_methods":q["test_sites_unknown_country_method_n"],
                    "upstream_improvement":q["MSE_changes"]["M1_upstream_reported_mass"]["MSE_improvement_vs_M0_percent"],
                    "downstream_placebo_improvement":q["MSE_changes"]["M2_downstream_negative_control"]["MSE_improvement_vs_M0_percent"],
                    "upstream_outperforms_both":q["MSE_changes"]["M1_upstream_reported_mass"]["predeclared_direction_upstream_outperforms_M0_and_placebo"]
                    },sort_keys=True),flush=True)
    report={
      "schema":"gesui_two_year_site_disjoint_temporal_model_transport_posthoc_v01",
      "date_utc":datetime.now(timezone.utc).isoformat(),
      "scientific_status":"POSTHOC_TRANSFER_STRESS_TEST__NOT_NEW_PREREGISTERED_REPLICATION",
      "original_2018_measured_source_and_EQR_run":42442266935,
      "original_2022_measured_source_and_EQR_run":42432384196,
      "original_monitoring_source_sha256":{
        f"{group}_{yr}":meta[yr+"sha"] for group,meta in FILES.items() for yr in ("2018","2022")},
      "actual_year_specific_WWTP_discharge_mass_used":True,
      "all_heldout_test_sites_exclude_source_training_year_EEA_monitoring_site_keys":True,
      "source_model_Ridge_alpha_fixed":ALPHA,
      "no_2018_2022_hyperparameter_search":True,
      "no_2018_2022_same_site_training_data_leakage":True,
      "predictive_model_not_predeclared_before_2018_outcomes":True,
      "not_sewage_treatment_causal_inference":True,
      "results":result,
      "limitations":[
        "Posthoc exploratory cross-year transport after seeing original 2018 and 2022 prediction signals; not a new confirmatory test",
        "All original EQR scores are from a single EEA Biology 2026 harmonized release and can share country classification methods",
        "Changing reporting countries/method classes causes country and one-hot source-covariate mismatch; explicitly reported",
        "2018 historical reporting of WWTP N/P measured loads differs in coverage from 2022, zero-registered not zero-sewage",
        "Other agricultural loads, urban runoff and topology/geometry errors are not causally controlled",
        "Year-specific source positive groups may have different sample/site populations; this is NOT a before-after treatment design",
      ]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_TWO_WAY_2018_2022_SITE_DISJOINT_TEMPORAL_TRANSPORT_QC_SUCCESS",
       json.dumps({"groups":{g:{"overlap":x["site_year_original_overlap"],
              "new_sites_test2022":x["site_disjoint_test2022_n"],
              "new_sites_test2018":x["site_disjoint_test2018_n"]}
              for g,x in result.items()},"not_confirmatory":True},sort_keys=True),flush=True)
if __name__=="__main__":main()
