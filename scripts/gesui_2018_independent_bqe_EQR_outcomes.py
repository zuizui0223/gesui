#!/usr/bin/env python3
"""Independent-year 2018 WWTP reported measured N/P → 2018 biological EQR.

Primary 2018 benthic and secondary phytobenthos river biological response
scores are opened ONLY AFTER genuine historic 2018 WWTP source and source-
only 2018 EQR site IDs/fold definitions were committed to this branch.
The 2022 model family, 3degree block construction, support gates, Ridge
alpha=10 and phys DOWNSTREAM negative controls were fixed in advance.
No 2022 nutrient-source or biological responses used in the model.
"""
from __future__ import annotations
import argparse,hashlib,io,json,math,zipfile
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler,OneHotEncoder

SOURCE_REGISTRY_SHA="5a074f5182711c4ab74ef9bdacfa3918e9fcb434dd0b75c67072a0f06459c3e1"
BIOLOGY_SOURCE_SHA="f11001a57172b618bee95f58fd0ec51c8e71be7b9c01a70927dde21b7ca6e2a9"
ZIP_MEMBER=(
  "eea_t_waterbase-biology-2026_p_1990-2024_v01_r00/"
  "Waterbase_v2025_1_T_WISE2_BiologyEQRData.csv"
)
HISTORICAL_SOURCE_RUN=42441999786
YEAR="2018"
BQE={"EEA_13-01-4":"benthic_invertebrates",
     "EEA_124-04-9":"phytobenthos",
     "EEA_14-02-8":"fish"}
SRC={
 "P":("2018_measured_upstream_P_t_per_year",
      "2018_measured_downstream_P_t_per_year_placebo"),
 "N":("2018_measured_upstream_N_t_per_year",
      "2018_measured_downstream_N_t_per_year_placebo")
}
MODEL_SPEC={
 "M0_country_method_hydrogeo":(False,False),
 "M1_plus_2018_upstream_measured_mass":(True,False),
 "M2_plus_2018_downstream_source_placebo":(False,True),
 "M3_upstream_and_downstream":(True,True)}
REQUIRED={
 "min_numeric_sites":100,
 "min_upstream_2018_measured_source_exposed":20,
 "min_reported_countries":5,
 "min_upstream_source_countries":3,
 "min_3degree_geographic_blocks":10,
}
ALPHA=10.
FOLDS=5
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    return h.hexdigest()

def read_2018_outcomes(path):
    with zipfile.ZipFile(path) as z:
        h=hashlib.sha256()
        with z.open(ZIP_MEMBER) as fd:
            for chunk in iter(lambda:fd.read(1024*1024),b""):h.update(chunk)
        if h.hexdigest()!=BIOLOGY_SOURCE_SHA:
            raise ValueError("Official EEA Biology original EQR file bytes changed")
        cols=[
          "countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
          "parameterWaterBodyCategory","observedPropertyDeterminandBiologyEQRCode",
          "phenomenonTimeReferenceYear","resultNormalisedEQRValue","resultEQRValue",
          "resultObservationStatus","resultEcologicalStatusClassValue"]
        data=pd.read_csv(z.open(ZIP_MEMBER),usecols=cols,dtype="string",low_memory=False)
    filt=(data.parameterWaterBodyCategory.eq("RW")&
          data.phenomenonTimeReferenceYear.eq(YEAR)&
          data.observedPropertyDeterminandBiologyEQRCode.isin(BQE))
    data=data.loc[filt].copy()
    keys=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme"]
    for k in keys:data[k]=data[k].str.strip()
    data["site_key"]=data[keys].fillna("").agg("|".join,axis=1)
    if data.duplicated(["site_key","observedPropertyDeterminandBiologyEQRCode"]).any():
        raise ValueError("Duplicated 2018 biological site/quality-element/year observed")
    observed=pd.to_numeric(data.resultNormalisedEQRValue,errors="coerce")
    raw=pd.to_numeric(data.resultEQRValue,errors="coerce")
    data["numeric_EQR_2018"]=observed.to_numpy(dtype=float,na_value=np.nan)
    data["raw_EQR_2018"]=raw.to_numpy(dtype=float,na_value=np.nan)
    return data,h.hexdigest()

def design(df,up,down,add_up,add_down):
    x=pd.DataFrame(index=df.index)
    x["lat"]=pd.to_numeric(df.lat,errors="raise")
    x["lon"]=pd.to_numeric(df.lon,errors="raise")
    x["log_Q"]=np.log1p(pd.to_numeric(df.DIS_AV_CMS,errors="raise"))
    x["log_area"]=np.log1p(pd.to_numeric(df.UPLAND_SKM,errors="raise"))
    nums=["lat","lon","log_Q","log_area"]
    if add_up:
        x["upstream_reported_mass_2018"]=np.log1p(pd.to_numeric(df[up],errors="raise"))
        nums.append("upstream_reported_mass_2018")
    if add_down:
        x["downstream_reported_mass_2018"]=np.log1p(pd.to_numeric(df[down],errors="raise"))
        nums.append("downstream_reported_mass_2018")
    x["country"]=df.countryCode.astype("string").fillna("UNKNOWN").astype(str)
    x["country_method"]=x.country+"|"+df.procedureClassificationSystem.astype("string").fillna("UNKNOWN").astype(str)
    if not np.isfinite(x[nums].to_numpy(float)).all():
        raise ValueError("Nonfinite numerical predictors in source frozen 2018 EQR model")
    return x,nums

def spatial_cv(frame,nutrient):
    groups=pd.to_numeric(frame.pre_2018_source_only_spatial_fold,errors="raise").to_numpy(int)
    if set(groups)!={1,2,3,4,5}:
        raise ValueError("Pre-frozen source-only 2018 geographic folds incomplete")
    up,down=SRC[nutrient]
    y=frame.numeric_EQR_2018.to_numpy(float)
    metrics={};residuals={}
    for name,(add_up,add_down) in MODEL_SPEC.items():
        x,nums=design(frame,up,down,add_up,add_down)
        predicted=np.full(len(frame),np.nan)
        folds=[]
        for k in range(1,FOLDS+1):
            test=np.flatnonzero(groups==k)
            train=np.flatnonzero(groups!=k)
            if not len(test) or not len(train):raise ValueError("Empty 2018 regional heldout fold")
            pipeline=Pipeline([
              ("preprocess",ColumnTransformer([
                ("numeric",StandardScaler(),nums),
                ("categorical",OneHotEncoder(handle_unknown="ignore"),
                 ["country","country_method"])
              ])),
              ("regression",Ridge(alpha=ALPHA))])
            pipeline.fit(x.iloc[train],y[train])
            pred=pipeline.predict(x.iloc[test])
            predicted[test]=pred
            folds.append({"fold":k,"test_site_n":len(test),
                 "test_original_3deg_blocks":int(frame.iloc[test].region_3deg_block.nunique()),
                 "MSE":float(np.mean((y[test]-pred)**2)),
                 "MAE":float(np.mean(np.abs(y[test]-pred)))})
        if not np.isfinite(predicted).all():
            raise ValueError("Missing original 2018 EQR predictions")
        err=y-predicted
        metrics[name]={"heldout_MSE":float(np.mean(err**2)),
                       "heldout_MAE":float(np.mean(np.abs(err))),
                       "regional_fold_results":folds}
        residuals[name]=err**2
    baseline=metrics["M0_country_method_hydrogeo"]["heldout_MSE"]
    if baseline<=0:raise ValueError("Degenerate 2018 biological EQR baseline")
    improvements={}
    for name in MODEL_SPEC:
        if name=="M0_country_method_hydrogeo":continue
        change=baseline-metrics[name]["heldout_MSE"]
        perblock=[]
        for block in sorted(frame.region_3deg_block.astype(str).unique()):
            f=frame.region_3deg_block.astype(str).eq(block).to_numpy()
            perblock.append(float((residuals["M0_country_method_hydrogeo"][f]-
                                  residuals[name][f]).mean()))
        improvements[name]={
           "MSE_absolute_improvement":float(change),
           "MSE_relative_improvement_pct":float(100*change/baseline),
           "positive_3deg_geographic_blocks":sum(x>0 for x in perblock),
           "total_3deg_geographic_blocks":len(perblock)}
    plus=metrics["M1_plus_2018_upstream_measured_mass"]["heldout_MSE"]
    placebo=metrics["M2_plus_2018_downstream_source_placebo"]["heldout_MSE"]
    temporal_replication=bool(plus<baseline and plus<placebo)
    return {"nutrient":nutrient,
      "same_station_denominator":True,
      "original_2018_preoutcome_fold_assignments_used":True,
      "original_2022_BQE_data_used_in_2018_model":False,
      "Ridge_fixed_alpha":ALPHA,
      "heldout_MSE_MAE_and_folds":metrics,
      "geographic_MSE_changes":improvements,
      "passes_predeclared_temporal_prediction_direction_rule":temporal_replication}

def one_group(original,source,code,label):
    x=source[source.observedPropertyDeterminandBiologyEQRCode.eq(code)].copy()
    r=original[original.observedPropertyDeterminandBiologyEQRCode.eq(code)].copy()
    if x.duplicated("site_key").any() or r.duplicated("site_key").any():
        raise ValueError("Duplicate original 2018 site identity")
    if not set(x.site_key)<=set(r.site_key):
        raise ValueError("A fixed 2018 source monitoring site missing from original archived EQR")
    merged=x.merge(
        r[["site_key","observedPropertyDeterminandBiologyEQRCode",
           "numeric_EQR_2018","raw_EQR_2018","resultObservationStatus",
           "resultEcologicalStatusClassValue"]],
        on=["site_key","observedPropertyDeterminandBiologyEQRCode"],
        how="left",suffixes=("_source","_original"),validate="one_to_one")
    values=merged.numeric_EQR_2018.to_numpy(float)
    valid=np.isfinite(values)&(values>=0)
    v=merged.loc[valid].copy().reset_index(drop=True)
    exposed_all=merged["2018_upstream_qualifying_WWTP_count"].astype(float).gt(0).to_numpy()
    exposed_valid=v["2018_upstream_qualifying_WWTP_count"].astype(float).gt(0)
    countries=int(v.countryCode.nunique())
    sourcecountries=int(v.loc[exposed_valid,"countryCode"].nunique())
    blocks=int(v.region_3deg_block.nunique())
    summary={
      "2018_original_BQE_source_only_river_sites":len(x),
      "2018_original_EQR_site_values_matched":int(merged.numeric_EQR_2018.notna().sum()),
      "2018_valid_original_numeric_normalized_EQR":len(v),
      "2018_missing_invalid_original_normalized_EQR":int((~valid).sum()),
      "2018_original_raw_EQR_numeric_among_normalized_missing":
          int(np.isfinite(merged.loc[~valid,"raw_EQR_2018"].to_numpy(float)).sum()),
      "2018_all_source_only_upstream_reported_measured_NP_sites":int(exposed_all.sum()),
      "2018_valid_EQR_upstream_reported_measured_NP_sites":int(exposed_valid.sum()),
      "2018_missing_upstream_reported_measured_NP_sites":
          int((exposed_all&~valid).sum()),
      "2018_country_n":countries,
      "2018_upstream_source_country_n":sourcecountries,
      "2018_spatial_3deg_blocks":blocks,
      "2018_numeric_EQR_country_counts":{
        str(k):int(n) for k,n in v.countryCode.value_counts().items()},
      "2018_missing_by_country":{
        str(k):int(n) for k,n in
        merged.loc[~valid,"countryCode"].value_counts().items()},
      "2018_missing_source_positive_by_country":{
        str(k):int(n) for k,n in
        merged.loc[exposed_all&~valid,"countryCode"].value_counts().items()},
      "2018_numeric_EQR_range":[float(v.numeric_EQR_2018.min()),
                                 float(v.numeric_EQR_2018.max())] if len(v) else None,
      "2018_original_EQR_above_1":int(v.numeric_EQR_2018.gt(1).sum()),
      "2018_mean_normalized_EQR":float(v.numeric_EQR_2018.mean()) if len(v) else None,
    }
    supported=bool(
       len(v)>=REQUIRED["min_numeric_sites"] and
       exposed_valid.sum()>=REQUIRED["min_upstream_2018_measured_source_exposed"] and
       countries>=REQUIRED["min_reported_countries"] and
       sourcecountries>=REQUIRED["min_upstream_source_countries"] and
       blocks>=REQUIRED["min_3degree_geographic_blocks"])
    summary["predeclared_2018_validation_source_and_EQR_support_pass"]=supported
    result={"element":label,"original_2018_EEA_EQR_code":code,
            "year":2018,"source_QC":summary}
    if supported:
        analyses={}
        for nutrient in ("P","N"):
            up=SRC[nutrient][0]
            rho,p=spearmanr(np.log1p(pd.to_numeric(v[up],errors="raise")),
                           v.numeric_EQR_2018.to_numpy(float))
            analyses[nutrient]={
             "descriptive_source_mass_vs_EQR_spearman_rho":float(rho),
             "descriptive_source_mass_vs_EQR_spearman_p_unadjusted":float(p),
             "fixed_source_year_geographic_holdout":spatial_cv(v,nutrient)}
            ms=analyses[nutrient]["fixed_source_year_geographic_holdout"]
            print("GESUI_INDEPENDENT_2018_BENTHIC_EQR_TEMPORAL_REPLICATION",
              json.dumps({
                "BQE":label,"nutrient":nutrient,"year":2018,
                "original_2018_valid_EQR_sites":len(v),
                "source_2018_upstream_qualifying_NP_sites":int(exposed_valid.sum()),
                "original_2018_reporting_countries":countries,
                "M0_MSE":ms["heldout_MSE_MAE_and_folds"]["M0_country_method_hydrogeo"]["heldout_MSE"],
                "M1_upstream_MSE":ms["heldout_MSE_MAE_and_folds"]["M1_plus_2018_upstream_measured_mass"]["heldout_MSE"],
                "M2_downstream_placebo_MSE":ms["heldout_MSE_MAE_and_folds"]["M2_plus_2018_downstream_source_placebo"]["heldout_MSE"],
                "upstream_MSE_improvement_pct":ms["geographic_MSE_changes"]["M1_plus_2018_upstream_measured_mass"]["MSE_relative_improvement_pct"],
                "source_vs_river_EQR_Spearman":float(rho),
                "temporal_direction_replication":ms["passes_predeclared_temporal_prediction_direction_rule"]
              },sort_keys=True),flush=True)
        result["status"]="2018_INDEPENDENT_YEAR_VALID_BQE_TEST_COMPLETE__OBSERVATIONAL"
        result["nutrients"]=analyses
    else:
        result["status"]="2018_NOT_TESTABLE_ORIGINAL_SOURCE_AND_BIOLOGY_SUPPORT_GATE_FAILED"
    print("GESUI_2018_BQE_OUTCOME_QUALITY_GATE",json.dumps({
       "BQE":label,"source":summary,
       "status":result["status"]},sort_keys=True),flush=True)
    return result,v

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--frozen-source",type=Path,required=True)
    p.add_argument("--source-receipt",type=Path,required=True)
    p.add_argument("--original-Biology-2026",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    if sha(args.frozen_source)!=SOURCE_REGISTRY_SHA:
        raise ValueError("Prior 2018 source-only station years/exposures/folds changed after freeze")
    receipt=json.loads(args.source_receipt.read_text())
    assert receipt["frozen_before_any_independent_2018_EQR_outcome_fitted"] is True
    assert receipt["original_2018_Biology_EQR_data_not_opened"] is True
    assert receipt["physical_2018_nutrient_to_BQE_station_year_source_only_run"]==HISTORICAL_SOURCE_RUN
    frame=pd.read_csv(args.frozen_source,dtype={
      "site_key":str,"countryCode":str,"monitoringSiteIdentifier":str,
      "monitoringSiteIdentifierScheme":str,"observedPropertyDeterminandBiologyEQRCode":str,
      "procedureClassificationSystem":str,
      "region_3deg_block":str},low_memory=False)
    if len(frame)!=6364 or frame[["site_key","observedPropertyDeterminandBiologyEQRCode"]].duplicated().any():
        raise ValueError("Original 2018 source-only site/BQE frame changed")
    observed,hash_e=read_2018_outcomes(args.original_Biology_2026)
    args.out.mkdir(parents=True,exist_ok=True)
    groups={}
    for code,label in BQE.items():
        outcome,eligible=one_group(observed,frame,code,label)
        groups[label]=outcome
        if outcome["status"].startswith("2018_INDEPENDENT_YEAR_VALID"):
            fields=["site_key","countryCode","observedPropertyDeterminandBiologyEQRCode",
                    "procedureClassificationSystem","lat","lon","waterBodyIdentifier",
                    "region_3deg_block","pre_2018_source_only_spatial_fold",
                    "HYRIV_ID","DIS_AV_CMS","UPLAND_SKM",
                    "numeric_EQR_2018","raw_EQR_2018",
                    "2018_measured_upstream_P_t_per_year","2018_measured_upstream_N_t_per_year",
                    "2018_upstream_qualifying_WWTP_count",
                    "2018_measured_downstream_P_t_per_year_placebo",
                    "2018_measured_downstream_N_t_per_year_placebo",
                    "resultObservationStatus_original"]
            keep=[x for x in fields if x in eligible]
            eligible[keep].to_csv(args.out/(label+"_2018_WWTP_source_real_EQR.csv"),
                                  index=False,float_format="%.12g")
    report={
      "schema":"gesui_2018_independent_temporal_source_year_real_NP_WWTP_BQE_EQR_v01",
      "date_utc":datetime.now(timezone.utc).isoformat(),
      "validation_year_2018_frozen_before_results":True,
      "original_2018_source_only_run":HISTORICAL_SOURCE_RUN,
      "2018_source_only_registry_sha256":SOURCE_REGISTRY_SHA,
      "official_EEA_Biology_original_uncompressed_SHA256":hash_e,
      "2018_EEA_UWWTD_source_archive_version":"eea_t_uwwtd-waterbase_p_2017-2018_v07_r01",
      "source_2018_measured_loads_are_not_2022_WWTP_loads":True,
      "biological_year_matches_WWTP_report_year":True,
      "spatial_site_source_values_and_folds_frozen_before_2018_EQR":True,
      "all_predeclared_biological_classes_preserved":list(BQE.values()),
      "original_2018_model_selection_not_tuned_to_2018_outcomes":True,
      "original_2022_benthic_prediction_MSE_gains":{"P_percent":1.103,"N_percent":1.397},
      "sewage_causal_benthic_effect_proven":False,
      "limitations":[
       "This is retrospective independently chosen year 2018, not a randomized WWTP emission change",
       "2018 EQR source country/method and monitor locations differ from 2022, and missing EQR must not be imputed",
       "Historical 2018 UWWTD measured annual load subset may differ from 2022 reporting completeness",
       "2018 measured WWTP load is not actual received N or P in a BQE sampling stream",
       "P and N source masses may collinear, inability to isolate P-vs-N cause",
       "Predefined source-only 3degree folds, country and national assessment method adjustment do not remove agriculture/urban/hydro confounds",
       "Official HydroRIVERS original not independently reproduced from the 3rd-party transformed v1 EU mirror",
       "Same EEA Biology2026 source release spans both 2018 and 2022, so temporal replication is not independent-institution replication"
      ],
      "results":groups}
    (args.out/"2018_ORIGINAL_TEMPORAL_BQE_REPLICATION.json").write_text(
      json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_INDEPENDENT_2018_TEMPORAL_REPLICATION_ALL_GROUPS",
      json.dumps({p:{"status":g["status"],
                     "site_n":g["source_QC"]["2018_valid_original_numeric_normalized_EQR"],
                     "source_exposed":g["source_QC"]["2018_valid_EQR_upstream_reported_measured_NP_sites"],
                     "country_n":g["source_QC"]["2018_country_n"]}
                   for p,g in groups.items()},sort_keys=True),flush=True)
if __name__=="__main__":
    main()
