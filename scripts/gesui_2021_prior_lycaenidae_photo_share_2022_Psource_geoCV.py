#!/usr/bin/env python3
"""Posthoc historical photo-composition falsifier: do 2021 Lycaenidae
photographic reporting patterns explain 2022 reported WWTP-P source signal?

Strictly use same original 43 Germany/Austria 2022 source/photo/worldcover cells,
source-lock 2021 3-family photographed observations before any stage2 models.
No new API requests or photo images; does NOT prove 2021 WWTP absent.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr,spearmanr

MIN_SUPPORT=30
MIN_DE=18
MIN_AT=9
MIN_BLOCKS=12
COUNTRIES=("DE","AT")
P_SOURCE="registered_P_t_year"

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for data in iter(lambda:f.read(1048576),b""):h.update(data)
    return h.hexdigest()
def key(r):
    return (r["countryCode"],int(float(r["lat_floor"])),int(float(r["lon_floor"])))
def grp(r):return r["countryCode"]+"|"+r["country_2deg_block"]

def family_logit(l,n):
    return math.log((float(l)+.5)/(float(n)-float(l)+.5))

def source_beta_xy(train,test,historical=False):
    def vector(rows,func):return np.asarray([func(r) for r in rows],float)
    cols={
      "country_AT":lambda r:1. if r["countryCode"]=="AT" else 0.,
      "log_Plantae_2021":lambda r:math.log1p(float(r["2021_Plantae_photo_n"])),
      "geographic_lat":lambda r:float(r["lat_floor"])+.5,
      "geographic_lon":lambda r:float(r["lon_floor"])+.5,
      "WC2021_built":lambda r:float(r["WC2021_512_built_land_fraction"]),
      "WC2021_crop":lambda r:float(r["WC2021_512_crop_land_fraction"]),
      "WC2021_log_countryLand":lambda r:math.log1p(float(r["WC2021_512_country_terrestrial_area_km2"])),
    }
    X=[np.ones(len(train))];T=[np.ones(len(test))]
    for name,fun in cols.items():
        a=vector(train,fun);b=vector(test,fun)
        if name=="country_AT":
            X.append(a);T.append(b)
        else:
            avg=float(np.mean(a));s=float(np.std(a))
            if s<=1e-12:raise ValueError("Degenerate source/covariate "+name)
            X.append((a-avg)/s);T.append((b-avg)/s)
    if historical:
        cols_history={
           "prior2021_Lycaenidae_photo_logit":lambda r:family_logit(r["2021_Lycaenidae_photo_n"],r["2021_three_family_photo_total"]),
           "log_prior2021_butterfly_targetgroup_photo_count":lambda r:math.log1p(float(r["2021_three_family_photo_total"]))
        }
        for name,fun in cols_history.items():
            a=vector(train,fun);b=vector(test,fun)
            m=float(np.mean(a));s=float(np.std(a))
            if s<=1e-12:raise ValueError("Unusable prior-year butterfly historical covariate "+name)
            X.append((a-m)/s);T.append((b-m)/s)
    X0=np.column_stack(X);T0=np.column_stack(T)
    source=lambda r:math.log1p(float(r[P_SOURCE]))
    sourcetrain=vector(train,source);sourcetest=vector(test,source)
    ztrain=np.zeros(len(train));ztest=np.zeros(len(test))
    for country in COUNTRIES:
        k=np.asarray([r["countryCode"]==country for r in train])
        l=np.asarray([r["countryCode"]==country for r in test])
        if k.sum()<3:raise ValueError("Insufficient train country source cases")
        m=float(np.mean(sourcetrain[k]));s=float(np.std(sourcetrain[k]))
        if s<=1e-12:raise ValueError("Within-country training WWTP-P original source constant")
        ztrain[k]=(sourcetrain[k]-m)/s
        ztest[l]=(sourcetest[l]-m)/s
    X1=np.column_stack([X0,ztrain])
    T1=np.column_stack([T0,ztest])
    return X0,T0,X1,T1

def cv(rows,response_year,historical=False):
    if historical and response_year==2021:
        raise ValueError("Prior-year butterfly share would be identical to its response")
    groups=sorted({grp(r) for r in rows})
    preds=[]
    for group in groups:
        train=[r for r in rows if grp(r)!=group]
        test=[r for r in rows if grp(r)==group]
        X0,T0,X1,T1=source_beta_xy(train,test,historical)
        if np.linalg.matrix_rank(X0)<X0.shape[1] or np.linalg.matrix_rank(X1)<X1.shape[1]:
            raise ValueError("Rank-deficient country blocked historical photo model")
        resp=lambda r:family_logit(
            r[f"{response_year}_Lycaenidae_photo_n"],
            r[f"{response_year}_three_family_photo_total"])
        ytrain=np.asarray([resp(x) for x in train],float)
        ytest=np.asarray([resp(x) for x in test],float)
        pred0=T0@np.linalg.lstsq(X0,ytrain,rcond=None)[0]
        pred1=T1@np.linalg.lstsq(X1,ytrain,rcond=None)[0]
        for row,y,a,b in zip(test,ytest,pred0,pred1):
            preds.append({"country":row["countryCode"],"group":group,"lat":row["lat_floor"],
               "lon":row["lon_floor"],"photo_share_logit":float(y),
               "baseline_pred":float(a),"P_source_pred":float(b),
               "baseline_sqerr":float((y-a)**2),"P_source_sqerr":float((y-b)**2)})
    mse0=float(np.mean([p["baseline_sqerr"] for p in preds]))
    mse1=float(np.mean([p["P_source_sqerr"] for p in preds]))
    z={
       "year":response_year,"uses_2021_Lycaenidae_share_and_2021_target_family_photo_effort_in_baseline":bool(historical),
       "cells":len(rows),"heldout_country2deg_blocks":len(groups),
       "baseline_MSE":mse0,"source_P_added_MSE":mse1,
       "EEA2022_reported_P_source_incremental_MSE_improvement_percent":100*(mse0-mse1)/mse0,
       "country":{},
       "source_P_is_not_verified_contemporaneous_2021_WWTP_exposure":True,
       "photos_are_reported_records_not_butterfly_individuals":True,
       "source_fitted_in_training_blocks_only":True,
       "per_country2deg_heldout_cell_predictions":preds,
    }
    for c in COUNTRIES:
        s=[p for p in preds if p["country"]==c]
        a=float(np.mean([p["baseline_sqerr"] for p in s]))
        b=float(np.mean([p["P_source_sqerr"] for p in s]))
        z["country"][c]={"cells":len(s),"baseline_MSE":a,"P_source_MSE":b,
            "gain_percent":100*(a-b)/a}
    return z

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--historical-2021-family-photo-source",type=Path,required=True)
    p.add_argument("--frozen-2022-worldcover-43cell-model",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    earlier=json.loads(a.frozen_2022_worldcover_43cell_model.read_text())
    hist=json.loads(a.historical_2021_family_photo_source.read_text())
    if hist["schema"]!="gesui_prior2021_butterfly_photo_composition_source_only_same43_original_2022_WorldCover_sample_v01":
        raise ValueError("Unexpected original 2021 historical source")
    if hist["original_2021_WorldCover2022_photo_source_full_SHA256"]!=sha(a.frozen_2022_worldcover_43cell_model):
        raise ValueError("2021 family history not source frozen to exact previous 43-cell cohort")
    if not hist["POSTHOC_after_positive_2022_source_photo_and_landcover_effect_seen"]:
        raise ValueError("History must carry posthoc classification")
    if not earlier["original_WorldCover_country_land_area_and_resolution_gate_pass"]:
        raise ValueError("Previously reported 2021 independent landcov source QC did not pass")
    prior={key(r):r for r in hist["source_2021_same43_cell_historical_butterfly_photo_records"]}
    original=earlier["all_original_same_sample_cells_2022_photo_counts_and_WorldCover_categorical_2021_covariates"]
    if len(original)!=43 or len(prior)!=43:raise ValueError("Original 2022 landcover sample or 2021 history has unexpected size")
    joined=[]
    for x in original:
        if key(x) not in prior:raise ValueError("Original 43 coordinates not fixed between 2021 and 2022")
        r=prior[key(x)]
        if not r["2021_photo_support_min15"]:continue
        y=dict(x)
        y["2021_Lycaenidae_photo_n"]=r["2021_Lycaenidae_photo_n"]
        y["2021_Pieridae_photo_n"]=r["2021_Pieridae_photo_n"]
        y["2021_Papilionidae_photo_n"]=r["2021_Papilionidae_photo_n"]
        y["2021_three_family_photo_total"]=r["2021_three_family_photo_total"]
        # Retrieve 2022 three-family known ratio, not original individual photorecord data.
        y["2022_three_family_photo_total"]=x["2022_total_3family_photos"]
        y["2022_Lycaenidae_photo_n"]=x["2022_Lycaenidae_photo_n"]
        y["2021_Plantae_photo_n"]=x["2021_Plantae_photo_n"]
        y["lat_floor"]=x["lat_floor"];y["lon_floor"]=x["lon_floor"]
        y["country_2deg_block"]=x["country_2deg_block"]
        y["registered_P_t_year"]=x["registered_P_t_year"]
        joined.append(y)
    c=Counter(x["countryCode"] for x in joined)
    groups=set(grp(x) for x in joined)
    ready=(hist["2021_photo_support_before_historical_comparison_model_gate"] and
          len(joined)>=MIN_SUPPORT and c.get("DE",0)>=MIN_DE and
          c.get("AT",0)>=MIN_AT and len(groups)>=MIN_BLOCKS)
    result={
      "schema":"gesui_2021_historical_butterfly_photo_composition_prior_share_vs_2022_WWTP_P_posthoc_QC_v01",
      "POSTHOC_after_existing_2022_WWTP_P_positive_landcover_adjusted_photo_reporting":True,
      "original_2021_frozen_prior_year_family_photo_source_SHA256":sha(a.historical_2021_family_photo_source),
      "original_2022_WorldCover_quality_gated_family_photo_source_SHA256":sha(a.frozen_2022_worldcover_43cell_model),
      "original_2022_cells":43,
      "same_original_2021_and_2022_complete_family_photo_supported_cells":len(joined),
      "same_supported_by_country":dict(c),"same_supported_2degree_blocks":len(groups),
      "predeclared_2021_photo_support_gate_before_temporal_models_pass":ready,
      "2021_original_photo_total_by_family":{f:sum(x[f"2021_{f}_photo_n"] for x in joined)
                        for f in ("Lycaenidae","Pieridae","Papilionidae")},
      "2022_original_photo_total_by_family_partial":{"Lycaenidae":sum(x["2022_Lycaenidae_photo_n"] for x in joined),
                       "all_three":sum(x["2022_three_family_photo_total"] for x in joined)},
      "never_butterfly_population_census_or_2021_WWTP_source_absence":True,
      "no_new_photo_or_observation_IDs_and_no_personal_data_downloaded":True,
      "date_utc":datetime.now(timezone.utc).isoformat()
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    if not ready:
        result.update({"status":"STOP_BEFORE_2021_HISTORICAL_BUTTERFLY_SHARE_MODEL",
                      "reason":"2021 family public photorecords insufficient by predeclared matched sample and country criteria"})
        a.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
        print("GESUI_2021_HISTORICAL_SOURCE_SUPPORT_STOP_NO_2022_WWTP_EFFECT_CLAIM",flush=True)
        return
    prevlog=np.asarray([family_logit(x["2021_Lycaenidae_photo_n"],x["2021_three_family_photo_total"]) for x in joined])
    nowlog=np.asarray([family_logit(x["2022_Lycaenidae_photo_n"],x["2022_three_family_photo_total"]) for x in joined])
    resp_corP,resp_corPP=pearsonr(prevlog,nowlog)
    resp_corS,resp_corSP=spearmanr(prevlog,nowlog)
    result["original_2021_vs_2022_same_geographic_cell_photo_logodds_correlations"]={
       "pearson_r":float(resp_corP),"spearman_rho":float(resp_corS),
       "nominal_pearson_p_exploratory":float(resp_corPP),
       "nominal_spearman_p_exploratory":float(resp_corSP),
       "mean_unweighted_2021_family_photo_logit":float(np.mean(prevlog)),
       "mean_unweighted_2022_family_photo_logit":float(np.mean(nowlog)),
       "mean_unweighted_2022_minus_2021_reported_logit":float(np.mean(nowlog-prevlog)),
       "NOT_true_butterfly_population_or_sewage_effect":True}
    for country in COUNTRIES:
        r=[x for x in joined if x["countryCode"]==country]
        y0=np.asarray([family_logit(x["2021_Lycaenidae_photo_n"],x["2021_three_family_photo_total"]) for x in r])
        y1=np.asarray([family_logit(x["2022_Lycaenidae_photo_n"],x["2022_three_family_photo_total"]) for x in r])
        s=np.log1p(np.asarray([float(x[P_SOURCE]) for x in r]))
        rhos,_=spearmanr(s,[float(x["2021_three_family_photo_total"]) for x in r])
        result.setdefault("2021_2022_photo_composition_by_country",{})[country]={
          "cells":len(r),"2021_total_3family_photo_n":sum(x["2021_three_family_photo_total"] for x in r),
          "2022_total_3family_photo_n":sum(x["2022_three_family_photo_total"] for x in r),
          "mean_2022_minus_2021_photo_logit":float(np.mean(y1-y0)),
          "rho_original_2022_reported_Psource_vs_2021_all3family_photo_reporting_total":float(rhos),
          "positive_recording_shift_cells":int(sum(y1>y0)),
          "negative_recording_shift_cells":int(sum(y1<y0)),
        }
    result["fixed_2degree_block_heldout_source_P_models"]={
       "2022_no_historical_butterfly_composition_baseline":cv(joined,2022,False),
       "2022_WITH_2021_butterfly_composition_and_targetgroup_photo_effort_baseline":cv(joined,2022,True),
       "2021_prior_year_photo_reporting_with_2022_future_reported_Psource":cv(joined,2021,False),
    }
    result["interpretation_firewall"]=[
      "Original 2022 EEA reported P source location likely predates 2022, thus 2021 photo response is not no-treatment group",
      "Photograph reporting share by 1degree cell and year is not butterfly abundance or local mud-puddling",
      "WorldCover2021 urban/crop fractions and 2021 plant photo effort do not represent all observer/host/habitat confounders",
      "Additional 2021 Lycaenidae composition can encode stable species habitat, observer bias and genuine ecology, not necessarily just detection",
      "Any positive source effect after prior-year butterfly share is predictor geography only and does not identify P physiology",
      "Prior chocho 22 of25 introduced-host potential resource-positive Lycaenidae and aquatic BQE studies are different estimands"
    ]
    result["status"]="POSTHOC_HISTORICAL_PHOTO_COMPOSITION_AND_SOURCE_PREDICTIVE_INCREMENT_COMPLETED"
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    for name,v in result["fixed_2degree_block_heldout_source_P_models"].items():
        print("GESUI_2021_2022_HISTORICAL_PHOTO_SHARE_SOURCE_P_GEOGRAPHIC_HELDOUT",
          json.dumps({"model":name,"n":v["cells"],"baseline_MSE":v["baseline_MSE"],
                     "source_model_MSE":v["source_P_added_MSE"],
                     "source_P_incremental_improvement_pct":v["EEA2022_reported_P_source_incremental_MSE_improvement_percent"],
                     "country":v["country"]},sort_keys=True),flush=True)
    print("GESUI_2021_2022_LYCAENIDAE_FAMILY_PHOTO_HISTORICAL_SOURCE_DIAGNOSTIC_COMPLETE",
      json.dumps({"same_supported":len(joined),"photo_logit_year_correlation_pearson":float(resp_corP),
                  "photo_logit_year_correlation_spearman":float(resp_corS)},sort_keys=True),flush=True)

if __name__=="__main__":main()
