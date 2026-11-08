#!/usr/bin/env python3
"""POSTHOC New-year 2023 vs original 2022 Lycaenidae family photo geography.

All 2023 photo counts were separately frozen before this script's geoCV.
Year 2022 EEA *reported wastewater source mass* is a spatial covariate
for 2023 only, NOT measured 2023 discharge or an experimentally manipulated
sewage exposure. No individual butterflies, adult puddling or host quality.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr, pearsonr

MIN_MATCHED=30
MIN_DE=18
MIN_AT=9
MIN_BLOCKS=12
SOURCE_P_KEY="registered_P_t_year"
COUNTRIES=("DE","AT")

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for data in iter(lambda:f.read(1048576),b""):h.update(data)
    return h.hexdigest()

def idkey(x):return (x["countryCode"],int(float(x["lat_floor"])),int(float(x["lon_floor"])))
def grp(x):return x["countryCode"]+"|"+str(x["country_2deg_block"])
def logit(l,n):return math.log((float(l)+.5)/(float(n)-float(l)+.5))
def get_y(x,year):
    return logit(x[f"{year}_Lycaenidae_photo_n"],x[f"{year}_three_family_photo_total"])

def arr(data,fn):return np.asarray([fn(r) for r in data],float)
def build(train,test,add_hist):
    fields={
      "country_is_AT":lambda r:1. if r["countryCode"]=="AT" else 0.,
      "2021_Plantae_photos":lambda r:math.log1p(float(r["2021_Plantae_photo_n"])),
      "latitude":lambda r:float(r["lat_floor"])+.5,
      "longitude":lambda r:float(r["lon_floor"])+.5,
      "WorldCover2021_built":lambda r:float(r["WC2021_512_built_land_fraction"]),
      "WorldCover2021_crop":lambda r:float(r["WC2021_512_crop_land_fraction"]),
      "WorldCover2021_country_land_area":lambda r:math.log1p(float(r["WC2021_512_country_terrestrial_area_km2"])),
    }
    if add_hist:
        fields.update({
          "previous_2022_Lycaenidae_logit":lambda r:get_y(r,2022),
          "previous_2022_butterfly_family_total_photos":lambda r:math.log1p(float(r["2022_three_family_photo_total"])),
        })
    X=[np.ones(len(train))];T=[np.ones(len(test))]
    for name,fn in fields.items():
        a=arr(train,fn);b=arr(test,fn)
        if name=="country_is_AT":
            X.append(a);T.append(b);continue
        m=float(np.mean(a));s=float(np.std(a))
        if s<=1e-12:raise ValueError("No training contrast in fixed feature "+name)
        X.append((a-m)/s);T.append((b-m)/s)
    X0=np.column_stack(X);T0=np.column_stack(T)
    # Within-country training-fold mean/sd; don't leak heldout source predictor
    trainP=arr(train,lambda r:math.log1p(float(r[SOURCE_P_KEY])))
    testP=arr(test,lambda r:math.log1p(float(r[SOURCE_P_KEY])))
    z=np.zeros(len(train),float);q=np.zeros(len(test),float)
    for c in COUNTRIES:
        a=np.asarray([r["countryCode"]==c for r in train])
        b=np.asarray([r["countryCode"]==c for r in test])
        if int(a.sum())<3:raise ValueError("Not enough training source tiles in country "+c)
        v=trainP[a];m=float(np.mean(v));sd=float(np.std(v))
        if sd<=1e-12:raise ValueError("Within-country P mass source has no train variance")
        z[a]=(v-m)/sd;q[b]=(testP[b]-m)/sd
    return X0,T0,np.column_stack([X0,z]),np.column_stack([T0,q])

def predict(X,T,y):
    if np.linalg.matrix_rank(X)<X.shape[1]:
        raise ValueError("Rank deficient 2degree heldout source design")
    return T@np.linalg.lstsq(X,y,rcond=None)[0]

def cv(rows,year,add_hist=False):
    if year==2022 and add_hist:raise ValueError("Cannot use 2022 predictor to predict itself")
    predicted=[]
    for group in sorted({grp(x) for x in rows}):
        tr=[x for x in rows if grp(x)!=group];te=[x for x in rows if grp(x)==group]
        X0,T0,X1,T1=build(tr,te,add_hist)
        yy=arr(tr,lambda r:get_y(r,year));yt=arr(te,lambda r:get_y(r,year))
        m0=predict(X0,T0,yy);m1=predict(X1,T1,yy)
        for q,y,p0,p1 in zip(te,yt,m0,m1):
            predicted.append({"country":q["countryCode"],"spatial_group":group,
             "country_1degree_source_latfloor":q["lat_floor"],
             "country_1degree_source_lonfloor":q["lon_floor"],
             "observed_lycaenid_3family_photo_logit":float(y),
             "baseline_prediction":float(p0),"2022_annual_reported_P_source_plus_prediction":float(p1),
             "baseline_sqerr":float((y-p0)**2),
             "source_sqerr":float((y-p1)**2)})
    mse0=float(np.mean([q["baseline_sqerr"] for q in predicted]))
    mse1=float(np.mean([q["source_sqerr"] for q in predicted]))
    out={
      "target_family_photo_reporting_year":year,
      "historical_2022_butterfly_photography_baseline_added":add_hist,
      "sample_n":len(rows),
      "heldout_region_2deg_blocks":len(set(q["spatial_group"] for q in predicted)),
      "MSE_baseline":mse0,
      "MSE_with_original_reported_2022_P_source":mse1,
      "source_P_relative_MSE_improvement_pct":100*(mse0-mse1)/mse0,
      "source_2022_year_is_not_same_as_2023_butterfly_behavior":year==2023,
      "per_cell_2deg_geographic_heldout_predictions":predicted,
      "bycountry":{},
      "not_actual_abundance_or_WWTP_causality":True,
    }
    for c in COUNTRIES:
        p=[q for q in predicted if q["country"]==c]
        aa=float(np.mean([q["baseline_sqerr"] for q in p]))
        bb=float(np.mean([q["source_sqerr"] for q in p]))
        out["bycountry"][c]={"n":len(p),"MSE_baseline":aa,"MSE_source":bb,
           "improvement_pct":100*(aa-bb)/aa}
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-2023-photography",type=Path,required=True)
    p.add_argument("--original-2022-landcover",type=Path,required=True)
    p.add_argument("--optional-2021-family-history",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    photo2023=json.loads(a.source_2023_photography.read_text())
    prior=json.loads(a.original_2022_landcover.read_text())
    hist2021=json.loads(a.optional_2021_family_history.read_text())
    if photo2023["schema"]!="gesui_2023_DE_AT_public_butterfly_family_photo_source_only_original43_WorldCover2021_v01":
        raise ValueError("Not original 2023 family source")
    if photo2023["original_2022_WorldCover_land_photo_source_sha256"]!=sha(a.original_2022_landcover):
        raise ValueError("2023 photo sample is not original immutable 2022 WorldCover study")
    if not photo2023["2023_outcome_year_source_gate_pass"]:
        raise ValueError("2023 photo support gate failed; don't model P source")
    if hist2021["original_2021_WorldCover2022_photo_source_full_SHA256"]!=sha(a.original_2022_landcover):
        raise ValueError("2021 photo source historical panel does not match exact original 43 cells")
    orig=prior["all_original_same_sample_cells_2022_photo_counts_and_WorldCover_categorical_2021_covariates"]
    target={idkey(r):r for r in photo2023["source_2023_all_original43_butterfly_photo_count_and_receipts"]}
    ancient={idkey(r):r for r in hist2021["source_2021_same43_cell_historical_butterfly_photo_records"]}
    if len(orig)!=43 or len(target)!=43 or len(ancient)!=43:
        raise ValueError("New year data cannot replace original panel cells")
    selected=[]
    for x in orig:
        key=idkey(x)
        p23=target[key]
        if not p23["2023_photo_support_ge20"]:continue
        r=dict(x)
        r["2022_three_family_photo_total"]=x["2022_total_3family_photos"]
        r["2022_Lycaenidae_photo_n"]=x["2022_Lycaenidae_photo_n"]
        r["2023_three_family_photo_total"]=p23["2023_three_family_photo_total"]
        r["2023_Lycaenidae_photo_n"]=p23["2023_Lycaenidae_photo_n"]
        anc=ancient[key]
        r["2021_three_family_photo_total"]=anc["2021_three_family_photo_total"]
        r["2021_Lycaenidae_photo_n"]=anc["2021_Lycaenidae_photo_n"]
        r["2021_historical_photo_support_min15"]=anc["2021_photo_support_min15"]
        selected.append(r)
    by=Counter(r["countryCode"] for r in selected)
    blocks=set(grp(r) for r in selected)
    ready=(len(selected)>=MIN_MATCHED and by["DE"]>=MIN_DE and by["AT"]>=MIN_AT and len(blocks)>=MIN_BLOCKS)
    result={
      "schema":"gesui_2023_new_year_vs_2022_EEA_reported_P_photographic_Lycaenidae_source_temporal_falsifier_v01",
      "POSTHOC_2023_falsifier_after_seeing_2022_butterfly_source_prediction":True,
      "original_2022_WorldCover43_source_SHA256":sha(a.original_2022_landcover),
      "2023_family_photo_source_original43_SHA256":sha(a.source_2023_photography),
      "2021_family_photo_source_original43_SHA256":sha(a.optional_2021_family_history),
      "original_WWTP_reported_source_cycle_2022_not_2023":True,
      "same_2023_photo_supported_1deg_cells_n":len(selected),
      "same_2023_photo_supported_by_country":dict(by),
      "same_2023_photo_supported_country_2deg_blocks":len(blocks),
      "fixed_2023_newyear_BMS_photo_composition_model_gate":ready,
      "original_2021_prior_year_family_source_supported":hist2021["2021_photo_support_before_historical_comparison_model_gate"],
      "original_all_2023_photo_source_query_source_only_count":photo2023["2023_photo_original_API_queries"],
      "no_new_photo_images_or_observer_personal_information_retrieved":True,
      "no_2023_reported_annual_P_mass_independently_verified":True,
      "photo_shares_not_population_censuses":True,
      "utc":datetime.now(timezone.utc).isoformat(),
    }
    if not ready:
        result.update({"status":"STOP_BEFORE_NEWYEAR_SOURCE_P_MODEL",
                       "reason":"original 2023 source support insufficient; no output fit"})
    else:
        a22=arr(selected,lambda r:get_y(r,2022))
        a23=arr(selected,lambda r:get_y(r,2023))
        pear,_=pearsonr(a22,a23);rho,_=spearmanr(a22,a23)
        result["2022_2023_photo_logit_stability"]={
         "pearson_r":float(pear),"spearman_rho":float(rho),
         "not_true_butterfly_population_stability":True}
        result["heldout_P_source_spatial_photo_reporting_models_fixed_same_sample"]={
         "original_2022_photo_reporting_no_history":cv(selected,2022,False),
         "new_2023_photo_reporting_no_history":cv(selected,2023,False),
         "new_2023_photo_reporting_with_2022_butterfly_history":cv(selected,2023,True),
        }
        both2021=[r for r in selected if r["2021_historical_photo_support_min15"]]
        result["three_year_common_2021_2022_2023_photo_source_cells"]={
         "n":len(both2021),"DE":sum(r["countryCode"]=="DE" for r in both2021),
         "AT":sum(r["countryCode"]=="AT" for r in both2021),
         "2deg_blocks":len(set(grp(r) for r in both2021)),
         "2021_original_total_photo_Lycaenidae":sum(int(r["2021_Lycaenidae_photo_n"]) for r in both2021),
         "2021_original_total_3family_photo":sum(int(r["2021_three_family_photo_total"]) for r in both2021),
         "2022_original_total_photo_Lycaenidae":sum(int(r["2022_Lycaenidae_photo_n"]) for r in both2021),
         "2022_original_total_3family_photo":sum(int(r["2022_three_family_photo_total"]) for r in both2021),
         "2023_original_total_photo_Lycaenidae":sum(int(r["2023_Lycaenidae_photo_n"]) for r in both2021),
         "2023_original_total_3family_photo":sum(int(r["2023_three_family_photo_total"]) for r in both2021),
        }
        result["status"]="2023_NEWYEAR_GEOGRAPHIC_PSOURCE_PHOTO_REPORTING_POSTHOC_FALSIFIER_COMPLETED"
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    if ready:
        models=result["heldout_P_source_spatial_photo_reporting_models_fixed_same_sample"]
        print("GESUI_2023_POSTHOC_LYCAENIDAE_NEWYEAR_P_SOURCE_HOLDOUT",json.dumps({
          k:{"n":v["sample_n"],"base_MSE":v["MSE_baseline"],
              "source_MSE":v["MSE_with_original_reported_2022_P_source"],
              "source_gain_percent":v["source_P_relative_MSE_improvement_pct"],
              "bycountry":v["bycountry"]} for k,v in models.items()},sort_keys=True),flush=True)
        print("GESUI_2023_PRIOR_YEAR_COMPOSITION_STABILITY",
          json.dumps(result["2022_2023_photo_logit_stability"],sort_keys=True),flush=True)
    else:
        print("GESUI_2023_SOURCE_SUPPORT_STOP_NO_P_OUTCOME_MODEL",json.dumps(result),flush=True)
if __name__=="__main__":main()
