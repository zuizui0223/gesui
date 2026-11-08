#!/usr/bin/env python3
"""POSTHOC: 2-degree region-held-out testing of published DE/AT 2022 Lycaenidae
photo reporting composition against measured-reported WWTP source covariates.

No new butterfly data/API queried. Reuses exact frozen 2021 Plantae effort
gate and all 2022 three-family photo counts from prior PR. Explicitly POSTHOC
because first developed after inspecting positive within-sample association.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import defaultdict,Counter
from pathlib import Path
from datetime import datetime,timezone
import numpy as np

SEED=20261108
COUNTRY_SHUFFLES=399
WITHIN_BLOCK_SHUFFLES=199
MIN_COUNTRY_N=10
FIELDS={
 "outfall_count":"public_EEA_both_measured_positive_unique_FW_outfalls",
 "reported_P_t_per_year":"public_EEA_total_reported_P_t_year",
 "reported_N_t_per_year":"public_EEA_total_reported_N_t_year",
}

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def getrows(j):
    old=j["all_fixed_cells_2022_family_photo_source_receipts"]
    return [x for x in old if isinstance(x.get("three_family_total_2022"),int)
            and x["three_family_total_2022"]>=j["min_2022_targetgroup_photos_per_eligible_cell"]]

def safe_z(v,t):
    v=np.asarray(v,dtype=float);t=np.asarray(t,dtype=float)
    m=float(np.mean(v));s=float(np.std(v))
    if s<1e-12:return np.zeros_like(v),np.zeros_like(t)
    return (v-m)/s,(t-m)/s

def build(train,test,field,override=None):
    # Model outcome is independently retained three-family photo composition.
    def col(arr,k,conv=float):return np.asarray([conv(r[k]) for r in arr],float)
    trsrc=col(train,field);tesrc=col(test,field)
    if override is not None:
        trsrc=np.asarray([override[(r["countryCode"],r["lat_floor"],r["lon_floor"])] for r in train],float)
        tesrc=np.asarray([override[(r["countryCode"],r["lat_floor"],r["lon_floor"])] for r in test],float)
    x_train=[np.ones(len(train)),(col(train,"countryCode",lambda x:1 if x=="AT" else 0))]
    x_test=[np.ones(len(test)),(col(test,"countryCode",lambda x:1 if x=="AT" else 0))]
    trplant=np.log1p(col(train,"2021_Plantae_photo_count"))
    teplant=np.log1p(col(test,"2021_Plantae_photo_count"))
    trlat=col(train,"lat_floor")+.5;telat=col(test,"lat_floor")+.5
    trlon=col(train,"lon_floor")+.5;telon=col(test,"lon_floor")+.5
    for a,b in ((trplant,teplant),(trlat,telat),(trlon,telon)):
        aa,bb=safe_z(a,b)
        x_train.append(aa);x_test.append(bb)
    X0=np.column_stack(x_train);T0=np.column_stack(x_test)
    # Training-only country-standardized source for model M1.
    src_a=np.log1p(trsrc);src_b=np.log1p(tesrc)
    country_train=[r["countryCode"] for r in train]
    country_test=[r["countryCode"] for r in test]
    sTrain=np.zeros(len(train));sTest=np.zeros(len(test))
    for c in ("DE","AT"):
        im=np.asarray([v==c for v in country_train])
        jm=np.asarray([v==c for v in country_test])
        vals=src_a[im]
        if len(vals)<2:raise ValueError("No training data for "+c)
        m=float(np.mean(vals));sd=float(np.std(vals))
        if sd<1e-12:raise ValueError("No source count variance in country "+c)
        sTrain[im]=(vals-m)/sd
        sTest[jm]=(src_b[jm]-m)/sd
    X1=np.column_stack([X0,sTrain])
    T1=np.column_stack([T0,sTest])
    return X0,T0,X1,T1

def response(rows):
    return np.asarray([math.log((float(x["Lycaenidae_2022_photo_n"])+.5)/
        (float(x["Pieridae_2022_photo_n"])+float(x["Papilionidae_2022_photo_n"])+.5)) for x in rows])

def pred(X,T,y):
    if np.linalg.matrix_rank(X)<X.shape[1]:
        raise ValueError("Rank-deficient fixed OLS CV design")
    return T@np.linalg.lstsq(X,y,rcond=None)[0]

def cv(rows,field,override=None,return_records=False):
    groups=sorted(set((x["countryCode"],x["geographic_2degree_block"]) for x in rows))
    errors0=[];errors1=[];out=[]
    for (c,b) in groups:
        tr=[x for x in rows if (x["countryCode"],x["geographic_2degree_block"])!=(c,b)]
        te=[x for x in rows if (x["countryCode"],x["geographic_2degree_block"])==(c,b)]
        X0,T0,X1,T1=build(tr,te,field,override=override)
        ytr=response(tr);yte=response(te)
        y0=pred(X0,T0,ytr);y1=pred(X1,T1,ytr)
        for i,r in enumerate(te):
            e0=float((yte[i]-y0[i])**2);e1=float((yte[i]-y1[i])**2)
            errors0.append(e0);errors1.append(e1)
            if return_records:
                out.append({
                    "country":r["countryCode"],
                    "2degree_block":r["geographic_2degree_block"],
                    "one_degree_cell_latfloor":r["lat_floor"],
                    "one_degree_cell_lonfloor":r["lon_floor"],
                    "Lycaenidae_photo_n":r["Lycaenidae_2022_photo_n"],
                    "three_family_photo_n":r["three_family_total_2022"],
                    "log_lycaenidae_to_other2_observed":float(yte[i]),
                    "baseline_pred":float(y0[i]),
                    "source_pred":float(y1[i]),
                    "baseline_sqerr":e0,
                    "source_sqerr":e1,
                })
    a=float(np.mean(errors0));b=float(np.mean(errors1))
    o={
       "heldout_groups_n":len(groups),
       "n_country_1degree_cells":len(rows),
       "baseline_MSE":a,
       "with_upstream_public_EEA_reported_source_MSE":b,
       "heldout_relative_MSE_reduction_percent":100*(a-b)/a if a>0 else None,
       "source_model_better_than_baseline":bool(b<a),
       "baseline_MAE":float(np.mean(np.sqrt(errors0))),
       "with_source_MAE":float(np.mean(np.sqrt(errors1))),
       "number_outcome_selected_hyperparameters":0,
    }
    if return_records:
        o["outcome_blind_fixed_2deg_holdout_predictions"]=out
        o["by_country"]= {}
        for c in ("DE","AT"):
            arr=[x for x in out if x["country"]==c]
            if len(arr)>=MIN_COUNTRY_N:
                mse0=float(np.mean([x["baseline_sqerr"] for x in arr]))
                mse1=float(np.mean([x["source_sqerr"] for x in arr]))
                o["by_country"][c]={"n":len(arr),"MSE_baseline":mse0,
                    "MSE_with_source":mse1,"relative_improvement_percent":100*(mse0-mse1)/mse0}
    return o

def get_override(rows,field,rng,blockwithin=False):
    v={}
    for country in ("DE","AT"):
        groups=defaultdict(list)
        for r in rows:
            if r["countryCode"]==country:
                b=r["geographic_2degree_block"] if blockwithin else "ENTIRE_COUNTRY"
                groups[b].append(r)
        for b,arr in groups.items():
            vals=[float(z[field]) for z in arr]
            vals=list(rng.permutation(vals))
            for r,val in zip(arr,vals):
                v[(r["countryCode"],r["lat_floor"],r["lon_floor"])]=val
    return v

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source-result",required=True,type=Path)
    parser.add_argument("--plant-gate",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    original=json.loads(args.source_result.read_text())
    original_prior=json.loads(args.plant_gate.read_text())
    if not original["2022_photo_share_model_gate"]:
        raise ValueError("Original photo-result did not pass the prespecified support gate")
    if original["source_only_2021_Plantae_gate_sha256"]!=sha(args.plant_gate):
        raise ValueError("Original 2021 botanical effort SHA mismatch: response leakage risk")
    if not original_prior["gate_2021_effort_eligibility_before_2022_butterfly_outcome"]:
        raise ValueError("Original preoutcome nonbutterfly gate absent")
    rows=getrows(original)
    if len(rows)!=44 or dict(Counter(x["countryCode"] for x in rows))!={"DE":29,"AT":15}:
        raise ValueError("No original outcome-dependent site replacement permitted")
    allkeys=[(x["countryCode"],x["lat_floor"],x["lon_floor"]) for x in rows]
    if len(set(allkeys))!=len(allkeys):raise ValueError("Duplicate country×1° cell")
    results={}
    for name,field in FIELDS.items():
        results[name]=cv(rows,field,return_records=True)
        print("GESUI_DE_AT_LYCAENIDAE_FIXED_2DEG_BLOCK_HELDOUT",
          json.dumps({"source":name,**{k:v for k,v in results[name].items()
                  if k in ("n_country_1degree_cells","heldout_groups_n",
                         "baseline_MSE","with_upstream_public_EEA_reported_source_MSE",
                         "heldout_relative_MSE_reduction_percent","source_model_better_than_baseline","by_country")}},sort_keys=True),flush=True)
    field=FIELDS["outfall_count"]
    actual=results["outfall_count"]["heldout_relative_MSE_reduction_percent"]
    rng=np.random.default_rng(SEED)
    x=[]
    for i in range(COUNTRY_SHUFFLES):
        override=get_override(rows,field,rng,False)
        try:
            res=cv(rows,field,override=override)
            x.append(float(res["heldout_relative_MSE_reduction_percent"]))
        except (ValueError,np.linalg.LinAlgError):
            continue
    within=[]
    for i in range(WITHIN_BLOCK_SHUFFLES):
        override=get_override(rows,field,rng,True)
        try:
            res=cv(rows,field,override=override)
            within.append(float(res["heldout_relative_MSE_reduction_percent"]))
        except (ValueError,np.linalg.LinAlgError):continue
    out={
      "schema":"gesui_2022_DE_AT_Lycaenidae_photo_reporting_posthoc_2deg_block_LOBO_validation_v01",
      "posthoc_NOT_preregistered_before_2022_photo_results":True,
      "executed_after_original_in_sample_beta_seen":True,
      "source_frozen_year_2022":True,
      "original_2022_photo_source_result_SHA256":sha(args.source_result),
      "original_2021_nonbutterfly_plant_effort_gate_SHA256":sha(args.plant_gate),
      "source_predictors_fixed_original_EEA2022":FIELDS,
      "model_baseline":"OLS cell log((L+.5)/(Pieridae+Papilionidae+.5)) ~ 2021 Plantae photo effort + lat + lon + country",
      "model_addition":"each separately add within-country standardized log1p(reported WWTP source count or annual measured mass), no tuning",
      "train_only_scaling_per_heldout_2deg_block":True,
      "within_country_source_standardization_train_only":True,
      "source_block_to_fold_geometry_original_prespecified":True,
      "CV_results":results,
      "country_permutation_null": {
         "source_predictor":"original source count only",
         "country_preserving_random_permutation_draws_requested":COUNTRY_SHUFFLES,
         "draws_valid":len(x),"seed":SEED,
         "empirical_percentage_of_random_permutation_gains_ge_observed":
              float(sum(v>=actual for v in x)/len(x)) if x else None,
         "interval_95pct_random_permutation_gains":[float(t) for t in np.quantile(x,[.025,.975])] if len(x)>20 else None,
         "null_is_NOT_randomized_WWTP_treatment":True},
      "within_2degree_block_permutation_stress": {
         "draws_requested":WITHIN_BLOCK_SHUFFLES,"valid":len(within),
         "empirical_fraction_random_gain_ge_observed":float(sum(v>=actual for v in within)/len(within)) if within else None,
         "original_source_2deg_pattern_preserved":True,
         "null_may_be_uninformative_due_to_singletons":True},
      "actual_2022_source_heldout_improvement_sign_predeclared_here_posthoc":bool(actual>0),
      "in_sample_positive_original_is_retained_separately":True,
      "not_butterfly_numbers_or_chemical_wastewater_effect":True,
      "remaining_confounders":["urban/cropland fractions","direct butterfly target-group observer effort",
              "varying partial land area of 1degree grid","2022 original source facilities reporting patterns",
              "host plant richness/quality","air/soil sodium/protein nutrient exposure","spatial photo geolocation errors"],
      "no_photo_or_individual_person_data_downloaded":True,
      "timestamp_utc":datetime.now(timezone.utc).isoformat(),
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_DE_AT_LYCAENIDAE_POSTHOC_GEOGRAPHIC_TRANSFER_VALIDATION",
      json.dumps({"primary_actual_outofblock_gain_pct":actual,
        "country_null_tail":out["country_permutation_null"]["empirical_percentage_of_random_permutation_gains_ge_observed"],
        "blocknull_tail":out["within_2degree_block_permutation_stress"]["empirical_fraction_random_gain_ge_observed"],
        "analytical_status":"POSTHOC"}),flush=True)

if __name__=="__main__":main()
