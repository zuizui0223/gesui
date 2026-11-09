#!/usr/bin/env python3
"""Measured 2022 WISE6 river TP ~ original strictly-upstream EEA2022 WWTP P,
negative control strictly DOWNSTREAM WWTP P, original HydroRIVERS directed
reaches, country/geography/river-flow/2021 ESA WorldCover adjustments.

POSTHOC ecological *environmental source-to-receptor observational test*.
Not a Lycaenidae ecology effect or causal nutrient delivery experiment.
Source-only 170 station/168 independent reaches was Git-frozen FIRST.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

EXPECTED_SOURCE_BLOB="420bab68988ca2586177993d916718f20fd1fa88"
EXPECTED_NUMERIC_BLOB="c18b076c74a9bd0ff9cd966c3960e8f35c9fc360"
EXPECTED_EEA_WISE6_ZIP="2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2"
MODELS={
 "strictly_upstream_registered_P_t_year":"original_EEA2022_upstream_registered_P_t_per_year",
 "strictly_upstream_registered_N_t_year":"original_EEA2022_upstream_registered_N_t_per_year",
 "strictly_upstream_qualified_WWTP_count":"original_EEA2022_upstream_qualified_source_n",
 "strictly_downstream_registered_P_t_year_PLACEBO":"original_EEA2022_downstream_placebo_P_t_per_year",
 "strictly_upstream_P_over_original_model_flow_mg_L_potential_proxy":"original_load_per_modeled_flow_proxy_mg_P_L",
}
BOOT_N=999
COUNTRY_PERM=199
WITHIN_BLOCK_PERM=199
SEED=20261109
MIN_POS=20
MIN_ZERO=20
MIN_POS_PER_COUNTRY=5
MIN_GROUPS=5

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def blob(path):
    b=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def sample_key(z):
    return z["original_country_scheme_site_SHA256"]

def blk(z):
    return z["countryCode"]+"|"+z["original_country_2degree_block"]

def mean(vals):
    return float(np.mean(np.asarray(vals,float)))

def make_reach_panel(stations):
    grouped=defaultdict(list)
    for x in stations:grouped[int(x["original_HydroRIVERS_HYRIV_ID"])].append(x)
    reaches=[]
    for rid,arr in sorted(grouped.items()):
        x=arr[0]
        for v in arr[1:]:
            for q in (
                 "countryCode","original_country_2degree_block",
                 "original_EEA2022_upstream_qualified_source_n",
                 "original_EEA2022_upstream_registered_P_t_per_year",
                 "original_EEA2022_downstream_placebo_P_t_per_year"):
                if v[q]!=x[q]:raise ValueError("Same original HydroRIVERS reach has conflicting location/source "+q)
        row=dict(x)
        # Only site-level TP medians are collapsed, not individual samples.
        row["original_receptor_site_count_same_reach"]=len(arr)
        row["measured_WISE6_2022_TP_reach_median_of_site_annual_medians_mgP_L"]=float(
            np.median([r["measured_WISE6_2022_TP_site_annual_median_mgP_L"] for r in arr]))
        for k in ("lat","lon","country_ESAWC2021_built_fraction",
                  "country_ESAWC2021_crop_fraction","country_ESAWC2021_terrestrial_area_km2"):
            row[k]=mean([r[k] for r in arr])
        reaches.append(row)
    return reaches

def y(arr):
    return np.log1p(np.asarray([float(r["measured_WISE6_2022_TP_reach_median_of_site_annual_medians_mgP_L"])
               for r in arr],float))

def geoblock_sources(rows):
    return sorted(set(blk(r) for r in rows))

def zscale(a,b,name):
    a=np.asarray(a,float);b=np.asarray(b,float)
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("Nonfinite train/test original environmental feature: "+name)
    sd=float(a.std())
    if sd<1e-12:raise ValueError("Train feature near zero variance: "+name)
    mu=float(a.mean())
    return (a-mu)/sd,(b-mu)/sd

def make_design(train,test,source_name,override=None):
    feats=[
      ("lat",lambda r:float(r["lat"])),
      ("lon",lambda r:float(r["lon"])),
      ("log1p_modeled_HydroRIVERS_flow_cms",
       lambda r:math.log1p(float(r["original_HydroRIVERS_modeled_flow_cms"]))),
      ("log1p_HydroRIVERS_upstream_basin_km2",
       lambda r:math.log1p(float(r["original_HydroRIVERS_modelled_upstream_catchment_km2"]))),
      ("independent_2021_WorldCover_built",
       lambda r:float(r["country_ESAWC2021_built_fraction"])),
      ("independent_2021_WorldCover_crop",
       lambda r:float(r["country_ESAWC2021_crop_fraction"])),
      ("log1p_country_terrestrial_area_km2",
       lambda r:math.log1p(float(r["country_ESAWC2021_terrestrial_area_km2"]))),
    ]
    X=[np.ones(len(train)),np.asarray([r["countryCode"]=="AT" for r in train],float)]
    T=[np.ones(len(test)),np.asarray([r["countryCode"]=="AT" for r in test],float)]
    for k,fn in feats:
        a=[fn(r) for r in train];b=[fn(r) for r in test]
        x,t=zscale(a,b,k);X.append(x);T.append(t)
    base=np.column_stack(X);base_test=np.column_stack(T)
    field=MODELS[source_name]
    def source(r):
        if override is None:
            return float(r[field]) if r[field] is not None else 0.
        return float(override[r["original_HydroRIVERS_HYRIV_ID"]])
    rawtr=np.asarray([math.log1p(source(r)) for r in train])
    rawte=np.asarray([math.log1p(source(r)) for r in test])
    ztr=np.zeros(len(train));zte=np.zeros(len(test))
    for c in ("DE","AT"):
        mtr=np.asarray([r["countryCode"]==c for r in train])
        mte=np.asarray([r["countryCode"]==c for r in test])
        if mtr.sum()<8:raise ValueError("Train country needs 8 or more independent river reaches")
        x=rawtr[mtr];s=float(x.std());mu=float(x.mean())
        if s<=1e-12:raise ValueError("No within-country training WWTP source contrast: "+source_name+" "+c)
        ztr[mtr]=(x-mu)/s;zte[mte]=(rawte[mte]-mu)/s
    return base,base_test,np.column_stack([base,ztr]),np.column_stack([base_test,zte])

def ols(X,Y,T):
    if np.linalg.matrix_rank(X)!=X.shape[1]:
        raise ValueError("Fixed geoblock training OLS rank deficient; no posthoc source replacement")
    return T@np.linalg.lstsq(X,Y,rcond=None)[0]

def blockcv(rows,name,override=None,return_predictions=False):
    blocks=geoblock_sources(rows)
    records=[]
    for b in blocks:
        train=[r for r in rows if blk(r)!=b]
        test=[r for r in rows if blk(r)==b]
        X0,T0,X1,T1=make_design(train,test,name,override)
        ytr=y(train);yt=y(test)
        pred0=ols(X0,ytr,T0);pred1=ols(X1,ytr,T1)
        for r,actual,base,src in zip(test,yt,pred0,pred1):
            records.append({
             "countryCode":r["countryCode"],
             "source_original_HydroRIVERS_reach":int(r["original_HydroRIVERS_HYRIV_ID"]),
             "original_country_2degree_block":blk(r),
             "source_registered_WWTP_upstream_n":int(r["original_EEA2022_upstream_qualified_source_n"]),
             "source_registered_WWTP_upstream_P_t_year":float(r["original_EEA2022_upstream_registered_P_t_per_year"]),
             "source_registered_WWTP_downstream_placebo_P_t_year":float(r["original_EEA2022_downstream_placebo_P_t_per_year"]),
             "original_WISE6_reach_2022_TP_mgP_L":float(r["measured_WISE6_2022_TP_reach_median_of_site_annual_medians_mgP_L"]),
             "actual_y_log1p_TP":float(actual),
             "baseline_pred":float(base),"added_source_pred":float(src),
             "baseline_squared_error":float((actual-base)**2),
             "source_squared_error":float((actual-src)**2),
             "baseline_absolute_error":float(abs(actual-base)),
             "source_absolute_error":float(abs(actual-src))
             })
    mse0=mean([r["baseline_squared_error"] for r in records])
    mse1=mean([r["source_squared_error"] for r in records])
    out={
     "original_2022_TP_independent_HydroRIVERS_reaches":len(rows),
     "country2degree_geographic_holdout_blocks":len(blocks),
     "baseline_MSE":mse0,
     "augmented_MSE":mse1,
     "incremental_2degree_block_heldout_MSE_improvement_percent":100*(mse0-mse1)/mse0,
     "baseline_MAE":mean([r["baseline_absolute_error"] for r in records]),
     "augmented_MAE":mean([r["source_absolute_error"] for r in records]),
     "source_predictor":name,
     "train_only_country_standardization":True,
     "no_outcome_based_source_selection_or_fold_tuning":True,
     "country_source_heterogeneity":{},
    }
    for country in ("DE","AT"):
        sub=[r for r in records if r["countryCode"]==country]
        m0=mean([r["baseline_squared_error"] for r in sub])
        m1=mean([r["source_squared_error"] for r in sub])
        out["country_source_heterogeneity"][country]={
         "n_independent_original_river_reaches":len(sub),
         "baseline_MSE":m0,"source_MSE":m1,
         "relative_MSE_improvement_percent":100*(m0-m1)/m0}
    if return_predictions:
        out["original_reach_by_reach_heldout_prediction_rows"]=records
    return out

def cluster_error_ci(rows):
    from collections import defaultdict
    grouped=defaultdict(list)
    for r in rows:grouped[r["original_country_2degree_block"]].append(r)
    rng=np.random.default_rng(SEED)
    blocks=list(grouped)
    nums=[]
    for _ in range(BOOT_N):
        keys=rng.choice(blocks,size=len(blocks),replace=True)
        sample=[x for b in keys for x in grouped[b]]
        b=mean([r["baseline_squared_error"] for r in sample])
        src=mean([r["source_squared_error"] for r in sample])
        if b>0:nums.append(100*(b-src)/b)
    return {"n":len(nums),"seed":SEED,"descriptive_95pct_prediction_error_gain_interval":
       [float(x) for x in np.quantile(nums,[0.025,0.975])],
       "not_randomized_intervention_uncertainty":True}

def shuffle(rows,seed,withinblock):
    rng=np.random.default_rng(seed)
    groups=defaultdict(list)
    for r in rows:
        group=(r["countryCode"],r["original_country_2degree_block"]) if withinblock else (r["countryCode"],"all")
        groups[group].append(r)
    out={}
    for k,z in groups.items():
        vs=rng.permutation([r["original_EEA2022_upstream_registered_P_t_per_year"] for r in z])
        for r,v in zip(z,vs):
            out[int(r["original_HydroRIVERS_HYRIV_ID"])]=float(v)
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-source-frozen-panel",required=True,type=Path)
    ap.add_argument("--original-numeric-TP",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    if blob(a.original_source_frozen_panel)!=EXPECTED_SOURCE_BLOB:
        raise ValueError("Pre-directed numerical source panel not Git blob frozen")
    if blob(a.original_numeric_TP)!=EXPECTED_NUMERIC_BLOB:
        raise ValueError("Original observed quantified WISE6 TP source file Git blob changed")
    source=json.loads(a.original_source_frozen_panel.read_text())
    numeric=json.loads(a.original_numeric_TP.read_text())
    if not source["strict_original_source_support_gate_20positive_20zero_5percountry_15reaches_5blocks"]:
        raise ValueError("Original source-only physical support failed, cannot access environmental outcome")
    if numeric["original_2022_EEA2026_WISE6_DisaggregatedData_ZIP_verified"]["original_archive_sha256"]!=EXPECTED_EEA_WISE6_ZIP:
        raise ValueError("Measured river TP original source version changed")
    sites=numeric["original_source_187_sites_all_numeric_and_censor_QC"]
    lookup={r["original_station_siteKey_scheme_country_SHA256"]:r for r in sites}
    if len(sites)!=187 or len(lookup)!=187:
        raise ValueError("Original source 187 site-year quality panel changed")
    station_rows=[]
    for r in source["all_qc_stations_with_original_directed_source_and_2021_WorldCover"]:
        x=lookup.get(r["original_country_scheme_site_SHA256"])
        if x is None or not x["numeric_site_support_ge4days_ge3months_LOQfraction_lehalf"]:
            raise ValueError("Original 170-source panel not aligned with unchanged numeric QC eligibility")
        tp=x["original_EEA2022_TP_station_median_of_daily_medians_mg_P_per_L"]
        if not isinstance(tp,(int,float)) or not math.isfinite(tp) or tp<0:
            raise ValueError("Original real river TP numeric value absent or invalid at source-selected station")
        station_rows.append({**r,"measured_WISE6_2022_TP_site_annual_median_mgP_L":float(tp)})
    if len(station_rows)!=170:raise ValueError("Original source-frozen 170 monitoring stations changed")
    # Collapse 2 duplicated Austrian same-reach observations; no organism/site
    # pseudoreplication. Source and flow are identical by reach.
    panel=make_reach_panel(station_rows)
    if len(panel)!=168:raise ValueError("Original 170 stations did not collapse to 168 river reaches")
    source_positive=[r for r in panel if r["original_EEA2022_upstream_qualified_source_n"]>0]
    source_zero=[r for r in panel if r["original_EEA2022_upstream_qualified_source_n"]==0]
    countries=Counter(r["countryCode"] for r in panel)
    if min(len(source_positive),len(source_zero))<20 or min(sum(r["countryCode"]==c for r in source_positive) for c in ("DE","AT"))<5:
        raise ValueError("Reach-wise eligible source/support gate failed")
    models={}
    for k in MODELS:
        models[k]=blockcv(panel,k,return_predictions=True)
        models[k]["country2degree_prediction_error_bootstrap"]=cluster_error_ci(
            models[k]["original_reach_by_reach_heldout_prediction_rows"])
        print("GESUI_2022_OBSERVED_RIVER_TP_STRICT_UPSTREAM_OR_DOWNSTREAM_NEGCONTROL",
          json.dumps({"source":k,"reaches":len(panel),
             "MSE_baseline":models[k]["baseline_MSE"],
             "MSE_plus_source":models[k]["augmented_MSE"],
             "MSE_percent_improvement":models[k]["incremental_2degree_block_heldout_MSE_improvement_percent"],
             "bycountry":models[k]["country_source_heterogeneity"]},sort_keys=True),flush=True)
    raw={}
    for c in ("DE","AT"):
        raw[c]={}
        for exposed,title in ((True,"registered_strict_upstream_WWTP_positive"),
                              (False,"zero_qualifying_registered_upstream_WWTP")):
            rr=[r for r in panel if (r["countryCode"]==c and (r["original_EEA2022_upstream_qualified_source_n"]>0)==exposed)]
            vals=[r["measured_WISE6_2022_TP_reach_median_of_site_annual_medians_mgP_L"] for r in rr]
            raw[c][title]={"n_original_reaches":len(rr),
               "median_actual_TP_mgP_L":float(np.median(vals)) if vals else None,
               "range_actual_TP_mgP_L":[float(min(vals)),float(max(vals))] if vals else None}
    sourceP=models["strictly_upstream_registered_P_t_year"]
    observed=sourceP["incremental_2degree_block_heldout_MSE_improvement_percent"]
    randomizations={}
    for groupmode,num,seed in (("within_country",COUNTRY_PERM,SEED+10),
                                ("within_country2degree_region",WITHIN_BLOCK_PERM,SEED+20)):
        rr=np.random.default_rng(seed)
        vals=[];errors=Counter()
        for i in range(num):
            groups=defaultdict(list)
            for r in panel:
                z=(r["countryCode"],r["original_country_2degree_block"]) if groupmode=="within_country2degree_region" else (r["countryCode"],"all")
                groups[z].append(r)
            override={}
            for _,v in groups.items():
                perm=rr.permutation([r["original_EEA2022_upstream_registered_P_t_per_year"] for r in v])
                for x,p in zip(v,perm):override[x["original_HydroRIVERS_HYRIV_ID"]]=float(p)
            try:
                gain=blockcv(panel,"strictly_upstream_registered_P_t_year",override=override)["incremental_2degree_block_heldout_MSE_improvement_percent"]
                vals.append(gain)
            except (ValueError,np.linalg.LinAlgError) as e:
                errors[type(e).__name__]+=1
        randomizations[groupmode]={
          "fixed_permutations":num,"valid":len(vals),"errors":dict(errors),
          "null_median_improvement_percent":float(np.median(vals)) if vals else None,
          "null_95pct_improvement_range":[float(v) for v in np.quantile(vals,[.025,.975])] if vals else None,
          "fraction_random_source_geography_gain_ge_original_with_plusone":(1+sum(x>=observed for x in vals))/(1+len(vals)) if vals else None,
          "not_causal_WWTP_assignment_permutation":True}
    result={
      "schema":"gesui_WISE6_2022_real_riverTP_168_independent_HydroRIVERS_reaches_strict_upstream_EEA_WWTP_P_spatial_source_test_v01",
      "timestamp_UTC":datetime.now(timezone.utc).isoformat(),
      "POSTHOC_source_receptor_test_after_photo_P_geographical_positive":True,
      "pre_numeric_source_only_panel_exact_Git_BLOB_SHA1":EXPECTED_SOURCE_BLOB,
      "original_WISE6_2022_187_numeric_TP_exact_Git_BLOB_SHA1":EXPECTED_NUMERIC_BLOB,
      "source_n_original_2022_TP_qualified_stations":170,
      "source_n_independent_HydroRIVERS_reaches":len(panel),
      "unique_original_country2degree_validation_blocks":len(geoblock_sources(panel)),
      "country_counts":dict(countries),
      "primary_directional_upstream_and_downstream_reported_source_Flow_land_geo_CV":models,
      "raw_observed_2022_riverTP_upstream_vs_no_qualifying_registered_source_descriptive":raw,
      "original_reported_upstream_P_source_country2degree_spatial_source_shuffles":randomizations,
      "all_original_reach_level_real_TP_2022_value_and_exposure_sources":[{
        "original_reach_HYRIV_ID":r["original_HydroRIVERS_HYRIV_ID"],
        "original_country":r["countryCode"],
        "original_country2deg_block":blk(r),
        "original_monitoring_sites_per_reach":r["original_receptor_site_count_same_reach"],
        "registered_strict_upstream_outfall_count_50km":r["original_EEA2022_upstream_qualified_source_n"],
        "registered_strict_upstream_P_t_year":r["original_EEA2022_upstream_registered_P_t_per_year"],
        "registered_strict_downstream_P_t_year":r["original_EEA2022_downstream_placebo_P_t_per_year"],
        "measured_reach_year_TP_mgP_L":r["measured_WISE6_2022_TP_reach_median_of_site_annual_medians_mgP_L"],
      } for r in panel],
      "original_source_same_reach_outfall_stations_relative_order_unknown_excluded":True,
      "original_hydro_from_thirdparty_GeoParquet_not_verified_against_official":True,
      "original_WWTP_source_load_reported_t_year_not_station_specific_chemical_dose":True,
      "no_site_source_selected_on_numeric_TP_magnitude":True,
      "no_iNaturalist_or_Lycaenidae_outcome_in_this_model":True,
      "no_proof_of_sewage_causality_or_adult_puddling_behavior":True,
      "unreviewed_independent_EEA_resultObservationStatus_expert_flags_remain_publication_blocker":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_2022_WISE6_TRUE_TP_STRICT_DIRECTED_WWTP_UPSTREAM_AND_DOWNSTREAM_CV_FINAL",
          json.dumps({"reaches":len(panel),"sites":len(station_rows),
           "country":dict(countries),"upstream":len(source_positive),"no_registered_upstream":len(source_zero),
           "upstream_P_MSE_gain_percent":models["strictly_upstream_registered_P_t_year"]["incremental_2degree_block_heldout_MSE_improvement_percent"],
           "downstream_PLACEBO_MSE_gain_percent":models["strictly_downstream_registered_P_t_year_PLACEBO"]["incremental_2degree_block_heldout_MSE_improvement_percent"],
           "country_permuted_P_tail":randomizations["within_country"]["fraction_random_source_geography_gain_ge_original_with_plusone"],
           "no_sewage_causality":True},sort_keys=True),flush=True)
if __name__=="__main__":main()
