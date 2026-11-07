#!/usr/bin/env python3
"""Portable 2022 EEA WWTP load ↔ river TN/TP audit with explicit fixed spatial folds.

Use exact source-only and receptor artifacts plus original 2026-11-07 fold
manifest, which makes GroupKFold heldout river blocks reproducible across
NumPy/scikit-learn versions. No live web/API or plant data in the verifier.
"""
from __future__ import annotations
from collections import Counter
from pathlib import Path
import argparse,hashlib,json,math,sys
import numpy as np
import pandas as pd
import sklearn
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from scipy.stats import spearmanr

SOURCE_SHA="79641d984b28fdbc11dbf02e32f32cda2818ec256f219439a1ccfae601da7dc3"
UP={"TP":"eea2022_upstream_P_measured_load_nativeunit",
    "TN":"eea2022_upstream_N_measured_load_nativeunit"}
DOWN={"TP":"eea2022_downstream_P_measured_load_nativeunit_placebo",
      "TN":"eea2022_downstream_N_measured_load_nativeunit_placebo"}
NAMES=("M0_geographic_hydrologic_baseline",
       "M1_plus_2022_upstream_measured_load",
       "M2_plus_2022_downstream_source_placebo",
       "M3_upstream_and_downstream")
DAYS=365.25
YEAR_SEC=DAYS*86400
N_MG_PER_TONNE=1e9
R_TOL=1e-7

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def alike(a,b):
    if not math.isclose(float(a),float(b),rel_tol=0,abs_tol=R_TOL):
        raise ValueError("Mismatch observed in portable CV replay: "+str((a,b)))

def design(data,par,rich=False):
    base=np.column_stack([
        data["WQ.lat"].to_numpy(float),
        data["WQ.lon"].to_numpy(float),
        np.log1p(data.DIS_AV_CMS.to_numpy(float)),
        np.log1p(data.UPLAND_SKM.to_numpy(float))
    ])
    if rich:
        base=np.column_stack((base,data.crop2.to_numpy(float),
                              np.log1p(data.pop_dens.to_numpy(float))))
    up=np.log1p(data[UP[par]].to_numpy(float))
    dn=np.log1p(data[DOWN[par]].to_numpy(float))
    return (base,np.column_stack((base,up)),
            np.column_stack((base,dn)),np.column_stack((base,up,dn)))

def folds_from_manifest(df,manifest):
    labels=df.spatial_3deg_block.astype(str).to_numpy()
    n=len(df)
    folds=[]
    allgroup=set(labels)
    visited=set()
    records=manifest["heldout_block_assignments"]
    # SHA verification freezes groups and fold number, not just total fold count.
    hash_value=hashlib.sha256(json.dumps(records,sort_keys=True,separators=(",",":")).encode()).hexdigest()
    if hash_value!=manifest["heldout_block_assignments_sha256"]:
        raise ValueError("Original fold assignment SHA changed")
    for record in records:
        chosen=set(record["geographic_block_labels"])
        if visited&chosen:raise ValueError("Geographic block appears in two test folds")
        visited|=chosen
        test=np.flatnonzero(np.isin(labels,list(chosen)))
        train=np.flatnonzero(~np.isin(labels,list(chosen)))
        if len(test)!=record["test_n"] or len(test)==0 or len(train)==0:
            raise ValueError("Saved fold size/sample grouping mismatch")
        folds.append((train,test))
    if len(folds)!=5 or visited!=allgroup or sum(len(x[1]) for x in folds)!=n:
        raise ValueError("Fold manifest fails exhaustive non-overlapping test")
    return folds

def replay(data,original,manifest,par,rich=False):
    if manifest["status"]=="SOURCE_SUPPORT_BELOW_CV_GATE":
        if not rich or par!="TN" or len(data.dropna(subset=["crop2","pop_dens"]))!=13:
            raise ValueError("Unexpected untestable rich TN condition")
        return {"status":"NOT_GENERALIZABLE_PREDECLARED_SOURCE_SUPPORT","n":13}
    if rich:data=data.dropna(subset=["crop2","pop_dens"]).copy()
    data=data.reset_index(drop=True)
    y=data.log_real_2022_median_mgL.to_numpy(float)
    if not np.allclose(y,np.log(data.real_2022_station_median_mgL.to_numpy(float)),
                       atol=1e-8,rtol=0):
        raise ValueError("2022 measured concentration log transform mismatch")
    groups=data.spatial_3deg_block.astype(str).to_numpy()
    folds=folds_from_manifest(data,manifest)
    matrices=design(data,par,rich)
    out={}
    for name,mat in zip(NAMES,matrices):
        if not np.isfinite(mat).all():raise ValueError("Missing source covariates")
        pred=np.full(len(y),np.nan)
        for train,test in folds:
            standard=StandardScaler()
            xx=standard.fit_transform(mat[train])
            reg=Ridge(alpha=10.0).fit(xx,y[train])
            pred[test]=reg.predict(standard.transform(mat[test]))
        if not np.isfinite(pred).all():raise ValueError("Missing replay predictions")
        mse=float(np.mean((y-pred)**2))
        orig=float(original["metrics"][name]["heldout_mse"])
        alike(mse,orig)
        out[name]=mse
    if len(data)!=original["n"]:
        raise ValueError("Source sample denominator differs from original report")
    return {"status":"MATCH_EXPLICIT_SAVED_SPATIAL_FOLDS",
            "n":len(data),"blocks":len(set(groups)),
            "fold_sha256":manifest["heldout_block_assignments_sha256"],
            "recalculated_MSE":out,
            "baseline_to_upstream_MSE_percent_improvement":
                 100*(out[NAMES[0]]-out[NAMES[1]])/out[NAMES[0]]}

def source_diagnostic(df,par):
    upstream=df[UP[par]].to_numpy(float)
    q=df.DIS_AV_CMS.to_numpy(float)
    observed=df.real_2022_station_median_mgL.to_numpy(float)
    if (q<=0).any() or (upstream<0).any() or (observed<=0).any():
        raise ValueError("Invalid source/river mass or chemistry input")
    potential=upstream*N_MG_PER_TONNE/(q*1000*YEAR_SEC)
    exposed=upstream>0
    rel=potential[exposed]/observed[exposed]
    if not np.isfinite(rel).all():raise ValueError("Invalid physically mixed ratio")
    return {
       "n_stations":len(df),"n_positive_upstream_reported_load":int(exposed.sum()),
       "reported_source_annual_load_units":"metric tonnes per year",
       "river_mean_Q_source_units":"cubic metres per second",
       "observed_2022_station_median_units":"mg/L",
       "annual_mass_dilution_is_hypothetical_not_measured_fraction":True,
       "fraction_of_all_stations_with_source":float(exposed.mean()),
       "source_annual_load_tonnes_median_among_exposed":
           float(np.median(upstream[exposed])) if exposed.any() else None,
       "hypothetical_fully_mixed_mgL_median_among_exposed":
           float(np.median(potential[exposed])) if exposed.any() else None,
       "hypothetical_fully_mixed_over_observed_station_concentration_ratio_median":
           float(np.median(rel)) if len(rel) else None,
       "hypothetical_ratio_10th_90th_percentiles":
           list(np.quantile(rel,[.1,.9]).astype(float)) if len(rel) else [],
       "stations_with_hypothetical_ratio_gt_10pct":int((rel>.1).sum()),
       "stations_with_hypothetical_ratio_gt_100pct":int((rel>1).sum()),
       "warning":"Not a bound on true seasonal load; uses source-reported subset, 2022 annual masses and river mean Q"
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--outcome",required=True,type=Path)
    p.add_argument("--original-folds",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    if sha(args.source)!=SOURCE_SHA:
        raise ValueError("EEA source-only original station table digest changed")
    prior=json.loads((args.outcome/"EEA2022_ORIGINAL_MEASURED_RECEPTOR_CV_AND_AUDIT.json").read_text())
    foldreport=json.loads(args.original_folds.read_text())
    if foldreport["status"]!="MATCH":
        raise ValueError("Original fold replay was unsuccessful")
    if prior["original_source_only_run"]!=42374105184 or not prior["2022_year_aligned_only"]:
        raise ValueError("Wrong original source/response frozen protocol")
    allresults={}
    table=pd.read_csv(args.source,dtype={"WQ.ID":str},low_memory=False)
    for par in ("TP","TN"):
        data=pd.read_csv(args.outcome/(par+"_2022_reported_wwtp_load_and_original_measured_river_station.csv"),
                         dtype={"WQ.ID":str},low_memory=False)
        if data["WQ.ID"].duplicated().any():raise ValueError("Duplicate river stations in response")
        subset=table[table.nutrient==par][["WQ.ID",UP[par],DOWN[par]]]
        joined=data[["WQ.ID",UP[par],DOWN[par]]].merge(
            subset,on="WQ.ID",suffixes=("_outcome","_source"),
            validate="one_to_one")
        if len(joined)!=len(data):raise ValueError("Receptor stations changed after source-only freeze")
        for field in (UP[par],DOWN[par]):
            if not np.allclose(joined[field+"_outcome"],joined[field+"_source"],atol=1e-5,rtol=0):
                raise ValueError("2022 EEA load source changed after response inspection")
        if len(data)!=prior["results"][par]["original_2022_station_n"]:
            raise ValueError("Original sample size changed")
        if not ((data.actual_2022_observation_count>=4)&
                 (data.distinct_2022_months>=3)).all():
            raise ValueError("Predeclared station quality threshold changed")
        orig=prior["results"][par]
        manifest=foldreport["results"][par]
        aud={
          "primary":replay(data,orig["primary_2022_source_CV"],manifest["primary"],par),
          "landuse":replay(data,orig["crop_population_2022_CV"],manifest["crop_pop"],par,rich=True),
          "physical_source_dilution":source_diagnostic(data,par)
        }
        rho=float(spearmanr(np.log1p(data[UP[par]]),data.log_real_2022_median_mgL).statistic)
        # Original 2022 source values were double precision during fitting.
        # Exported station panels use 10 significant decimal digits; near-tied
        # annual loads can merge, changing rank-correlation by ~1.6e-5 even
        # while full fitted MSEs agree to <1e-7 on frozen exact folds.
        original_rho=orig["raw_source_vs_river_observed_rho"]["spearman_rho"]
        if abs(rho-original_rho)>1e-4:
            raise ValueError("Nontrivial source rank correlation mismatch")
        aud["raw_source_vs_receptor_spearman_rho_verified"]=rho
        aud["rank_rho_difference_due_to_export_precision"]=abs(rho-original_rho)
        allresults[par]=aud
    out={"status":"PASS_CROSS_VERSION_EXPLICIT_SPATIAL_FOLD_REPLAY",
         "source_only_SHA256":SOURCE_SHA,
         "original_fold_audit_sha256":sha(args.original_folds),
         "original_2022_outcome_JSON_SHA256":sha(args.outcome/"EEA2022_ORIGINAL_MEASURED_RECEPTOR_CV_AND_AUDIT.json"),
         "original_2022_run_id":42375542151,
         "original_fold_audit_run_id":42376191822,
         "runtime_version":sklearn.__version__,
         "no_new_source_queries":True,
         "scientific_scope":"source-reported WWTP N/P annual masses and measured 2022 river chemistry, not causally attributable wastewater",
         "results":allresults}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_EEA2022_EXPLICIT_FOLD_PORTABLE_REPLAY",
          json.dumps({"status":out["status"],"sklearn_version":out["runtime_version"],
                      "results":{par:{"n":r["primary"]["n"],
                                       "dilution_potential_median":r["physical_source_dilution"]["hypothetical_fully_mixed_over_observed_station_concentration_ratio_median"],
                                       "upstream_model_MSE_delta_pct":r["primary"]["baseline_to_upstream_MSE_percent_improvement"],
                                       "fold_sha256":r["primary"]["fold_sha256"]}
                                 for par,r in allresults.items()}},
                     sort_keys=True),flush=True)

if __name__=="__main__":
    main()
