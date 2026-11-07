#!/usr/bin/env python3
"""Independent zero-network audit of the original EEA 2022 measured TN/TP CV.

A source-run and a separate real 2022 receptor-run are read as immutable GH
artifacts. Recompute fitted Ridge CV from station-level data without importing
the original model module. Record exact held-out river-geographic groups for
software-version portability. The original model uses GroupKFold with tie
handling that may change between package versions. Do NOT hide disagreement.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np,pandas as pd,scipy,sklearn
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
import sys

SOURCE_SHA="79641d984b28fdbc11dbf02e32f32cda2818ec256f219439a1ccfae601da7dc3"
SOURCE_RUN=42374105184
OUTCOME_RUN=42375542151
FILES={"TP":"TP_2022_reported_wwtp_load_and_original_measured_river_station.csv",
       "TN":"TN_2022_reported_wwtp_load_and_original_measured_river_station.csv"}
UP={"TP":"eea2022_upstream_P_measured_load_nativeunit",
    "TN":"eea2022_upstream_N_measured_load_nativeunit"}
DOWN={"TP":"eea2022_downstream_P_measured_load_nativeunit_placebo",
      "TN":"eea2022_downstream_N_measured_load_nativeunit_placebo"}
NAMES=("M0_geographic_hydrologic_baseline",
       "M1_plus_2022_upstream_measured_load",
       "M2_plus_2022_downstream_source_placebo",
       "M3_upstream_and_downstream")

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def matrix(data,rich=False):
    x=np.column_stack([
        data["WQ.lat"].to_numpy(dtype=float),
        data["WQ.lon"].to_numpy(dtype=float),
        np.log1p(data.DIS_AV_CMS.to_numpy(dtype=float)),
        np.log1p(data.UPLAND_SKM.to_numpy(dtype=float))
    ])
    if rich:
        x=np.column_stack([x,data.crop2.to_numpy(dtype=float),
                            np.log1p(data.pop_dens.to_numpy(dtype=float))])
    return x

def replay(source, saved, par,rich=False):
    if rich:
        source=source.dropna(subset=["crop2","pop_dens"]).copy()
    if (par=="TN" and rich):
        assert len(source)==13
        return {"status":"SOURCE_SUPPORT_BELOW_CV_GATE","n":len(source)}
    source=source.reset_index(drop=True)
    y=source.log_real_2022_median_mgL.to_numpy(float)
    assert np.allclose(y,np.log(source.real_2022_station_median_mgL.to_numpy(float)),
                       atol=1e-8,rtol=0)
    groups=source.spatial_3deg_block.astype(str).to_numpy()
    upstream=np.log1p(source[UP[par]].to_numpy(float))
    downstream=np.log1p(source[DOWN[par]].to_numpy(float))
    base=matrix(source,rich=rich)
    variants=(base,np.column_stack((base,upstream)),
              np.column_stack((base,downstream)),
              np.column_stack((base,upstream,downstream)))
    kfold=list(GroupKFold(n_splits=5).split(base,y,groups))
    assignments=[]
    for i,(train,test) in enumerate(kfold,1):
        assignments.append({"fold":i,"test_n":len(test),
            "geographic_block_labels":sorted(set(groups[test]))})
    lookup=saved["crop_population_2022_CV" if rich else "primary_2022_source_CV"]
    assert lookup["status"]=="FIXED_2022_SPATIAL_BLOCK_CV_COMPLETE"
    estimates={}
    for name,features in zip(NAMES,variants):
        if not np.isfinite(features).all():raise ValueError("Invalid predictor rows")
        predicted=np.full(len(y),np.nan)
        foldm=[]
        for i,(train,test) in enumerate(kfold,1):
            scale=StandardScaler()
            z_train=scale.fit_transform(features[train])
            est=Ridge(alpha=10).fit(z_train,y[train])
            predicted[test]=est.predict(scale.transform(features[test]))
            foldm.append({"fold":i,"n":len(test),
                          "mse":float(np.mean((y[test]-predicted[test])**2))})
        mse=float(np.mean((y-predicted)**2))
        target=float(lookup["metrics"][name]["heldout_mse"])
        estimates[name]={"recomputed_MSE":mse,"original_MSE":target,
                         "absolute_difference":abs(mse-target),
                         "folds":foldm}
    maxerr=max(x["absolute_difference"] for x in estimates.values())
    blocks=sorted(set(groups))
    foldstring=json.dumps(assignments,sort_keys=True,separators=(",",":"))
    return {"status":"MATCH" if maxerr<1e-7 else "MISMATCH",
            "n":len(y),"n_original":lookup["n"],
            "geographic_block_n":len(blocks),
            "heldout_block_assignments_sha256":hashlib.sha256(foldstring.encode()).hexdigest(),
            "heldout_block_assignments":assignments,
            "model_MSE_comparison":estimates,
            "maximum_MSE_disagreement":maxerr}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--original",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    if sha(args.source)!=(SOURCE_SHA):
        raise ValueError("Original source-only station data have changed")
    saved=json.loads((args.original/"EEA2022_ORIGINAL_MEASURED_RECEPTOR_CV_AND_AUDIT.json").read_text())
    if saved["original_source_only_run"]!=SOURCE_RUN or saved["2022_year_aligned_only"] is not True:
        raise ValueError("Unexpected original response protocol")
    outcomes={}
    for par in ("TP","TN"):
        df=pd.read_csv(args.original/FILES[par],dtype={"WQ.ID":str},low_memory=False)
        if df["WQ.ID"].duplicated().any():raise ValueError("Station duplicates")
        if not df["actual_2022_observation_count"].ge(4).all() or not df["distinct_2022_months"].ge(3).all():
            raise ValueError("2022 source quality gate changed")
        assert len(df)==saved["results"][par]["original_2022_station_n"]
        assert int(df[UP[par]].gt(0).sum())==saved["results"][par]["source_2022_upstream_measured_NP_load_positive"]
        subset=df.dropna(subset=["crop2","pop_dens"])
        assert len(subset)==saved["results"][par]["crop_population_complete_2022_stations"]
        outcomes[par]={"primary":replay(df,saved["results"][par],par),
                       "crop_pop":replay(df,saved["results"][par],par,rich=True)}
    report={"status":"MATCH" if all(v["status"]!="MISMATCH" for p in outcomes.values() for v in p.values()) else "MISMATCH",
            "source_run":SOURCE_RUN,"measured_2022_run":OUTCOME_RUN,
            "source_SHA256":SOURCE_SHA,
            "original_files_SHA256":{p:sha(args.original/n) for p,n in FILES.items()},
            "original_results_JSON_SHA256":sha(args.original/"EEA2022_ORIGINAL_MEASURED_RECEPTOR_CV_AND_AUDIT.json"),
            "scientific_scope":"observational source N/P tons per year vs 2022 measured river mg/L; no causal attribution",
            "library_versions":{"python":sys.version.split()[0],"numpy":np.__version__,
                                "pandas":pd.__version__,"scipy":scipy.__version__,"scikit_learn":sklearn.__version__},
            "no_live_network_calls":True,
            "results":outcomes}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_INDEPENDENT_2022_TNTP_CV_REPLAY",json.dumps({
        "audit_status":report["status"],"versions":report["library_versions"],
        "results":{k:{t:{"n":v["n"],"status":v["status"],
                          "max_abs_mse_difference":v.get("maximum_MSE_disagreement"),
                          "fold_manifest_SHA":v.get("heldout_block_assignments_sha256")}
                      for t,v in z.items()} for k,z in outcomes.items()}},
        sort_keys=True),flush=True)
    if report["status"]!="MATCH":
        raise SystemExit(2)

if __name__=="__main__":
    main()
