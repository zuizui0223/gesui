#!/usr/bin/env python3
"""Independent *offline* audit of real EU nutrient/source CV from source artifacts.

Recomputes metrics from frozen source-only and measured-station outputs.
No external nutrient, vegetation, weather, or Figshare requests.
"""
from __future__ import annotations
from pathlib import Path
import argparse,hashlib,json,math
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
ROOT_EXPECTED={
    "nutrient_monitor_station_source_features.csv":"47fb5bb92726602170a98473d02d2c810f98e49a51ef867d017190a85ae2f10b",
    "MODEL_RESULTS_AND_AUDIT.json":"912ae7a912522ff9836203e4a47c8a9696564cfcaeb3cd7a4f8990ce23c5ac34",
    "TP_monitored_river_wwtp_derived_station_panel.csv":"3e284ba3e4a3ee6b0dfc5cace1384d9d4a23b4331b31694d37fa4288fc4edbe5",
    "TN_monitored_river_wwtp_derived_station_panel.csv":"0d85bd1ba51e263684a690cc0cb147f1569eb0727600e9716701f8c267d6df04"
}
def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def validate_equal(a,b,tol=2e-5):
    if not math.isclose(float(a),float(b),rel_tol=0,abs_tol=tol):
        raise ValueError(f"Independent numerical mismatch: {a} vs {b}")
def cv(df,features):
    y=np.log(df.measured_median_mgL.to_numpy(dtype=float))
    groups=df.spatial_3deg_block.astype(str).to_numpy()
    x=np.column_stack([np.asarray(f,dtype=float) for f in features])
    if not np.isfinite(x).all():raise ValueError("NaN/inf in audit predictors")
    preds=np.full(len(df),np.nan)
    splits=GroupKFold(n_splits=5)
    for train,test in splits.split(x,y,groups):
        standard=StandardScaler()
        xx=standard.fit_transform(x[train])
        x_test=standard.transform(x[test])
        fitted=Ridge(alpha=10.0).fit(xx,y[train])
        preds[test]=fitted.predict(x_test)
    return float(np.mean((y-preds)**2))
def audit_one(df,record):
    assert not set(("C_mean","C_median","L","Y"))&set(df.columns)
    assert not df["WQ.ID"].duplicated().any()
    df=df.reset_index(drop=True)
    validate_equal(float(df.measured_median_mgL.median()),record["median_observed_station_concentration_mgL"])
    assert len(df)==record["source_total_sites_with_measured_2011_2024"]
    assert int(df.upstream_50km_ww_m3day.gt(0).sum())==record["original_upstream_WWTP_positive_sites"]
    lat=df["WQ.lat"].to_numpy(float)
    lon=df["WQ.lon"].to_numpy(float)
    q=np.log1p(df.DIS_AV_CMS.to_numpy(float))
    area=np.log1p(df.UPLAND_SKM.to_numpy(float))
    upstream=np.log1p(df.upstream_50km_ww_m3day.to_numpy(float))
    downstream=np.log1p(df.downstream_50km_ww_m3day.to_numpy(float))
    base=[lat,lon,q,area]
    tests={
       "M0_geography_hydrology":base,
       "M1_plus_upstream_ww":base+[upstream],
       "M2_plus_downstream_placebo":base+[downstream],
       "M3_upstream_and_downstream":base+[upstream,downstream]}
    primary={}
    for k,v in tests.items():
        observed=cv(df,v)
        published=record["primary_source_group_cv"]["metrics"][k]["heldout_mse"]
        validate_equal(observed,published)
        primary[k]=observed
    rho,_=spearmanr(upstream,df.log_measured_median_mgL)
    validate_equal(rho,record["source_associations"]["source_upstream_WWTP_volume"]["rho"])
    source_extra=record["secondary_covariates_source_available"]
    extra=df.dropna(subset=source_extra).copy().reset_index(drop=True)
    assert len(extra)==record["secondary_landcover_population_deposition_sites"]
    assert extra.upstream_50km_ww_m3day.gt(0).sum()==record["secondary_WWTP_positive_sites"]
    secondary={}
    if len(extra)>=40 and extra.spatial_3deg_block.nunique()>=5:
        lat=extra["WQ.lat"].to_numpy(float);lon=extra["WQ.lon"].to_numpy(float)
        q=np.log1p(extra.DIS_AV_CMS.to_numpy(float))
        area=np.log1p(extra.UPLAND_SKM.to_numpy(float))
        bs=[lat,lon,q,area,extra.crop2.to_numpy(float),
            np.log1p(extra.pop_dens.to_numpy(float))]
        if "N_deposition" in source_extra:bs.append(np.log1p(extra.N_deposition.to_numpy(float)))
        up=np.log1p(extra.upstream_50km_ww_m3day.to_numpy(float))
        down=np.log1p(extra.downstream_50km_ww_m3day.to_numpy(float))
        sx={
            "M0_geography_hydrology":bs,
            "M1_plus_upstream_ww":bs+[up],
            "M2_plus_downstream_placebo":bs+[down],
            "M3_upstream_and_downstream":bs+[up,down]}
        for name,cols in sx.items():
            calc=cv(extra,cols)
            expected=record["landuse_complete_source_group_cv"]["metrics"][name]["heldout_mse"]
            validate_equal(calc,expected)
            secondary[name]=calc
    else:
        assert record["landuse_complete_source_group_cv"]["status"]=="INSUFFICIENT_SOURCE_SUPPORT"
    return {"n":len(df),"exposed":int(df.upstream_50km_ww_m3day.gt(0).sum()),
            "independent_source_covariate_n":len(extra),
            "independent_source_covariate_exposed_n":int(extra.upstream_50km_ww_m3day.gt(0).sum()),
            "all_panel_holdout_mse_verified":primary,
            "landuse_panel_holdout_mse_verified":secondary,
            "spearman_upstream_rho_verified":float(rho),
            "unique_3degree_blocks":int(df.spatial_3deg_block.nunique())}
def run(args):
    checked={}
    for name,digest in ROOT_EXPECTED.items():
        where=(args.source_dir if "source_features.csv" in name else args.outcome_dir)/name
        checked[name]=sha(where)
        if checked[name]!=digest:raise ValueError(f"Source artifact checksum changed: {name}")
    source=pd.read_csv(args.source_dir/"nutrient_monitor_station_source_features.csv",
                       dtype={"WQ.ID":str},low_memory=False)
    saved=json.loads((args.outcome_dir/"MODEL_RESULTS_AND_AUDIT.json").read_text())
    assert saved["source_only_input_sha256"]==ROOT_EXPECTED["nutrient_monitor_station_source_features.csv"]
    result={}
    for par in ("TP","TN"):
        frame=pd.read_csv(args.outcome_dir/f"{par}_monitored_river_wwtp_derived_station_panel.csv",
                          dtype={"WQ.ID":str},low_memory=False)
        orig=source[source.nutrient==par]
        assert set(frame["WQ.ID"]).issubset(set(orig["WQ.ID"]))
        matched=frame[["WQ.ID","upstream_50km_ww_m3day","downstream_50km_ww_m3day"]].merge(
            orig[["WQ.ID","upstream_50km_ww_m3day","downstream_50km_ww_m3day"]],
            on="WQ.ID",validate="one_to_one",suffixes=("_new","_source"))
        assert len(matched)==len(frame)
        for col in ("upstream_50km_ww_m3day","downstream_50km_ww_m3day"):
            if not np.allclose(matched[col+"_new"],matched[col+"_source"],atol=.01,rtol=0):
                raise ValueError("Modelled upstream exposure values changed")
        result[par]=audit_one(frame,saved["results"][par])
    summary={
      "schema":"gesui_independent_offline_real_river_measured_tp_tn_audit_v01",
      "status":"PASS_SOURCE_SHA_AND_SPATIAL_CV_RECONSTRUCTION",
      "raw_nutrient_source_redownloaded":False,
      "source_ARTIFACT_run_id":42371629424,
      "outcome_ARTIFACT_run_id":42372235775,
      "verified_input_sha256":checked,
      "claim_scope":"observational modelled wastewater water flows and monitored river chemistry, NOT sewage-specific nutrient load",
      "results":result}
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print("GESUI_INDEPENDENT_EU_RIVER_MACRO_AUDIT",json.dumps(summary,sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",required=True,type=Path)
    p.add_argument("--outcome-dir",required=True,type=Path)
    p.add_argument("--report",required=True,type=Path)
    run(p.parse_args())
