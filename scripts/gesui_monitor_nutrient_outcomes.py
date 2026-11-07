#!/usr/bin/env python3
"""European river TP/TN observed concentration vs directed modeled WWTP discharge.

Predeclared protocol docs/EU_TNTP_WWTP_MACRO_PREDECLARED_20261107.md.
Read *only original Peters v5 harmonised date,value measurements*, no WRTDS
or RF estimates. Single monitoring station = one response, no pseudoreplicated
sample dates, group 3-degree spatial blocks for all predictive comparisons.
Nutrient outcome does NOT establish wastewater-specific N/P chemistry loading.
"""
from __future__ import annotations
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
import argparse,hashlib,io,json,math,zipfile
import numpy as np,pandas as pd,requests

ARTICLE="https://api.figshare.com/v2/articles/30546371"
DOI="10.25400/lincolnuninz.30546371.v5"
OBS_FILE_IDS={"TP":70896288,"TN":70896279}
ORIGINAL_SOURCE_FEATURE_SHA="47fb5bb92726602170a98473d02d2c810f98e49a51ef867d017190a85ae2f10b"
START=pd.Timestamp("2011-01-01")
END=pd.Timestamp("2024-12-31")
MIN_OBS=12
MIN_YEARS=2
MAX_CV_FOLDS=5
MIN_N=40
MIN_BLOCKS=5
ALPHA=10.0
def get(url,limit):
    r=requests.get(url,stream=True,timeout=(10,100),
         headers={"User-Agent":"gesui-EU-monitored-nutrient-statistical-test/0.1"})
    r.raise_for_status();buf=bytearray()
    for block in r.iter_content(131072):
        buf.extend(block)
        if len(buf)>limit:
            r.close()
            raise ValueError("Source downloaded larger than fixed allowed file size")
    r.close();return bytes(buf)
def load_measured(source:pd.DataFrame,nutrient:str,index:dict)->tuple[pd.DataFrame,dict]:
    key=f"{nutrient}.2_-_Harmonised.zip"
    f=index[key]
    if f["id"]!=OBS_FILE_IDS[nutrient]:
        raise ValueError("Original harmonised archive file ID changed")
    raw=get(f["download_url"],32_000_000)
    digest=hashlib.md5(raw).hexdigest()
    if digest!=f["supplied_md5"]:
        raise ValueError("Source observed nutrient archive MD5 changed")
    keep=[];notfound=0;missing=0;excluded=0;dupe_dates=0
    z=zipfile.ZipFile(io.BytesIO(raw))
    entries=set(z.namelist())
    for row in source.itertuples(index=False):
        name=f"{nutrient}.2_-_Harmonised/{row[1]}.csv"
        # robust to pandas namedtuple renaming of WQ.ID -> _1
        if name not in entries:
            notfound+=1;continue
        data=pd.read_csv(io.BytesIO(z.read(name)))
        if not {"date","value"}<=set(data.columns):
            raise ValueError("Original source time series columns changed")
        data["date"]=pd.to_datetime(data["date"],errors="coerce")
        data["value"]=pd.to_numeric(data["value"],errors="coerce")
        invalid=int((data.date.isna()|data.value.isna()).sum())
        missing+=invalid
        data=data.dropna(subset=["date","value"])
        years=data.date.dt.year
        sub=data.loc[data.date.between(START,END)&data.value.gt(0)].copy()
        dupe_dates+=int(sub.date.duplicated().sum())
        if len(sub)<MIN_OBS or sub.date.dt.year.nunique()<MIN_YEARS:
            excluded+=1;continue
        median=float(sub.value.median())
        if not math.isfinite(median) or median<=0:
            excluded+=1;continue
        keep.append({
            "WQ.ID":str(row[1]),"measured_median_mgL":median,
            "log_measured_median_mgL":math.log(median),
            "measured_positive_observations_2011_2024":len(sub),
            "measured_unique_years_2011_2024":int(sub.date.dt.year.nunique()),
            "first_observed_date":str(sub.date.min().date()),
            "last_observed_date":str(sub.date.max().date()),
            "measured_months_2011_2024":int(sub.date.dt.to_period("M").nunique())
        })
    out=pd.DataFrame(keep)
    return out,{
        "original_archive_name":key,"original_figshare_id":f["id"],"archive_original_md5":digest,
        "archive_bytes":len(raw),
        "source_selected_monitored_stations":len(source),
        "source_missing_station_files":notfound,
        "source_missing_date_or_value_rows":missing,
        "positive_measured_2011_2024_low_depth_stations_excluded":excluded,
        "duplicate_date_records_in_eligible_raw_windows":dupe_dates,
        "measured_sites_passing_12_samples_and_two_years":len(out)
    }
def feature_matrix(df:pd.DataFrame,extras=()):
    frame=pd.DataFrame(index=df.index)
    frame["lat"]=pd.to_numeric(df["WQ.lat"],errors="raise")
    frame["lon"]=pd.to_numeric(df["WQ.lon"],errors="raise")
    frame["log_flow"]=np.log1p(pd.to_numeric(df.DIS_AV_CMS,errors="raise").clip(lower=0))
    frame["log_upstream_area"]=np.log1p(pd.to_numeric(df.UPLAND_SKM,errors="raise").clip(lower=0))
    up=np.log1p(pd.to_numeric(df.upstream_50km_ww_m3day,errors="raise").clip(lower=0))
    down=np.log1p(pd.to_numeric(df.downstream_50km_ww_m3day,errors="raise").clip(lower=0))
    for col in extras:
        if col=="crop2":
            frame[col]=pd.to_numeric(df.crop2,errors="raise")
        elif col=="pop_dens":
            frame[col]=np.log1p(pd.to_numeric(df.pop_dens,errors="raise").clip(lower=0))
        elif col=="N_deposition":
            frame[col]=np.log1p(pd.to_numeric(df.N_deposition,errors="raise").clip(lower=0))
        else:raise ValueError("Unfrozen secondary covariate "+col)
    return frame,up,down
def fixed_group_cv(df,extras=()):
    n=len(df);blocks=df.spatial_3deg_block.nunique()
    if n<MIN_N or blocks<MIN_BLOCKS:
        return {"status":"INSUFFICIENT_SOURCE_SUPPORT","n":n,"spatial_blocks":blocks}
    baseline,up,down=feature_matrix(df,extras)
    variants={
       "M0_geography_hydrology":baseline,
       "M1_plus_upstream_ww":baseline.assign(log_upstream_ww=up),
       "M2_plus_downstream_placebo":baseline.assign(log_downstream_ww=down),
       "M3_upstream_and_downstream":baseline.assign(log_upstream_ww=up,log_downstream_ww=down)
    }
    target=df.log_measured_median_mgL.to_numpy(dtype=float)
    groups=df.spatial_3deg_block.astype(str).to_numpy()
    splits=list(GroupKFold(n_splits=MAX_CV_FOLDS).split(baseline,target,groups))
    metrics={}
    perfold=[]
    residuals={}
    for model_name,x in variants.items():
        x=x.to_numpy(dtype=float)
        if not np.isfinite(x).all():raise ValueError("NaN or infinite covariate in fixed model "+model_name)
        predictions=np.full(n,np.nan)
        for fold,(train,test) in enumerate(splits,1):
            m=Pipeline([("standardizer",StandardScaler()),("ridge",Ridge(alpha=ALPHA))])
            m.fit(x[train],target[train])
            pred=m.predict(x[test])
            predictions[test]=pred
            perfold.append({"model":model_name,"fold":fold,"test_sites":len(test),
                  "test_blocks":len(set(groups[test])),
                  "test_mse":float(np.mean((target[test]-pred)**2)),
                  "test_mae":float(np.mean(np.abs(target[test]-pred)))})
        if not np.isfinite(predictions).all():raise ValueError("Predictions incomplete")
        errors=target-predictions
        metrics[model_name]={
            "n":n,"geographic_3degree_blocks":blocks,
            "heldout_mse":float(np.mean(errors**2)),
            "heldout_mae":float(np.mean(np.abs(errors)))}
        residuals[model_name]=errors
    base_mse=metrics["M0_geography_hydrology"]["heldout_mse"]
    metrics["delta_vs_M0"]={model:{
        "mse_pct_improvement":float(100*(base_mse-v["heldout_mse"])/base_mse),
        "mse_absolute_reduction":float(base_mse-v["heldout_mse"])}
        for model,v in metrics.items() if model in variants and model!="M0_geography_hydrology"}
    # A spatial block is the unit for paired model diagnostics, not a sampling occurrence.
    block_sums={}
    for name in ("M1_plus_upstream_ww","M2_plus_downstream_placebo","M3_upstream_and_downstream"):
        improvements=[]
        for group in sorted(set(groups)):
            mask=groups==group
            improvements.append({"block":group,"n":int(mask.sum()),
                "mean_squared_error_difference_baseline_minus_added":
                  float(np.mean(residuals["M0_geography_hydrology"][mask]**2-
                                residuals[name][mask]**2))})
        block_sums[name]={
          "blocks_with_positive_delta":sum(x["mean_squared_error_difference_baseline_minus_added"]>0 for x in improvements),
          "blocks_total":len(improvements)}
    return {"status":"FIXED_SPATIAL_GROUP_CV_COMPLETE",
            "n":n,"spatial_3degree_blocks":blocks,"alpha_fixed":ALPHA,
            "metrics":metrics,"fold_metrics":perfold,"spatial_block_signs":block_sums}
def simple_associations(df):
    exposure=np.log1p(df.upstream_50km_ww_m3day.to_numpy(float))
    downstream=np.log1p(df.downstream_50km_ww_m3day.to_numpy(float))
    y=df.log_measured_median_mgL.to_numpy(float)
    def r(x):
        if len(x)<10 or len(set(x))<2:return {"status":"INSUFFICIENT_VARIATION","n":len(x)}
        rho,p=spearmanr(x,y)
        return {"n":len(x),"rho":float(rho),"p_exploratory":float(p)}
    return {"source_upstream_WWTP_volume":r(exposure),
            "source_downstream_WWTP_placebo":r(downstream)}
def analyze(df,nutrient):
    base=len(df)
    if base==0:return {"status":"NO_MEASUREMENTS_PASS_PRESPECIFIED_RULE"}
    exposed=int(df.upstream_50km_ww_m3day.gt(0).sum())
    secondary=df.loc[df[["crop2","pop_dens","N_deposition"]].notna().all(axis=1)].copy()
    secondary_exposed=int(secondary.upstream_50km_ww_m3day.gt(0).sum())
    out={
      "nutrient":nutrient,"status":"OBSERVATIONAL_SOURCE_ALIGNED_PREDICTION_NOT_CAUSAL",
      "source_total_sites_with_measured_2011_2024":base,
      "original_upstream_WWTP_positive_sites":exposed,
      "secondary_landcover_population_deposition_sites":len(secondary),
      "secondary_WWTP_positive_sites":secondary_exposed,
      "mean_observed_concentration_mgL":float(df.measured_median_mgL.mean()),
      "median_observed_station_concentration_mgL":float(df.measured_median_mgL.median()),
      "source_associations":simple_associations(df),
      "primary_source_group_cv":fixed_group_cv(df),
      "landuse_complete_source_group_cv":fixed_group_cv(secondary,("crop2","pop_dens","N_deposition"))}
    if exposed<20:
        out["exposure_support_warning"]="Under 20 stations have nonzero observed upstream WWTP proxy; nutrient associations unstable"
    if secondary_exposed<20:
        out["landuse_exposure_support_warning"]="Under 20 land-complete stations have nonzero upstream WWTP proxy; rich-model result particularly fragile"
    return out
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-features",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    b=args.source_features.read_bytes()
    if hashlib.sha256(b).hexdigest()!=ORIGINAL_SOURCE_FEATURE_SHA:
        raise ValueError("Original source-only station exposure table digest mismatch")
    features=pd.read_csv(io.BytesIO(b),dtype={"WQ.ID":str},low_memory=False)
    assert set(features.nutrient)=={"TP","TN"} and not features[["nutrient","WQ.ID"]].duplicated().any()
    assert not set(("C_mean","C_median","L","Y")) & set(features.columns)
    article=json.loads(get(ARTICLE,500000))
    if article["doi"]!=DOI:raise ValueError("Peters v5 DOI changed")
    idx={f["name"]:f for f in article["files"]}
    args.out.mkdir(parents=True,exist_ok=True)
    results={};sources=[]
    # Chemical measurements are only read after this source design and frozen protocol.
    for nutrient in ("TP","TN"):
        source=features[features.nutrient==nutrient].copy()
        measured,provenance=load_measured(source,nutrient,idx)
        sources.append({"nutrient":nutrient,**provenance})
        if len(measured):
            data=source.merge(measured,how="inner",on="WQ.ID",validate="one_to_one")
            if data.WQ_ID.duplicated().any() if "WQ_ID" in data else False:
                raise ValueError("Unexpected numeric taxon identifier conflict")
            data.to_csv(args.out/f"{nutrient}_monitored_river_wwtp_derived_station_panel.csv",
                        index=False,float_format="%.8g")
            result=analyze(data,nutrient)
        else:result={"status":"NO_MEASURED_SOURCE_STATIONS"}
        results[nutrient]=result
        print("GESUI_EU_RIVER_TNTP_FIXED_MODEL",json.dumps(result,sort_keys=True),flush=True)
    out={
       "schema":"gesui_EU_monitor_river_TNTP_measured_vs_sourceWWTP_v01",
       "source_only_input_sha256":ORIGINAL_SOURCE_FEATURE_SHA,
       "Peters_doi":DOI,
       "source_measurement_receipts":sources,
       "source_run_id":42371629424,
       "date_utc":datetime.now(timezone.utc).isoformat(),
       "original_plant_response_examined":False,
       "measurements_original_not_WRTDS_or_RF":True,
       "matched_observation_years":"2011-2024",
       "min_positive_measurements_per_station":MIN_OBS,
       "min_unique_observation_years":MIN_YEARS,
       "landuse_census_selects_Peters_original_model_input_subset":True,
       "Peters_published_C_mean_C_median_L_Y_leakage_columns_excluded":True,
       "river_network_is_derived_third_party_not_independent_official":True,
       "causal_claims":False,
       "constraints":[
          "Estimated WWTP water volume is not nutrient mass load or local measured discharged TP",
          "Urbanization, upstream agriculture, sediment sources, seasonality and monitoring selection can confound",
          "Peters published environmental covariate subset selected on original model inputs",
          "Spatial block CV is observational source prediction not causal mechanism",
          "Zero registered WWTP upstream50km does not mean no sewer discharge",
          "Station medians are individual response units; sampling detection frequencies vary"],
       "results":results}
    (args.out/"MODEL_RESULTS_AND_AUDIT.json").write_text(json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_EU_NUTRIENT_SOURCE_CV_COMPLETE",json.dumps({
      "source_only_input_sha256":ORIGINAL_SOURCE_FEATURE_SHA,
      "TP":{k:results["TP"].get(k) for k in
         ("source_total_sites_with_measured_2011_2024","original_upstream_WWTP_positive_sites",
          "secondary_landcover_population_deposition_sites","secondary_WWTP_positive_sites")},
      "TN":{k:results["TN"].get(k) for k in
         ("source_total_sites_with_measured_2011_2024","original_upstream_WWTP_positive_sites",
          "secondary_landcover_population_deposition_sites","secondary_WWTP_positive_sites")}},sort_keys=True),flush=True)
if __name__=="__main__":
    main()
