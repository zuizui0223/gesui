#!/usr/bin/env python3
"""2022 EEA measured WWTP N/P annual source loads versus Peters 2022 river chemistry.

The observational source-only EEA discharge/river crosswalk was frozen BEFORE
access to any Peters receptor concentration values. Sources: official 2022
EEA UWWTD N/P loads in t/year at unique FW discharge outlets; Peters v5
original 2022 *measured* TN/TP station observations in mg/L. No imputed WRTDS/RF.
No plant responses or causal wastewater attribution.
"""
from __future__ import annotations
import argparse,hashlib,io,json,math,zipfile
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
import numpy as np,pandas as pd,requests
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

FIGSHARE="https://api.figshare.com/v2/articles/30546371"
DOI="10.25400/lincolnuninz.30546371.v5"
EEA_SOURCE_FILE_SHA="79641d984b28fdbc11dbf02e32f32cda2818ec256f219439a1ccfae601da7dc3"
ORIGINAL_SOURCE_ONLY_GH_RUN=42374105184
HARMONIZED_IDS={"TP":70896288,"TN":70896279}
COLUMNS={
  "TP":"eea2022_upstream_P_measured_load_nativeunit",
  "TN":"eea2022_upstream_N_measured_load_nativeunit"
}
COLUMNS_DOWN={
  "TP":"eea2022_downstream_P_measured_load_nativeunit_placebo",
  "TN":"eea2022_downstream_N_measured_load_nativeunit_placebo"
}
CALENDAR_YEAR=2022
MIN_SAMPLES=4
MIN_DISTINCT_MONTHS=3
MIN_CV_STATIONS=40
MIN_CV_GROUPS=5
MIN_CV_EXPOSED=10
ALPHA=10.0
S_PER_YEAR=365.25*86400

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def get(url,cap=32_000_000):
    r=requests.get(url,stream=True,timeout=(12,90),
       headers={"User-Agent":"gesui-EEA-real-2022-river-chemistry-study/0.1"})
    r.raise_for_status()
    b=bytearray()
    for x in r.iter_content(131072):
        b.extend(x)
        if len(b)>cap:
            r.close()
            raise ValueError("Original source file exceeded bounded bytes")
    r.close()
    return bytes(b)

def load_2022_observed(original,nutrient,files):
    item=files[f"{nutrient}.2_-_Harmonised.zip"]
    assert int(item["id"])==HARMONIZED_IDS[nutrient]
    raw=get(item["download_url"])
    digest=hashlib.md5(raw).hexdigest()
    if digest!=item["supplied_md5"]:
        raise ValueError("Original Peters observed-source MD5 changed")
    z=zipfile.ZipFile(io.BytesIO(raw))
    names=set(z.namelist())
    result=[]
    nmissing=0
    ninvalid=0
    nperiod=0
    nzero=0
    no_support=0
    dupedates=0
    allsourceobs=0
    for station in original["WQ.ID"].astype(str):
        fn=f"{nutrient}.2_-_Harmonised/{station}.csv"
        if fn not in names:
            nmissing+=1
            continue
        rows=pd.read_csv(io.BytesIO(z.read(fn)))
        if set(("date","value"))-set(rows.columns):
            raise ValueError("Peters original site schema changed")
        allsourceobs+=len(rows)
        date=pd.to_datetime(rows["date"],errors="coerce")
        value=pd.to_numeric(rows["value"],errors="coerce")
        inyear=date.dt.year.eq(CALENDAR_YEAR)
        nperiod+=int(inyear.sum())
        ninvalid+=int(inyear & (date.isna()|value.isna()|~np.isfinite(value)))
        nzero+=int(inyear & value.le(0))
        sub=pd.DataFrame({"date":date,"value":value})
        sub=sub.loc[inyear & sub.value.gt(0) & np.isfinite(sub.value)].copy()
        sub=sub.dropna()
        dupedates+=int(sub.date.duplicated().sum())
        months=int(sub.date.dt.month.nunique())
        if len(sub)<MIN_SAMPLES or months<MIN_DISTINCT_MONTHS:
            no_support+=1;continue
        observed_median=float(sub.value.median())
        if observed_median<=0 or not math.isfinite(observed_median):
            no_support+=1;continue
        result.append({"WQ.ID":station,
            "real_2022_station_median_mgL":observed_median,
            "log_real_2022_median_mgL":math.log(observed_median),
            "actual_2022_observation_count":len(sub),
            "distinct_2022_months":months,
            "distinct_2022_dates":int(sub.date.nunique()),
            "measured_min_2022_mgL":float(sub.value.min()),
            "measured_max_2022_mgL":float(sub.value.max())})
    result=pd.DataFrame(result)
    audit={"original_measurement_file":item["name"],
       "figshare_file_id":int(item["id"]),"original_measured_zip_md5":digest,
       "source_file_bytes":len(raw),"original_sites_requested":len(original),
       "all_original_source_rows_in_station_files":allsourceobs,
       "source_station_files_missing":nmissing,
       "source_2022_rows_including_invalids":nperiod,
       "source_2022_rows_invalid":ninvalid,
       "source_2022_rows_zero_or_negative":nzero,
       "2022_duplicate_sample_date_entries":dupedates,
       "stations_failing_2022_min4_observations_and_3months":no_support,
       "source_stations_measured_2022_qualified":len(result)}
    return result,audit

def covariates(df,rich=False):
    data=pd.DataFrame({
      "latitude":df["WQ.lat"].to_numpy(float),
      "longitude":df["WQ.lon"].to_numpy(float),
      "log_river_mean_Q_m3sec":np.log1p(df.DIS_AV_CMS.to_numpy(float)),
      "log_river_upstream_area_km2":np.log1p(df.UPLAND_SKM.to_numpy(float))})
    if rich:
        data["published_crop2"]=df.crop2.to_numpy(float)
        data["published_log_pop_density"]=np.log1p(df.pop_dens.to_numpy(float))
    return data

def evaluate_cv(df,nutrient,rich=False):
    groups=df["spatial_3deg_block"].astype(str).to_numpy()
    n=len(df)
    ng=len(set(groups))
    n_exposed=int(df[COLUMNS[nutrient]].gt(0).sum())
    gate={"source_n":n,"source_3degree_blocks":ng,
          "source_upstream_measured_load_2022_exposed_n":n_exposed,
          "min_station_requirement":MIN_CV_STATIONS,
          "min_block_requirement":MIN_CV_GROUPS,
          "min_exposed_requirement":MIN_CV_EXPOSED}
    if n<MIN_CV_STATIONS or ng<MIN_CV_GROUPS or n_exposed<MIN_CV_EXPOSED:
        return {"status":"SOURCE_SUPPORT_INSUFFICIENT__NO_GENERALIZABLE_CV",**gate}
    predictors=covariates(df,rich=rich)
    nlog=np.log1p(df[COLUMNS[nutrient]].to_numpy(float))
    dlog=np.log1p(df[COLUMNS_DOWN[nutrient]].to_numpy(float))
    variants={
      "M0_geographic_hydrologic_baseline":predictors,
      "M1_plus_2022_upstream_measured_load":predictors.assign(log_upstream_measured_load_t_per_yr=nlog),
      "M2_plus_2022_downstream_source_placebo":predictors.assign(log_downstream_source_load_t_per_yr=dlog),
      "M3_upstream_and_downstream":predictors.assign(
          log_upstream_measured_load_t_per_yr=nlog,
          log_downstream_source_load_t_per_yr=dlog)}
    y=df.log_real_2022_median_mgL.to_numpy(float)
    cv=GroupKFold(n_splits=MIN_CV_GROUPS)
    folds=list(cv.split(predictors,y,groups))
    out={}
    stored={}
    for model,X in variants.items():
        xx=X.to_numpy(float)
        if not np.isfinite(xx).all():raise ValueError("NaN or inf in fixed candidate predictor design")
        pred=np.full(n,np.nan)
        testgroups=[]
        for index,(train,test) in enumerate(folds,1):
            m=make_pipeline(StandardScaler(),Ridge(alpha=ALPHA))
            m.fit(xx[train],y[train])
            pred[test]=m.predict(xx[test])
            testgroups.append({"fold":index,"n":len(test),
                               "heldout_3degree_blocks":len(set(groups[test])),
                               "mse":float(np.mean((y[test]-pred[test])**2))})
        if not np.isfinite(pred).all():raise ValueError("Model failed to produce every held-out prediction")
        err=y-pred
        mse=float(np.mean(err**2))
        out[model]={"heldout_mse":mse,"heldout_mae":float(np.mean(np.abs(err))),
                    "folds":testgroups}
        stored[model]=err
    base=out["M0_geographic_hydrologic_baseline"]["heldout_mse"]
    changes={}
    for model,stats in out.items():
        if model=="M0_geographic_hydrologic_baseline":continue
        diff=base-stats["heldout_mse"]
        signs={}
        for group in sorted(set(groups)):
            mask=groups==group
            lossdiff=float(np.mean(stored["M0_geographic_hydrologic_baseline"][mask]**2 -
                                  stored[model][mask]**2))
            signs[group]=lossdiff
        changes[model]={"mse_absolute_reduction":diff,
                        "mse_percent_improvement":100*diff/base,
                        "blocks_better":sum(x>0 for x in signs.values()),
                        "blocks_total":len(signs)}
    return {"status":"FIXED_2022_SPATIAL_BLOCK_CV_COMPLETE",
            "n":n,"3degree_blocks":ng,"2022_load_exposed_stations":n_exposed,
            "fixed_Ridge_alpha":ALPHA,
            "metrics":out,"changes_from_geographic_baseline":changes}

def evaluate(df,nutrient):
    if not len(df):return {"status":"NO_2022_ORIGINAL_MEASURED_STATIONS_PASS"}
    col=COLUMNS[nutrient]
    upstream=df[col].to_numpy(float)
    y=df.log_real_2022_median_mgL.to_numpy(float)
    mask=df[col].gt(0)
    src_rho=None
    if len(set(upstream))>1 and len(df)>=20:
        rho,p=spearmanr(np.log1p(upstream),y)
        src_rho={"spearman_rho":float(rho),"p_descriptive":float(p),"n":len(df)}
    # Native EEA 'tons/year' is confirmed from the official Eionet Art15
    # data element definitions for measured N/P source discharge.
    q=df.DIS_AV_CMS.to_numpy(float)
    # 1 metric tonne = 1e9 mg. Annual Q at receiving river = Q * 1e3 L/m3 *365.25*86400s.
    potential=(upstream*1e9)/(q*1e3*S_PER_YEAR)
    if not np.isfinite(potential).all():raise ValueError("Invalid modeled annual mixed concentration")
    complete=df.dropna(subset=["crop2","pop_dens"]).copy()
    return {
      "original_2022_station_n":len(df),
      "source_2022_upstream_measured_NP_load_positive":int(mask.sum()),
      "source_2022_downstream_NP_load_positive":int(df[COLUMNS_DOWN[nutrient]].gt(0).sum()),
      "measured_station_2022_concentration_median_mgL":float(df.real_2022_station_median_mgL.median()),
      "measured_station_2022_concentration_range_mgL":[float(df.real_2022_station_median_mgL.min()),
                                                          float(df.real_2022_station_median_mgL.max())],
      "source_reported_annual_tonnes_upstream_median_among_exposed":float(np.median(upstream[mask])) if mask.any() else None,
      "source_hypothetical_fully_mixed_mgL_median_among_exposed":float(np.median(potential[mask])) if mask.any() else None,
      "source_hypothetical_fully_mixed_mgL_max":float(np.max(potential)) if len(df) else None,
      "raw_source_vs_river_observed_rho":src_rho,
      "primary_2022_source_CV":evaluate_cv(df,nutrient),
      "crop_population_complete_2022_stations":len(complete),
      "crop_population_complete_exposed_stations":int(complete[col].gt(0).sum()),
      "crop_population_2022_CV":evaluate_cv(complete,nutrient,rich=True),
      "notes":["Source reported mass is measured annual WWTP load, NOT observed receptor attribution",
        "Annual Q is HydroRIVERS modeled mean flow, not measured 2022 discharge",
        "Station 2022 medians are original Peters measured sample concentrations only",
        "No original 2011-2024 GBIF plant effect is fit in this gate"]
    }

def main():
    cli=argparse.ArgumentParser()
    cli.add_argument("--source-only",type=Path,required=True)
    cli.add_argument("--out",type=Path,required=True)
    args=cli.parse_args()
    if sha256(args.source_only)!=EEA_SOURCE_FILE_SHA:
        raise ValueError("2022 official EEA source-only measured load input changed from frozen SHA")
    table=pd.read_csv(args.source_only,dtype={"WQ.ID":str},low_memory=False)
    if set(table.nutrient)!={"TP","TN"} or len(table)!=1712:
        raise ValueError("Fixed EEA source-only TP+TN source support table changed")
    if table[["nutrient","WQ.ID"]].duplicated().any():
        raise ValueError("Duplicate station+nutrient would pseudoreplicate")
    article=json.loads(get(FIGSHARE,500000))
    if article.get("doi")!=DOI:raise ValueError("Original Peters v5 DOI version changed")
    files={item["name"]:item for item in article["files"]}
    args.out.mkdir(parents=True,exist_ok=True)
    outcomes={}
    receipts=[]
    for nutrient in ("TP","TN"):
        source=table[table.nutrient==nutrient].copy()
        measured,measure_qc=load_2022_observed(source,nutrient,files)
        receipts.append({"parameter":nutrient,**measure_qc})
        if len(measured):
            joined=source.merge(measured,on="WQ.ID",validate="one_to_one")
            if not len(joined)==len(measured):
                raise ValueError("Measured station outcome failed exact original source join")
            if (joined.DIS_AV_CMS<=0).any():
                raise ValueError("Invalid river discharge used in physical mixing estimate")
            cols=["nutrient","WQ.ID","WQ.lat","WQ.lon","HYRIV_ID","spatial_3deg_block",
                  "DIS_AV_CMS","UPLAND_SKM","crop2","pop_dens",
                  COLUMNS[nutrient],COLUMNS_DOWN[nutrient],
                  "eea2022_upstream_wwtp_both_positive_count",
                  "real_2022_station_median_mgL","log_real_2022_median_mgL",
                  "actual_2022_observation_count","distinct_2022_months","distinct_2022_dates"]
            joined[cols].to_csv(args.out/(nutrient+"_2022_reported_wwtp_load_and_original_measured_river_station.csv"),
                                index=False,float_format="%.10g")
            outcomes[nutrient]=evaluate(joined,nutrient)
        else:
            outcomes[nutrient]={"status":"NO_2022_ORIGINAL_MEASURED_STATIONS_PASS"}
        print("GESUI_EEA2022_ORIGINAL_RECEPTOR_CHEMISTRY",
              json.dumps({"nutrient":nutrient,"original_source":measure_qc,
                  "evaluation":outcomes[nutrient]},sort_keys=True),flush=True)
    report={
      "schema":"gesui_2022_EEA_measured_WWTP_mass_vs_Peters_original_measured_receptor_TNTP_v01",
      "status":"SOURCE_VERIFIED_OBSERVATIONAL_NO_CAUSAL_WASTEWATER_ATTRIBUTION",
      "date_utc":datetime.now(timezone.utc).isoformat(),
      "original_source_only_run":ORIGINAL_SOURCE_ONLY_GH_RUN,
      "original_EEA_source_only_sha256":EEA_SOURCE_FILE_SHA,
      "original_Peters_figshare_v5":DOI,
      "source_measured_N_P_annual_load_units":"tons/year",
      "source_units_confirmed_from_Eionet_official_measured_field_dictionary":True,
      "source_2022_EEA_loads_are_reported_from_measured_concentration_and_annual_volume":True,
      "Peters_measured_2022_original_zip_receipts":receipts,
      "original_plant_GBIF_responses_read":False,
      "2022_year_aligned_only":True,
      "2022_station_eligibility":{"minimum_positive_observations":MIN_SAMPLES,
                                  "minimum_distinct_months":MIN_DISTINCT_MONTHS},
      "fixed_5_group_cv":{"min_stations":MIN_CV_STATIONS,"min_geographic_blocks":MIN_CV_GROUPS,
                          "min_source_exposed_stations":MIN_CV_EXPOSED,
                          "ridge_alpha":ALPHA},
      "effects_causal":False,
      "cautions":["EEA reported facility annual N/P mass is not a uniquely identified driver of ambient station concentration",
                 "No 2022 catchment-specific precipitation/river flow or contemporaneous agricultural emissions",
                 "Only positive BOTH-measured EEA WWTP source sites with exactly one FW discharge point considered",
                 "HydroRIVERS Europe GeoParquet mirror was not independently reproduced against official dataset",
                 "Strong station selection and historical reporting differences; nonuniform spatial coverage",
                 "One 2022 source year cannot establish long-term vegetation effects"],
      "results":outcomes
    }
    (args.out/"EEA2022_ORIGINAL_MEASURED_RECEPTOR_CV_AND_AUDIT.json").write_text(
        json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_EEA2022_NUTRIENT_SOURCE_VS_RECEPTOR_COMPLETE",
         json.dumps({p:{"qualified_sites":outcomes[p].get("original_2022_station_n"),
                         "reported_source_exposed_sites":outcomes[p].get("source_2022_upstream_measured_NP_load_positive"),
                         "primary_cv_status":outcomes[p].get("primary_2022_source_CV",{}).get("status")}
                    for p in ("TP","TN")},sort_keys=True),flush=True)

if __name__=="__main__":main()
