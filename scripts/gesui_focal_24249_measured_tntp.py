#!/usr/bin/env python3
"""Predeclared 24249 TN/TP source-availability check and TP-only exploratory gradient.

One source-first selected outfall; original measured/harmonized site files from
Peters2026 figshare v5. No plant response or wastewater causal inference.
Missing TN source sites cannot be imputed from modeled fields.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
from pathlib import Path
import hashlib
import io
import json
import math
from collections import Counter
import zipfile
import numpy as np
import pandas as pd
import requests

ARTICLE="https://api.figshare.com/v2/articles/30546371"
DOI="10.25400/lincolnuninz.30546371.v5"
SOURCES=("TN.2_-_Harmonised.zip","TP.2_-_Harmonised.zip",
         "TN.2_meta.csv","TP.2_meta.csv","TN.2_unavailable_ID.csv")
STATIONS={
    "TN":{"upstream":"WQ016702","downstream":"WQ016695"},
    "TP":{"upstream":"WQ976406","downstream":"WQ977059",
          "alternative_upstream":"WQ977062"}
}
PRIMARY_YEAR_START=2012
PRIMARY_YEAR_END=2013
TP_SENS_YEAR_START=2008
TP_SENS_YEAR_END=2022
OUTFALL_ID=24249
FRACTION=1.7102364894632606e-06

def get(url,cap):
    r=requests.get(url,stream=True,timeout=(15,100),
                   headers={"User-Agent":"gesui-NTP-locked-focal-source-audit/0.1"})
    r.raise_for_status()
    b=bytearray()
    for chunk in r.iter_content(131072):
        b.extend(chunk)
        if len(b)>cap:
            r.close()
            raise ValueError("Bounded source ZIP exceeded cap")
    r.close()
    return bytes(b)

def safe_csv(blob):
    return pd.read_csv(io.BytesIO(blob),low_memory=False,dtype={"WQ.ID":str})

def paired_months(up,down,start,end):
    def reshape(g):
        a=g.loc[g["date"].dt.year.between(start,end)&g["value"].gt(0)].copy()
        a["month"]=a["date"].dt.to_period("M").astype(str)
        return a.groupby("month")["value"].median()
    u=reshape(up);d=reshape(down)
    panel=pd.concat([u.rename("upstream_mg_L"),d.rename("downstream_mg_L")],axis=1).dropna()
    panel=panel.reset_index().sort_values("month")
    if not panel.empty:
        panel["downstream_to_upstream_ratio"]=panel.downstream_mg_L/panel.upstream_mg_L
        panel["log_downstream_minus_upstream"]=np.log(panel.downstream_to_upstream_ratio)
    def nwithin(df):
        return int(df.date.dt.year.between(start,end).sum())
    report={
        "period":f"{start}-{end}",
        "upstream_original_measurement_count_in_window":nwithin(up),
        "downstream_original_measurement_count_in_window":nwithin(down),
        "valid_same_calendar_month_pairs":len(panel),
        "sufficient_six_month_pairs":len(panel)>=6,
        "downstream_over_upstream_months":int((panel.downstream_mg_L>panel.upstream_mg_L).sum())
            if len(panel) else 0,
        "equal_months":int((panel.downstream_mg_L==panel.upstream_mg_L).sum())
            if len(panel) else 0,
        "median_monthly_downstream_upstream_ratio":float(panel.downstream_to_upstream_ratio.median())
            if len(panel)>=6 else None,
        "mean_log_ratio":float(panel.log_downstream_minus_upstream.mean())
            if len(panel)>=6 else None,
        "median_log_ratio":float(panel.log_downstream_minus_upstream.median())
            if len(panel)>=6 else None,
        "monthly_values_are_independent_observations":False,
    }
    return panel,report

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    article=json.loads(get(ARTICLE,500000))
    if article["doi"]!=DOI:raise ValueError("Figshare article DOI changed")
    info={f["name"]:f for f in article["files"]}
    raw={}
    receipt=[]
    for name in SOURCES:
        f=info[name]
        data=get(f["download_url"],32_000_000)
        digest=hashlib.md5(data).hexdigest()
        if digest!=f["supplied_md5"]:
            raise ValueError("Peters source md5 differs: "+name)
        raw[name]=data
        receipt.append({"name":name,"original_id":f["id"],"md5":digest,
                        "actual_bytes":len(data),"license":"CC BY 4.0"})
    meta={}
    for nutrient in ("TN","TP"):
        meta[nutrient]=safe_csv(raw[f"{nutrient}.2_meta.csv"]).set_index("WQ.ID")
    unavailable=raw["TN.2_unavailable_ID.csv"].decode("utf8","replace")
    measurements={}
    source_qc={}
    for nutrient in ("TN","TP"):
        z=zipfile.ZipFile(io.BytesIO(raw[f"{nutrient}.2_-_Harmonised.zip"]))
        members=set(z.namelist())
        m={}
        for role,station in STATIONS[nutrient].items():
            path=f"{nutrient}.2_-_Harmonised/{station}.csv"
            present=path in members
            meta_row=meta[nutrient].loc[station] if station in meta[nutrient].index else None
            q={"station_id":station,"source_zip_site_file_exists":present,
               "listed_in_TN_unavailable_ID_file":station in unavailable if nutrient=="TN" else None,
               "metadata_available_flag":str(meta_row.get("WQ.available")) if meta_row is not None else None,
               "metadata_date_first":str(meta_row.get("start")) if meta_row is not None else None,
               "metadata_date_last":str(meta_row.get("end")) if meta_row is not None else None}
            if present:
                b=z.read(path)
                df=pd.read_csv(io.BytesIO(b))
                if set(df.columns)!={"date","value"}:
                    raise ValueError("Unexpected measured input schema for "+path)
                df["date"]=pd.to_datetime(df["date"],errors="coerce")
                df["value"]=pd.to_numeric(df["value"],errors="coerce")
                q["source_rows"]=len(df)
                q["missing_date_or_value"]=int(df[["date","value"]].isna().any(axis=1).sum())
                q["nonpositive_values"]=int((df.value<=0).sum())
                df=df.dropna(subset=["date","value"])
                df=df[df.value>0].copy()
                m[role]=df
                q["usable_positive_obs"]=len(df)
                q["earliest_obs"]=str(df.date.min().date()) if len(df) else None
                q["latest_obs"]=str(df.date.max().date()) if len(df) else None
            source_qc[f"{nutrient}_{role}"]=q
            print("GESUI_FOCAL_NUTRIENT_SOURCE_AVAILABILITY",json.dumps(
                {"nutrient":nutrient,"role":role,**q},sort_keys=True),flush=True)
        measurements[nutrient]=m
    primary_tn_files=all(x in measurements["TN"] for x in ("upstream","downstream"))
    primary_tp_files=all(x in measurements["TP"] for x in ("upstream","downstream"))
    primary_both_available=primary_tn_files and primary_tp_files
    args.out.mkdir(parents=True,exist_ok=True)
    comparisons={}
    if primary_tp_files:
        for a,b,label in [(2012,2013,"common_2012_2013"),
                          (2008,2022,"TP_premarked_2008_2022_sensitivity")]:
            panel,analysis=paired_months(measurements["TP"]["upstream"],
                                         measurements["TP"]["downstream"],a,b)
            panel.to_csv(args.out/("TP_"+label+"_monthly_pairs.csv"),index=False,float_format="%.8g")
            comparisons[label]=analysis
            print("GESUI_FOCAL_SOURCE_MEASURED_TP_MONTHS",json.dumps(analysis,sort_keys=True),flush=True)
    # Keep the TN-specific preregistered outcome as untestable if either site is absent.
    if not primary_tn_files:
        comparisons["TN_common_2012_2013"]={
           "status":"NOT_TESTABLE_SOURCE_ID_MISSING_FROM_PUBLIC_HARMONISED_ZIP",
           "no_imputation_or_modelled_values_substituted":True}
    result={
        "schema":"gesui_focal_24249_Peters2026_measured_NP_source_gate_v01",
        "retrieved_utc":datetime.now(timezone.utc).isoformat(),
        "source_doi":DOI,"source_file_checksums":receipt,
        "outfall_id":OUTFALL_ID,"estimated_ww_flow_fraction":FRACTION,
        "data_source":"harmonised_observed_concentration_mg_L_no_WRTDS",
        "predeclared_up_down_stations":STATIONS,
        "per_station_source_qc":source_qc,
        "TN_both_primary_site_files_available":primary_tn_files,
        "TP_both_primary_site_files_available":primary_tp_files,
        "combined_TN_TP_primary_estimand_testable":primary_both_available,
        "analyses":comparisons,
        "primary_result_status":"NOT_TESTABLE" if not primary_both_available else "SOURCE_AVAILABLE_ONLY",
        "cautions":[
          "Observations are in-river, NOT WWTP effluent concentrations or loads",
          "One outfall, with extreme dilution and asymmetric network station distances",
          "Upstream and downstream locations may receive other pollution from multiple sources",
          "Source metadata overlap is not proof of concurrent concentration measures",
          "No plant trait-response fitted or general wastewater mechanism identified",
          "No unavailable TN source station records reconstructed from estimates"
        ]}
    (args.out/"SOURCE_AND_MEASURED_TP_REPORT.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_FOCAL_NUTRIENT_ANALYSIS_STATUS",json.dumps({
        "combined_primary_TNTP_testable":primary_both_available,
        "TP_2012_2013":comparisons.get("common_2012_2013"),
        "TP_2008_2022":comparisons.get("TP_premarked_2008_2022_sensitivity"),
        "TN_status":comparisons.get("TN_common_2012_2013")},
        sort_keys=True),flush=True)
if __name__=="__main__":main()
