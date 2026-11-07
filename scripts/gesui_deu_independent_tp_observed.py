#!/usr/bin/env python3
"""Measured TP for six PRESELECTED, station-independent German WWTP outfalls.

No source response is used in selection; source-only frozen JSON is input.
Not a test of sewage-P causality and not a TN substitution.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile
import pandas as pd
import numpy as np
import requests

ARTICLE="https://api.figshare.com/v2/articles/30546371"
ARTICLE_DOI="10.25400/lincolnuninz.30546371.v5"
FILENAME="TP.2_-_Harmonised.zip"
FILE_ID=70896288
MD5="b596dbe9833c49426147a99f10f3fe31"
YEARS=(2011,2022)

def fetch(url,limit):
    r=requests.get(url,stream=True,timeout=(15,100),
      headers={"User-Agent":"gesui-locked-deu-independent-tp-measurements/0.1"})
    r.raise_for_status()
    b=bytearray()
    for block in r.iter_content(131072):
        b.extend(block)
        if len(b)>limit:
            r.close()
            raise ValueError("Download exceeds source-fixed byte cap")
    r.close()
    return bytes(b)

def original_source_cohort(path):
    root=json.loads(path.read_text())
    assert root["selection_made_without_plant_responses_or_nutrient_concentrations"] is True
    assert root["strict_independence_audit"]["TP"]["all_station_ids_unique"] is True
    a=root["independent_source_only_cohorts"]["TP"]
    assert len(a)==6 and len(root["independent_source_only_cohorts"]["TN"])==0
    assert len({x["wwtp_id"] for x in a})==6
    assert len({x["hydrobasin_l12"] for x in a})==6
    assert {x["country"] for x in a}=={"DEU"}
    assert {x["independent_source_fraction_stratum"] for x in a}=={"high","low"}
    assert sum(x["independent_source_fraction_stratum"]=="high" for x in a)==3
    ids=[v["station_id"] for x in a
        for v in x["measured_source_stations"]["TP"].values()]
    assert len(ids)==len(set(ids))==12
    return a

def month_medians(t:pd.DataFrame):
    z=t.loc[t.date.dt.year.between(*YEARS) & t.value.gt(0)].copy()
    z["month"]=z.date.dt.to_period("M").astype(str)
    return z.groupby("month").value.median()

def run(args):
    cohort=original_source_cohort(args.source)
    source=json.loads(fetch(ARTICLE,500_000))
    if source.get("doi")!=ARTICLE_DOI:
        raise ValueError("Figshare data source changed")
    file=next((x for x in source["files"] if x["name"]==FILENAME),None)
    if file is None or file["id"]!=FILE_ID or file.get("supplied_md5")!=MD5:
        raise ValueError("Original TP archive identity changed")
    zip_bytes=fetch(file["download_url"],32_000_000)
    if hashlib.md5(zip_bytes).hexdigest()!=MD5:
        raise ValueError("Measured TP archive original MD5 mismatch")
    z=zipfile.ZipFile(io.BytesIO(zip_bytes))
    station_cache={}
    original=[]
    result=[]
    paired=[]
    for e in cohort:
        wid=int(e["wwtp_id"])
        meta=e["measured_source_stations"]["TP"]
        timeseries={}
        stats={}
        for side in ("UPSTREAM","DOWNSTREAM"):
            station=meta[side]["station_id"]
            name=f"TP.2_-_Harmonised/{station}.csv"
            if name not in z.namelist():
                raise ValueError("Source-preselected station missing from original TP archive: "+station)
            t=pd.read_csv(io.BytesIO(z.read(name)))
            if set(t.columns)!={"date","value"}:
                raise ValueError("Unexpected original measured TP schema")
            total=len(t)
            t["date"]=pd.to_datetime(t["date"],errors="coerce")
            t["value"]=pd.to_numeric(t["value"],errors="coerce")
            missing=int(t[["date","value"]].isna().any(axis=1).sum())
            nonpositive=int((t.value<=0).sum())
            t=t.dropna(subset=["date","value"])
            t=t[t.value>0]
            t=t[t.date.dt.year.between(*YEARS)].copy()
            stats[side]={"station_id":station,"all_original_rows":total,
                 "source_missing_rows":missing,"nonpositive_original_rows":nonpositive,
                 "valid_positive_rows_in_2011_2022":len(t),
                 "source_network_distance_km":float(meta[side]["network_distance_km"])}
            for row in t.itertuples(index=False):
                original.append({"outfall_id":wid,"role":side,"station_id":station,
                                 "date":row.date.date().isoformat(),
                                 "observed_TP_mg_L":float(row.value)})
            timeseries[side]=month_medians(t)
        join=pd.concat([timeseries["UPSTREAM"].rename("upstream_TP_mg_L"),
                        timeseries["DOWNSTREAM"].rename("downstream_TP_mg_L")],
                        axis=1).dropna()
        join=join.sort_index()
        if len(join):
            join["log_down_over_up"]=np.log(join.downstream_TP_mg_L/join.upstream_TP_mg_L)
        for month,row in join.iterrows():
            paired.append({"outfall_id":wid,
                "source_fraction_stratum":e["independent_source_fraction_stratum"],
                "yearmonth":month,
                "upstream_TP_mg_L":float(row.upstream_TP_mg_L),
                "downstream_TP_mg_L":float(row.downstream_TP_mg_L),
                "log_down_over_up":float(row.log_down_over_up)})
        eligible=len(join)>=12
        result.append({
           "outfall_id":wid,"country":e["country"],
           "hydrobasin_l12":e["hydrobasin_l12"],
           "source_modeled_effluent_water_fraction":float(e["effluent_fraction_proxy"]),
           "source_fraction_stratum":e["independent_source_fraction_stratum"],
           "source_unique_monitoring_pair":True,
           "station_coverage":stats,
           "calendar_months_with_measured_up_down":len(join),
           "primary_source_months_gate_met":eligible,
           "median_monthly_down_up_ratio":float(
               np.exp(join.log_down_over_up.median())) if eligible else None,
           "median_monthly_log_down_up":float(join.log_down_over_up.median()) if eligible else None,
           "fraction_months_downstream_higher":float((join.log_down_over_up>0).mean())
               if eligible else None,
           "actual_TP_is_WWTP_effluent_load":False,
        })
        print("GESUI_INDEPENDENT_TP_OUTFALL",json.dumps({
            "wwtp_id":wid,"stratum":e["independent_source_fraction_stratum"],
            "matched_calendar_months":len(join),"source_gate_met":eligible,
            "median_ratio":result[-1]["median_monthly_down_up_ratio"]},sort_keys=True),flush=True)
    args.out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(original).sort_values(["outfall_id","role","date"]).to_csv(
        args.out/"observed_TP_station_dates_2011_2022.csv",index=False,float_format="%.8g")
    pd.DataFrame(paired).sort_values(["outfall_id","yearmonth"]).to_csv(
        args.out/"paired_monthly_TP_mg_L.csv",index=False,float_format="%.8g")
    report={
        "status":"SIX_SOURCE_BLIND_INDEPENDENT_TP_STATIONS_REAL_OBSERVED",
        "source_doi":ARTICLE_DOI,
        "source_zip_file_id":FILE_ID,
        "source_original_zip_MD5":MD5,
        "original_preoutcome_source_selection_run":42371319376,
        "source_selected_outfall_count":6,
        "metadata_selected_independent_station_ids":12,
        "observed_TP_positive_values_2011_2022":len(original),
        "paired_calendar_month_rows":len(paired),
        "source_units_with_12_or_more_paired_months":sum(x["primary_source_months_gate_met"] for x in result),
        "source_year_window":YEARS,
        "all_results":result,
        "limitations":[
          "Measured source TP concentrations represent receiving river, not the WWTP effluent",
          "Within river station distance 7-25km, unaccounted farming, tributaries, previous releases",
          "Only six source-selected Germany outfalls, no across-country statistical generalization",
          "Monthly medians are serially correlated and not independent experimental replicates",
          "Measured TP source station availability not proof of simultaneous TN",
          "No plant results were read when picking source stations or stratification"
        ],
    }
    (args.out/"TP_SOURCE_RESULTS.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("GESUI_INDEPENDENT_DEU_TP_SOURCE_RESULT",json.dumps({k:v for k,v in report.items()
        if k not in ("all_results","limitations")},sort_keys=True),flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())
if __name__=="__main__":
    main()
