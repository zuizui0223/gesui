#!/usr/bin/env python3
"""AESN official 2022 station-month performance indicators: source schema only.

This is distinct from daily/monthly 100MB raw autosurveillance measurements.
Retrieve official data.gouv.fr catalog, confirm one exact 2022 CSV filename,
then GET ONLY first 128KiB to inspect field names. No plant measurements read,
no treated-effluent N/P mass inferred, no ecological outcomes.
"""
import argparse,csv,hashlib,json,re
from datetime import datetime,timezone
from io import StringIO
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request,urlopen

CATALOG="https://www.data.gouv.fr/api/1/datasets/assainissement-stations-depuration-indicateurs-de-performance/"
FILENAME="stepcol-indicateurs-mensuels-2022.csv"
MAX_CATALOG=2_000_000
PREFIX=131072

def approved_resource(metadata):
    if not isinstance(metadata,dict) or not isinstance(metadata.get("resources"),list):
        raise ValueError("Invalid official source catalog")
    resources=metadata["resources"]
    names=[(r.get("title"),r.get("filename")) for r in resources]
    matches=[r for r in resources if FILENAME in (r.get("title"),r.get("filename"))]
    if len(matches)!=1:raise ValueError("Source 2022 resource not uniquely listed")
    r=matches[0]
    url=str(r.get("url") or "")
    p=urlparse(url)
    if p.scheme!="https" or p.netloc!="static.data.gouv.fr" or not p.path.endswith("/"+FILENAME):
        raise ValueError("Not the officially published 2022 performance CSV URL")
    sz=int(r.get("filesize") or 0)
    if not 25_000_000<sz<90_000_000:
        raise ValueError("Official 2022 CSV size unexpectedly different")
    return {"url":url,"expected_bytes":sz,"id":r.get("id"),
            "official_catalog_titles":[{"title":str(x or ""),"filename":str(y or "")} for x,y in names],
            "resource_format":r.get("format")}

def header(prefix):
    try:
        raw=prefix.decode("utf-8-sig")
        encoding="utf-8-sig"
    except UnicodeDecodeError:
        raw=prefix.decode("cp1252")
        encoding="cp1252"
    line=raw.splitlines()[0]
    delim=max((";",",","\t","|"),key=lambda d:line.count(d))
    cols=next(csv.reader(StringIO(line),delimiter=delim))
    if not 5<=len(cols)<=150 or len(set(cols))!=len(cols) or any(not x.strip() for x in cols):
        raise ValueError("Unrecognized original performance data schema")
    patterns={
      "station_id":r"station|sandre|steu|ouvrage|code.*step",
      "date_and_reporting_period":r"mois|an.*e|year|date|periode|p.*riode",
      "parameter_name":r"param|azote|phosph|ngl|pt|polluant",
      "outlet_inlet":r"sortie|entree|entrée|point|rejet|a3|a4|lieu",
      "unit_mass_flow":r"unit|unité|unite|kg|tonne|flux|masse|debit|débit|volume|concent|charge|valeur",
      "assessment_and_method":r"qual|methode|méthode|statut|rendement|calcul|estim"
    }
    return {"encoding":encoding,"delimiter":delim,"field_count":len(cols),
      "field_names":cols,
      "possible_semantic_columns":{key:[c for c in cols if re.search(p,c,re.I)] for key,p in patterns.items()},
      "year_2022_records_NOT_yet_examined":True,
      "Seine_Aval_exact_station_match_NOT_yet_checked":True,
      "outlet_NGL_PT_flow_weighted_flux_NOT_yet_confirmed":True,
      "N_P_annual_mass_not_measured_from_header":True}

def probe():
    obj={"schema":"gesui_AESN_monthly_station_performance_2022_OFFICIAL_HEADER_ONLY_v1",
         "utc":datetime.now(timezone.utc).isoformat(),
         "official_data_gouv_catalog":CATALOG,
         "source_status":"ACCESS_STOP",
         "original_2022_measurement_rows_opened":False,
         "no_EQR_NUTRI_HNS_used":True,
         "paired_water_land_nutrient_balance_GO":False}
    try:
        req=Request(CATALOG,headers={"Accept":"application/json","User-Agent":"gesui-source-schema-AESN/1.0"},method="GET")
        with urlopen(req,timeout=40) as res:
            body=res.read(MAX_CATALOG+1)
            if res.status!=200 or len(body)>MAX_CATALOG:raise ValueError("Catalog unavailable or oversized")
        metadata=json.loads(body)
        d=approved_resource(metadata)
        obj["catalog_sha256"]=hashlib.sha256(body).hexdigest()
        obj["official_resource"]={k:v for k,v in d.items() if k!="official_catalog_titles"}
        obj["other_original_resource_names"]=d["official_catalog_titles"]
        req=Request(d["url"],headers={"Range":f"bytes=0-{PREFIX-1}","Accept":"text/csv,application/octet-stream",
                          "User-Agent":"gesui-source-schema-AESN/1.0"},method="GET")
        with urlopen(req,timeout=55) as res:
            if urlparse(res.url).netloc!="static.data.gouv.fr":
                raise ValueError("Redirect to non-source host")
            raw=res.read(PREFIX)
            receipt={"status":res.status,"content_range":res.headers.get("Content-Range"),
                     "read_prefix_bytes":len(raw),"sha256_prefix_only":hashlib.sha256(raw).hexdigest()}
        obj["first_bytes_source_receipt"]=receipt
        obj["original_CSV_header"]=header(raw)
        obj["source_status"]="2022_MONTHLY_PERFORMANCE_SOURCE_HEADER_CONFIRMED"
    except Exception as exc:
        obj["source_error"]=type(exc).__name__+": "+str(exc)[:250]
    return obj

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    result=probe()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_AESN_PERFORMANCE_2022_HEADER_SOURCE",json.dumps({
       "status":result["source_status"],"error":result.get("source_error"),
       "header":result.get("original_CSV_header"),
       "url":result.get("official_resource",{}).get("url")},ensure_ascii=False))
