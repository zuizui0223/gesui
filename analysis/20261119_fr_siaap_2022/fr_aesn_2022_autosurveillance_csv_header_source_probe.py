#!/usr/bin/env python3
"""Read only header bytes from one official AESN monthly autosurveillance CSV.

SOURCE-only pilot, no numeric measurement values processed, no ecological
responses and no full ~100MB download. The CSV URL comes from the official
public data.gouv.fr dataset JSON, not a guessed URL.
"""
import argparse
import csv
import hashlib
import json
import re
from datetime import datetime,timezone
from io import StringIO
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request,urlopen

DATASET_ID="64833d77507be41507260862"
CATALOG=f"https://www.data.gouv.fr/api/1/datasets/{DATASET_ID}/"
MONTH_NAME="autosurveillance-donnees-2022-01.csv"
MAX_CATALOG_BYTES=5_000_000
MAX_HEADER_BYTES=131072

def original_resource_url(metadata):
    if not isinstance(metadata,dict):raise ValueError("Invalid official catalog")
    matches=[r for r in metadata.get("resources",[])
             if r.get("title")==MONTH_NAME or r.get("filename")==MONTH_NAME]
    if len(matches)!=1:raise ValueError("Exactly one original 2022-01 CSV resource needed")
    r=matches[0]
    if r.get("format","").lower()!="csv":
        raise ValueError("Official resource is not CSV")
    url=str(r.get("url") or "")
    parsed=urlparse(url)
    if parsed.scheme!="https" or parsed.netloc!="static.data.gouv.fr" or not parsed.path.endswith("/"+MONTH_NAME):
        raise ValueError("URL origin/path not the original official CSV")
    return url,{"resource_id":r.get("id"),"declared_source_size_bytes":r.get("filesize"),
                "official_title":MONTH_NAME}

def parse_first_header(buffer):
    if not buffer:raise ValueError("No first bytes from official CSV")
    try:
        raw=buffer.decode("utf-8-sig")
        encoding="utf-8-sig"
    except UnicodeDecodeError:
        raw=buffer.decode("cp1252")
        encoding="cp1252"
    first=raw.splitlines()[0]
    delimiter=max((";",",","\t","|"),key=lambda x:first.count(x))
    columns=next(csv.reader(StringIO(first),delimiter=delimiter))
    if not (3<=len(columns)<=150):
        raise ValueError("Unexpected column count in CSV header")
    if any(not c.strip() for c in columns):
        raise ValueError("Blank column names")
    if len(set(columns))!=len(columns):
        raise ValueError("Duplicated column names")
    # Names identify schema candidates only, not actual measured TN/PT data.
    categories={
     "potential_STEU_ID":r"station|steu|ouvrage|sys.?traitement|code.*step|cd.*station|code.?sandre",
     "potential_sampling_date":r"date|jour|mois|annee|année|dateprelev",
     "potential_sandRE_point":r"point|type.?point|code.?point|prelev",
     "potential_NGL_PT_parameter":r"param|polluant|azote|phosphore|ngl|pt",
     "potential_value_unit_and_status":r"valeur|result|resultat|unite|unité|statut|qualite|qualité|flux|debit|débit|concentration"}
    return {
      "columns":columns,"column_count":len(columns),
      "header_encoding":encoding,"detected_delimiter":delimiter,
      "matching_candidate_column_names":{
        k:[c for c in columns if re.search(pattern,c,re.I)]
        for k,pattern in categories.items()},
      "actual_outlet_A4_not_yet_verified":True,
      "2022_SeineAval_specific_records_not_read":True,
      "individual_sample_and_flux_values_not_read":True,
      "annual_measured_NGL_PT_not_calculated":True
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    report={"schema":"gesui_AESN_original_2022_Jan_autosurveillance_bounded_CSV_header_source_v1",
            "official_catalog":CATALOG,"utc":datetime.now(timezone.utc).isoformat(),
            "status":"ACCESS_STOP","measured_annual_SeineAval_N_P_mass_GO":False,
            "biology_EQR_NUTRI_HNS_values_read":False}
    try:
        with urlopen(Request(CATALOG,headers={"Accept":"application/json","User-Agent":"gesui-source-only-aesn-CSV-header/1.0"},method="GET"),timeout=40) as response:
            catalog_bytes=response.read(MAX_CATALOG_BYTES+1)
            if len(catalog_bytes)>MAX_CATALOG_BYTES or response.status!=200:
                raise ValueError("Catalog GET or size error")
        catalog=json.loads(catalog_bytes)
        url,meta=original_resource_url(catalog)
        report["catalog_SHA256"]=hashlib.sha256(catalog_bytes).hexdigest()
        report["original_resource"]=meta
        req=Request(url,headers={
            "Range":f"bytes=0-{MAX_HEADER_BYTES-1}",
            "Accept":"text/csv,application/octet-stream",
            "User-Agent":"gesui-source-only-aesn-CSV-header/1.0"},method="GET")
        with urlopen(req,timeout=50) as response:
            if urlparse(response.url).netloc!="static.data.gouv.fr":
                raise ValueError("Unexpected redirected origin")
            buf=response.read(MAX_HEADER_BYTES)
            receipt={"HTTP_status":response.status,
               "content_type":response.headers.get("Content-Type",""),
               "content_range":response.headers.get("Content-Range"),
               "returned_prefix_bytes":len(buf),
               "returned_prefix_SHA256":hashlib.sha256(buf).hexdigest(),
               "partial_prefix_only_not_original_full_file_SHA":True}
        report["original_CSV_url"]=url
        report["byte_prefix_receipt"]=receipt
        report["header"]=parse_first_header(buf)
        report["status"]="OFFICIAL_RESOURCE_CSV_HEADER_CONFIRMED"
    except Exception as e:
        report["source_error"]=type(e).__name__+": "+str(e)[:250]
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf8")
    print("GESUI_AESN_2022_JAN_CSV_SOURCE_HEADER",json.dumps({
        "status":report["status"],"source_error":report.get("source_error"),
        "receipt":report.get("byte_prefix_receipt"),
        "header":report.get("header")},ensure_ascii=False,sort_keys=True))

if __name__=="__main__":
    main()
