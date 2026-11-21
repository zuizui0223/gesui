#!/usr/bin/env python3
"""Official UK PRTR 2017/2018 original CSV source-only schema/identity inventory.

This is distinct from EA WQE treated-effluent sample concentrations and from
EEA UWWTD2018 strict positive P/N outlet cohort. It is an operator-reported
FACILITY/PERMIT annual release register, with potentially diffuse releases,
method tags and pollutant-specific reporting thresholds. It is NOT an
observed river outfall mass without permit/outfall crosswalk.

Original official 2026-publication CSV URLs from GOV.UK:
https://www.gov.uk/government/publications/download-industrial-emissions-and-waste-transfer-uk-prtr-data
The official 2017/18 source CSV may be accessible while separate EA permit
register or original DEFRA ODS downloads are 403. Only original source.
No facility-level source rows are republished to Git; small aggregate manifest
only. Missing pollutant rows are NOT zeros.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,urllib.error,urllib.request
from pathlib import Path
import tempfile

SOURCES={
 "2017":"https://assets.publishing.service.gov.uk/media/6a63303291d54bc7d5127f24/uk_prtr_dataset_2017.csv",
 "2018":"https://assets.publishing.service.gov.uk/media/6a633032243ee000f7127f31/uk_prtr_dataset_2018.csv"
}
MAX_SOURCE_BYTES=12*1024*1024
MAX_ROWS=300000
HEADERS={"User-Agent":"gesui-official-UK-PRTR-2017-2018-original-source/20261121",
         "Accept":"text/csv,application/octet-stream,*/*"}
CANDIDATE_TOKENS=("year","report","plant","site","facility","permit","operator",
    "postcode","region","grid","east","north","lat","lon",
    "activity","sector","industry","waste","sewage","treatment",
    "medium","releas","pollutan","substan","determinand","method",
    "amount","quantity","unit","water","mass","measurement")
SAFE_CATEGORICAL=("media","medium","activity","sector","industry","release",
    "method","pollutant","substance","unit","compartment","waste")

def download(url,dest):
    x={"http":None,"content_type":None,"status":"SOURCE_ACCESS_STOP"}
    try:
        req=urllib.request.Request(url,headers=HEADERS)
        sha=hashlib.sha256();n=0
        with urllib.request.urlopen(req,timeout=75) as r,dest.open("wb") as f:
            x["http"]=int(r.status)
            x["content_type"]=r.headers.get("Content-Type","")
            if r.status!=200:raise RuntimeError("Official UK PRTR original source not HTTP200")
            while True:
                chunk=r.read(65536)
                if not chunk:break
                n+=len(chunk)
                if n>MAX_SOURCE_BYTES:
                    raise RuntimeError("Official source exceeded bounded original CSV download")
                sha.update(chunk);f.write(chunk)
        if n<1000:raise RuntimeError("Official 2017/18 CSV suspiciously small")
        x.update({"status":"OFFICIAL_ORIGINAL_CSV_BYTES_OBTAINED",
                  "source_bytes":n,"source_sha256":sha.hexdigest()})
        return x
    except urllib.error.HTTPError as e:
        x.update({"http":e.code,"error":"HTTPError: "+str(e)[:230]})
    except Exception as e:
        x["error"]=type(e).__name__+": "+str(e)[:400]
    return x

def inspect(path,year):
    out={"status":"SCHEMA_STOP"}
    try:
        with path.open("rb") as f:head=f.read(8192)
        if b"\x00" in head or head[:20].lower().lstrip().startswith(b"<html"):
            raise RuntimeError("Official source not plain CSV")
        try:
            text=head.decode("utf-8-sig")
            encoding="utf-8-sig"
        except UnicodeDecodeError:
            text=head.decode("cp1252")
            encoding="cp1252"
        try:
            dialect=csv.Sniffer().sniff(text,delimiters=",;\t")
        except csv.Error:
            dialect=csv.excel
        fieldnames=None;rows_n=0
        unique=collections.defaultdict(collections.Counter)
        distinct_source_year=collections.Counter()
        numeric_likely=collections.Counter()
        original_fields=[]
        with path.open("r",encoding=encoding,errors="replace",newline="") as f:
            reader=csv.DictReader(f,dialect=dialect)
            fieldnames=list(reader.fieldnames or [])
            if len(fieldnames)<4 or len(fieldnames)>150 or len(set(fieldnames))!=len(fieldnames):
                raise RuntimeError("Official CSV has missing/ambiguous original field names")
            visible=[c for c in fieldnames if any(t in c.lower() for t in CANDIDATE_TOKENS)]
            for row in reader:
                rows_n+=1
                if rows_n>MAX_ROWS:raise RuntimeError("Original source exceeded bounded row census")
                for col in fieldnames:
                    lc=col.lower()
                    v=str(row.get(col,"") or "").strip()
                    if not v:continue
                    if ("year" in lc and len(v)<=20):
                        distinct_source_year[col+"|"+v]+=1
                    if (any(token in lc for token in SAFE_CATEGORICAL) and len(v)<65):
                        unique[col][v]+=1
                    if ("amount" in lc or "quantity" in lc or "mass" in lc):
                        try:
                            fval=float(v.replace(",",""))
                            if fval>=0 and fval<1e20:numeric_likely[col]+=1
                        except (ValueError,TypeError):pass
        categ={}
        for col,ct in unique.items():
            categ[col]={"distinct_labels":len(ct),
               "top15":ct.most_common(15) if len(ct)<=1000 else []}
        out.update({
          "status":"OFFICIAL_SOURCE_CSV_SCHEMA_AND_ORIGINAL_CATEGORICAL_LABELS_VERIFIED",
          "csv_encoding":encoding,"csv_delimiter":dialect.delimiter,
          "original_csv_header":fieldnames,
          "original_row_count":rows_n,
          "source_fields_relevant_to_facility_PN_water":visible,
          "source_year_candidate_counts":dict(distinct_source_year),
          "limited_source_category_census":categ,
          "nonnegative_numeric_mass_candidate_counts":dict(numeric_likely),
          "no_plant_level_source_mass_or_identity_republished":True,
          "does_not_assume_all_water_release_rows_are_WWTP":True,
          "publisher_reported_mass_not_independent_measured_outfall_mass":True
        })
    except Exception as e:
        out.update({"status":"SOURCE_SCHEMA_STOP","error":type(e).__name__+": "+str(e)[:450]})
    return out

def main(a):
    result={"schema":"GESUI_UK_PRTR_2017_18_OFFICIAL_CSV_SOURCE_SCHEMA_ONLY_20261121_V1",
      "publisher":"UK DEFRA PRTR annual 2017/18 (published 2026)",
      "official_metadata":"https://www.gov.uk/government/publications/download-industrial-emissions-and-waste-transfer-uk-prtr-data",
      "original_facility_reporting_threshold_high_WWTW_GE_100000_PE":True,
      "facility_permit_annual_mass_can_include_fugitive_or_nonpoint":True,
      "source_only":True,"not_WWTP_outfall_river_routed":True,
      "not_independent_proof_of_treatment_plant_effect_on_ChemPop":True,
      "no_original_rows_or_facility_identifiers_republished":True,
      "source_by_year":{}}
    with tempfile.TemporaryDirectory(prefix="gesui_prtr_original_source_") as temp:
        for yr,url in SOURCES.items():
            p=Path(temp)/(yr+".csv")
            receipt=download(url,p)
            if receipt["status"]=="OFFICIAL_ORIGINAL_CSV_BYTES_OBTAINED":
                receipt["source_csv_schema"]=inspect(p,yr)
            result["source_by_year"][yr]={"official_url":url,**receipt}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_UK_PRTR_2017_18_SOURCE_SCHEMA="+json.dumps({
      yr:{"status":d["status"],"http":d.get("http"),"error":d.get("error"),
          "source_bytes":d.get("source_bytes"),"sha256":d.get("source_sha256"),
          "source_schema":d.get("source_csv_schema",{})}
      for yr,d in result["source_by_year"].items()
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    main(ap.parse_args())
