#!/usr/bin/env python3
"""Authentic Wessex Water 2010–2019 STW effluent chemistry source-only gate.

No facility-to-river attribution, no treatment effect, no simulated source
flow, no biological regression. Reads two PUBLIC Wessex Water effluent CSV
archives from official website; verifies date, site/determinand identifiers
and N/P source record opportunity. Geo XLSX is downloaded ONLY for source
integrity/availability, not parsed as an ecological site/point layer.

Original published data/website may block automated access; in that case
write a transparent source access STOP rather than impute zeros, scrape from
third parties, or assume published titles imply availability.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import io
import json
import math
import re
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

BASE="https://marketplace.wessexwater.co.uk"
FILES={
 "effluent_2010_17":{
     "url":BASE+"/media/d1unvji4/wrcs-effluent-sample-data-2010-2017.csv",
     "format":"csv", "max_bytes":62*1024*1024,"expected_first_year":2010,"expected_last_year":2017},
 "effluent_2018_19":{
     "url":BASE+"/media/ykkldpbz/effluent_2018-19_formatted.csv",
     "format":"csv", "max_bytes":14*1024*1024,"expected_first_year":2018,"expected_last_year":2019},
 "wrc_location_2014_19":{
     "url":BASE+"/media/zh5bthwp/wrc-location-data-1419.xlsx",
     "format":"xlsx_manifest_only","max_bytes":2*1024*1024},
 "small_wrc_phosphorus_removal":{
     "url":BASE+"/media/meukdxb2/small-sites-requiring-p-removal.xlsx",
     "format":"xlsx_manifest_only","max_bytes":2*1024*1024}
}
UA="gesui-scientific-public-wastewater-source-audit/2026-11-20"
NUTRIENT_PATTERN=re.compile(r"phosph|phospha|ortho|nitra|nitrit|ammon|nitrogen|total\s*p\b|total\s*n\b",re.I)
COL_DATE=re.compile(r"(?:sample|sampling|test|collection)?[ _-]*(?:date|datetime|dateandtime|time)",re.I)
COL_SITE=re.compile(r"(?:sample|sampl|wrc|works|location|site).*?(?:id|name|description|code)",re.I)
COL_ANALYTE=re.compile(r"determ|analy|parameter|substance|pollutan",re.I)
COL_RESULT=re.compile(r"result|value|measurement|concentration",re.I)
COL_UNIT=re.compile(r"unit|measunit",re.I)
COL_QUAL=re.compile(r"qualif|less|detect|flag",re.I)
MAX_ROWS=2000000
SAMPLE_PREFIX_LEN=1024*24


def probe_source(name, source, target):
    info={
        "name":name,"url":source["url"],"format":source["format"],
        "status":"SOURCE_ACCESS_STOP","http_code":None,
        "byte_size":None,"sha256":None,"content_type":None,
        "source_page":"https://marketplace.wessexwater.co.uk/dataset/water-recycling-centre-influent-and-effluent-data",
        "error":None,
    }
    try:
        req=urllib.request.Request(source["url"],headers={
            "User-Agent":UA,
            "Accept":"text/csv,application/octet-stream,application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,*/*"
        },method="GET")
        h=hashlib.sha256()
        n=0
        with urllib.request.urlopen(req,timeout=65) as r, target.open("wb") as f:
            info["http_code"]=int(r.status)
            info["content_type"]=r.headers.get("Content-Type")
            if info["http_code"]!=200:
                raise RuntimeError("Original public source did not return HTTP 200")
            length=r.headers.get("Content-Length")
            if length and int(length)>source["max_bytes"]:
                raise RuntimeError("Authentic source exceeds original expected safety cap")
            while True:
                block=r.read(1024*128)
                if not block:break
                n+=len(block)
                if n>source["max_bytes"]:
                    raise RuntimeError("Source exceeded hard download cap")
                h.update(block)
                f.write(block)
        info.update({"byte_size":n,"sha256":h.hexdigest(),
                     "status":"SOURCE_BYTES_OBTAINED"})
        with target.open("rb") as f:
            magic=f.read(8)
        if source["format"]=="csv" and magic.startswith(b"PK"):
            raise RuntimeError("Original CSV URL returned ZIP rather than CSV")
        if source["format"]=="xlsx_manifest_only":
            if not magic.startswith(b"PK\x03\x04"):
                raise RuntimeError("Original XLSX URL did not return OOXML zip")
            info["status"]="XLSX_ORIGINAL_BYTES_VERIFIED_NOT_READ_AS_GEOGRAPHY"
            info["original_geography_still_requires_public_precise_outfall_coords"]=True
        return info
    except urllib.error.HTTPError as e:
        info.update({"http_code":e.code,"error":"HTTPError: "+str(e)[:200]})
    except Exception as e:
        info["error"]=type(e).__name__+": "+str(e)[:320]
    return info


def fieldchoose(headers,pattern):
    matches=[c for c in headers if pattern.search(c)]
    return {"candidates":matches,"first":matches[0] if matches else None}


def date_year(raw):
    v=str(raw or "").strip()
    if not v:return None
    for p in (r"\b(20(?:0[0-9]|1[0-9]|2[0-6]))[-/]",r"[-/](20(?:0[0-9]|1[0-9]|2[0-6]))\b"):
        q=re.search(p,v)
        if q:return int(q.group(1))
    q=re.search(r"\b([0-3]?\d)[-/][0-1]?\d[-/](\d\d)\b",v)
    if q:
        year=2000+int(q.group(2))
        if year<=2026:return year
    return None


def profile_csv(path,label):
    output={"status":"SCHEMA_READ_STOP","source_file":label}
    try:
        raw=path.open("rb").read(SAMPLE_PREFIX_LEN)
        if b"\x00" in raw:
            raise RuntimeError("Original response looks binary, not CSV")
        encoding="utf-8-sig"
        try:decoded=raw.decode(encoding)
        except UnicodeDecodeError:
            encoding="cp1252"
            decoded=raw.decode(encoding)
        # Proper CSV records can have quoted newlines; sniff but do not
        # guess a delimiter from one site id or output values.
        sample=decoded[:16000]
        try: dialect=csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error: dialect=csv.excel
        with path.open("r",encoding=encoding,newline="",errors="strict") as f:
            rd=csv.DictReader(f,dialect=dialect)
            header=rd.fieldnames or []
            if len(header)<3 or len(set(header))!=len(header):
                raise RuntimeError("Unexpected header duplicates / schema")
            columns={
                "site":fieldchoose(header,COL_SITE),
                "sample_date":fieldchoose(header,COL_DATE),
                "analyte":fieldchoose(header,COL_ANALYTE),
                "value":fieldchoose(header,COL_RESULT),
                "units":fieldchoose(header,COL_UNIT),
                "qualifier":fieldchoose(header,COL_QUAL)
            }
            n=0;dates=collections.Counter();analytes=collections.Counter()
            units=collections.Counter();np_analytes=collections.Counter()
            distinct_sites=set();date_parsefail=0;np_numeric_positive=0
            np_numeric_missing=0;unique_chem_site=set()
            for row in rd:
                n+=1
                if n>MAX_ROWS:raise RuntimeError("Source records exceed hard row cap")
                key=columns["site"]["first"]
                if key and row.get(key):
                    distinct_sites.add(str(row[key]).strip())
                dk=columns["sample_date"]["first"]
                if dk:
                    yr=date_year(row.get(dk,""))
                    if yr:dates[str(yr)]+=1
                    else:date_parsefail+=1
                an=columns["analyte"]["first"]
                if an:
                    q=str(row.get(an,"")).strip()
                    if q:analytes[q]+=1
                    if NUTRIENT_PATTERN.search(q):
                        np_analytes[q]+=1
                        if key and row.get(key):
                            unique_chem_site.add(str(row[key]).strip())
                        rk=columns["value"]["first"]
                        if rk:
                            try:
                                val=float(str(row.get(rk,"")).strip().replace(",",""))
                                if math.isfinite(val) and val>=0:np_numeric_positive+=1
                                else:np_numeric_missing+=1
                            except (TypeError,ValueError):np_numeric_missing+=1
                uk=columns["units"]["first"]
                if uk and row.get(uk):
                    units[str(row[uk]).strip()]+=1
        output.update({
            "status":"ORIGINAL_CSV_HEADER_AND_SITE_ANALYTE_RECORD_CENSUS",
            "source_columns_original":header,
            "column_mapping_source_only":columns,
            "original_row_n":n,
            "n_original_distinct_site_labels":len(distinct_sites),
            "sample_year_count":dict(sorted(dates.items())),
            "unparsed_sample_date_rows":date_parsefail,
            "n_distinct_analyte_labels":len(analytes),
            "top_all_analytes":analytes.most_common(18),
            "source_all_nutrient_candidate_labels":dict(sorted(np_analytes.items())),
            "n_distinct_source_chemistry_sites":len(unique_chem_site),
            "numeric_nonnegative_nutrient_record_rows":np_numeric_positive,
            "missing_nonnumeric_nutrient_record_rows":np_numeric_missing,
            "unit_source_values":dict(sorted(units.items())),
            "no_source_discharge_flow_or_mass_calculated":True,
            "does_not_match_to_ChemPop_or_EA_Biosys_sites":True,
        })
    except Exception as e:
        output.update({"status":"SOURCE_SCHEMA_STOP",
                       "error":type(e).__name__+": "+str(e)[:450]})
    return output


def main(a):
    out={
        "schema":"GESUI_WESSEX_WRC_FINAL_EFFLUENT_2010_19_OFFICIAL_SOURCE_V1",
        "source_only_no_ecological_outcomes_opened":True,
        "not_original_2010_19_discharge_mass":True,
        "not_historical_source_resolved_river_receiving_route":True,
        "supplier_dataset_CC_BY_published":True,
        "site_geo_original_xlsx_manifest_only_no_columns_inferred":True,
        "daily_flow_dataset_publisher_currently_removed_for_update":True,
        "links": {
            "source_page":"https://marketplace.wessexwater.co.uk/dataset/water-recycling-centre-influent-and-effluent-data",
            "geography_page":"https://marketplace.wessexwater.co.uk/dataset/water-recycling-centre-geographic-data",
            "flow_page":"https://marketplace.wessexwater.co.uk/dataset/water-recycling-centre-flow-data",
            "riverhub_original_discharge_location_information":"https://riverhub.co.uk/monitoring",
            "EA_authoritative_water_quality_explorer":"https://environment.data.gov.uk/dataset/8034d47a-ba8a-4978-aca0-bd5d6e870536",
        },
        "source_receipts":{},
        "csv_source_profiles":{},
        "biology_stw_exposure_model_authorized":False,
        "source_gates":{},
    }
    with tempfile.TemporaryDirectory(prefix="gesui_wessex_source_") as temp:
        for name,setting in FILES.items():
            f=Path(temp)/(name+"."+(setting["format"].split("_")[0]))
            receipt=probe_source(name,setting,f)
            out["source_receipts"][name]=receipt
            if setting["format"]=="csv" and receipt["status"]=="SOURCE_BYTES_OBTAINED":
                out["csv_source_profiles"][name]=profile_csv(f,name)
    observed={k:v["status"] for k,v in out["source_receipts"].items()}
    by_csv=[out["csv_source_profiles"].get(k,{}) for k in ("effluent_2010_17","effluent_2018_19")]
    enough=all(x.get("status")=="ORIGINAL_CSV_HEADER_AND_SITE_ANALYTE_RECORD_CENSUS" for x in by_csv)
    out["source_gates"]={
        "both_original_effluent_2010_2019_CSV_verified":enough,
        "wrc_geography_original_bytes_verified":
            out["source_receipts"]["wrc_location_2014_19"]["status"]=="XLSX_ORIGINAL_BYTES_VERIFIED_NOT_READ_AS_GEOGRAPHY",
        "original_outfall_point_and_downstream_network_route_verified":False,
        "native_WRC_SampleSiteId_to_outfall_and_discharge_plant_year_verified":False,
        "measured_treatment_effluent_nutrient_concentrations_origin_confirmed":enough,
        "original_measured_daily_final_outfall_flow_validated":False,
        "WWTP_P_or_N_kg_year_causal_source_exposure_authorized":False,
        "Wessex_region_source_can_be_compared_with_ChemPop_only_after_original_genuine_upstream_routes":False,
    }
    out["status"]="ORIGINAL_SOURCE_ACCESS_AUDITED_NO_BIOLOGICAL_SLOPES"
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_WESSEX_WRC_SOURCE_ONLY="+json.dumps({
        "status":out["status"],"source_status":observed,
        "source_errors":{k:v.get("error") for k,v in out["source_receipts"].items()},
        "CSV_profiles":{k:{m:v.get(m) for m in ("status","source_columns_original",
             "original_row_n","n_original_distinct_site_labels",
             "sample_year_count","source_all_nutrient_candidate_labels",
             "n_distinct_source_chemistry_sites","numeric_nonnegative_nutrient_record_rows",
             "missing_nonnumeric_nutrient_record_rows","unit_source_values","error")}
            for k,v in out["csv_source_profiles"].items()},
        "gates":out["source_gates"]
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
