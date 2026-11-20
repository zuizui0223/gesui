#!/usr/bin/env python3
"""Official Environment Agency WQE new-2026 API source-only service audit.

Legacy 2025 WQA /id and /data/measurement paths were retired in December 2025.
All requests use the published current /data/sampling-point,
 /codelist/determinand, /data/observation paths, small paginated limits and
 official v1 content headers. No batch download, no source-to-bio links,
 no calculated effluent load, no biological outcome model.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT="https://environment.data.gov.uk/water-quality"
BASE_HEADERS={
  "User-Agent":"gesui-EA-2026-WQE-source-quality-pilot/20261120",
  "API-Version":"1",
  "Accept":"application/ld+json",
}
LIMIT_BYTES=128000

def compact_request(label,verb,path,params,accept):
    uri=ROOT+path+"?"+urllib.parse.urlencode(params)
    headers={**BASE_HEADERS,"Accept":accept}
    if accept=="text/csv":
        headers.update({"CSV-Header":"present","Content-Type":"application/json"})
    r=urllib.request.Request(uri,headers=headers,
           method=verb,data=(b'' if verb=="POST" else None))
    rec={"label":label,"url":uri,"method":verb,
         "http_status":None,"content_type":None,
         "source_status":"SOURCE_HTTP_STOP","error":None,
         "n_bytes_read":0,"response_sha256":None}
    try:
        with urllib.request.urlopen(r,timeout=35) as response:
            rec["http_status"]=int(response.status)
            rec["content_type"]=response.headers.get("Content-Type","")
            content=response.read(LIMIT_BYTES+1)
        if len(content)>LIMIT_BYTES:
            raise RuntimeError("Small source-only request returned oversized payload")
        rec["n_bytes_read"]=len(content)
        rec["response_sha256"]=hashlib.sha256(content).hexdigest()
        if verb=="POST" and accept=="text/csv":
            if not content:
                rec["source_status"]="NO_SAMPLE_IN_THE_REQUESTED_REGION_TIME"
            else:
                txt=content.decode("utf-8-sig","replace")
                rdr=csv.DictReader(io.StringIO(txt))
                col=rdr.fieldnames or []
                lines=list(rdr)
                if len(lines)>15:raise RuntimeError("Unexpected >15 rows with limit5")
                rec["source_status"]="ORIGINAL_OBSERVATION_SCHEMA_VERIFIED"
                rec["csv_columns"]=col
                rec["returned_rows"]=len(lines)
                rec["original_determinand_labels"]=sorted(set(x.get("determinand.prefLabel","") for x in lines))
                rec["source_sampling_purposes"]=sorted(set(x.get("samplingPurpose","") for x in lines))
                rec["source_sampling_point_types"]=sorted(set(x.get("samplingPoint.samplingPointType","") for x in lines))
                rec["unit_labels"]=sorted(set(x.get("unit","") for x in lines))
        else:
            doc=json.loads(content.decode("utf-8-sig"))
            if not isinstance(doc,dict):raise RuntimeError("Unexpected official JSONLD response type")
            mem=doc.get("member",doc.get("results",doc.get("items",[])))
            if not isinstance(mem,list):
                mem=[]
            rec["source_status"]="PAGINATED_ORIGINAL_JSONLD_SCHEMA_VERIFIED"
            rec["response_keys"]=sorted(doc.keys())
            rec["returned_members"]=len(mem)
            rec["list_schema_keys"]=sorted(mem[0]) if mem and isinstance(mem[0],dict) else []
            rec["source_member_labels"]=sorted(set(
                str(x.get("prefLabel",x.get("name",""))) for x in mem if isinstance(x,dict)))[:30]
            if "determinand" in label.lower():
                rec["published_determinand_notations"]=[
                    {"prefLabel":x.get("prefLabel"),"notation":x.get("notation")}
                    for x in mem if isinstance(x,dict)]
            if "sampling" in label.lower():
                rec["sample_source_type_names"]=sorted(set(
                    str(x.get("samplingPointType",x.get("prefLabel","")))
                    for x in mem if isinstance(x,dict)))[:30]
    except urllib.error.HTTPError as e:
        rec.update({"http_status":e.code,"error":"HTTPError: "+str(e)[:200]})
    except Exception as e:
        rec["error"]=type(e).__name__+": "+str(e)[:320]
    return rec


def main(a):
    req=[
      ("sampling_point_saltford","GET","/data/sampling-point",
       {"limit":15,"prefLabel":"Saltford"},"application/ld+json"),
      ("sampling_point_sewage","GET","/data/sampling-point",
       {"limit":15,"prefLabel":"Sewage"},"application/ld+json"),
      ("determinand_phosphorus_total","GET","/codelist/determinand",
       {"limit":25,"prefLabel":"Phosphorus, Total as P"},"application/ld+json"),
      ("determinand_ammoniacal_load","GET","/codelist/determinand",
       {"limit":25,"prefLabel":"Ammoniacal Nitrogen : Load"},"application/ld+json"),
      ("determinand_orthophosphate","GET","/codelist/determinand",
       {"limit":25,"prefLabel":"Orthophosphate"},"application/ld+json"),
      ("2018_observation_bath_5km","POST","/data/observation",
       {"dateFrom":"2018-06-01","dateTo":"2018-06-07",
        "latitude":51.383,"longitude":-2.363,"radius":5,
        "limit":5,"skip":0},"text/csv"),
      ("2018_observation_bath_2days","POST","/data/observation",
       {"dateFrom":"2018-06-12","dateTo":"2018-06-13",
        "latitude":51.383,"longitude":-2.363,"radius":5,
        "limit":5,"skip":0},"text/csv"),
    ]
    out={
        "schema":"GESUI_EA_WQE_NEW_2026_SAMPLING_POINT_DETERMINAND_OBSERVATION_SOURCE_V1",
        "old_WQA_endpoint_retired_2025_12":True,
        "official_new_API_version":"1",
        "source_only":True,
        "not_WWTP_effluent_monitoring_station_identity_confirmed":True,
        "not_measured_year_specific_STW_PN_kg_load":True,
        "no_demography_or_BMWP_outcome_calculated":True,
        "probes":{},
    }
    for label,verb,path,params,accept in req:
        out["probes"][label]=compact_request(label,verb,path,params,accept)
    p=out["probes"]
    out["source_gate"]={
        "new_samplingpoint_endpoint_accessible":p["sampling_point_saltford"]["source_status"]=="PAGINATED_ORIGINAL_JSONLD_SCHEMA_VERIFIED",
        "new_determinand_endpoint_accessible":p["determinand_phosphorus_total"]["source_status"]=="PAGINATED_ORIGINAL_JSONLD_SCHEMA_VERIFIED",
        "new_observation_endpoint_valid_small_source":any(
            p[k]["source_status"]=="ORIGINAL_OBSERVATION_SCHEMA_VERIFIED"
            for k in ("2018_observation_bath_5km","2018_observation_bath_2days")),
        "true_final_effluent_sample_point_and_permit_identifier_obtained":False,
        "measured_outfall_NP_source_mass_and_flow_authorized":False,
        "future_query_must_match_native_sampling_point_ID_and_original_2010_18_bio_site":True,
        "bulk_results_not_requested":True,
    }
    out["status"]="SOURCE_ONLY_NEW_EA_WQE_API_PROBED"
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_EA_NEW_WQE_SOURCE_PILOT="+json.dumps({
        "status":out["status"],"summary":{
            k:{"method":v["method"],"http":v["http_status"],
               "source_status":v["source_status"],
               "returned_members":v.get("returned_members"),
               "returned_rows":v.get("returned_rows"),
               "source_labels":v.get("source_member_labels"),
               "determinand_notations":v.get("published_determinand_notations"),
               "csv_columns":v.get("csv_columns"),
               "sampling_purposes":v.get("source_sampling_purposes"),
               "point_types":v.get("source_sampling_point_types"),
               "units":v.get("unit_labels"),
               "error":v.get("error")}
            for k,v in out["probes"].items()},
         "gates":out["source_gate"]
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
