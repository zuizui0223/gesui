#!/usr/bin/env python3
"""Test whether EA WQE 2026 sampling-point API actually FILTERS WaterCompany SA.

The API is paginated and source type is nested JSON-LD samplingPointType with
notation SA. Probe filters with SA and negative-control SB; assume NOTHING
when unknown request parameters are silently ignored. Only schema/type
counts are emitted, no original station addresses or IDs.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,urllib.error,urllib.parse,urllib.request
from pathlib import Path

HOST="https://environment.data.gov.uk/water-quality"
LIMIT=100
UA="gesui-source-only-EA-WQE-samplingPointType-contract-pilot-20261120"
PARAMS={
  "baseline":{},
  "sa_type_plain":{"samplingPointType":"SA"},
  "sb_type_plain_negative":{"samplingPointType":"SB"},
  "sa_type_nested":{"samplingPointType.notation":"SA"},
  "sa_type_prefLabel":{"samplingPointType":"SEWAGE DISCHARGES - FINAL/TREATED EFFLUENT - WATER COMPANY"},
  "sewage_label":{"prefLabel":"Sewage"}
}
EXPECTED_SA="SA"
EXPECTED_SB="SB"


def scan(key,extra):
    query={"limit":LIMIT,"skip":0,**extra}
    url=HOST+"/data/sampling-point?"+urllib.parse.urlencode(query)
    req=urllib.request.Request(url,method="POST",data=b"",headers={
       "User-Agent":UA,"API-Version":"1","Accept":"application/ld+json",
       "Content-Type":"application/json"})
    out={"status":"SOURCE_STOP","url":url,"HTTP":None,"error":None}
    try:
        with urllib.request.urlopen(req,timeout=45) as resp:
            out["HTTP"]=int(resp.status)
            raw=resp.read(250000)
        if len(raw)>=250000:raise RuntimeError("Too many source metadata bytes")
        data=json.loads(raw)
        members=data.get("member",[])
        if not isinstance(members,list):
            raise RuntimeError("Expected JSON-LD member array")
        hist=collections.Counter()
        samples={}
        for m in members:
            k=m.get("samplingPointType",{})
            if isinstance(k,dict):
                code=str(k.get("notation",""))
            else:code=str(k)
            hist[code]+=1
            if code not in samples:
                samples[code]={"type":k.get("prefLabel") if isinstance(k,dict) else str(k)}
        out.update({
          "status":"PUBLISHED_METADATA_PAGE_VERIFIED",
          "sampled_member_count":len(members),
          "totalItems":data.get("totalItems"),
          "type_notation_histogram":dict(hist),
          "sampled_type_labels":samples,
          "responseSHA256":hashlib.sha256(raw).hexdigest(),
          "no_original_sampling_point_ids_rehosted":True
        })
    except urllib.error.HTTPError as e:
        out.update({"HTTP":e.code,"error":"HTTPError "+str(e)[:220]})
    except Exception as e:
        out["error"]=type(e).__name__+": "+str(e)[:260]
    return out


def main(a):
    results={k:scan(k,v) for k,v in PARAMS.items()}
    base=results["baseline"]
    sa=results["sa_type_plain"]
    sb=results["sb_type_plain_negative"]
    nested=results["sa_type_nested"]
    sa_pass=(
      base.get("status")=="PUBLISHED_METADATA_PAGE_VERIFIED"
      and sa.get("status")=="PUBLISHED_METADATA_PAGE_VERIFIED"
      and sb.get("status")=="PUBLISHED_METADATA_PAGE_VERIFIED"
      and sa.get("sampled_member_count",0)>0
      and sb.get("sampled_member_count",0)>0
      and list(sa.get("type_notation_histogram",{}))==[EXPECTED_SA]
      and list(sb.get("type_notation_histogram",{}))==[EXPECTED_SB]
      and base.get("totalItems")!=sa.get("totalItems")
      and sa.get("totalItems")!=sb.get("totalItems")
    )
    nested_pass=(nested.get("sampled_member_count",0)>0 and
                 list(nested.get("type_notation_histogram",{}))==[EXPECTED_SA] and
                 nested.get("totalItems")!=base.get("totalItems"))
    out={
       "schema":"GESUI_WQE_2026_SA_SAMPLINGPOINT_FILTER_API_CONTRACT_V1",
       "source_only":True,
       "results":results,
       "filter_gates":{
          "plain_samplingPointType_SA_vs_SB_actual_filter_verified":bool(sa_pass),
          "nested_samplingPointType_notation_SA_verified":bool(nested_pass),
          "is_full_nationwide_original_SA_point_registry_authorized":bool(sa_pass or nested_pass),
          "unrecognized_query_filter_parameter_must_not_be_treated_as_applied":True,
          "no_permit_link_or_outfall_chemistry_route_verified":True
       }
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_WQE_SA_TYPE_FILTER_CONTRACT="+json.dumps({
       "filters":out["filter_gates"],
       "probes":{k:{x:v.get(x) for x in ("HTTP","status","sampled_member_count","totalItems","type_notation_histogram","error")}
                 for k,v in results.items()}
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
