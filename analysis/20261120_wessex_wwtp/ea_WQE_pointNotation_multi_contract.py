#!/usr/bin/env python3
"""Verify official EA WQE 2026 multi pointNotation syntax and strict filtering.

Two independent, genuine 2018 Water Company final-treated effluent SA
sampling-point codes are recovered from the already source-verified Bath
2018 TP sample. Query their union both as repeated pointNotation query
parameters and as a comma-separated value; retain whichever official API
actually selects ONLY those two and returns both, matching exact original
single-site observations. Never publish original sampling point IDs, raw
coordinates or chemistry values; outputs are QA counts/hashed source only.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,urllib.parse,urllib.error,urllib.request
from pathlib import Path
import ea_wqe_2017_2018_SA_tp_siteyear_census as ea

BASE="https://environment.data.gov.uk/water-quality/data/observation"
HEADERS={"User-Agent":"gesui-EA-WQE2026-verified-multi-STW-notation-test",
    "API-Version":"1","Accept":"text/csv","CSV-Header":"present",
    "Content-Type":"application/json"}
SOURCE_LIMIT=250
SOURCE_MAX_BYTES=180000

def get_notations():
    codes=set()
    for skip in range(0,2000,ea.LIMIT):
        rows,receipt=ea.fetch(2018,skip)
        if rows is None:raise RuntimeError("Original verified EA year source page unavailable")
        for row in rows:
            if (row["samplingPoint.samplingPointType"]==ea.SA
                and row["determinand.notation"]=="0348"
                and row["sampleMaterialType"]=="FINAL SEWAGE EFFLUENT"):
                code=row["samplingPoint.notation"].strip()
                if code:codes.add(code)
        if len(rows)<ea.LIMIT:break
    if len(codes)!=49:
        raise RuntimeError("2018 Bath original 49 EA native SA TP source point count changed")
    ordered=sorted(codes)
    return [ordered[0], ordered[-1]]

def call(params):
    endpoint=BASE+"?"+urllib.parse.urlencode(params,doseq=True)
    rec={"status":"SOURCE_STOP","url_query_mode":
         "repeated" if isinstance(params.get("pointNotation"),list) else "single_or_comma",
         "HTTP":None,"reason":None}
    try:
        req=urllib.request.Request(endpoint,headers=HEADERS,method="POST",data=b"")
        with urllib.request.urlopen(req,timeout=40) as response:
            rec["HTTP"]=response.status
            data=response.read(SOURCE_MAX_BYTES+1)
        if len(data)>SOURCE_MAX_BYTES:raise RuntimeError("Source page exceeds bound")
        rec["response_sha256"]=hashlib.sha256(data).hexdigest()
        if not data:
            rec.update({"status":"EMPTY_SOURCE_OBSERVATION","n_source_rows":0})
            return rec,set(),{}
        rows=list(csv.DictReader(io.StringIO(data.decode("utf-8-sig","replace"))))
        if len(rows)>SOURCE_LIMIT:raise RuntimeError("Returned above fixed 250 rows")
        ids=collections.Counter(row.get("samplingPoint.notation","").strip() for row in rows)
        for row in rows:
            if row.get("determinand.notation")!="0348" or not row.get("phenomenonTime","").startswith("2018"):
                raise RuntimeError("Selected wrong original year or material determinand")
        rec.update({"status":"ORIGINAL_SOURCE_OBSERVATIONS_RETURNED",
             "n_source_rows":len(rows),"unique_native_point_codes":len(ids),
             "rows_per_native_point_hashed":hashlib.sha256(
                repr(sorted(ids.values())).encode()).hexdigest()})
        return rec,set(ids),dict(ids)
    except urllib.error.HTTPError as e:
        rec.update({"HTTP":e.code,"reason":"HTTPError "+str(e)[:220]})
    except Exception as e:
        rec["reason"]=type(e).__name__+": "+str(e)[:350]
    return rec,set(),{}

def main(a):
    codes=get_notations()
    parameters={"dateFrom":"2018-01-01","dateTo":"2018-12-31",
                "determinand":"0348","limit":SOURCE_LIMIT,"skip":0}
    singles=[]
    for code in codes:
        r,returned,counts=call({**parameters,"pointNotation":code})
        if r["status"]!="ORIGINAL_SOURCE_OBSERVATIONS_RETURNED" or returned!={code}:
            raise RuntimeError("Original pointNotation single STW source not exclusive")
        singles.append((r,counts))
    wanted=set(codes)
    n_expected=sum(sum(c.values()) for _,c in singles)
    mode_result={}
    for label,filtervals in (
        ("repeated_parameters",codes),
        ("comma_joined_parameter",",".join(codes)),
    ):
        rec,ids,counts=call({**parameters,"pointNotation":filtervals})
        correct=(ids==wanted and len(counts)==2
                 and sum(counts.values())==n_expected
                 and all(counts[k]==sum(s_count.get(k,0) for _,s_count in singles)
                         for k in wanted))
        mode_result[label]={**rec,
            "exact_original_singletons_union_selected":bool(correct),
            "no_other_site_code_returned":bool(ids.issubset(wanted)),
            "n_expected_from_two_individual_points":n_expected}
    data={
        "schema":"GESUI_EA_WQE2026_MULTI_POINTNOTATION_ORIGINAL_STW_FILTER_CONTRACT_V1",
        "source_only":True,
        "original_Bath2018_TP_SA_source_native_points":49,
        "n_verified_input_points":len(codes),
        "input_original_points_SHA256":hashlib.sha256(
            "\n".join(sorted(codes)).encode()).hexdigest(),
        "single_points":{
            "each_selected_only_one_original_native_point":True,
            "original_counts":[v["n_source_rows"] for v,_ in singles],
        },
        "multiple_point_filter_modes":mode_result,
        "valid_batch_mode":[k for k,v in mode_result.items()
                            if v["exact_original_singletons_union_selected"]],
        "authorization_to_filter_unobserved_STW_and_true_year_source":any(
            v["exact_original_singletons_union_selected"] for v in mode_result.values()),
        "not_confirmation_of_2017_18_activity_of_all_current_5722_SA_registry":True,
        "not_source_to_receptor_flow_network":True
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(data,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("GESUI_EA_WQE_MULTI_POINTNOTATION_CONTRACT="+json.dumps(
        data,sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
