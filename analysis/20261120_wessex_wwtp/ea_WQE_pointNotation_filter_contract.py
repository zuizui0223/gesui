#!/usr/bin/env python3
"""Validate official WQE 2026 observation pointNotation filtering in 2017/18.

Use one *real* Bath-region Water Company FINAL treated effluent source point
obtained from official EA 2018 public TP observations. Query pointNotation
for that native point, and for an impossible negative-control notation. Do
NOT export station IDs, names, coordinates, or measured concentration values.
A new national original-site/year census is not authorized unless API
point-identity filtering is proven to be strict.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,urllib.error,urllib.parse,urllib.request
from pathlib import Path
import ea_wqe_2017_2018_SA_tp_siteyear_census as src

BASE="https://environment.data.gov.uk/water-quality/data/observation"
HEADERS={"User-Agent":"gesui-WQE-2026-native-pointNotation-source-contract",
  "API-Version":"1","Accept":"text/csv","CSV-Header":"present",
  "Content-Type":"application/json"}
PARAM_BASE={"dateFrom":"2018-01-01","dateTo":"2018-12-31",
            "determinand":"0348","limit":200,"skip":0}

def query(notation):
    params={**PARAM_BASE,"pointNotation":notation}
    request=urllib.request.Request(
      BASE+"?"+urllib.parse.urlencode(params),
      data=b"",method="POST",headers=HEADERS)
    result={"status":"SOURCE_STOP","http":None,"error":None}
    try:
        with urllib.request.urlopen(request,timeout=45) as response:
            result["http"]=int(response.status)
            raw=response.read(180000)
        if result["http"]==204 or not raw:
            result.update({"status":"NO_SOURCE_MATCH","n":0,"matching":0,"other_notation":0})
            return result
        rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace"))))
        if len(rows)>200:raise RuntimeError("API response exceeds fixed page limit")
        matches=sum(x.get("samplingPoint.notation")==notation for x in rows)
        result.update({"status":"SOURCE_RETURNED_RECORDS",
           "n":len(rows),"matching":matches,"other_notation":len(rows)-matches,
           "sampling_point_type_codes":dict(__import__("collections").Counter(
             x.get("samplingPoint.samplingPointType","") for x in rows)),
           "unit_labels":sorted(set(x.get("unit","") for x in rows)),
           "response_SHA256":hashlib.sha256(raw).hexdigest()})
    except urllib.error.HTTPError as ex:
        result.update({"http":ex.code,"error":"HTTPError "+str(ex)[:220]})
    except Exception as ex:
        result["error"]=type(ex).__name__+": "+str(ex)[:350]
    return result

def main(a):
    official_source=[]
    for offset in range(0,1750,250):
        page,receipt=src.fetch(2018,offset)
        if page is None:raise RuntimeError("Original local EA official 2018 TP page failed")
        for record in page:
            if (record.get("samplingPoint.samplingPointType")==src.SA
                and record.get("determinand.notation")=="0348"
                and record.get("sampleMaterialType")=="FINAL SEWAGE EFFLUENT"):
                notation=record.get("samplingPoint.notation","").strip()
                if notation and notation not in official_source:
                    official_source.append(notation)
        if len(page)<250:break
    if len(official_source)<30:raise RuntimeError("Original source Water Company SA count unexpectedly <30")
    positive=query(sorted(official_source)[0])
    second=query(sorted(official_source)[-1])
    invalid=query("GESUI_SOURCE_NOT_A_REAL_SAMPLING_POINT_999999999")
    valid=(positive.get("n",0)>0 and positive.get("other_notation") == 0
           and second.get("n",0)>0 and second.get("other_notation")==0
           and positive.get("response_SHA256")!=second.get("response_SHA256")
           and invalid.get("n",0)==0
           and invalid.get("status")!="SOURCE_RETURNED_RECORDS"
           and invalid.get("http") in (200,204,400))
    result={"schema":"GESUI_EA_WQE_POST_OBSERVATION_NATIVE_POINTNOTATION_20261120_V1",
            "source_only":True,
            "original_2018_watercompany_SA_tp_station_labels_local":len(official_source),
            "positive_station1":positive,"positive_station2":second,
            "impossible_station_negative_control":invalid,
            "true_pointNotation_filter_verified":bool(valid),
            "no_original_station_codes_values_or_raw_observations_republished":True,
            "cannot_claim_year_specific_STW_at_ChemPop_until_nearby_native_STW_sample_ID_verified":True}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_WQE_POINTNOTATION_ID_FILTER_CONTRACT="+json.dumps(result,sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    main(p.parse_args())
