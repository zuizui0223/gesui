#!/usr/bin/env python3
"""New EA WQE 2026 API: actual observed effluent *load* source gate, 2017–18.

Small purposefully bounded official API queries. Validate native determinand
notation, sample type, reporting purpose and units of samples for:
  total P concentration 0348, orthophosphate load as P 6008,
  ammoniacal nitrogen load as N 6003.

This IS NOT a source census of every 2018 WWTP, and contains no site names,
coordinates or chemical values in the public GitHub artifact. Each query
requests at most 20 observations, with source coverage intentionally UNKNOWN
until a true pagination/sample-point design and official permit join exists.
"""
from __future__ import annotations

import argparse,csv,hashlib,io,json,urllib.error,urllib.parse,urllib.request
from pathlib import Path

BASE="https://environment.data.gov.uk/water-quality"
HEADERS={
  "User-Agent":"gesui-WQE-public-bounded-historical-permit-effluent-source/20261120",
  "API-Version":"1","Accept":"text/csv",
  "CSV-Header":"present","Content-Type":"application/json"}
DETERMINANDS={
 "total_phosphorus_as_P_0348":"0348",
 "orthophosphate_load_as_P_6008":"6008",
 "ammoniacal_N_load_as_N_6003":"6003",
 "ammoniacal_N_load_3561":"3561",
}
YEARS=("2017","2018")
LOCATIONS={
 "SouthWest_Bath_50km":{"latitude":"51.38","longitude":"-2.36","radius":"50"},
 "England_first_records":{}
}
MAXROWS=20
MAXBYTES=80000


def query_2018(label,params):
    endpoint=BASE+"/data/observation"
    u=endpoint+"?"+urllib.parse.urlencode(params)
    req=urllib.request.Request(
        u,method="POST",data=b"",
        headers=HEADERS)
    output={
      "label":label,"url":u,
      "status":"SOURCE_REQUEST_STOP",
      "http":None,"content_type":None,"source_rows_requested":MAXROWS,
      "n_source_rows_returned":0,"response_SHA256":None,
      "error":None,
    }
    try:
        with urllib.request.urlopen(req,timeout=30) as resp:
            output["http"]=int(resp.status)
            output["content_type"]=resp.headers.get("Content-Type","")
            data=resp.read(MAXBYTES+1)
            if len(data)>MAXBYTES:raise RuntimeError("Exceeded bounded source probe bytes")
        output["response_SHA256"]=hashlib.sha256(data).hexdigest()
        if output["http"]==204 or not data:
            output["status"]="NO_OBSERVATIONS_IN_SMALL_SOURCE_REQUEST"
            return output
        if output["http"]!=200:raise RuntimeError("Unexpected successful response status")
        txt=data.decode("utf-8-sig","replace")
        rd=csv.DictReader(io.StringIO(txt))
        cols=rd.fieldnames or []
        must={"samplingPoint.notation","samplingPoint.samplingPointType",
              "samplingPurpose","sampleMaterialType",
              "determinand.notation","result","unit"}
        if not must.issubset(set(cols)):
            raise RuntimeError("Original EA observation schema changed")
        rows=list(rd)
        if len(rows)>MAXROWS:raise RuntimeError("Source response exceeded fixed 20-row limit")
        unique_point=set(x["samplingPoint.notation"] for x in rows)
        output.update({
           "status":"ORIGINAL_EA_WQE_LOAD_OBSERVATION_SOURCE_SCHEMA_VERIFIED",
           "source_columns":cols,
           "n_source_rows_returned":len(rows),
           "n_original_distinct_sampling_point_codes_in_sample":len(unique_point),
           "source_determinand_codes":sorted(set(x["determinand.notation"] for x in rows)),
           "source_determinand_names":sorted(set(x.get("determinand.prefLabel","") for x in rows)),
           "source_sampling_purposes":sorted(set(x["samplingPurpose"] for x in rows)),
           "source_sampling_point_type_labels":sorted(set(x["samplingPoint.samplingPointType"] for x in rows)),
           "sample_material_types":sorted(set(x["sampleMaterialType"] for x in rows)),
           "source_unit_labels":sorted(set(x["unit"] for x in rows)),
           "censored_result_count_in_sample":sum(
             str(x["result"]).strip().startswith(("<",">")) for x in rows),
           "this_is_a_small_page_not_source_population_coverage":True,
           "no_original_site_IDs_or_measurements_republished":True,
        })
    except urllib.error.HTTPError as e:
        output.update({"http":e.code,"error":"HTTPError: "+str(e)[:220]})
    except Exception as e:
        output["error"]=type(e).__name__+": "+str(e)[:340]
    return output


def samplingpoint_post():
    url=BASE+"/data/sampling-point?"+urllib.parse.urlencode({"limit":20,"prefLabel":"Sewage"})
    r=urllib.request.Request(
        url,method="POST",data=b"",headers={
         "User-Agent":HEADERS["User-Agent"],
         "API-Version":"1","Accept":"application/ld+json",
         "Content-Type":"application/json"})
    out={"url":url,"http":None,"source_status":"SOURCE_HTTP_STOP"}
    try:
        with urllib.request.urlopen(r,timeout=30) as resp:
            out["http"]=resp.status
            raw=resp.read(MAXBYTES+1)
            if len(raw)>MAXBYTES:raise RuntimeError("Unexpected response too large")
            d=json.loads(raw)
            rows=d.get("member",[])
            out.update({"source_status":"OFFICIAL_SAMPLE_POINT_PAGINATED_SOURCE_ACCESS",
               "n_retrieved_source_members":len(rows) if isinstance(rows,list) else None,
               "response_keys":sorted(d),
               "first_member_keys":sorted(rows[0]) if rows and isinstance(rows[0],dict) else [],
               "source_member_type_values":sorted(set(
                   str(x.get("samplingPointType","")) for x in rows if isinstance(x,dict))),
               "response_SHA256":hashlib.sha256(raw).hexdigest()})
    except urllib.error.HTTPError as e:
        out.update({"http":e.code,"error":"HTTPError: "+str(e)[:220]})
    except Exception as e:
        out["error"]=type(e).__name__+": "+str(e)[:220]
    return out


def main(a):
    out={
      "schema":"GESUI_EA_WQE_2026_HISTORICAL_PN_LOAD_DETERMINAND_PILOT_V1",
      "source_only":True,
      "no_bulk_source_census_or_BMWP_regression":True,
      "source_not_WWTP_to_ChemPop_network_routing":True,
      "pre_search_locked_original_determinand_notations":DETERMINANDS,
      "original_requested_years":list(YEARS),
      "original_source_size_limit_per_query":MAXROWS,
      "probes":{},"point_post_probe":samplingpoint_post(),
    }
    for name,code in DETERMINANDS.items():
      for year in YEARS:
        for region,coord in LOCATIONS.items():
          if region=="England_first_records" and year!="2018":continue
          p={"dateFrom":year+"-01-01","dateTo":year+"-12-31",
             "determinand":code,"limit":MAXROWS,"skip":0}
          p.update(coord)
          label=name+"|"+year+"|"+region
          out["probes"][label]=query_2018(label,p)
    out["status"]="BOUNDED_OFFICIAL_HISTORICAL_LOAD_MONITORING_SOURCE_PILOT"
    out["source_gates"]={
      "actual_data_with_P_load_value_units_and_treated_effluent_sample_type":
        any(v.get("status")=="ORIGINAL_EA_WQE_LOAD_OBSERVATION_SOURCE_SCHEMA_VERIFIED"
            and any("EFFLUENT" in s.upper() for s in v.get("source_sampling_point_type_labels",[]))
            and any("TONNE" in s.upper() or "KILOGRAM" in s.upper()
                    for s in v.get("source_unit_labels",[]))
            for key,v in out["probes"].items() if "load" in key),
      "underlying_national_true_mass_source_support_estimated":False,
      "full_historical_sampling_point_permit_crosswalk_verified":False,
      "river_network_upstream_direction_verified":False,
      "source_resolved_outfall_BMWP_models_authorized":False,
      "restricted_EA_permit_reuse_dataset_not_ingested":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_EA_WQE_PN_LOAD_PILOT="+json.dumps({
      "status":out["status"],
      "point_post":out["point_post_probe"],
      "probes":{
        k:{m:v.get(m) for m in ("http","status","n_source_rows_returned",
           "n_original_distinct_sampling_point_codes_in_sample",
           "source_sampling_point_type_labels","source_sampling_purposes",
           "source_unit_labels","source_determinand_codes",
           "sample_material_types","error")}
        for k,v in out["probes"].items()
      },
      "source_gates":out["source_gates"],
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    main(p.parse_args())
