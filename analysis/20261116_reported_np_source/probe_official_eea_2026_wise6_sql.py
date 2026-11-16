#!/usr/bin/env python3
"""Bounded public EEA WISE6 Discodata 2026 station TN/TP access probe.

Source-only. At most three tiny (TOP 1/2) SQL GETs; JSON response capped 120 KB;
never fetch the 4.3 GB full ICM archive or rehost observation values. Official
public web documentation describes Discodata SQL JSON API; WISE6 data freshness
is separately warned about. A failed query is UNKNOWN, not absence of water N/P.
"""
from __future__ import annotations
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import requests

BASE="https://discodata.eea.europa.eu/sql"
OFFICIAL_DATAHUB="https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1"
OFFICIAL_2026_DIRECT="https://sdi.eea.europa.eu/data/82e194cc-2506-411c-b341-9eb50b6ad5e0"
MAX_BYTES=120000
AGG="[WISE_SOE].[latest].[Waterbase_T_WISE6_AggregatedData]"
AGG_OLD="[WISE_SOE].[v1r1].[Waterbase_T_WISE6_AggregatedData]"
DISAGG="[WISE_SOE].[latest].[Waterbase_T_WISE6_DisaggregatedData]"
QUERIES=[
("aggregated_latest_schema",f"SELECT TOP 1 * FROM {AGG}"),
("aggregated_latest_2018_2022_TNTP",
 f"SELECT TOP 2 monitoringSiteIdentifier, phenomenonTimeReferenceYear, observedPropertyDeterminandLabel, resultMeanValue FROM {AGG} WHERE phenomenonTimeReferenceYear IN (2018,2022) AND observedPropertyDeterminandLabel IN ('Total nitrogen','Total phosphorus')"),
("disaggregated_latest_schema",f"SELECT TOP 1 * FROM {DISAGG}"),
]
HEADERS={"User-Agent":"gesui-20261116-EEA-WISE6-data-availability-probe/0.1",
         "Accept":"application/json,*/*"}
def source_sql_probe(name,sql):
    item={"id":name,
          "SQL_source_sha256":hashlib.sha256(sql.encode()).hexdigest(),
          "records_requested_max":2,
          "source_URL_origin":BASE,
          "not_exhaustive_observations":True}
    try:
        resp=requests.get(BASE,params={"query":sql,"p":1,"nrOfHits":2},
              headers=HEADERS,timeout=(10,24),stream=True)
        item["http_status"]=int(resp.status_code)
        item["content_type"]=resp.headers.get("Content-Type")
        payload=bytearray()
        for chunk in resp.iter_content(chunk_size=32768):
            payload.extend(chunk)
            if len(payload)>MAX_BYTES:raise RuntimeError("API response exceeds source probe cap")
        item["response_bytes"]=len(payload)
        item["response_sha256"]=hashlib.sha256(payload).hexdigest()
        try:
            obj=json.loads(payload)
        except Exception:
            item["result_kind"]="not_json"
            item["possible_html_or_portal_response"]=True
            return item
        if not isinstance(obj,dict):
            item["result_kind"]="unexpected_json_type";return item
        rows=obj.get("results")
        errors=obj.get("errors")
        if isinstance(rows,list):
            item["result_kind"]="result_rows"
            item["returned_rows_at_most_2"]=len(rows)
            item["first_row_field_names"]=list(rows[0]) if rows and isinstance(rows[0],dict) else []
            # Deliberately omit real individual station ID, dates, N and P values.
            item["original_observation_values_exported"]=False
        elif isinstance(errors,list):
            item["result_kind"]="SQL_API_errors"
            item["error_codes"]=[x.get("errorcode") for x in errors if isinstance(x,dict)]
            item["error_names"]=[str(x.get("error",""))[:180] for x in errors if isinstance(x,dict)]
        else:
            item["result_kind"]="unexpected_JSON_schema"
            item["top_level_keys"]=list(obj)[:25]
        return item
    except Exception as e:
        item["result_kind"]="request_failed"
        item["error"]=type(e).__name__+": "+str(e)[:180]
        return item

def original_download_metadata():
    out={"url":OFFICIAL_2026_DIRECT}
    try:
        q=requests.head(OFFICIAL_2026_DIRECT,
              headers={"User-Agent":HEADERS["User-Agent"]},timeout=(10,24),
              allow_redirects=True)
        out.update(http_status=int(q.status_code),
                   final_url=q.url.split("?")[0],
                   content_type=q.headers.get("Content-Type"),
                   size_header=q.headers.get("Content-Length"))
        out["HEAD_metadata_only_no_4GB_download"]=True
    except Exception as e:
        out["status"]="HEAD_FAILED"
        out["error"]=type(e).__name__+": "+str(e)[:180]
    return out

def main(a):
    outputs=[source_sql_probe(name,q) for name,q in QUERIES]
    has_schema=any(x.get("result_kind")=="result_rows" and x.get("returned_rows_at_most_2",0)>0 for x in outputs)
    has_year_TNTP=any(x["id"]=="aggregated_latest_2018_2022_TNTP" and
          x.get("result_kind")=="result_rows" and x.get("returned_rows_at_most_2",0)>0
          for x in outputs)
    out={"schema":"GESUI_EEA_2026_WISE6_ICM_BOUNDED_SQL_AVAILABILITY_SOURCE_GATE_V1",
         "EEA_official_2026_ICM_page":OFFICIAL_DATAHUB,
         "EEA_direct_dataset_overall_size_reported_public_page":"4.3 GB",
         "public_discordata_SQL_API_HELP":"https://discodata.eea.europa.eu/Help.html",
         "original_2026_icm_metadata":original_download_metadata(),
         "bounded_tests":outputs,
         "at_least_one_published_sql_view_readable":has_schema,
         "target_year_2018_or_2022_TNTP_aggregated_row_probe_success":has_year_TNTP,
         "original_source_station_identifiers_FROM_183H7_NOT_YET_MATCHED":True,
         "direct_matching_of_EEA_biology_station_and_chemistry_NOT_DONE":True,
         "API_SQL_failure_NOT_absence_of_measured_TNTP":True,
         "source_probe_does_not_classify_positive_or_negative_N_pollution":True,
         "WISE6_dataset_broad_release_NOT_proof_of_2026_DISCODATA_freshness":True,
         "no_large_archive_downloaded":True,
         "no_biological_response_or_causal_effect_fitted":True}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"EEA_2026_WISE6_TNTP_SQL_ACCESS_SOURCE_ONLY.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_EEA_2026_WISE6_SQL_SOURCE_PROBE="+json.dumps(out,ensure_ascii=False,sort_keys=True),flush=True)
    print("SITE_YEAR_TNTP_LINK_GATE=STOP_NOT_YET_SOURCE_IDENTIFIED",flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
