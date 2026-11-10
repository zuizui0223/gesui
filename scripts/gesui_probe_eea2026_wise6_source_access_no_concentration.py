#!/usr/bin/env python3
"""Official EEA 2026 WISE6 source ACCESS/SCHEMA audit, not chemistry values.

Use only the official public Discodata SQL API; query SELECT TOP 1 metadata
columns and official bulk-download HEAD. No TP concentrations, station outcome
values, butterfly observations, or protected monitoring-site data are saved.
If source is not accessible, record error as DATA ACCESS limitation, not zero.
"""
from __future__ import annotations
import argparse,hashlib,json,time,urllib.parse,urllib.request,urllib.error
from datetime import datetime,timezone
from pathlib import Path

BASE="https://discodata.eea.europa.eu"
TABLES={
 "WISE6_disaggregated_latest":"[WISE_SOE].[latest].[Waterbase_T_WISE6_DisaggregatedData]",
 "WISE6_aggregated_latest":"[WISE_SOE].[latest].[Waterbase_T_WISE6_AggregatedData]",
 "WISE6_spatial_latest":"[WISE_SOE].[latest].[Waterbase_S_WISE_SpatialObject_DerivedData]",
 "WISE6_disaggregated_v2r1":"[WISE_SOE].[v2r1].[Waterbase_T_WISE6_DisaggregatedData]",
}
QUERIES={
 "WISE6_disaggregated_latest":"SELECT TOP 1 countryCode, monitoringSiteIdentifier, monitoringSiteIdentifierScheme, parameterWaterBodyCategory, phenomenonTimeSamplingDate FROM [WISE_SOE].[latest].[Waterbase_T_WISE6_DisaggregatedData]",
 "WISE6_aggregated_latest":"SELECT TOP 1 countryCode, monitoringSiteIdentifier, monitoringSiteIdentifierScheme, parameterWaterBodyCategory, phenomenonTimeReferenceYear FROM [WISE_SOE].[latest].[Waterbase_T_WISE6_AggregatedData]",
 "WISE6_spatial_latest":"SELECT TOP 1 countryCode, monitoringSiteIdentifier, monitoringSiteIdentifierScheme, lat, lon FROM [WISE_SOE].[latest].[Waterbase_S_WISE_SpatialObject_DerivedData]",
 "WISE6_disaggregated_v2r1":"SELECT TOP 1 countryCode, monitoringSiteIdentifier, monitoringSiteIdentifierScheme, parameterWaterBodyCategory, phenomenonTimeSamplingDate FROM [WISE_SOE].[v2r1].[Waterbase_T_WISE6_DisaggregatedData]",
}
DOWNLOADS={
 "original_WISE6_disaggregated":"https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_T_WISE6_DisaggregatedData",
 "original_WISE6_aggregated":"https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_T_WISE6_AggregatedData",
 "original_WISE_soE_spatial":"https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_S_WISE_SpatialObject_DerivedData",
}
def get(url,method="GET"):
    source={"original_url":url,"method":method,"http":None}
    for attempt in range(2):
        try:
            req=urllib.request.Request(url,method=method,
                headers={"User-Agent":"gesui-2026-official-EEA-WISE6-public-source-only-2022-TP-QC/1.0",
                    "Accept":"application/json,*/*"})
            with urllib.request.urlopen(req,timeout=32) as r:
                meta={"content_type":r.headers.get("Content-Type"),
                      "content_length_header":r.headers.get("Content-Length"),
                      "last_modified":r.headers.get("Last-Modified"),
                      "ETag":r.headers.get("ETag")}
                raw=r.read(220000) if method=="GET" else b""
                source.update({"http":r.status,"response_metadata":meta,
                               "response_bytes_read":len(raw)})
                if raw:
                    source["source_SHA256_first_read_bytes"]=hashlib.sha256(raw).hexdigest()
                    try:
                        parsed=json.loads(raw)
                        source["valid_JSON_response"]=True
                        source["JSON_root_type"]=type(parsed).__name__
                        source["JSON_root_fields"]=list(parsed.keys())[:20] if isinstance(parsed,dict) else []
                        if isinstance(parsed,list):
                            source["record_count_returned"]=len(parsed)
                            source["first_record_fields"]=list(parsed[0].keys()) if parsed and isinstance(parsed[0],dict) else []
                        if isinstance(parsed,dict):
                            source["first_data_array_field_names"]={
                                k:(list(v[0].keys()) if v and isinstance(v[0],dict) else [])
                                for k,v in parsed.items() if isinstance(v,list)
                                } 
                            source["error_summary"]=str(parsed.get("message") or parsed.get("error") or "")[:250]
                    except Exception:
                        source["valid_JSON_response"]=False
                        source["first_180_response_non_numeric_diagnostic"]=raw[:180].decode("utf8",errors="replace")
                return source
        except urllib.error.HTTPError as e:
            source.update({"http":e.code,"error":"HTTP "+str(e.code)})
            if attempt==1 or e.code not in (429,500,502,503):return source
        except Exception as e:
            source.update({"error":type(e).__name__+": "+str(e)})
            if attempt==1:return source
        time.sleep(1.2)
    return source

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    results={}
    for n,q in QUERIES.items():
        u=BASE+"/sql?"+urllib.parse.urlencode({"query":q,"p":1,"nrOfHits":1})
        x=get(u)
        x["source_SQL_SELECT_only_monitoring_site_metadata"]=q
        results[n]=x
        print("GESUI_2026_EEA_WISE6_SOURCE_METADATA_READ",
              json.dumps({"name":n,"http":x.get("http"),
                 "JSON_root_fields":x.get("JSON_root_fields"),
                 "first_data_array_field_names":x.get("first_data_array_field_names"),
                 "error":x.get("error"),"message":x.get("error_summary")},sort_keys=True),flush=True)
    for n,url in DOWNLOADS.items():
        x=get(url,"HEAD")
        results[n]=x
        print("GESUI_2026_EEA_WISE6_SOURCE_DOWNLOAD_DISCOVERY",
              json.dumps({"name":n,"http":x.get("http"),
                     "content_length_header":x.get("response_metadata",{}).get("content_length_header"),
                     "error":x.get("error")},sort_keys=True),flush=True)
    out={
      "schema":"gesui_EEA_2026_WISE6_2022_TP_original_waterquality_source_access_QC_no_concentrations_v01",
      "retrieved_UTC":datetime.now(timezone.utc).isoformat(),
      "EEA_official_Waterbase2026_landing":"https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1",
      "2026_release_published":"2026-07-07",
      "official_wise_soE_note":"WISE6 Discodata tables may not be updated in current portal according to official freshwater WISE website",
      "source_only_SQL_TOP1_no_values_queried":True,
      "original_source_URL_and_response_SHA_receipts":results,
      "2022_TP_country_station_counts":None,
      "no_numeric_TP_or_other_chemistry_concentrations_accessed":True,
      "no_butterfly_occurrences_or_observer_data_accessed":True,
      "data_access_failure_not_equal_to_no_2022_waterquality_sites":True,
      "do_not_confuse_total_phosphorus_with_phosphate_or_source_reported_P_tonnes":True,
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_2026_EEA_WISE6_PUBLIC_SOURCE_ACCESS_GATE_COMPLETE",flush=True)

if __name__=="__main__":main()
