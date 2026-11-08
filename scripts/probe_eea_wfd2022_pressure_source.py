#!/usr/bin/env python3
"""Source metadata probe for WFD-2022 reported agricultural/urban/hydromorphologic pressures.

Read only public EEA official schema, bounded preview and total counts.
No numerical 2022 biological EQR, WWTP-effect residuals or new covariate
selection based on ecological signs. Pin source layer URL and response hashes.
"""
from __future__ import annotations
import argparse,hashlib,json
from datetime import datetime,timezone
from pathlib import Path
import requests

SQL="https://discodata.eea.europa.eu/sql"
PRESSURE_TABLE="[WISE_WFD].[latest].[SWB_SurfaceWaterBody_swSignificantPressureType]"
ARCGIS=(
 "https://water.discomap.eea.europa.eu/arcgis/rest/services/WISE_WFD/"
 "WFD2022_SurfaceWaterBody_swSignificantPressureType_WM/MapServer"
)
HEADERS={"User-Agent":"gesui-EEA-Biology2022-independent-WFD-pressure-source-schema/0.1",
         "Accept":"application/json"}
MAX_BYTES=1_500_000

def req(url,params):
    try:
        r=requests.get(url,params=params,headers=HEADERS,timeout=(12,70))
        raw=r.content[:MAX_BYTES+1]
        if len(raw)>MAX_BYTES:raise ValueError("Exceeded metadata response bound")
        response={"request_url":r.url,"HTTP_status":r.status_code,
           "response_bytes":len(raw),"response_SHA256":hashlib.sha256(raw).hexdigest(),
           "content_type":r.headers.get("Content-Type")}
        try:obj=json.loads(raw)
        except Exception:
            response["response_preview"]=raw[:350].decode("utf8","replace")
            return None,response
        if isinstance(obj,dict) and "error" in obj:
            response["remote_error"]=str(obj["error"])[:400]
        return obj,response
    except Exception as e:
        return None,{"request_url":url,"error":type(e).__name__+": "+str(e)[:300]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    sources={}
    table_queries={
       "WFD_latest_pressure_first_rows":f"SELECT TOP 4 * FROM {PRESSURE_TABLE}",
       "WFD_latest_pressure_total":f"SELECT COUNT_BIG(*) AS source_record_count FROM {PRESSURE_TABLE}",
    }
    for tag,q in table_queries.items():
        v,receipt=req(SQL,{"query":q,"p":1,"nrOfHits":5})
        if isinstance(v,dict):
            body=v.get("results",[])
            receipt["result_keys"]=list(v.keys())[:20]
            receipt["record_count_returned"]=len(body) if isinstance(body,list) else None
            receipt["sample_field_names"]=list(body[0].keys()) if isinstance(body,list) and body and isinstance(body[0],dict) else []
            receipt["sample_values"]=[{k:(str(x)[:110] if x is not None else None) for k,x in d.items()}
                 for d in body[:2]] if isinstance(body,list) else []
        sources[tag]=receipt
        print("GESUI_WFD2022_PUBLIC_SOURCE_PROBE",json.dumps({"name":tag,**receipt},sort_keys=True),flush=True)
    m,rec=req(ARCGIS,{"f":"pjson"})
    layers=m.get("layers",[]) if isinstance(m,dict) else []
    rec["layers"]=[{"id":x.get("id"),"name":x.get("name"),"subLayerIds":x.get("subLayerIds")} for x in layers[:25]]
    sources["WFD_2022_ArcGIS_pressure_layer_list"]=rec
    print("GESUI_WFD2022_PUBLIC_SOURCE_PROBE",json.dumps({"name":"ArcGIS_layers",**rec},sort_keys=True),flush=True)
    if layers:
        for layer in [x for x in layers if x.get("id") in (1,2,4,5,7,8,16,17)]:
            id=layer.get("id")
            sub,info=req(ARCGIS+"/"+str(id),{"f":"pjson"})
            if isinstance(sub,dict):
                info["name"]=sub.get("name")
                info["geometryType"]=sub.get("geometryType")
                info["maxRecordCount"]=sub.get("maxRecordCount")
                info["field_names"]=[t.get("name") for t in (sub.get("fields") or [])[:90]]
                info["capabilities"]=sub.get("capabilities")
            sources["WFD_2022_ArcGIS_layer_"+str(id)]=info
            print("GESUI_WFD2022_PUBLIC_SOURCE_LAYER",json.dumps(info,sort_keys=True),flush=True)
    report={
      "schema":"gesui_2022_EEA_WFD_significant_pressure_source_metadata_probe_v01",
      "retrieved_utc":datetime.now(timezone.utc).isoformat(),
      "official_WFD_third_cycle_significant_pressure_source":True,
      "source_documentation":"https://water.europa.eu/freshwater/resources/metadata/wfd-dashboards/surface-water-bodies-significant-pressures-table",
      "data_endpoints":sources,
      "preoutcome_agricultural_and_urban_pressure_categories_not_yet_validated":True,
      "actual_WWTP_BQE_EQR_or_effects_read":False,
      "no_model_fitted_and_no_biological_scores_used":True,
      "no_site_specific_pressure_values_redistributed":True,
      "cautions":[
        "WFD significant pressure is a regulatory classification, not continuous fertilizer mass",
        "Pressure classes are method/country/reporting-cycle dependent, may be outcome-assessment derived",
        "A biology EQR waterBodyIdentifier requires a source quality / exact key crosswalk; no geometry-only nearest catchment proxy",
        "Pressure variables correlated with biological EQR may create adjustment/collider bias if assigned downstream of ecological condition",
        "The 2022 WFD third-cycle waterbody assessment may not equal EEA WWTP operational/biological monitoring dates"
      ]}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
if __name__=="__main__":
    main()
