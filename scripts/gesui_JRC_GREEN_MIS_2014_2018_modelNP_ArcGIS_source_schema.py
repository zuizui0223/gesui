#!/usr/bin/env python3
"""JRC OFFICIAL MAP service (not the inaccessible CDF FileGDB): layer SOURCE QC.

The JRC Water Hub inmap/Losses_to_water_MIS ArcGIS MapServer states that
GREEN N/P point sources were aggregated to country, NUTS2 and river basin,
averaged 2014–2018. Query published schema/attributes availability only.
NEVER equate them to 2016 observed WWTP discharge, wetland nutrient dosing,
site-level ecological effects, or 2018 EEA biological EQR scores.
"""
from __future__ import annotations
import argparse,hashlib,json,time,urllib.request,urllib.parse,urllib.error
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

BASE="https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/Losses_to_water_MIS/MapServer"
HEAD={"User-Agent":"gesui-JRC-GREEN-2014-2018-pointsource-arcgis-metadata/1.0","Accept":"application/json"}
MAX_LAYER=100
def fetch(path,opts):
    url=BASE+path+"?"+urllib.parse.urlencode({"f":"json",**opts})
    err=None
    for i in range(3):
        try:
            req=urllib.request.Request(url,headers=HEAD)
            with urllib.request.urlopen(req,timeout=55) as r:
                raw=r.read(850_000)
                j=json.loads(raw)
                if not isinstance(j,dict):raise ValueError("Official JRC MapServer non-JSON")
                return {"url":url,"status":r.status,"SHA256":hashlib.sha256(raw).hexdigest(),
                        "JSON":j}
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,ValueError) as e:
            err=(type(e).__name__+": "+str(e))[:180]
            if isinstance(e,urllib.error.HTTPError) and e.code in (400,401,403,404):
                break
            if i<2:time.sleep(i+1)
    return {"url":url,"error":err}
def tidy_layer(l):
    return {k:l.get(k) for k in ("id","name","type","parentLayerId","subLayerIds","defaultVisibility","minScale","maxScale")}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    service=fetch("",{})
    meta=service.get("JSON",{})
    ll=meta.get("layers",[])
    if len(ll)>MAX_LAYER:raise ValueError("JRC original official MapServer unexpected layer count")
    logs=[]
    for item in ll:
        layerid=item.get("id")
        if not isinstance(layerid,int):continue
        full=fetch("/"+str(layerid),{})
        d=full.get("JSON",{})
        fields=[{"name":v.get("name"),"alias":v.get("alias"),"type":v.get("type"),
                 "length":v.get("length")} for v in d.get("fields",[])[:300]]
        pointfield=[v for v in fields if any(k in str(v["name"]).lower()+" "+str(v["alias"]).lower()
            for k in ("point","green","nutr","emission","waste","sewage","phosph","nitrog",
                      "1990","2016","2014","2018","load","agri","industry","rbd","nuts","country"))]
        logs.append({
          "layer":tidy_layer(item),
          "http":full.get("status"),"error":full.get("error") or d.get("error"),
          "official_source_metadata_SHA256":full.get("SHA256"),
          "original_service_layer_type":d.get("type"),
          "original_geometryType":d.get("geometryType"),
          "service_fields":fields,
          "source_point_and_nutrient_candidate_columns":pointfield,
          "original_source_maxRecordCount":d.get("maxRecordCount"),
          "original_supportsStatistics":d.get("advancedQueryCapabilities",{}).get("supportsStatistics"),
          "reported_copyright":str(d.get("copyrightText",""))[:700],
          "layer_description_excerpt":str(d.get("description",""))[:700],
          "layer_name":d.get("name"),
          "field_count":len(fields),
          "supportsGeometry":bool(d.get("geometryType")),
        })
        if len(logs)%15==0:
            print("JRC_GREEN_MIS_LAYER_SCHEMA_READ",len(logs),"/",len(ll),flush=True)
    text=" "+str(meta.get("serviceDescription",""))+" "+str(meta.get("description",""))
    actual_names=[x["layer_name"] or x["layer"]["name"] for x in logs]
    pointLayers=[x for x in logs if any(
        q in str(x["layer_name"] or x["layer"]["name"]).lower()
        for q in ("point source","point_source","pointsource","p_sou","n_sou","point src"))]
    results={
       "schema":"gesui_OFFICIAL_JRC_GREEN_MIS_2014to2018_model_NP_pointsource_MapServer_layer_metadata_v01",
       "utc":datetime.now(timezone.utc).isoformat(),
       "original_JRC_official_WaterHub_MapServer":BASE,
       "official_original_service_http":service.get("status"),
       "original_service_SHA256":service.get("SHA256"),
       "original_service_error":service.get("error") or meta.get("error"),
       "original_JRC_serviceDescription":str(meta.get("serviceDescription",""))[:5500],
       "source_description":str(meta.get("description",""))[:5500],
       "source_reported_pointsource_NP_country_NUTS2_basin_averaged_2014to2018":True,
       "source_modelled_not_measured_NP_from_sewage_plus_industry":True,
       "original_layer_count":len(ll),
       "original_layer_names":actual_names,
       "original_candidate_point_source_layers":len(pointLayers),
       "layer_schema_and_query_support":logs,
       "source_numeric_2016_WWTP_specific_tonnes_available":False,
       "2016_EEA_benthic_EQR_outcome_scores_read":False,
       "source_2018_MIS_basins_NOT_direct_NHDplus_EEA_microsite_exposures":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(results,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("JRC_MIS_ACTUAL_OFFICIAL_SOURCE",json.dumps({
       "http":results["official_original_service_http"],"count":len(ll),
       "point_NP_layers":len(pointLayers),
       "candidate_layer_names":[x["layer_name"] for x in pointLayers],
       "all_layer_names":actual_names[:90],
       "errors":sum(x["error"] is not None for x in logs),
    },ensure_ascii=False),flush=True)
if __name__=="__main__":main()
