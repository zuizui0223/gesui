#!/usr/bin/env python3
"""Find OFFICIAL JRC GREEN 1990-2021 basin-load hosted FeatureServer alternative.

The identical JRC geodatabase archive contains CDF-compressed FileGDB tables:
pyogrio OpenFileGDB cannot read it. This probes *only* the same JRC portal's
official item relationships, search index and published services for a
source-equivalent feature layer or table. Does NOT invent feature counts,
download chemistry data, read EEA 2016 biology/EQR or attribute modelled
sea-outlet nutrients to individual WWTPs.
"""
from __future__ import annotations
import argparse,hashlib,json,time,urllib.request,urllib.parse,urllib.error
from datetime import datetime,timezone
from pathlib import Path

ROOT="https://arcgis-maps.jrc.ec.europa.eu/portal/sharing/rest"
ITEMID="8b105c72109243798913853ad9d05072"
REQMAX=1_500_000
UA="gesui-JRC-GREEN-2016-hosted-layer-discovery/1.0"
def request(path,params):
    q=urllib.parse.urlencode({"f":"json",**params})
    url=ROOT+"/"+path.lstrip("/")+"?"+q
    request=urllib.request.Request(url,headers={
        "Accept":"application/json","User-Agent":UA})
    for i in range(3):
        try:
            with urllib.request.urlopen(request,timeout=55) as r:
                raw=r.read(REQMAX)
                j=json.loads(raw)
                if not isinstance(j,dict):raise ValueError("Official JRC Portal non-object")
                return {"request_url":url,"http":r.status,
                    "response_bytes_truncated_at_1point5MB":len(raw)>=REQMAX,
                    "source_response_sha256":hashlib.sha256(raw).hexdigest(),
                    "response":j}
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,ValueError) as e:
            if isinstance(e,urllib.error.HTTPError) and e.code in (400,401,403,404):
                return {"request_url":url,"http":e.code,
                        "source_error":str(e)[:170]}
            if i==2:return {"request_url":url,"source_error":type(e).__name__+": "+str(e)[:200]}
            time.sleep(i+1.0)
    return {"request_url":url,"source_error":"Unknown"}
def shrink(item):
    d=item.get("response",{})
    return {
        "http":item.get("http"),"url":item.get("request_url"),
        "SHA256":item.get("source_response_sha256"),"error":item.get("source_error"),
        "portal_error":d.get("error"),
        "total":d.get("total"),"start":d.get("start"),"num":d.get("num"),
        "source_item_id":d.get("id"),"source_item_title":d.get("title"),
        "source_item_type":d.get("type"),"source_url":d.get("url"),
        "item_access":d.get("access"),"item_owner":d.get("owner"),
        "matched_official_portal_items":[
           {"id":x.get("id"),"title":x.get("title"),"type":x.get("type"),
            "url":x.get("url"),"size":x.get("size"),"owner":x.get("owner"),
            "access":x.get("access"),
            "modified":x.get("modified")}
           for x in d.get("results",[])[:100]
        ],
        "related_official_portal_items":[
            {"id":x.get("id"),"title":x.get("title"),"type":x.get("type"),
             "url":x.get("url"),"access":x.get("access")}
            for x in d.get("relatedItems",[])[:100]
        ],
        "data_names":sorted([k for k in d if k.lower() in
                       ("features","layers","tables","name","types","copyrighttext")]),
    }
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    queries={
      "original_archive_item":("content/items/"+ITEMID,{}),
      "related_service2data_reverse":("content/items/"+ITEMID+"/relatedItems",{
          "relationshipType":"Service2Data","direction":"reverse"}),
      "related_service2data_forward":("content/items/"+ITEMID+"/relatedItems",{
          "relationshipType":"Service2Data","direction":"forward"}),
      "related_service2service_reverse":("content/items/"+ITEMID+"/relatedItems",{
          "relationshipType":"Service2Service","direction":"reverse"}),
      "search_OLD_filename":("search",{"q":"greenBasinLoad_vrsApr2024","num":100}),
      "search_GREEN_basin":("search",{"q":"greenBasinLoad","num":100}),
      "search_green_model":("search",{"q":"GREEN nutrient","num":100}),
      "search_owner_JRC_water":("search",{"q":"owner:water_job_f","num":100})
    }
    data={}
    candidates={}
    for key,(route,params) in queries.items():
        res=request(route,params)
        clean=shrink(res)
        data[key]=clean
        for entry in clean["matched_official_portal_items"]+clean["related_official_portal_items"]:
            if entry.get("id"):candidates[entry["id"]]=entry
        print("JRC_GREEN_HOSTED_SOURCE_PROBE",key,json.dumps({
          k:clean[k] for k in ("http","portal_error","total","source_item_type","source_item_title","source_url","error")
        },ensure_ascii=False),flush=True)
    # Metadata only for candidate service IDs; never assume ArcGIS GDB=FeatureService.
    service=[x for x in candidates.values() if
       any(q in str(x.get("type","")).lower() for q in ("feature service","map service","scene service","geodata service"))]
    for item in service[:40]:
        url=item.get("url")
        if url and url.startswith("https://"):
            # Probe public MapServer / FeatureServer service metadata, not source features.
            req=urllib.request.Request(url+"?f=pjson",headers={"Accept":"application/json","User-Agent":UA})
            try:
                with urllib.request.urlopen(req,timeout=45) as r:
                    raw=r.read(200_000)
                    m=json.loads(raw)
                    item["service_source_metadata"]={
                        "http":r.status,"SHA256":hashlib.sha256(raw).hexdigest(),
                        "service_description":str(m.get("serviceDescription",""))[:400],
                        "maxRecordCount":m.get("maxRecordCount"),
                        "layers":[{"id":v.get("id"),"name":v.get("name"),"type":v.get("type")}
                                  for v in m.get("layers",[])[:40]],
                        "tables":[{"id":v.get("id"),"name":v.get("name"),"type":v.get("type")}
                                  for v in m.get("tables",[])[:40]],
                        "service_error":m.get("error")
                    }
            except Exception as e:
                item["service_source_metadata"]={"error":type(e).__name__+": "+str(e)[:160]}
    result={
       "schema":"gesui_original_JRC_GREEN_CDF_GDB_source_equivalent_official_REST_service_discovery_v01",
       "UTC":datetime.now(timezone.utc).isoformat(),
       "original_JRC_GDB_archive_item":ITEMID,
       "original_JRC_reported_GDB_GEO_archive_bytes":47068808,
       "original_open_GDAL_OpenFileGDB_failed_CDF_compressed_data":True,
       "original_GDB_CDF_requires_Esri_FileGDB_SDK_not_GDAL_OpenFileGDB":True,
       "only_official_JRC_Portal_search_and_item_relationships_queried":True,
       "official_JRC_rest_source_probe":data,
       "all_official_service_candidate_metadata":list(candidates.values())[:140],
       "hosted_service_candidate_count":len(service),
       "no_EEA_2016_BQE_EQR_values_or_measured_UWWTD_POTW_sources_read":True,
       "do_not_claim_hosted_equivalence_without_same_source_year_and_GDB_field_validation":True,
       "JRC_GREEN1990_2021_is_sea_basin_MODELED_NOT_2016_measured_WWTP_effluent":True
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("JRC_GREEN_HOSTED_OFFICIAL_CANDIDATE_SUMMARY",json.dumps({
        "n_official_candidate_ids":len(candidates),"n_service_items":len(service),
        "service_candidates":[{"id":x["id"],"title":x["title"],"url":x["url"]}
                              for x in service[:20]]},ensure_ascii=False),flush=True)
if __name__=="__main__":main()
