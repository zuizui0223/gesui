#!/usr/bin/env python3
"""Source-only iNaturalist 2018/2022 European butterfly puddling FEASIBILITY.

Requests only public, aggregate query counts and at most a few representative
observation IDs for inspecting API/schema correctness. No photo is classified
as feeding, no sewage exposure is attached, no coordinates/users are archived.

Status: candidate-annotation audit, NOT photo-validation or WWTP effect.
"""
from __future__ import annotations
import argparse, datetime, hashlib, json, time, urllib.error
import urllib.parse, urllib.request
from pathlib import Path

ROOT="https://api.inaturalist.org/v1"
LEGACY="https://www.inaturalist.org"
FIELD_ID=17428
BBOX={"swlat":"44","swlng":"-12","nelat":"64","nelng":"26"}
TAXA={"Lycaenidae":"Lycaenidae","Pieridae":"Pieridae","Papilionidae":"Papilionidae"}
YEARS=(2018,2022)
UA="gesui-source-only-puddling-feasibility/0.1 (public aggregate records; academic research)"
MAX_CALLS=40
COUNT=0

def get_json(base,path,params=None):
    global COUNT
    COUNT+=1
    if COUNT>MAX_CALLS:raise RuntimeError("Safety ceiling on public API queries")
    q=urllib.parse.urlencode(params or {})
    url=base+path+("?" + q if q else "")
    out={"url":url,"api_domain":urllib.parse.urlparse(url).hostname}
    for attempt in range(2):
        try:
            request=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":"application/json"})
            with urllib.request.urlopen(request,timeout=28) as r:
                raw=r.read(400000)
                try:j=json.loads(raw.decode("utf-8"))
                except Exception:
                    out.update({"status":r.status,"body_non_json":True,"sha256":hashlib.sha256(raw).hexdigest()})
                    return out,None
                out.update({"status":r.status,"body_sha256":hashlib.sha256(raw).hexdigest()})
                return out,j
        except urllib.error.HTTPError as e:
            out.update({"status":e.code,"error":str(e),"retry_after":e.headers.get("Retry-After")})
            if e.code not in (429,500,502,503) or attempt:break
        except Exception as e:
            out.update({"error":type(e).__name__+": "+str(e)})
            if attempt:break
        time.sleep(1.8)
    return out,None

def rootparams(taxon,year):
    return {**BBOX,"taxon_name":taxon,"d1":f"{year}-01-01","d2":f"{year}-12-31",
            "photos":"true","geo":"true","geoprivacy":"open","per_page":"2"}

def filter_counts(taxon,year,mods):
    p=rootparams(taxon,year);p.update(mods)
    receipt,j=get_json(ROOT,"/observations",p)
    info={"request":receipt,
      "returned_total":j.get("total_results") if isinstance(j,dict) else None,
      "returned_n":len(j.get("results",[])) if isinstance(j,dict) else None,
      "sample_for_schema":[],
      "candidate_not_photo_confirmed":True}
    if isinstance(j,dict):
        for obs in j.get("results",[])[:2]:
            fields=obs.get("ofvs") or obs.get("observation_field_values") or []
            relevant=[{"field_id":f.get("observation_field_id"),
                      "field_name":(f.get("name") or (f.get("observation_field") or {}).get("name")),
                      "value":f.get("value")}
                      for f in fields]
            info["sample_for_schema"].append({
                "id":obs.get("id"),"taxon":(obs.get("taxon") or {}).get("name"),
                "has_photos":bool(obs.get("photos")),
                "coordinates_obscured":obs.get("coordinates_obscured"),
                "geoprivacy":obs.get("geoprivacy"),
                "position_accuracy_m":obs.get("positional_accuracy"),
                "observation_fields":relevant[:12]})
    return info

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    metadata={}
    for path in (f"/observation_fields/{FIELD_ID}.json",):
        req,x=get_json(LEGACY,path)
        metadata["legacy_field_details"]={"request":req}
        if isinstance(x,dict):
            metadata["legacy_field_details"].update({
                "id":x.get("id"),"name":x.get("name"),"description":x.get("description"),
                "allowed_values":x.get("allowed_values"),"datatype":x.get("datatype"),
                "raw_fields":list(x)[:30]})
    req,x=get_json(ROOT,f"/observation_fields/{FIELD_ID}")
    metadata["v1_field_details"]={"request":req}
    if isinstance(x,dict):
        metadata["v1_field_details"].update({"total_results":x.get("total_results"),
           "name":x.get("name"),"results":x.get("results",[])[:2]})
    print("GESUI_PUDDLING_FIELD_METADATA",
          json.dumps(metadata,ensure_ascii=False,default=str)[:1800],flush=True)

    tests={}
    # Field-name search is provisional until exact public metadata confirms it.
    # Query only bounded historical years; result count is a candidate pool,
    # not number of real puddling photos or observations near wastewater.
    for family in TAXA:
        tests[family]={}
        for year in YEARS:
            results={
              "all_photo_georeferenced":filter_counts(family,year,{}),
              "text_puddling_candidate":filter_counts(family,year,{"q":"puddling"}),
              "field_id_candidate":filter_counts(family,year,{f"field:{FIELD_ID}":"Yes"}),
            }
            tests[family][str(year)]=results
            print("GESUI_PUBLIC_INAT_PUDDLING_SOURCE_ONLY",json.dumps({
                "taxon":family,"year":year,
                "counts":{k:v["returned_total"] for k,v in results.items()},
                "statuses":{k:v["request"].get("status") for k,v in results.items()},
            },sort_keys=True),flush=True)
            time.sleep(0.7)

    status={
       "schema":"gesui_iNaturalist_Europe_butterfly_2018_2022_puddling_source_candidate_gate_v01",
       "run_date_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat(),
       "source":"iNaturalist official public v1 observations",
       "field_id_target":FIELD_ID,"metadata":metadata,
       "query_date_is_observation_date_not_upload_date":True,
       "years":YEARS,"region_bbox":BBOX,"searches":tests,
       "no_photo_images_downloaded":True,
       "photo_confirmed_puddling_n":None,
       "no_real_sewage_association_or_ecological_prediction_fitted":True,
       "no_coarse_observation_counts_interpreted_as_absence":True,
       "critically_uncertain":[
          "field 17428 has mixed aggregation/roosting/puddling categories",
          "field-by-ID API filtering may be unsupported and therefore zero must be schema-verified",
          "text search matches are only candidate labels, not confirmed behavior",
          "presence of a photograph alone does not show substrate contact",
          "geoprivacy, taxon geoprivacy, coordinate accuracy and observer effort are nonrandom",
          "modern iNaturalist public API may reject queries, producing NA rather than zeros",
       ]
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(status,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_PUDDLING_SOURCE_CANDIDATE_QC_COMPLETED",
      json.dumps({"api_calls":COUNT,"no_puddle_confirmed":True,"no_sewage_result":True}),flush=True)
if __name__=="__main__":
    main()
