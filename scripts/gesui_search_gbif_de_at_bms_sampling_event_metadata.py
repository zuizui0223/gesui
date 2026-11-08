#!/usr/bin/env python3
"""Read official GBIF DATASET registry metadata only for German/Austrian BMS sources.

No GBIF occurrence, observation, species-name/individual count, event source
data, or ecological outcome touched. This is dataset *discovery*, not an eBMS
data-access bypass. Detect publisher original Event/DwC-A endpoints only.
"""
from __future__ import annotations
import argparse,hashlib,json,time,urllib.error,urllib.parse,urllib.request
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path

BASE="https://api.gbif.org/v1"
SEARCH_TERMS=(
  "Tagfalter-Monitoring Deutschland","Tagfaltermonitoring","Viel-Falter",
  "Butterfly Monitoring Austria","Schmetterlingsmonitoring",
  "Lepidoptera monitoring Germany","Lepidoptera monitoring Austria",
  "BMS Germany","BMS Austria","Austrian butterflies",
  "German butterflies monitoring","butterfly transect"
)
MAX_PAGE=100
MAX_CANDIDATE_DETAIL=18
RETRIES=3

def api(path,params=None):
    url=BASE+path+("?" + urllib.parse.urlencode(params) if params else "")
    receipt={"url":url}
    for attempt in range(RETRIES):
        try:
            request=urllib.request.Request(url,headers={
             "User-Agent":"gesui-public-2022-DE-AT-butterfly-event-metadata-discovery/1.0",
             "Accept":"application/json"})
            with urllib.request.urlopen(request,timeout=42) as response:
                raw=response.read(4_000_000)
                j=json.loads(raw.decode("utf-8"))
                receipt.update({"http":response.status,
                      "retrieved_body_sha256":hashlib.sha256(raw).hexdigest(),
                      "bytes":len(raw)})
                return receipt,j
        except urllib.error.HTTPError as e:
            receipt.update({"http":e.code,"error":str(e)})
            if e.code not in (429,500,502,503) or attempt==RETRIES-1:break
        except Exception as e:
            receipt["error"]=type(e).__name__+": "+str(e)
            if attempt==RETRIES-1:break
        time.sleep(2*(attempt+1))
    return receipt,None

def is_eligible_metadata(record):
    title=str(record.get("title") or record.get("name") or "")
    description=str(record.get("description") or "")
    p=str(record.get("publishingCountry") or "")
    # Rank only non-sensitive, *dataset metadata*. Catch DE/AT country even
    # when English title lacks local country strings.
    terms=(title+" "+description[:700]).casefold()
    country_hit=any(x in terms for x in
          ("germany","deutschland","tagfalter","schmetterling",
           "austria","österreich","viel-falter","tirol","vorarlberg"))
    bms_hit=any(x in terms for x in
          ("butterfly","butterflies","lepidoptera","tagfalter",
           "schmetterling","bms","viel-falter","transect"))
    return (country_hit and bms_hit) or (p in ("DE","AT") and bms_hit)

def slim(record):
    return {k:record.get(k) for k in (
      "key","title","type","license","publishingCountry",
      "doi","created","modified","pubDate","homepage","hostingOrganizationKey",
      "publishingOrganizationKey","publishingOrganizationTitle",
      "sourceId","external","numOccurrences",
      "installationKey") if k in record}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    receipts=[];found={}
    for q in SEARCH_TERMS:
        rec,j=api("/dataset/search",{"q":q,"limit":MAX_PAGE})
        if not isinstance(j,dict) or not isinstance(j.get("results"),list):
            receipts.append({"query":q,"request":rec,"error":"No valid result list"})
            continue
        arr=j["results"]
        receipts.append({
          "query":q,"request":rec,
          "original_search_count_total":j.get("count"),
          "records_returned_first100":len(arr),
          "not_exhaustive_if_more_than_limit":int(j.get("count") or 0)>MAX_PAGE,
          "matching_title_metadata":[slim(x) for x in arr[:10]],
        })
        for x in arr:
            if not x.get("key"):continue
            key=x["key"]
            if key not in found:found[key]={"metadata":slim(x),"matched_queries":[]}
            found[key]["matched_queries"].append(q)
        print("GESUI_GBIF_DE_AT_BMS_PUBLIC_DATASET_QUERY",
          json.dumps({"query":q,"total":j.get("count"),
            "returned":len(arr),"status":rec.get("http")},ensure_ascii=False),flush=True)
        time.sleep(0.35)
    candidate=[
        {**v,"dataset_key":k} for k,v in found.items()
        if is_eligible_metadata(v["metadata"]) or
           (v["metadata"].get("type")=="SAMPLING_EVENT" and
            v["metadata"].get("publishingCountry") in ("DE","AT"))
    ]
    # DWC source type and likely geography first, query frequency second.
    candidate.sort(key=lambda v:(
      0 if v["metadata"].get("type")=="SAMPLING_EVENT" else 1,
      -len(v["matched_queries"]),str(v["metadata"].get("title"))))
    details=[]
    for c in candidate[:MAX_CANDIDATE_DETAIL]:
        key=c["dataset_key"]
        rec,j=api("/dataset/"+key)
        info={"dataset_key":key,"metadata_request":rec,"catalogue_preview":c}
        if isinstance(j,dict):
            info["official_gbif_metadata"]=slim(j)
            info["dataset_type"]=j.get("type")
            info["citation"]=str(j.get("citation") or "")[:600]
            info["original_dataset_endpoints"]=[
                {"type":e.get("type"),"url":e.get("url")}
                for e in j.get("endpoints",[]) if isinstance(e,dict)
            ][:20]
            info["homepage"]=j.get("homepage")
            info["description_preview"]=str(j.get("description") or "")[:650]
            info["publisher_contacts_present"]=bool(j.get("contacts"))
            info["event_core_confirmed_from_original_archive"]=False
            info["2022_event_records_examined"]=False
            info["eligibility_status"]=(
                "CANDIDATE_SAMPLING_EVENT_METADATA_ONLY" if j.get("type")=="SAMPLING_EVENT"
                else "NOT_EVENT_CORE_METADATA" )
        details.append(info)
        print("GESUI_GBIF_DE_AT_BMS_METADATA_CANDIDATE",
              json.dumps({"key":key,"title":c["metadata"].get("title"),
                    "type":info.get("dataset_type"),"endpoint_count":len(info.get("original_dataset_endpoints",[])),
                    "2022_event_opened":False},ensure_ascii=False),flush=True)
        time.sleep(0.3)
    selected_event=[x for x in details if x.get("dataset_type")=="SAMPLING_EVENT"]
    out={
      "schema":"gesui_GBIF_public_DE_AT_butterfly_monitoring_dataset_metadata_source_search_20261108_v01",
      "retrieved_UTC":datetime.now(timezone.utc).isoformat(),
      "GBIF_official_registry":"https://api.gbif.org/v1/dataset/search",
      "query_terms_fixed":list(SEARCH_TERMS),
      "search_limit_per_query":MAX_PAGE,
      "query_receipts":receipts,
      "number_distinct_query_response_datasets":len(found),
      "number_possible_de_at_butterfly_metadata":len(candidate),
      "candidate_inventory":[{**x,"detail_inspected":i<MAX_CANDIDATE_DETAIL}
                 for i,x in enumerate(candidate)],
      "candidate_detail_inspections":details,
      "number_detail_inspected_sampling_event_candidates":len(selected_event),
      "broad_DE_AT_event_candidates_can_be_nonbutterfly":True,
      "source_discovery_search_is_NOT_exhaustive":True,
      "verified_open_2022_DwC_Event_core_n":0,
      "explicit_UNVERIFIED_original_Event_core":True,
      "no_download_of_butterfly_counts_or_species_records":True,
      "no_equating_occurrence_dataset_with_visit_sampling":True,
      "does_not_imply_EEA_outfall_matching":True,
      "public_registry_query_does_not_bypass_coordinator_access":True,
      "next_source_only_action":"For credible SAMPLING_EVENT candidates, verify public publisher ZIP and original Event core; if none, request coarse near/far site counts through eBMS/UFZ/Viel-Falter coordinators.",
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_DE_AT_OPEN_BMS_DATASET_REGISTRY_DISCOVERY_COMPLETE",
      json.dumps({"searched_queries":len(SEARCH_TERMS),
          "all_metadata":len(found),"candidate_de_at_butterfly_metadata":len(candidate),
          "detail_event_candidates":len(selected_event),"2022_event_core_verified":0,
          "NO_BUTTERFLY_OUTCOME_OPENED":True},sort_keys=True),flush=True)
if __name__=="__main__":
    main()
