#!/usr/bin/env python3
"""2026-11-15 original GBIF registry DATASET-METADATA ONLY route to 2022 butterfly events.

No API call to occurrence endpoints; no event/occurrence rows, coordinates,
species, HNS, effluent conditions or biological responses opened. No access
entitlement is assumed from public registry metadata. Frozen PR95 primary
eligibility criteria, year, endpoints and STOP remain unchanged.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

FROZEN_PROTOCOL_GIT_BLOB="b9ce415c36d398eee77d8318a20aef4873e0bf44"
API="https://api.gbif.org/v1/dataset/"
# Pre-specified original dataset identities; visit GBIF registry only.
CANDIDATES=[
 ("SE","Swedish Butterfly Monitoring Scheme","be77e203-486c-4651-91b9-8347968b728c","previously audited 2022 264 point-site event source"),
 ("BE","Flanders all-butterfly monitoring","81b909c7-5612-4148-8e70-99fea2cdb572","published 5km generalized positions; catalog temporal scope 1991–2020"),
 ("BE","Flanders priority species monitoring","63938753-1fec-4c08-ae39-e9f8a6576521","reported 1/5/10km generalized UTM positions; precise section positions require approval"),
 ("DK","Denmark EBMS transects","caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c","sampling started 2023"),
 ("DK","Denmark EBMS fifteen-minute counts","99f9976a-dfdd-4359-bf4a-5538f29458cf","sampling started 2023"),
 ("ES","Sierra Nevada butterfly monitoring","d7a3df60-e856-438b-851d-74653ef6c0d2","2021 archived version not necessarily including 2022"),
 ("FI","Finnish agricultural landscape butterfly monitoring","eacd5f7d-60aa-49e9-a7cd-36d1c46b22ac","GBIF occurrence dataset not sampling-event core; dated transect counts may not preserve visits/zeros"),
 ("NL","NDFF butterfly-dragonfly monitoring programme","fda5ee54-cd11-4da7-9314-0d328f1ccfda","metadata-only registry; original sections/events require NDFF request"),
 ("NL","Dutch Butterfly Monitoring Indexes","1897eef5-1e10-4842-97f7-cad98925b59d","published indexes, not 2022 visits/section coordinates"),
 ("NL","NDFF Butterfly Data","421e6128-74cb-48ac-b76e-e419d39e2527","2022 publication date; occurrence data not verified as 2022 survey events"),
 ("GB","UK Butterfly Monitoring Scheme","1e266c3d-92ef-4d5a-8e4a-c04742c772c3","2026 publisher archive but precise section event access/year requires direct audit"),
]
EVENT_TOKENS={"SAMPLING_EVENT","SAMPLINGEVENT","EVENT"}
def sha_blob(data:bytes)->str:
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()

def query_json(url):
    if not url.startswith(API) or ("occurrence" in url.lower()):
        raise ValueError("Only GBIF dataset registry API allowed, NOT species or occurrence")
    last=None
    for retry in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"GESUI-sourceonly-20261115/1.0 no butterfly records",
                "Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=45) as r:
                raw=r.read(1024*1024+1)
                if len(raw)>1024*1024:raise RuntimeError("Unexpected registry response >1 MB")
            return json.loads(raw),hashlib.sha256(raw).hexdigest()
        except Exception as e:
            last=e
            if retry<2:time.sleep(retry+1)
    return None,str(last)

def extract(item):
    # Be conservative: metadata summaries are not exact annual event availability.
    contacts=[]
    for p in item.get("contacts") or []:
        mail=p.get("email")
        if isinstance(mail,list):mail=[] # never publish personal/irrelevant addresses
        contacts.append({"type":p.get("type"),"role":p.get("roles")})
    return {
      "key":item.get("key"),"title":item.get("title"),"type":item.get("type"),
      "subtype":item.get("subtype"),"license":item.get("license"),
      "publishingOrganizationKey":item.get("publishingOrganizationKey"),
      "published":item.get("pubDate"),"modified":item.get("modified"),
      "created":item.get("created"),"homepage":item.get("homepage"),
      "doi":item.get("doi"),
      "temporalCoverages":item.get("temporalCoverages"),
      "geographicCoverages":item.get("geographicCoverages"),
      "endpoints":[{"type":e.get("type"),"url":e.get("url")} for e in (item.get("endpoints") or [])],
      "n_metadata_contacts":len(contacts),
      "registry_only_no_observation_count":True,
      "archive_field_2022_precision_NOT_VERIFIED":True,
    }

def run(args):
    original=args.freeze.read_bytes()
    if sha_blob(original)!=FROZEN_PROTOCOL_GIT_BLOB:
        raise RuntimeError("Frozen primary v1.0 protocol modified; STOP")
    results=[]
    for country,label,key,known_constraint in CANDIDATES:
        data,digest=query_json(API+key)
        if not data or data.get("key")!=key:
            results.append({"country":country,"label":label,"dataset_key":key,
                "status":"REGISTRY_UNAVAILABLE","error":str(digest)[:180],
                "known_from_original_source":known_constraint})
            continue
        record=extract(data)
        results.append({"country":country,"label":label,"dataset_key":key,
          "status":"REGISTRY_VERIFIED_METADATA_ONLY",
          "registry_receipt_sha256":digest,
          "known_from_original_source":known_constraint,
          "registry":record,
          "source_2022_event_0to100m_verified":False,
          "qualified_8_cell_source_reference_pairs_verified":0})
        print("SOURCE_REGISTRY_ROW="+json.dumps({"country":country,"label":label,
          "type":record["type"],"published":record["published"],
          "key":key},ensure_ascii=False),flush=True)
    # Strictly *NO* GBIF occurrence/endpoints retrieved; no source matches computed.
    out={
      "schema":"GESUI_20261115_2022_3COUNTRY_BUTTERFLY_DATASET_REGISTRY_METADATA_ONLY_V1",
      "frozen_protocol_git_blob":FROZEN_PROTOCOL_GIT_BLOB,
      "n_registry_candidates_preselected":len(CANDIDATES),
      "n_registry_metadata_available":sum(x["status"]=="REGISTRY_VERIFIED_METADATA_ONLY" for x in results),
      "n_2022_full_8cell_matched_countries_verified":0,
      "GBIF_event_occurrence_or_species_results_queried":False,
      "outfall_site_matching_or_HNS_performed":False,
      "2022_2022_crosscountry_data_granted":False,
      "frozen_biological_primary_status":"STOP_SOURCE_OPPORTUNITY_UNPROVEN",
      "additional_country_precision_must_be_authorized_and_verified":True,
      "results":results,
      "official_eBMS_2022_2024_data_access":"https://butterfly-monitoring.net/ebms-data%20access",
      "eBMS_data_request_requires_signed_licence_and_coordinator_approval":True,
      "no_request_was_sent":True,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_DATASET_REGISTRY_SOURCEONLY="+json.dumps({
       "metadata_available":out["n_registry_metadata_available"],
       "candidates":len(CANDIDATES),
       "outcome_test_ran":False,
       "main_status":out["frozen_biological_primary_status"]}),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--freeze",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
