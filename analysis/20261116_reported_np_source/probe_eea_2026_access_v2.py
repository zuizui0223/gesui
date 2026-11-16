#!/usr/bin/env python3
"""Official EEA Waterbase 2026 metadata vs older Discodata availability.

Investigate exact publisher Part 2 source entry and minimal SQL predicate
failure. All HTTP requests capped at 200 KiB body, SQL max one row; no
individual N/P observation value is retained. No model or ecological output.
The 'latest' SQL database MUST NOT be called 2026 publisher source without
verified source-version correspondence.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import urlparse
import xml.etree.ElementTree as ET
import requests

DATASET_UUID="82e194cc-2506-411c-b341-9eb50b6ad5e0"
EEA_HOST="https://sdi.eea.europa.eu"
EEA_SHARE=EEA_HOST+"/data/"+DATASET_UUID
CATALOG=EEA_HOST+"/catalogue/datahub/api/records/"+DATASET_UUID
DISCO="https://discodata.eea.europa.eu/sql"
SQL_TABLE="[WISE_SOE].[latest].[Waterbase_T_WISE6_AggregatedData]"
MAX_BYTES=200_000
TIMEOUT=(8,18)
UA="gesui-EEA-2026-ICM-public-metadata-predicate-QA/20261116"

def bounded_request(url,method="GET",params=None,headers=None,body=None,auth=None):
    kw=dict(params=params,headers={"User-Agent":UA,"Accept":"application/json,application/xml,text/html,*/*",**(headers or {})},
       timeout=TIMEOUT,stream=True,allow_redirects=False)
    if auth is not None:kw["auth"]=auth
    try:
        resp=requests.request(method,url,data=body,**kw)
        payload=bytearray()
        for part in resp.iter_content(chunk_size=16384):
            payload.extend(part)
            if len(payload)>MAX_BYTES:
                raise RuntimeError("Bounded metadata response exceeded 200 kB, not ingested")
        out={"http_status":int(resp.status_code),"content_type":resp.headers.get("Content-Type"),
            "content_length_header":resp.headers.get("Content-Length"),
            "location_hostname":urlparse(resp.headers.get("Location","")).hostname,
            "body_bytes":len(payload),
            "body_sha256":hashlib.sha256(payload).hexdigest()}
        return out,bytes(payload)
    except Exception as e:
        return {"http_status":None,"failure":type(e).__name__+": "+str(e)[:140]},b""

def sqltest(tag,query):
    meta,data=bounded_request(DISCO,params={"query":query,"p":1,"nrOfHits":1})
    out={"tag":tag,"SQL_SHA256":hashlib.sha256(query.encode()).hexdigest(),
         **meta,"maximum_rows_requested":1,"full_chemical_values_not_recorded":True}
    if data and "json" in str(meta.get("content_type","")).lower():
        try:
            payload=json.loads(data)
            if isinstance(payload,dict) and isinstance(payload.get("results"),list):
                rows=payload["results"]
                out["status"]="RESULTS_JSON"
                out["n_actual_rows_0_or_1"]=len(rows)
                if rows and isinstance(rows[0],dict):
                    out["returned_column_names"]=list(rows[0])
                    if "phenomenonTimeReferenceYear" in rows[0]:
                        out["returned_year"]=rows[0]["phenomenonTimeReferenceYear"]
                    if "observedPropertyDeterminandCode" in rows[0]:
                        out["returned_determinand_code"]=rows[0]["observedPropertyDeterminandCode"]
                    if "resultUom" in rows[0]:
                        out["reported_unit"]=rows[0]["resultUom"]
            elif isinstance(payload,dict) and isinstance(payload.get("errors"),list):
                out["status"]="SQL_ERROR_JSON"
                out["error_codes"]=[e.get("errorcode") for e in payload["errors"] if isinstance(e,dict)]
                out["error_category"]=[str(e.get("error",""))[:120] for e in payload["errors"] if isinstance(e,dict)]
            else:out["status"]="JSON_OTHER"
        except Exception as e:
            out["status"]="INVALID_JSON"
            out["parsing_exception"]=type(e).__name__
    else:
        out["status"]="BLOCK_OR_HTML_OR_EMPTY"
    return out

def parse_catalog_info():
    meta,raw=bounded_request(CATALOG,headers={"Accept":"application/json"})
    out={"endpoint":CATALOG,**meta,"type":"EEA_2026_official_catalogue"}
    if "json" not in str(meta.get("content_type","")).lower() or not raw:
        out["catalogue_parse"]="NO_JSON"
        return out
    try:
        j=json.loads(raw)
        out["catalogue_parse"]="JSON"
        out["keys"]=list(j)[:25] if isinstance(j,dict) else []
        urls=set()
        def walk(v):
            if isinstance(v,dict):
                for q in v.values():walk(q)
            elif isinstance(v,list):
                for q in v:walk(q)
            elif isinstance(v,str) and v.startswith("https://") and (
                "sdi.eea.europa.eu" in v or "eea.europa.eu/data/" in v):
                urls.add(v.split("?")[0][:400])
        walk(j)
        out["public_eea_source_resource_links"]=sorted(urls)[:30]
        out["original_catalogue_values_not_republished"]=True
    except Exception as ex:
        out["catalogue_parse"]="ERROR"
        out["error"]=type(ex).__name__+": "+str(ex)[:100]
    return out

def inspect_public_share():
    meta,html=bounded_request(EEA_SHARE,headers={"Accept":"text/html"})
    out={"share_url":EEA_SHARE,**meta}
    if not html or "html" not in str(meta.get("content_type","")).lower():
        return out,[]
    s=html.decode("utf-8",errors="replace")
    patterns=[
        r'data-sharing-token=["\x27]([^"\x27]{3,100})',
        r'data-share-token=["\x27]([^"\x27]{3,100})',
        r'"shareToken":\s*"([^"]{3,100})"',
        r'data-token=["\x27]([^"\x27]{3,100})'
    ]
    out["is_nextcloud"]=bool(re.search(r'nextcloud|ownCloud|files_sharing|public\.php',s,re.I))
    out["html_contains_public_filelist_markup"]=bool(re.search(r'fileslist|filelist|data-file',s,re.I))
    # Only public share identifiers; do not log non-share auth/session tokens.
    public_tokens=[]
    for pattern in patterns[:3]:
        public_tokens.extend(re.findall(pattern,s)[:2])
    if not public_tokens:public_tokens=[DATASET_UUID]
    out["share_identifier_candidates_used_n"]=len(set(public_tokens))
    return out,sorted(set(public_tokens))[:2]

def test_public_dav(token):
    # Nextcloud's documented public-share DAV. No private login or cookies.
    paths=[
        ("nextcloud_public_dav_v2",EEA_HOST+"/public.php/dav/files/"+token+"/"),
        ("nextcloud_public_dav_legacy",EEA_HOST+"/public.php/webdav/"),
    ]
    ans=[]
    for name,url in paths:
        # legacy uses publicly shared token only as HTTP Basic username; blank password.
        auth=(token,"") if name=="nextcloud_public_dav_legacy" else None
        meta,raw=bounded_request(url,method="PROPFIND",
            headers={"Depth":"1","Content-Type":"application/xml"},
            body=b'<?xml version="1.0"?><d:propfind xmlns:d="DAV:"><d:prop><d:displayname/><d:getcontentlength/><d:resourcetype/></d:prop></d:propfind>',
            auth=auth)
        item={"test":name,**meta,"only_public_share":True,
              "not_a_zip_or_scientific_measurement_download":True}
        if raw and "xml" in str(meta.get("content_type","")).lower():
            try:
                root=ET.fromstring(raw)
                ns={"d":"DAV:"}
                resp=root.findall(".//d:response",ns)
                entries=[]
                for v in resp:
                    href=v.findtext("d:href",default="",namespaces=ns)
                    n=v.findtext(".//d:displayname",default="",namespaces=ns)
                    s=v.findtext(".//d:getcontentlength",default="",namespaces=ns)
                    if href or n:
                        entries.append({"name":n[:90],"size_bytes":s or None,
                                        "is_directory":v.find(".//d:collection",ns) is not None})
                item["webdav_listing_entries"]=entries[:100]
                item["n_items"]=len(resp)
            except Exception as e:
                item["xml_parse_error"]=str(e)[:100]
        ans.append(item)
    return ans

def main(args):
    prefixes=[
      ("unfiltered_year",f"SELECT TOP 1 phenomenonTimeReferenceYear FROM {SQL_TABLE}"),
      ("filter_year_2018",f"SELECT TOP 1 phenomenonTimeReferenceYear FROM {SQL_TABLE} WHERE phenomenonTimeReferenceYear = 2018"),
      ("filter_year_2022",f"SELECT TOP 1 phenomenonTimeReferenceYear FROM {SQL_TABLE} WHERE phenomenonTimeReferenceYear = 2022"),
      ("filter_N_code",f"SELECT TOP 1 observedPropertyDeterminandCode FROM {SQL_TABLE} WHERE observedPropertyDeterminandCode = 'EEA_31615-01-7'"),
      ("filter_P_code",f"SELECT TOP 1 observedPropertyDeterminandCode FROM {SQL_TABLE} WHERE observedPropertyDeterminandCode = 'CAS_7723-14-0'"),
      ("filter_2022_N_code",f"SELECT TOP 1 phenomenonTimeReferenceYear,observedPropertyDeterminandCode FROM {SQL_TABLE} WHERE phenomenonTimeReferenceYear = 2022 AND observedPropertyDeterminandCode = 'EEA_31615-01-7'"),
      ("filter_2022_P_code",f"SELECT TOP 1 phenomenonTimeReferenceYear,observedPropertyDeterminandCode FROM {SQL_TABLE} WHERE phenomenonTimeReferenceYear = 2022 AND observedPropertyDeterminandCode = 'CAS_7723-14-0'"),
      ("filter_river_RW",f"SELECT TOP 1 parameterWaterBodyCategory FROM {SQL_TABLE} WHERE parameterWaterBodyCategory = 'RW'")
    ]
    sql=[sqltest(*p) for p in prefixes]
    cat=parse_catalog_info()
    share,tokens=inspect_public_share()
    listings=[]
    for token in tokens:
        listings+=test_public_dav(token)
    successes=[x for x in sql if x.get("status")=="RESULTS_JSON" and x.get("n_actual_rows_0_or_1",0)>0]
    year_ok=[x for x in successes if x["tag"] in ("filter_year_2018","filter_year_2022","filter_2022_N_code","filter_2022_P_code")]
    out={
       "schema":"GESUI_20261116_OFFICIAL_EEA_ICM_2026_SOURCE_ACCESS_NEXTCLOUD_AND_SIMPLE_SQL_V2",
       "source_2026_publish_date":"2026-07-07",
       "source_2026_datahub_url":"https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1",
       "WISE6_public_warning_disco_tables_not_up_to_date":True,
       "latest_disco_is_NOT_VERIFIED_equal_to_2026_publisher":True,
       "sql_public_bounded_predicate_probes":sql,
       "EEA_official_metadata_catalogue":cat,
       "EEA_official_download_share_page":share,
       "EEA_public_nextcloud_webdav_metadata_only":listings,
       "single_predicate_year_rows_found_in_disco":bool(year_ok),
       "official_publisher_2026_filelist_accessible":any(x.get("n_items",0)>1 for x in listings),
       "measured_2018_2022_same_EEA_station_TN_and_TP_NOT_FOUND_OR_TESTED":True,
       "no_actual_site_measurements_published_in_git_output":True,
       "no_EQR_HNS_EIV_or_causal_model_run":True}
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"EEA_2026_ACCESS_V2_SOURCE_ONLY.json").write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_EEA_ICM_2026_SOURCE_ACCESS_V2="+json.dumps({
        "SQL":[{"tag":x["tag"],"status":x["status"],"http":x["http_status"],"result_n":x.get("n_actual_rows_0_or_1"),
                "returned_year":x.get("returned_year"),"determinand":x.get("returned_determinand_code"),
                "errorcode":x.get("error_codes")} for x in sql],
        "catalogue":{k:cat.get(k) for k in ["http_status","catalogue_parse","public_eea_source_resource_links"]},
        "share":{k:share.get(k) for k in ["http_status","is_nextcloud","html_contains_public_filelist_markup"]},
        "webdav":[{"endpoint":i["test"],"http":i.get("http_status"),"n_items":i.get("n_items"),
                   "entries":i.get("webdav_listing_entries")} for i in listings],
        "discodata_still_not_2026_verified":True,
        "publisher_2026_filelist_accessible":out["official_publisher_2026_filelist_accessible"]
    },ensure_ascii=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
