#!/usr/bin/env python3
"""Source-only European macroinvertebrate sample×taxon availability inventory.

Do NOT download the raw multihundred-MB biological datasets at this phase.
This tests *only* official source identity, licensing, manifest, HTTP HEAD and
small byte-range response characteristics. No EQR-to-trait or nutrient model
may be fit from an EQR-only EEA aggregate or from a metadata listing.
"""
from __future__ import annotations

import argparse, hashlib, html.parser, json, time, urllib.error, urllib.parse, urllib.request
from pathlib import Path

CHEMPOP_ID="23f194b0-1e21-4f14-be9c-6eabdcc1feb7"
PAN_EU_ID=15268354
URLS={
    "eea_biology_landing":"https://www.eea.europa.eu/en/datahub/datahubitem-view/7bd881f0-59ca-4f0a-a43c-7e0cd0c66d6c",
    "chempop_v2_metadata":"https://catalogue.ceh.ac.uk/documents/"+CHEMPOP_ID+"?format=json",
    "chempop_v2_data":"https://data-package.ceh.ac.uk/data/"+CHEMPOP_ID,
    "chempop_v2_support":"https://data-package.ceh.ac.uk/sd/"+CHEMPOP_ID+".zip",
    "pan_eu_2010_2015_api":"https://zenodo.org/api/records/"+str(PAN_EU_ID),
    "ea_biosys_surveys":"https://environment.data.gov.uk/ecology/explorer/downloads/INV_OPEN_DATA.zip",
    "ea_biosys_taxon_info":"https://environment.data.gov.uk/ecology/explorer/downloads/OPEN_DATA_TAXON_INFO.zip",
}
ALLOWED_DOMAINS={
    "www.eea.europa.eu","catalogue.ceh.ac.uk","data-package.ceh.ac.uk",
    "zenodo.org","environment.data.gov.uk",
}
UA="gesui-source-first-taxonomic-data-audit/20261119 (only metadata and tiny headers)"
LIMIT=131072

class SafeRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        host=urllib.parse.urlsplit(newurl).hostname
        if host not in ALLOWED_DOMAINS:raise ValueError("Redirect outside official source domains: "+str(host))
        return super().redirect_request(req,fp,code,msg,headers,newurl)


opener=urllib.request.build_opener(SafeRedirect())

def request(name, url, method, limit=4096):
    hdr={"User-Agent":UA,"Accept":"application/json,text/html,text/plain,*/*"}
    if method=="RANGE":hdr["Range"]="bytes=0-4095"
    r=urllib.request.Request(url,headers=hdr,method="GET" if method=="RANGE" else method)
    record={"source":name,"method":method,"url":url,"status":None,"content_type":None,
            "content_length":None,"final_url":None,"read_bytes":0,
            "content_sample_sha256":None,"error":None}
    t0=time.time()
    try:
        with opener.open(r,timeout=18) as response:
            record["status"]=int(response.status)
            record["final_url"]=response.geturl()
            record["content_type"]=response.headers.get("Content-Type","")
            record["content_length"]=response.headers.get("Content-Length","")
            record["accept_ranges"]=response.headers.get("Accept-Ranges","")
            record["content_range"]=response.headers.get("Content-Range","")
            data=b"" if method=="HEAD" else response.read(limit)
            record["read_bytes"]=len(data)
            if data:record["content_sample_sha256"]=hashlib.sha256(data).hexdigest()
    except urllib.error.HTTPError as ex:
        record["status"]=ex.code
        record["error"]=str(ex)[:260]
    except Exception as ex:
        record["error"]=type(ex).__name__+": "+str(ex)[:260]
    record["elapsed_seconds"]=round(time.time()-t0,2)
    return record

def official_metadata():
    url=URLS["chempop_v2_metadata"]
    r=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":"application/json"})
    info={"source":"chempop_v2_metadata","status":"SOURCE_UNVERIFIED"}
    try:
        with opener.open(r,timeout=18) as response:
            raw=response.read(LIMIT)
            if not raw:raise RuntimeError("Empty metadata")
            d=json.loads(raw)
        if (d.get("id")!=CHEMPOP_ID or str(d.get("version")) not in ("2","2.0")
            or "macroinvertebrate" not in d.get("title","").lower()
            or "2025" not in d.get("publicationDate","")):
            raise RuntimeError("Official v2 / title / publication identity mismatch")
        uri=[z.get("url") for z in d.get("onlineResources",[]) if z.get("function")=="download"]
        info.update({"status":"ORIGINAL_SOURCE_METADATA_VERIFIED",
            "source_id":d.get("id"),"version":d.get("version"),
            "publication_date":d.get("publicationDate"),
            "metadata_sha256":hashlib.sha256(raw).hexdigest(),
            "download_uris":uri,"licence":d.get("licences",[]),
            "title":d.get("title"),
            "endpoint_is_original_1519_England_sites_not_EEA_BQE":True})
    except Exception as ex:
        info["status"]="SOURCE_METADATA_STOP"
        info["error"]=type(ex).__name__+": "+str(ex)[:320]
    return info


def zenodo_metadata():
    url=URLS["pan_eu_2010_2015_api"]
    r=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":"application/json"})
    info={"source":"pan_eu_2010_2015_api","status":"SOURCE_UNVERIFIED"}
    try:
        with opener.open(r,timeout=18) as response:
            raw=response.read(LIMIT)
            d=json.loads(raw)
        files=d.get("files",[])
        if d.get("id")!=PAN_EU_ID or not files:
            raise RuntimeError("Zenodo record id/files mismatch")
        info.update({"status":"SOURCE_METADATA_VERIFIED","record_id":d.get("id"),
             "original_created":d.get("created"),
             "license":d.get("metadata",{}).get("license",{}),
             "files":[{"key":f.get("key"),"size_bytes":f.get("size"),
                       "checksum":f.get("checksum"),
                       "links":f.get("links",{}).get("self")} for f in files],
             "not_collocated_EEA_2018_2022_WWTP_sources":True})
    except Exception as ex:
        info["status"]="SOURCE_METADATA_STOP"
        info["error"]=type(ex).__name__+": "+str(ex)[:320]
    return info

def main(a):
    result={
      "schema":"GESUI_EU_BENTHIC_TAXON_SOURCE_ONLY_20261119_V1",
      "date":"2026-11-19",
      "source_only":True,"do_not_fit_biological_response":True,
      "eea_waterbase2026_has_normalized_EQR_not_site_taxon_counts":True,
      "alternative_data_are_independent_populations_not_same_EEA_source_sample":True,
      "not_measured_WWTP_routed_receiving_river_NP":True,
      "metadata":{},
      "http_probes":{},
    }
    result["metadata"]["chempop"]=official_metadata()
    result["metadata"]["pan_eu"]=zenodo_metadata()
    for key,url in URLS.items():
        if key.endswith("metadata") or key.endswith("api"):continue
        head=request(key,url,"HEAD")
        result["http_probes"][key]={"HEAD":head}
        if key!="eea_biology_landing":
            # Capture only first 4KiB even when server ignores Range.
            # Never persist the binary probe, only hash and headers.
            result["http_probes"][key]["RANGE"]=request(key,url,"RANGE",limit=4096)
    result["analysis_gate"]={
      "multi_country_taxon_analysis":"SOURCE_ONLY_NO_SOURCE_SHAPE_OR_TRAIT_LINK_VERIFIED",
      "England_same_sample_sewage_trait":"SOURCE_ONLY_NO_JOINABLE_TAXON_ABUNDANCE_YET",
      "reusing_EEA2018_2022_EQR_as_taxon_occurrences":"STOP_EEA_ONLY_EQR",
      "need": [
          "Verified actual taxon×sample original tables and licencing",
          "Paired genuinely measured/comparable nutrient load or modelled effluent exposure",
          "Native site+sample+year identity with nonleaking geographic source groups",
          "Taxon identifiers and official reproducible tolerance-score lookup",
          "Treatment-vs-urban+agriculture+flow controls and negative controls",
          "Pre-outcome support thresholds before any trait/environment model"
      ]
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_BENTHIC_SOURCE_ONLY="+json.dumps({
        "chempop_source_status":result["metadata"]["chempop"]["status"],
        "pan_eu_source_status":result["metadata"]["pan_eu"]["status"],
        "http":{key:{m:{"status":v["status"],"content_type":v["content_type"],
                         "content_length":v["content_length"],"error":v["error"],
                         "read_bytes":v["read_bytes"]}
                     for m,v in values.items()} for key,values in result["http_probes"].items()},
        "analysis_gate":result["analysis_gate"],
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
