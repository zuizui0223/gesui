#!/usr/bin/env python3
"""Source-only identification of original national GZÜV freshwater station
dataset's actual official public distributions and machine-readability.

Queries the previously authenticated documented EU Open Data catalog search.
No private API probing; no biologic outcomes. Public WFS/WMS allowed only if
metadata officially advertises endpoint (no hidden website internals).
"""
from __future__ import annotations
import argparse,hashlib,json,re,time
from pathlib import Path
from urllib.parse import urlparse
import requests

SOURCE_ENDPOINT="https://data.europa.eu/api/hub/search/search"
EXPECTED_ID="b054ce75-f0b1-449e-b11c-9603a629c9c5"
SCHEMA="GESUI_20261116_GZUEV_FLOWING_WATER_ORIGINAL_DISTRIBUTION_SOURCE_V1"
QUERY={"q":"GZÜV Messnetz Fließgewässer","filters":"dataset","limit":25}
OFFICIAL_CATALOG="data-gv-at"

def _get(session,url,**kw):
    try:
        r=session.get(url,timeout=(12,55),allow_redirects=True,**kw)
        return r
    except requests.RequestException as e:
        return {"error_type":type(e).__name__,"error":str(e)[:400]}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    session=requests.Session()
    session.headers.update({"User-Agent":"gesui-public-original-GZUEV-source-academic/1.0","Accept":"application/json"})
    r=_get(session,SOURCE_ENDPOINT,params=QUERY)
    if isinstance(r,dict):raise RuntimeError("Official documented EU catalogue original source unavailable "+str(r))
    if r.status_code!=200:raise RuntimeError("Original official EU catalog status "+str(r.status_code))
    js=r.json()
    docs=js.get("result",js).get("results",[])
    hits=[z for z in docs if z.get("id")==EXPECTED_ID]
    if len(hits)!=1:raise RuntimeError("Original Austria GZUEV flowing-water dataset not uniquely found in official catalog")
    d=hits[0]
    if d.get("catalog",{}).get("id")!=OFFICIAL_CATALOG:
        raise RuntimeError("Original GZUEV station dataset source catalog differs from Austria data.gv.at")
    names=d.get("title") or {}
    if not isinstance(names,dict) or "GZÜV" not in (names.get("de","")+names.get("en","")):
        raise RuntimeError("Original BMLUK GZUEV water source title unexpected")
    if d.get("publisher",{}).get("name")!="BMLUK":
        raise RuntimeError("Original Austria publisher unexpected")
    dist=d.get("distributions",[])
    if not isinstance(dist,list) or not dist:
        raise RuntimeError("Original GZUEV public metadata has no distributions")
    allowed=[]
    # Distribution public metadata is the authority for each URL. We do not
    # construct a URL by guessing a file name or a private application path.
    for item in dist:
        res={"id":item.get("id"),"title":item.get("title",{}),
             "format":item.get("format",{}),"license":item.get("license"),
             "access_url":item.get("access_url",[]),
             "download_url":item.get("download_url",[]),
             "access_service":item.get("access_service",[])}
        link=[]
        for key in ("download_url","access_url"):
            vals=res.get(key) or []
            if isinstance(vals,str):vals=[vals]
            link.extend((key,v) for v in vals if isinstance(v,str) and v.startswith("https://"))
        for service in res["access_service"]:
            for u in service.get("endpoint_url",[]) or []:
                if isinstance(u,str) and u.startswith("https://"):link.append(("public_catalogue_advertised_endpoint",u))
        probes=[]
        for role,url in link:
            q=_get(session,url,headers={"Accept":"*/*"},stream=True)
            if isinstance(q,dict):
                probes.append({"link_role":role,"public_catalogued_URL":url,**q})
            else:
                probes.append({"link_role":role,"public_catalogued_URL":url,
                    "http_status":q.status_code,"final_url":q.url,
                    "media_type":q.headers.get("Content-Type"),
                    "content_length":q.headers.get("Content-Length"),
                    "response_content_disposition":q.headers.get("Content-Disposition"),
                    "is_public_2xx":q.status_code in (200,206)})
                q.close()
            time.sleep(.15)
        res["official_catalogued_public_links_probed"]=probes
        allowed.append(res)
    result={
        "schema":SCHEMA,
        "official_catalogue_request":r.url,
        "response_source_SHA256":hashlib.sha256(r.content).hexdigest(),
        "official_original_dataset_ID":EXPECTED_ID,
        "title":names,
        "original_source_publisher":d.get("publisher"),
        "original_catalog":d.get("catalog",{}).get("id"),
        "original_dataset_metadata_modified":d.get("modified"),
        "original_original_catalogue_description":d.get("description"),
        "official_original_distributions_complete":allowed,
        "source_GZUEV_geographic_station_metadata_not_2022_observation_values":True,
        "source_public_2022_TNTP_and_MZB_same_station_YEAR_NOT_YET_DOWNLOADED":True,
        "no_2022_benthic_EQR_or_sewage_effect_computed":True}
    (a.output/"GZUEV_ORIGINAL_STATION_PUBLIC_DISTRIBUTION_AND_ACCESS.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_GZUEV_SOURCE_DISTRIBUTION="+json.dumps({
        "original_dataset_id":EXPECTED_ID,
        "de_title":names.get("de"),"en_title":names.get("en"),
        "distributions":len(allowed),
        "public_urls":[{"role":p["link_role"],"url":p["public_catalogued_URL"],
                        "http":p.get("http_status"),"content_type":p.get("media_type"),
                        "content_length":p.get("content_length")}
                   for d in allowed for p in d["official_catalogued_public_links_probed"]],
        "only_station_geodata_not_2022_TNTP_MZB":True},ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
