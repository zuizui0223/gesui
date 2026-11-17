#!/usr/bin/env python3
"""Original 67 Austria same-water 2022 source cohort vs public HR-VPP TPROD discovery.

SOURCE-ONLY: restores official original source field sites, scans Copernicus
public STAC 10m LAEA Total Productivity collections and obtains *metadata*
and short HTTP range header only. No ecological endpoint opened, no vegetation
raster pixel values, no fitted model, no source-gate relaxation.
"""
from __future__ import annotations
import argparse
import collections
from datetime import datetime
import hashlib
import json
import re
import sys
from pathlib import Path
from urllib.parse import urlsplit
import requests

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_austria_wisa_official_source"))
from AT68_2022_real_WWTP_P_to_ambient_TP_to_benthic_EQR import cohort  # noqa: E402

STAC="https://stac.dataspace.copernicus.eu/v1"
EXPECTED={"15_3":9,"15_4":14,"15_5":21,"16_4":13,"16_5":10}
EXCLUDE="FW54110117"
SCHEMA="GESUI_20261117_ORIGINAL_AT67_2022_COPERNICUS_HRVPP_TPROD_10M_LAEA_PUBLIC_SOURCE_GATE_V1"
UA="gesui-20261117-official-HRVPP-2022-10m-source-only-no-biology"
TARGET_YEAR=2022
MAX_COLLECTION_PAGES=20
REGIONS=list(sorted(EXPECTED))

def req_json(session,url,**kwargs):
    r=session.get(url,timeout=(12,60),**kwargs)
    r.raise_for_status()
    q=r.json()
    if not isinstance(q,dict):raise RuntimeError("Original STAC service no JSON object")
    return q

def source_site_frame():
    source=cohort()
    sites=[dict(z) for z in source if z["station"]!=EXCLUDE]
    if (len(sites)!=67 or len({z["station"] for z in sites})!=67 or
        dict(collections.Counter(z["block"] for z in sites))!=EXPECTED):
        raise RuntimeError("Original physically sampled source 67 river sites 5 blocks changed")
    pick=[]
    for block in REGIONS:
        xx=[z for z in sites if z["block"]==block]
        s=sorted(xx,key=lambda z:hashlib.sha256(z["station"].encode()).hexdigest())[0]
        pick.append({"source_station":s["station"],"block":block,
                     "lat":float(s["lat"]),"lon":float(s["lon"])})
    return pick,{"source_genuine_2022_water_monitor_sites":67,
                 "source_region_station_counts":EXPECTED,
                 "source_site_exact_EEA_and_official_BMLUK_2022_TP_and_EQR":True,
                 "one_known_954m_spatial_mismatch_site_excluded_BEFORE_raster":EXCLUDE}

def collection_scan(session):
    url=STAC+"/collections?limit=100"
    seen=set()
    items=[]
    urls=[]
    for page in range(MAX_COLLECTION_PAGES):
        if url in seen:raise RuntimeError("STAC collections next-link cycle")
        seen.add(url)
        q=req_json(session,url)
        if not isinstance(q.get("collections"),list):
            raise RuntimeError("STAC collections page missing list")
        items.extend(q["collections"])
        urls.append(urlsplit(url).path)
        nexturls=[x.get("href") for x in q.get("links",[])
                  if x.get("rel")=="next" and isinstance(x.get("href"),str)]
        if not nexturls:
            break
        url=nexturls[0]
    return items,urls,len(nexturls)>0

def source_matches(items):
    all_prod=[]
    for c in items:
        name=(str(c.get("id",""))+" "+str(c.get("title",""))).lower()
        if "tprod" in name or "total productivity" in name:
            all_prod.append(c)
    exact=[]
    for c in all_prod:
        name=(str(c.get("id",""))+" "+str(c.get("title",""))).lower()
        if "europe" in name and "laea" in name and re.search(r"\b10\s*m\b|10m",name):
            exact.append(c)
    return exact,all_prod

def check_asset(session,item):
    assets=item.get("assets",{})
    short={"original_STAC_item_id":item.get("id"),
           "source_asset_keys":list(assets),
           "source_signed_or_restricted_href_not_exposed":True,
           "public_cog_header_status":"STOP_NO_PUBLIC_TPROD_TIFF_ASSET"}
    found=[]
    for key,val in assets.items():
        href=val.get("href") or ""
        low=(str(key)+" "+str(val.get("type"))+" "+str(href)).lower()
        if ("tprod" in low or any(s in low for s in ("geotiff","image/tiff",".tif"))) and href.startswith("https://"):
            found.append((key,href))
    short["eligible_tiff_asset_count"]=len(found)
    for key,href in found[:2]:
        try:
            r=session.get(href,headers={"Range":"bytes=0-8191"},
                timeout=(12,45),stream=True)
            body=next(r.iter_content(chunk_size=8192),b"")
            is_tiff=body[:4] in (b"II*\x00",b"MM\x00*",b"II+\x00",b"MM\x00+")
            public=r.status_code in (200,206) and is_tiff
            short["public_cog_test"]={
                "asset_key":key,"http_status":r.status_code,
                "received_bytes_first_chunk":len(body),
                "actually_tiff_header":is_tiff,
                "public_read_without_credentials":public,
                "asset_provider":urlsplit(href).netloc
            }
            r.close()
            if public:
                short["public_cog_header_status"]="PASS_PUBLIC_TIFF_HEADER_SOURCE_ONLY"
                return short
        except (requests.RequestException,StopIteration) as exc:
            short["public_cog_test"]={"asset_key":key,
                     "error_type":type(exc).__name__,"public_read_without_credentials":False}
    if not found and assets:
        short["public_cog_header_status"]="STOP_NO_PUBLIC_TPROD_TIFF_ASSET"
    else:
        short["public_cog_header_status"]="STOP_PUBLIC_RASTER_REQUIRES_AUTH_OR_NO_RANGE"
    return short

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    out={"schema":SCHEMA,"status":"SOURCE_ONLY_RUNNING",
         "prior_2022_S2_same_station_10site_natural_mask_smoke_5_of_10_STOP_UNCHANGED":True,
         "source_new_2022_VPP_TPROD_metadata_not_riparian_biomass":True,
         "no_2022_terrestrial_ecology_regression_or_WWTP_causality":True,
         "source_reconstruction":None}
    try:
        sites,source=source_site_frame()
        out["source_reconstruction"]=source
        out["sha_selected_stations_per_original_region"]=[{
          "site":x["source_station"],"region":x["block"]} for x in sites]
        print("GESUI_HRVPP_ORIGINAL_AT67_SOURCE=PASS n67 five original blocks",flush=True)
    except Exception as exc:
        out["status"]="STOP_ORIGINAL_AUSTRIA67_PROVENANCE"
        out["source_error"]=type(exc).__name__+": "+str(exc)[:240]
        finish(a,out)
        return

    session=requests.Session()
    session.headers.update({"User-Agent":UA,"Accept":"application/json"})
    try:
        cols,queries,pagelimit=collection_scan(session)
    except Exception as exc:
        out["status"]="STOP_CATALOG_REACHABILITY"
        out["source_error"]=type(exc).__name__+": "+str(exc)[:240]
        finish(a,out)
        return
    candidates,allproducts=source_matches(cols)
    out["source_official_STAC_endpoint"]=STAC
    out["source_total_collection_records_scanned"]=len(cols)
    out["source_collection_pages_scanned"]=len(queries)
    out["source_collection_scan_pagelimit_not_exhaustive"]=pagelimit
    out["source_any_TPROD_collection_titles"]=[
        {"id":p.get("id"),"title":p.get("title")}
        for p in allproducts]
    out["source_predeclared_2022_Europe_LAEA_10m_TPROD_collection_matches"]=[
        {"id":p.get("id"),"title":p.get("title")}
        for p in candidates]
    if pagelimit or len(candidates)!=1:
        out["status"]="STOP_NO_UNIQUE_2022_10M_TPROD_COLLECTION"
        finish(a,out)
        return
    c=candidates[0];cid=c["id"]
    out["matched_public_STAC_collection_id"]=cid
    out["source_collection_temporal_extent"]=c.get("extent",{}).get("temporal")
    checks=[]
    for site in sites:
        lon,lat=site["lon"],site["lat"]
        params={"collections":cid,
                "bbox":f"{lon-.005},{lat-.005},{lon+.005},{lat+.005}",
                "datetime":"2022-01-01T00:00:00Z/2022-12-31T23:59:59Z",
                "limit":3}
        q={"site":site["source_station"],"region":site["block"],
           "source_2022_metadata_items":0}
        try:
            itemq=req_json(session,STAC+"/search",params=params)
            feats=itemq.get("features") or []
            q["source_2022_metadata_items"]=len(feats)
            q["source_item_ids"]=[x.get("id") for x in feats]
            q["sample_item_source"]=[{"id":x.get("id"),
               "datetime":x.get("properties",{}).get("datetime"),
               "assets":list(x.get("assets",{}))}
               for x in feats[:2]]
            if feats:
                q.update(check_asset(session,feats[0]))
        except Exception as exc:
            q["source_error"]=type(exc).__name__+": "+str(exc)[:200]
        checks.append(q)
        print("GESUI_HRVPP_METADATA_SITE_REGION="+str(q["region"])+
              " n="+str(q["source_2022_metadata_items"])+
              " tiff="+str(q.get("public_cog_header_status")),flush=True)
    out["source_sha_selected_one_per_region_metadata_probe"]=checks
    n_meta=sum(x.get("source_2022_metadata_items",0)>0 for x in checks)
    n_tiff=sum(x.get("public_cog_header_status")=="PASS_PUBLIC_TIFF_HEADER_SOURCE_ONLY" for x in checks)
    out["source_metadata_regions_with_2022_items"]=n_meta
    out["source_public_TIFF_header_regions"]=n_tiff
    if n_meta<5:out["status"]="STOP_2022_10M_TPROD_FIVE_REGION_ITEM_SUPPORT"
    elif n_tiff<5:out["status"]="STOP_PUBLIC_2022_10M_TPROD_TIFF_ACCESS"
    else:out["status"]="PASS_PUBLIC_SOURCE_HEADERS_5REGIONS_NO_NUMERIC_RIPARIAN_MODEL"
    finish(a,out)

def finish(a,out):
    name=a.output/"AT67_HRVPP_2022_10M_TPROD_STAC_SOURCE_ONLY_GATE.json"
    name.write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_AT67_2022_HRVPP_10M_STAC_SOURCE_ONLY="+json.dumps({
        k:out.get(k) for k in ("status","source_reconstruction",
                "source_total_collection_records_scanned",
                "source_any_TPROD_collection_titles",
                "matched_public_STAC_collection_id",
                "source_metadata_regions_with_2022_items",
                "source_public_TIFF_header_regions")
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
