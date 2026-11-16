#!/usr/bin/env python3
"""Bounded check of the real OFFICIAL Austrian GZÜV river station public share.
The URL is copied unchanged from verified original BMLUK dataset distribution
b054ce75-f0b1-449e-b11c-9603a629c9c5, NOT guessed.
Does not bypass authentication or open 2022 chemistry or benthos values.
"""
from __future__ import annotations
import argparse,base64,hashlib,json,re
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin,urlparse
import requests

SHARE="https://docs.umweltbundesamt.at/s/SNJqq7TEprxJtpJ"
SCHEMA="GESUI_20261116_ORIGINAL_BMLUK_GZUEV_RIVER_OFFICIAL_PUBLIC_SHARE_FILE_DISCOVERY_V1"

class HTML(HTMLParser):
    def __init__(self):super().__init__();self.a=[];self.scripts=[];self.inputs=[];self.forms=[]
    def handle_starttag(self,tag,attrs):
        attrs=dict(attrs)
        if tag=="a":self.a.append(attrs)
        elif tag=="script":self.scripts.append(attrs)
        elif tag=="input":self.inputs.append(attrs)
        elif tag=="form":self.forms.append(attrs)

def sha(s):return hashlib.sha256(s).hexdigest()

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    ss=requests.Session()
    ss.headers.update({"User-Agent":"gesui-original-public-GZUEV-stations-academic/1.0",
                       "Accept":"text/html,application/xhtml+xml"})
    r=ss.get(SHARE,timeout=(12,55),allow_redirects=True)
    r.raise_for_status()
    if urlparse(r.url).hostname!="docs.umweltbundesamt.at":
        raise RuntimeError("Original official BMLUK shared geodata redirected outside Austrian UBA source host")
    h=r.text
    p=HTML();p.feed(h)
    linked=[]
    for x in p.a:
        href=x.get("href","")
        if not href:continue
        src=urljoin(r.url,href)
        if urlparse(src).hostname!="docs.umweltbundesamt.at":continue
        linked.append({"url":src,"title":x.get("title",""),"download":x.get("download")})
    possible_downloads=[l for l in linked if (
        "download" in l["url"].lower() or l.get("download") is not None or
        "zip" in l["url"].lower() or "gpkg" in l["url"].lower())]
    initial=[]
    for x in p.inputs:
        if "initial-state" in (x.get("id","")+x.get("name","")).lower():
            val=x.get("value","")
            item={"public_element_id":x.get("id",""),"public_state_base64_chars":len(val)}
            if val:
                try:
                    text=base64.b64decode(val).decode("utf-8")
                    item["public_decoded_original_state_excerpt"]=text[:1100]
                    item["possible_original_public_dataset_filenames"]=[
                        t for t in re.findall(r'[^"\s<>]{1,130}\.(?:zip|gpkg|shp|geojson|csv|gml)',text,re.I)[:18]]
                except Exception as e:
                    item["source_state_decoding_error_type"]=type(e).__name__
            initial.append(item)
    snippets=[]
    for m in list(re.finditer(r'(?:download|filename|filesize|publicshare|sharetoken|webdav|fileName)',h,re.I))[:20]:
        excerpt=re.sub(r'[\s\n]+',' ',h[max(0,m.start()-60):min(len(h),m.end()+140)])
        snippets.append(excerpt[:230])
    result={
       "schema":SCHEMA,
       "original_official_catalog_dataset_id":"b054ce75-f0b1-449e-b11c-9603a629c9c5",
       "original_public_share_URL":SHARE,
       "original_site_http_status":r.status_code,
       "original_public_response_SHA256":sha(r.content),
       "original_public_html_bytes":len(r.content),
       "final_source_url":r.url,
       "public_HTML_title_match":re.findall(r"<title[^>]*>(.*?)</title>",h,re.I|re.S)[:3],
       "public_HTML_form_actions":[{k:v for k,v in x.items() if k in ("action","method","id")} for x in p.forms[:12]],
       "source_site_static_scripting_filepaths":[x.get("src") for x in p.scripts if x.get("src")][:18],
       "original_public_static_link_count":len(linked),
       "public_page_downloads_explicit_in_anchor_tags":possible_downloads[:25],
       "public_page_other_same_origin_links":linked[:20],
       "public_initial_state_excerpts":initial[:18],
       "public_download_related_markup_snippets":snippets,
       "original_GZUEV_geo_file_bytes_NOT_YET_VERIFIED":True,
       "source_not_2022_TNTP_or_benthos_measurements":True,
       "no_unlisted_download_path_attempted":True,
       "no_private_auth_or_NDVI_benthos_ecology":True}
    (a.output/"ORIGINAL_BMLUK_GZUEV_PUBLIC_SHARE_HTML_LINK_SOURCE_AUDIT.json").write_text(
       json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_GZUEV_ORIGINAL_SHARE_SOURCE="+json.dumps({
       "http":r.status_code,"response_bytes":len(r.content),"title":result["public_HTML_title_match"],
       "same_host_links":len(linked),
       "explicit_public_download_anchors":possible_downloads[:8],
       "initial_state_blocks":len(initial),
       "decoded_metadata_filenames":[v.get("possible_original_public_dataset_filenames") for v in initial],
       "public_share_is_not_yet_2022_water_benthos_values":True},ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
