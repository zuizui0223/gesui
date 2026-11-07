#!/usr/bin/env python3
"""Read-only source-availability probe for Peters 2026 figshare and EEA ICM 2026.

No empirical biology or nutrient effect calculated; no large bulk data pulled.
Only request source metadata with strict byte limit and recording of response
code, redirect destination, content-type, short sample, and selected file lists.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import urlparse
import requests

DOI="10.25400/lincolnuninz.30546371"
LIMIT=600000
SOURCES=[
("figshare_global","https://api.figshare.com/v2/articles/30546371"),
("lincoln_institutional","https://data.lincoln.ac.nz/api/articles/30546371"),
("lincoln_article","https://data.lincoln.ac.nz/articles/dataset/Global_Nutrient_Yields/30546371"),
("lincoln_article_underscore","https://data.lincoln.ac.nz/articles/dataset/_/30546371"),
("figshare_doi","https://doi.org/"+DOI),
("eea_waterbase_landing","https://www.eea.europa.eu/en/datahub/datahubitem-view/fbf3717c-cd7b-4785-933a-d0cf510542e1"),
("eea_webdav","https://sdi.eea.europa.eu/webdav/datastore/public/eea_t_waterbase-water-quality-icm-2026_p_1900-2025_v01_r00/"),
("eea_discodata","https://discodata.eea.europa.eu/"),
]

def request_one(label,url):
    try:
        r=requests.get(url,headers={"User-Agent":"gesui-nutrient-coverage-source-gate/0.1",
                                    "Accept":"application/json,text/html;q=0.8,*/*;q=0.6",
                                    "Range":"bytes=0-599999"},
                       timeout=(10,28),stream=True,allow_redirects=True)
        raw=b""
        for chunk in r.iter_content(65536):
            raw+=chunk
            if len(raw)>LIMIT:
                raw=raw[:LIMIT]
                break
        r.close()
        resp={"source":label,"initial_url":url,"final_url":r.url,
              "status":r.status_code,
              "headers":{x:r.headers.get(x) for x in
                ("Content-Type","Content-Length","Content-Range","Last-Modified","ETag")
                 if r.headers.get(x) is not None},
              "bytes_sampled":len(raw),
              "sample_sha256":hashlib.sha256(raw).hexdigest(),
              "truncated_or_larger_than_byte_limit":bool(len(raw)==LIMIT)}
        parsed=None
        if raw.lstrip().startswith((b"{",b"[")):
            try:parsed=json.loads(raw)
            except json.JSONDecodeError:pass
        if isinstance(parsed,dict):
            fields=("id","title","doi","published_date","license","license_id","url_public_html","url",
                    "is_embargoed","total_size")
            resp["json_metadata"]={key:parsed.get(key) for key in fields if key in parsed}
            if isinstance(parsed.get("files"),list):
                resp["file_count"]=len(parsed["files"])
                resp["files"]=[{k:v.get(k) for k in
                    ("id","name","size","is_link_only","download_url","supplied_md5","computed_md5")}
                    for v in parsed["files"]][:150]
            if "error" in parsed:resp["response_error"]=str(parsed["error"])[:500]
            if "message" in parsed:resp["response_message"]=str(parsed["message"])[:500]
        snippet=raw.decode("utf-8","replace")
        if "<html" in snippet.lower():
            filenames=re.findall(r'(?:href=["\'])([^"\'?]+)(?:["\'])',snippet,re.I)
            resp["html_links_sample"]=[x for x in filenames if any(w in x.lower()
                for w in ("nutri","water","monitor","data","csv","sqlite","zip","download","article"))][:60]
        resp["text_sample"]=re.sub(r"\s+"," ",snippet[:750])
        return resp
    except Exception as ex:
        return {"source":label,"initial_url":url,"request_error":type(ex).__name__+": "+str(ex)[:350]}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    records=[request_one(name,url) for name,url in SOURCES]
    args.out.parent.mkdir(parents=True,exist_ok=True)
    payload={"schema":"gesui_peters2026_EEA2026_source_probe_v01",
             "retrieved_utc":datetime.now(timezone.utc).isoformat(),
             "data_not_analyzed":True,"no_bulk_data_download":True,
             "source_doi":DOI,"responses":records}
    args.out.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+"\n")
    for x in records:
        print("GESUI_NUTRIENT_SOURCE_PROBE",
              json.dumps({k:v for k,v in x.items() if k not in ("text_sample","files","html_links_sample")},
                         ensure_ascii=False,sort_keys=True),flush=True)
        if "files" in x:
            print("GESUI_NUTRIENT_SOURCE_FILES",
                  json.dumps({"source":x["source"],"files":x["files"][:70]},ensure_ascii=False),flush=True)
    return payload

if __name__=="__main__":main()
