#!/usr/bin/env python3
"""Source-file format gate for actual measured TN/TP; no outcome inference."""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import io
import json
from pathlib import Path
import re
import zipfile
import requests
ART="https://api.figshare.com/v2/articles/30546371"
DOI="10.25400/lincolnuninz.30546371.v5"
FILES=("TN.2_-_Harmonised.zip","TP.2_-_Harmonised.zip")
IDS={"TN":("WQ016702","WQ016695"),"TP":("WQ976406","WQ977059","WQ977062")}
def download(url,maxbytes):
    r=requests.get(url,stream=True,timeout=(15,100),
      headers={"User-Agent":"gesui-focal-TNTP-measured-source-schema/0.1"})
    r.raise_for_status()
    b=bytearray()
    for x in r.iter_content(131072):
        b.extend(x)
        if len(b)>maxbytes:raise ValueError("source ZIP larger than bound")
    r.close()
    return bytes(b)
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    article=requests.get(ART,timeout=30).json()
    if article["doi"]!=DOI:raise ValueError("source v5 DOI changed")
    idx={x["name"]:x for x in article["files"]}
    rows=[]
    for name in FILES:
        info=idx[name]
        raw=download(info["download_url"],32_000_000)
        md5=hashlib.md5(raw).hexdigest()
        if md5!=info["supplied_md5"]:raise ValueError("original ZIP MD5 mismatch "+name)
        z=zipfile.ZipFile(io.BytesIO(raw))
        names=z.namelist()
        taxon=name.split(".")[0]
        target=IDS[taxon]
        matches=[x for x in names if any(k.lower() in x.lower() for k in target)]
        preview={}
        for path in (matches[:8] if matches else names[:4]):
            if path.endswith("/"):continue
            with z.open(path) as f: sample=f.read(2300).decode("utf-8","replace")
            preview[path]=sample[:1500]
        row={"name":name,"id":info["id"],"source_bytes":len(raw),
             "source_zip_md5":md5,"archive_entry_count":len(names),
             "archive_first_entries":names[:15],
             "target_station_file_matches":matches[:30],
             "files_sample_text":preview,
             "all_target_names":list(target),
             "source_category":"measured_harmonised_observations_not_WRTDS"}
        rows.append(row)
        print("GESUI_FOCAL_TNTP_ZIP_SCHEMA",
              json.dumps({k:v for k,v in row.items() if k not in ("files_sample_text","archive_first_entries")},
                ensure_ascii=False),flush=True)
        print("GESUI_FOCAL_TNTP_SOURCE_SAMPLE",json.dumps(preview,ensure_ascii=False),flush=True)
    out={"source_doi":DOI,"read_time_utc":datetime.now(timezone.utc).isoformat(),
         "source_licenses":article.get("license"),"sources":rows,
         "data_only_no_effects_analyzed":True}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
if __name__=="__main__":main()
