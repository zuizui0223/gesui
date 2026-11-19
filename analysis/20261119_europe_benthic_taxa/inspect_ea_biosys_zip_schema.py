#!/usr/bin/env python3
"""Strict source-only EA Biosys ZIP/Taxon Info schema inspection + ChemPop manifest.

Does not read biological outcomes, do joins, or model wastewater. Downloads
original EA public archives temporarily on GitHub Actions to inspect ZIP
members, CSV HEADERS and minimal structural source information; raw archives
and biological records are never uploaded or committed.
"""
from __future__ import annotations
import argparse, csv, hashlib, html.parser, io, json, os, shutil, tempfile, time, urllib.request, zipfile
from pathlib import Path

CHEMPOP="https://data-package.ceh.ac.uk/data/23f194b0-1e21-4f14-be9c-6eabdcc1feb7"
BENTHOS="https://environment.data.gov.uk/ecology/explorer/downloads/INV_OPEN_DATA.zip"
TAXON="https://environment.data.gov.uk/ecology/explorer/downloads/OPEN_DATA_TAXON_INFO.zip"
PUBYEAR=2026
HEADERS={"User-Agent":"gesui-scientific-original-EA-source-schema-audit/20261119",
         "Accept":"*/*"}

class Links(html.parser.HTMLParser):
    def __init__(self):
        super().__init__();self.hrefs=[];self.title=""
    def handle_starttag(self,tag,attrs):
        if tag.lower()=="a":
            attrs=dict(attrs);href=attrs.get("href")
            if href:self.hrefs.append(href)

def chempop_manifest():
    info={"source_url":CHEMPOP,"status":"UNRESOLVED"}
    try:
        req=urllib.request.Request(CHEMPOP,headers=HEADERS)
        with urllib.request.urlopen(req,timeout=35) as response:
            if response.status!=200:raise RuntimeError("HTTP "+str(response.status))
            raw=response.read(30000)
            ct=response.headers.get("Content-Type","")
        if "text/html" not in ct:raise RuntimeError("Data URL not HTML directory")
        txt=raw.decode("utf-8","replace")
        links=Links();links.feed(txt)
        info.update({"status":"HTML_INDEX_OBTAINED",
            "content_length_read":len(raw),
            "sha256":hashlib.sha256(raw).hexdigest(),
            "hyperlinks":links.hrefs,
            "source_title_text_snippet":txt[:2200]})
    except Exception as e:
        info.update({"status":"SOURCE_STOP","error":type(e).__name__+": "+str(e)[:250]})
    return info

def download_temp(url,target,maxbytes):
    req=urllib.request.Request(url,headers=HEADERS)
    s=hashlib.sha256();size=0
    with urllib.request.urlopen(req,timeout=40) as src,target.open("wb") as dst:
        if src.status!=200:raise RuntimeError("No original ZIP success")
        clen=src.headers.get("Content-Length")
        if clen and int(clen)>maxbytes:raise RuntimeError("Source ZIP exceeds fixed audit-only size ceiling")
        while True:
            chunk=src.read(1024*1024)
            if not chunk:break
            size+=len(chunk)
            if size>maxbytes:raise RuntimeError("Source ZIP exceeded fixed audit-only size ceiling mid-download")
            s.update(chunk);dst.write(chunk)
    return {"source_bytes":size,"source_sha256":s.hexdigest()}

def inspect_zip(url,folder,name,maxbytes):
    info={"source":name,"url":url,"status":"SOURCE_STOP"}
    with tempfile.TemporaryDirectory(prefix="gesui_benthos_zip_") as tmp:
        target=Path(tmp)/"official.zip"
        try:
            info.update(download_temp(url,target,maxbytes=maxbytes))
            if not zipfile.is_zipfile(target):
                raise RuntimeError("Downloaded original source is not a valid ZIP")
            with zipfile.ZipFile(target) as z:
                bad=z.testzip()
                if bad is not None:raise RuntimeError("ZIP crc failed")
                members=z.infolist()
                if len(members)>100:raise RuntimeError("Unexpected >100 files in official ZIP; stop")
                output=[]
                for mem in members:
                    r={"member":mem.filename,"uncompressed_size":mem.file_size,
                       "compressed_size":mem.compress_size,"zip_crc32":format(mem.CRC,"08x"),
                       "is_dir":mem.is_dir()}
                    if not mem.is_dir() and mem.filename.lower().endswith((".csv",".txt")) and mem.file_size:
                        with z.open(mem) as handle:
                            head=handle.readline(65536)
                            if len(head)>60000:raise RuntimeError("Unreasonably long CSV header")
                        decoded=head.decode("utf-8-sig","replace").strip()
                        try:r["csv_header"]=next(csv.reader([decoded]))
                        except Exception:r["header_raw"]=decoded[:500]
                    output.append(r)
                info["status"]="OFFICIAL_ORIGINAL_ZIP_CRC_MEMBERS_SCHEMA_VERIFIED"
                info["members"]=output
                info["member_count"]=len(members)
        except Exception as e:
            info.update({"status":"SOURCE_STOP","error":type(e).__name__+": "+str(e)[:350]})
    return info

def main(a):
    results={"schema":"GESUI_20261119_EA_BIOSYS_CHEMPOP_SOURCE_SCHEMA_V1",
             "source_only_no_biotic_outcomes_fitted":True,
             "not_EEA_normalized_BQE_taxon_counts":True,
             "publisher_year":PUBYEAR,
             "chempop":chempop_manifest(),
             "EA_BIOSYS_river_macroinvertebrate":inspect_zip(BENTHOS,None,"EA_BIOSYS_river_macroinvertebrate",100*1024*1024),
             "EA_BIOSYS_TaxonInfo":inspect_zip(TAXON,None,"EA_BIOSYS_TaxonInfo",5*1024*1024)}
    # Distinguish identifiers from taxonomic score system via literal source schema;
    # never infer that a taxon is tolerant from a sample-level BMWP/ASPT index.
    results["biology_go_no_go"]={
        "phase1_model_authorized":False,
        "no_source_matched_WWTP_load_test":True,
        "no_source_verified_taxon_tolerance_reference":True,
        "next":"Verify original site_sample_analysisID joins and published taxon-score values; confirm ChemPop sewage-exposure variable and year before outcome modelling"
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(results,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_EA_BIOSYS_SCHEMA_SOURCE_ONLY="+json.dumps({
        "chempop":results["chempop"],
        "biosys":{k:{"status":results[k]["status"],
                    "archive_bytes":results[k].get("source_bytes"),
                    "member_count":results[k].get("member_count"),
                    "members":[{"filename":m["member"],"columns":m.get("csv_header",[]),
                                "n_uncompressed_bytes":m["uncompressed_size"]}
                               for m in results[k].get("members",[])]}
                  for k in ("EA_BIOSYS_river_macroinvertebrate","EA_BIOSYS_TaxonInfo")},
        "biology_go_no_go":results["biology_go_no_go"]
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    main(p.parse_args())
