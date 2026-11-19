#!/usr/bin/env python3
"""Original England ChemPop v2 ZIP schema-only source inspection.

Requires exact EIDC v2 UUID. Bounded original bulk download, SHA256, ZIP
members and CSV headers only, no taxon outcomes, environmental regression,
or exposed/controls selection. No full source rehosting in Git history/artifact.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,tempfile,urllib.request,zipfile
from pathlib import Path

UUID="23f194b0-1e21-4f14-be9c-6eabdcc1feb7"
URL="https://data-package.ceh.ac.uk/data/"+UUID+".zip"
MAXBYTES=400*1024*1024
HEADERS={"User-Agent":"gesui-official-chempop-v2-source-only/20261119",
         "Accept":"application/zip,application/octet-stream,*/*"}

def get_head():
    rec={"url":URL,"status":None,"content_type":None,"content_length":None,"error":None}
    try:
        req=urllib.request.Request(URL,method="HEAD",headers=HEADERS)
        with urllib.request.urlopen(req,timeout=25) as r:
            rec["status"]=int(r.status)
            rec["content_type"]=r.headers.get("Content-Type","")
            s=r.headers.get("Content-Length")
            rec["content_length"]=int(s) if s else None
    except Exception as e:
        rec["error"]=type(e).__name__+": "+str(e)[:350]
    return rec

def make(a):
    result={"schema":"GESUI_CHEMPOP_V2_AUTHENTIC_ZIP_FILES_SCHEMA_V1",
            "original_EIDC_UUID":UUID,"source_licence":"Open Government Licence",
            "source_only":True,"no_outcomes_models_fitted":True,
            "not_same_site_sample_as_original_2018_2022_EEA_BQE":True,
            "HEAD":get_head(),"status":"SOURCE_STOP","members":[]}
    head=result["HEAD"]
    if head["status"]!=200:
        result["reason"]="OFFICIAL_ZIP_HEAD_NOT_200"
    elif head["content_length"] is None or head["content_length"]>MAXBYTES:
        result["reason"]="SOURCE_FILE_TOO_LARGE_OR_UNKNOWN_SIZE_FOR_THIS_BOUNDED_PILOT"
    else:
        try:
            with tempfile.TemporaryDirectory(prefix="gesui_chem_v2_") as temp:
                file=Path(temp)/"v2_original.zip"
                h=hashlib.sha256();n=0
                req=urllib.request.Request(URL,headers=HEADERS)
                with urllib.request.urlopen(req,timeout=45) as r, file.open("wb") as fd:
                    if r.status!=200:raise RuntimeError("Original archive stream failed")
                    while True:
                        chunk=r.read(2**20)
                        if not chunk:break
                        n+=len(chunk)
                        if n>MAXBYTES:raise RuntimeError("Original archive exceeds maximum")
                        fd.write(chunk);h.update(chunk)
                if n!=head["content_length"]:
                    raise RuntimeError("Original ZIP HEAD size != actual received size")
                if not zipfile.is_zipfile(file):
                    raise RuntimeError("Official v2 downloaded file is not ZIP")
                with zipfile.ZipFile(file) as z:
                    files=z.infolist()
                    if len(files)>300:raise RuntimeError("Unexpected huge number of files")
                    manifest=[]
                    for f in files:
                        meta={"name":f.filename,"size_uncompressed_bytes":f.file_size,
                              "size_compressed_bytes":f.compress_size,"CRC32":format(f.CRC,"08x")}
                        if not f.is_dir() and f.filename.lower().endswith(".csv"):
                            with z.open(f) as fi:
                                first=fi.readline(80000)
                            if len(first)>65536:raise RuntimeError("A source CSV header exceeds bound")
                            meta["columns"]=next(csv.reader([first.decode("utf-8-sig","replace").strip()]))
                        manifest.append(meta)
                    result["status"]="SOURCE_AUTHENTIC_ZIP_MEMBERS_AND_CSV_HEADERS_VERIFIED"
                    result["ZIP_source_SHA256"]=h.hexdigest()
                    result["ZIP_size_bytes"]=n
                    result["members"]=manifest
        except Exception as e:
            result["reason"]=type(e).__name__+": "+str(e)[:420]
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_CHEMPOP_V2_ZIP_SCHEMA="+json.dumps({
        "status":result["status"],"reason":result.get("reason"),
        "HEAD":result["HEAD"],"n_files":len(result["members"]),
        "files":[{"name":x["name"],"columns":x.get("columns",[]),
                  "size":x["size_uncompressed_bytes"]} for x in result["members"]],
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    make(p.parse_args())
