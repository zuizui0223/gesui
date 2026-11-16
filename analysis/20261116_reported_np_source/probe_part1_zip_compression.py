#!/usr/bin/env python3
"""Safely read only 256-byte ZIP local header for original EEA ICM Part1.

Discover compression method that Python 3.11 zipfile cannot extract.
If server ignores HTTP byte range, stop after 256 bytes and close.
No full archive, no individual N/P samples.
"""
from __future__ import annotations
import argparse,hashlib,json,subprocess
from pathlib import Path
import requests
URL="https://sdi.eea.europa.eu/datashare/public.php/dav/files/ZcJ37WK6oaWonFJ/WISE6_DisaggregatedData-csv.zip"
EXPECTED=1639231716
METHOD_NAMES={0:"STORE",8:"DEFLATE",9:"DEFLATE64",12:"BZIP2",14:"LZMA",93:"ZSTD",20:"ZSTD_OLD",99:"AES"}
def cmd(cmd):
    p=subprocess.run(cmd,capture_output=True,text=True,timeout=8,check=False)
    return {"returncode":p.returncode,"stdout_first":p.stdout.strip()[:280],"stderr_first":p.stderr.strip()[:150]}
def main(a):
    with requests.get(URL,headers={"Range":"bytes=0-255","Accept":"application/zip",
        "User-Agent":"gesui-Part1-original-only-header-compression-diagnostic/20261116"},
        timeout=(8,30),stream=True) as r:
        r.raise_for_status()
        data=next(r.iter_content(chunk_size=256),b"")
        meta={"http_status":r.status_code,"content_range":r.headers.get("Content-Range"),
              "content_type":r.headers.get("Content-Type"),
              "content_length":r.headers.get("Content-Length")}
    if len(data)<30 or data[:4]!=b"PK\x03\x04":raise RuntimeError("Original EEA first byte-range is not ZIP local entry header")
    method=int.from_bytes(data[8:10],"little")
    name_len=int.from_bytes(data[26:28],"little")
    source_name=data[30:30+name_len].decode("utf-8","replace") if len(data)>=30+name_len else "not_wholly_in_256_byte_range"
    out={"schema":"GESUI_EEA_2026_PART1_ORIGINAL_ZIP_COMPRESSION_PREFLIGHT_V1",
         "original_publisher_url":URL,"original_listed_ZIP_bytes":EXPECTED,
         "file_local_header_magic":"PK0304",
         "first_entry_compression_method":method,"first_entry_method_name":METHOD_NAMES.get(method,"OTHER_UNSUPPORTED"),
         "first_entry_uncompressed_filename_or_partial":source_name,
         "original_ZIP_full_SHA256_from_independent_run_43798990690":"ff779521749b40b952cce35a4671240c3b6c47d87f32a8fef9240a4ac503a082",
         "response":meta,
         "bsdtar_availability":cmd(["bash","-lc","command -v bsdtar || true"]),
         "sevenzip_availability":cmd(["bash","-lc","command -v 7z || command -v 7zz || true"]),
         "unzip_availability":cmd(["bash","-lc","command -v unzip || true"]),
         "no_all_data_or_chemistry_extracted":True}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"EEA2026_PART1_ZIP_COMPRESSION_SOURCE_PREFLIGHT.json").write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_EEA_PART1_ZIP_COMPRESSION="+json.dumps(out,ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path);main(p.parse_args())
