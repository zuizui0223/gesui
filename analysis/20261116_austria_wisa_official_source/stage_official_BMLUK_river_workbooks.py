#!/usr/bin/env python3
"""Byte-verify and temporarily stage *only* 2022–2024 official BMLUK
river results workbooks for local analysis with dedicated spreadsheet tools.
This stages original publicly published workbooks only into the user's
short-lived Actions artifact; never commits raw original Excel to Git.
No ecological effects fitted. Preserve source URL and SHA of each file.
"""
from __future__ import annotations
import argparse,hashlib,json,zipfile
from pathlib import Path
import requests
URL="https://www.bmluk.gv.at/dam/jcr%3A1532ed46-4444-4430-bf89-c3507729d70f/Wasserg%C3%BCte-Jahresbericht%202022-2024%20Ergebnistabellen.zip"
EXPECTED_SHA="e7cd1118cee231f73cd79596446b9a4f9aaa7ed56f57627867dc86b8ee1eb8b9"
FILES=("FW1_Fließgewässer_Chemie_2023_2024.xlsx",
       "FW2_Fließgewässer_Biologie_Überblicksmessnetz_2019.xlsx",
       "FW3_Fließgewässer_operativeMST_2022_2024.xlsx")
def shab(b):return hashlib.sha256(b).hexdigest()
def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    r=requests.get(URL,timeout=(15,75),headers={"User-Agent":"gesui-BMLUK-2022-to2024-results-science/1.0"})
    if r.status_code!=200 or not r.url.startswith("https://www.bmluk.gv.at/"):
        raise RuntimeError("Original BMLUK national quality report data not delivered by publisher")
    b=r.content
    if shab(b)!=EXPECTED_SHA:
        raise RuntimeError("Original official publisher 2022–2024 result ZIP source digest changed")
    path=a.output/"TEMP_original_BMLUK_report_source.zip"
    path.write_bytes(b)
    manifest=[]
    with zipfile.ZipFile(path) as z:
        for filename in FILES:
            if filename not in z.namelist():
                raise RuntimeError("Expected exact original Austrian national published workbook member missing")
            payload=z.read(filename)
            if len(payload)<15000 or payload[:2]!=b"PK":
                raise RuntimeError("Original public national XLSX not genuine Excel OpenXML ZIP")
            dest=a.output/filename
            dest.write_bytes(payload)
            manifest.append({"publisher":"BMLUK Austria","filename":filename,
                "source_orig_zip_SHA256":EXPECTED_SHA,
                "original_workbook_sha256":shab(payload),
                "original_workbook_bytes":len(payload),
                "source_public_url":URL,
                "no_ecological_coefficients_analyzed":True})
    path.unlink()
    (a.output/"ORIGINAL_BMLUK_RIVER_REPORT_WORKBOOK_SOURCE_MANIFEST.json").write_text(
        json.dumps(manifest,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_ORIGINAL_BMLUK_RIVER_REPORT_XLSX_STAGE="+json.dumps(manifest,sort_keys=True,ensure_ascii=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
