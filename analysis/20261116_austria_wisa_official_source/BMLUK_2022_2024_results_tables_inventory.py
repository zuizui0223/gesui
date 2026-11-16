#!/usr/bin/env python3
"""Source-only publisher original BMLUK 2022–2024 water-quality report result
tables official ZIP. Inventory genuine files and source type/grain only;
do not read nutrient/biotic responses or treat report years as 2022 station.
"""
from __future__ import annotations
import argparse,hashlib,json,zipfile
from pathlib import Path
import requests

URL="https://www.bmluk.gv.at/dam/jcr%3A1532ed46-4444-4430-bf89-c3507729d70f/Wasserg%C3%BCte-Jahresbericht%202022-2024%20Ergebnistabellen.zip"
SCHEMA="GESUI_20261116_OFFICIAL_BMLUK_2022_2024_REPORT_RESULT_TABLES_PUBLIC_ZIP_SOURCE_INVENTORY_V1"
MAX_SOURCE_BYTES=4_000_000
def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    r=requests.get(URL,timeout=(14,65),headers={"User-Agent":"gesui-WISA-BMLUK-2022-report-sourceonly/1.0"},
                   stream=True)
    if r.status_code!=200:
        raise RuntimeError("BMLUK official 2022-24 results tables source HTTP "+str(r.status_code))
    if not r.url.startswith("https://www.bmluk.gv.at/"):
        raise RuntimeError("Original official BMLUK source redirected away from public publisher host")
    dest=a.output/"TEMP_ORIGINAL_BMLUK_2022_2024_REPORT_TABLES.zip"
    digest=hashlib.sha256();n=0
    with dest.open("wb") as f:
        for b in r.iter_content(64*1024):
            if not b:continue
            n+=len(b)
            if n>MAX_SOURCE_BYTES:raise RuntimeError("Official published 631k source unexpectedly >4 MB")
            digest.update(b);f.write(b)
    if not zipfile.is_zipfile(dest):
        raise RuntimeError("Official report results ZIP not genuine ZIP bytes")
    with zipfile.ZipFile(dest) as z:
        rows=[]
        for inf in z.infolist():
            if inf.is_dir():continue
            suffix=Path(inf.filename).suffix.lower()
            rows.append({"original_file":inf.filename,
                  "original_member_uncompressed_bytes":inf.file_size,
                  "source_extension":suffix,
                  "original_member_is_xlsx":suffix==".xlsx",
                  "original_member_is_CSV_or_TSV":suffix in (".csv",".tsv"),
                  "only_names_grain_inventory_no_data_cells_read":True})
    if not rows:raise RuntimeError("Official report result tables archive empty")
    result={"schema":SCHEMA,
        "original_publisher_url":URL,
        "original_public_results_years_as_named":"2022–2024",
        "source_original_zip_bytes":n,"source_original_zip_SHA256":digest.hexdigest(),
        "original_public_archive_members":rows,
        "original_reports_may_be_aggregate_statistics_not_original_2022_sites":True,
        "source_no_station_2022_TNTP_MZB_values_opened":True,
        "no_wastewater_land_effect_fitted":True}
    (a.output/"BMLUK_2022_24_PUBLIC_WATER_REPORT_RESULTS_TABLES_SOURCE_INVENTORY.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    dest.unlink()
    print("GESUI_BMLUK_2022_24_REPORT_TABLE_FILES="+json.dumps(result,ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True);main(p.parse_args())
