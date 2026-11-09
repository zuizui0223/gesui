#!/usr/bin/env python3
"""Official EEA 2026 Waterbase WISE6 2022 TP archive schema SOURCE ONLY.

Discodata SQL returned 403, but official bulk published Waterbase2026 ZIP
URLs were HEAD-verified 200. Downloads original aggregated annual WQ and
spatial reference ZIP, hashes bytes, inspects ZIP member filenames/headers
ONLY. Never reads 2022 TP concentration values, images, butterfly outcomes.
No original source file redistributed through GitHub commits/artifacts.
"""
from __future__ import annotations
import argparse,hashlib,io,json,zipfile,urllib.request,time
from datetime import datetime,timezone
from pathlib import Path

SOURCES={
 "2026_WISE6_aggregated_original":
   "https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_T_WISE6_AggregatedData",
 "2026_WISE_SOe_spatial_original":
   "https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_S_WISE_SpatialObject_DerivedData",
}
MAX_BYTES={"2026_WISE6_aggregated_original":260_000_000,
           "2026_WISE_SOe_spatial_original":9_000_000}

def download(url,out,limit):
    h=hashlib.sha256();total=0
    req=urllib.request.Request(url,headers={
        "User-Agent":"gesui-2026-original-EEA-WISE6-public-WQ-source-only-ZIP-header-audit/1.0",
        "Accept":"application/zip,*/*"})
    with urllib.request.urlopen(req,timeout=135) as r,out.open("wb") as f:
        if r.status!=200:raise ValueError("Official Waterbase source download not 200")
        content_type=r.headers.get("Content-Type")
        for b in iter(lambda:r.read(1024*1024),b""):
            total+=len(b)
            if total>limit:raise ValueError("Published official WQ zip exceeded allowed pilot download max")
            f.write(b);h.update(b)
    if not zipfile.is_zipfile(out):raise ValueError("Original EEA source is not ZIP")
    return {"download_sha256":h.hexdigest(),"download_bytes":total,
            "original_content_type":content_type}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--tmpdir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    a.tmpdir.mkdir(parents=True,exist_ok=True);reports={}
    for key,url in SOURCES.items():
        path=a.tmpdir/(key+".zip")
        record={"original_EEA_public_URL":url,"contains_only_original_Waterbase2026":True}
        try:
            source=download(url,path,MAX_BYTES[key])
            record["source_download"]=source
            with zipfile.ZipFile(path) as z:
                entries=[]
                for f in z.infolist():
                    if f.is_dir():continue
                    item={"original_filename":f.filename,
                          "original_uncompressed_bytes":f.file_size,
                          "original_ZIP_compressed_bytes":f.compress_size,
                          "original_zip_crc32":f.CRC}
                    if f.filename.lower().endswith((".csv",".txt")):
                        with z.open(f) as stream:
                            line=stream.readline(15000).decode("utf-8-sig",errors="replace").strip()
                            delimiter="\t" if line.count("\t")>line.count(",") else ","
                            header=line.split(delimiter)
                            item["source_first_header_names"]=header[:120]
                            item["source_column_number"]=len(header)
                            item["one_original_header_line_only_no_data_rows_read"]=True
                    entries.append(item)
                record["original_ZIP_members_source_only"]=entries[:80]
        except Exception as e:
            record["source_access_error"]=type(e).__name__+": "+str(e)
        reports[key]=record
        print("GESUI_EEA2026_WISE6_OFFICIAL_SOURCE_ZIP_METADATA",
              json.dumps({"source":key,"download_bytes":record.get("source_download",{}).get("download_bytes"),
               "members":record.get("original_ZIP_members_source_only",[])[:5],
               "access_error":record.get("source_access_error")},sort_keys=True),flush=True)
    out={
      "schema":"gesui_EEA2026_WISE6_original_bulk_archives_2022_TP_metadata_header_only_v01",
      "retrieved_utc":datetime.now(timezone.utc).isoformat(),
      "original_bulks":reports,
      "original_disaggregated_2026_2GB_not_downloaded":True,
      "original_EEA_WISE6_SQL_403_not_absence_of_data":True,
      "no_river_TP_numeric_outcomes_or_nutrient_concentrations_read":True,
      "no_butterfly_photo_individual_count_or_EQR_outcomes_read":True,
      "only_original_published_ZIP_archives_inspected":True,
      "no_original_200MB_zip_stored_in_GitHub_repository":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_EEA2026_WISE6_ORIGINAL_2022_TP_ACCESSIBLE_BULK_ZIP_SCHEMA_ONLY_COMPLETE",flush=True)
if __name__=="__main__":main()
