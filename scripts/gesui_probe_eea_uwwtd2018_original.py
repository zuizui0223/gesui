#!/usr/bin/env python3
"""Source-only inspect the public EEA 2017-18 UWWTD historical archive.

No biological EQR outcomes accessed. Record entire original 2017-18 ZIP
SHA256, inner file names/sizes/CRCs, bounded CSV heads and nutrient/
treatment/outfall source schema with exact primary 2018 reference-year flags.
The official public source files are NOT rehosted here.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,zipfile
from datetime import datetime,timezone
from pathlib import Path

SOURCE_PAGE="https://sdi.eea.europa.eu/data/dbbccc77-ec41-4b19-8ec5-093bba1ba7a8"
DOWNLOAD="https://sdi.eea.europa.eu/datashare/s/mGmXg3zcEqKQLTG/download"
SOURCE_VERSION="eea_t_uwwtd-waterbase_p_2017-2018_v07_r01"
MAX_PREVIEW=8192
MAX_MEMBERS=350

def digest_file(path):
    sha=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):sha.update(b)
    return sha.hexdigest()

def header_info(buf):
    decoded=buf.decode("utf-8-sig","replace")
    lines=decoded.splitlines()
    if not lines:return {}
    first=lines[0].replace("\ufeff","")
    candidates=[",",";","\t","|"]
    sep=max(candidates,key=lambda x:first.count(x))
    try:fields=next(csv.reader([first],delimiter=sep))
    except Exception:fields=[first[:1000]]
    return {"preview_delimiter":sep,"preview_column_names":fields[:200],
            "row1_preview_redacted":True}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--archive",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    path=args.archive
    if not zipfile.is_zipfile(path):
        raise ValueError("Official 2017–18 public archive is not ZIP: do not assume source access")
    sha=digest_file(path)
    records=[]
    with zipfile.ZipFile(path) as z:
        members=[x for x in z.infolist() if not x.is_dir()]
        if len(members)>MAX_MEMBERS:
            raise ValueError("Unexpected more than 350 entries in official UWWTD archive")
        for item in members:
            row={"relative_path":item.filename,"bytes_uncompressed":item.file_size,
                 "bytes_compressed":item.compress_size,
                 "source_member_CRC32":format(item.CRC,"08x")}
            name=item.filename.lower()
            if name.endswith((".csv",".txt",".tsv",".sql",".md")):
                with z.open(item) as f:buf=f.read(MAX_PREVIEW)
                row.update(header_info(buf))
                row["first_8KB_SHA256"]=hashlib.sha256(buf).hexdigest()
                row["source_binary_not_redistributed"]=True
            if name.endswith(".zip"):
                # Check nested inventories without extracting or downloading
                # additional copies; inspect source path hierarchy.
                if item.file_size<120_000_000:
                    raw=z.read(item)
                    if zipfile.is_zipfile(io.BytesIO(raw)):
                        with zipfile.ZipFile(io.BytesIO(raw)) as nested:
                            sub=[x.filename for x in nested.infolist() if not x.is_dir()]
                            row["nested_member_count"]=len(sub)
                            row["nested_member_first_names"]=sub[:50]
                else:
                    row["nested_archive_skipped_too_large"]=True
            records.append(row)
    candidates=[x for x in records if any(v in x["relative_path"].casefold()
               for v in ("uwwtp","discharge","reportperiod","readme","metadata",
                         "tables","coded","schema","plant","waterbase"))]
    potential_fields=sorted(set(field for x in candidates for field in x.get("preview_column_names",[])
                if any(k in field.lower() for k in ("nremov","premov","ndischarge","pdischarge",
                                                    "measur","totaln","totalp","year",
                                                    "load","flow","code","latitude",
                                                    "longitude","treatment","volume"))))
    report={
       "schema":"gesui_EEA_2017_2018_original_UWWTD_source_inventory_v01",
       "retrieved_utc":datetime.now(timezone.utc).isoformat(),
       "official_archive_source_page":SOURCE_PAGE,
       "official_direct_download":DOWNLOAD,
       "source_version":SOURCE_VERSION,
       "original_downloaded_archive_bytes":path.stat().st_size,
       "original_archive_sha256":sha,
       "original_archive_members":len(records),
       "all_public_original_archive_members":records,
       "potential_UWWTP_outfall_metadata":candidates,
       "possible_measured_NP_source_columns_2018":potential_fields,
       "actual_2018_NP_measured_loads_verified":False,
       "reported_source_year_units_and_field_definitions_not_yet_validated":True,
       "2018_biological_EQR_records_read":False,
       "large_original_source_ZIP_redistributed":False,
       "next_step":"Check source 2018 UWWTP and DischargePoints exact ID join, year and method, before biological response",
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_UWWTD_2018_SOURCE_ARCHIVE_INSPECTED",json.dumps({
       "file_bytes":path.stat().st_size,
       "SHA256":sha,"original_members":len(records),
       "member_names":[x["relative_path"] for x in records][:120],
       "potential_load_fields":potential_fields[:100],
       "NO_2018_BIOLOGICAL_EQR_READ":True},sort_keys=True),flush=True)
    for row in candidates:
        print("GESUI_UWWTD_2018_SOURCE_MEMBER_SCHEMA",
           json.dumps({"name":row["relative_path"],
                        "original_bytes":row["bytes_uncompressed"],
                        "original_CRC32":row["source_member_CRC32"],
                        "columns":row.get("preview_column_names"),
                        "nested_names":row.get("nested_member_first_names")},sort_keys=True),flush=True)

if __name__=="__main__":main()
