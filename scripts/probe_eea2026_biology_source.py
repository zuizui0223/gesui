#!/usr/bin/env python3
"""Source-only EEA Waterbase Biology 2026 EQR archive inspection.

No ecological response is modeled. Download official dataset and record byte
checksum, uncompressed file inventory and original CSV column schemas.
Avoid redistributing confidential coordinates or full licensed raw observations.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,zipfile
from datetime import datetime,timezone
from pathlib import Path
URL="https://sdi.eea.europa.eu/datashare/s/HCxTFz3Ac8xbHX6/download"
DOI_PAGE="https://www.eea.europa.eu/en/datahub/datahubitem-view/7bd881f0-59ca-4f0a-a43c-7e0cd0c66d6c"
VERSION="eea_t_waterbase-biology-2026_p_1990-2024_v01_r00"
MAX_PREVIEW=131072
MAX_FILES=400

def sha_file(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def get_schema(arc,item):
    out={"relative_path":item.filename,"original_uncompressed_bytes":item.file_size,
         "zip_crc32":format(item.CRC,"08x"),"zip_compressed_bytes":item.compress_size}
    if item.filename.lower().endswith((".csv",".txt",".tsv",".sql",".json",".xml")):
        with arc.open(item) as f:
            part=f.read(MAX_PREVIEW)
        out["first_preview_sha256"]=hashlib.sha256(part).hexdigest()
        txt=part.decode("utf-8-sig","replace")
        lines=txt.splitlines()
        if lines:
            # Try CSV or TSV header. Avoid sampling individual biological records:
            # the archive catalogue/column names are all that is needed here.
            delimiter="\t" if "\t" in lines[0] else ";"
            if lines[0].count(",")>lines[0].count(";"):delimiter=","
            try:
                header=next(csv.reader([lines[0]],delimiter=delimiter))
                out["column_header"]=header[:200]
            except Exception:out["column_header_parse_error"]=True
            out["first_line_size"]=len(lines[0])
        out["preview_truncated"]=item.file_size>MAX_PREVIEW
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--archive",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    if not args.archive.exists() or args.archive.stat().st_size<50_000:
        raise ValueError("EEA official Biology 2026 dataset missing or too small")
    if not zipfile.is_zipfile(args.archive):
        raise ValueError("Source is not a ZIP archive")
    source_sha=sha_file(args.archive)
    with zipfile.ZipFile(args.archive) as arc:
        members=[x for x in arc.infolist() if not x.is_dir()]
        if len(members)>MAX_FILES:raise ValueError("Unexpectedly large EEA archive file inventory")
        inventory=[get_schema(arc,x) for x in members]
    candidate=[x for x in inventory if any(t in x["relative_path"].lower()
        for t in ("eqr","biology","spatial","monitoring","waterbody","site","readme"))]
    report={
       "schema":"gesui_EEA_Waterbase_Biology_2026_public_source_inventory_v01",
       "source_official_page":DOI_PAGE,
       "source_version":VERSION,"download_url":URL,
       "original_archive_bytes":args.archive.stat().st_size,
       "original_archive_sha256":source_sha,
       "timestamp_utc":datetime.now(timezone.utc).isoformat(),
       "source_dataset_claimed_biology_year_coverage":"2012-2024 for EQR records",
       "source_full_archive_file_count":len(inventory),
       "candidate_BQE_or_site_source_files":candidate,
       "all_file_inventory":inventory,
       "numeric_BQE_EQR_records_analysed":False,
       "EEA_2022_WWTP_loads_crosslinked_to_EQR":False,
       "data_limitations":[
         "EEA official source is biological quality indicator EQR, not actual per-species abundance",
         "Different countries and biological quality elements require norm and classification harmonization",
         "EQR and WWTP load monitoring year and freshwater setting must be aligned",
         "Monitoring station locations and waterbody codes need verified river graph matches, not mere distance",
         "No environmental effect, causal association or species dispersion hypothesis tested at source-inventory stage"
       ]
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_EEA2026_BIOLOGY_SOURCE_INVENTORY",json.dumps({
         "official_source_sha256":source_sha,"archive_bytes":args.archive.stat().st_size,
         "source_files":len(inventory),
         "potential_biology_site_files":len(candidate),
         "potential_source_file_names":[v["relative_path"] for v in candidate][:60],
         "no_biological_outcome_analysed":True},sort_keys=True),flush=True)
    for v in candidate[:25]:
        print("GESUI_EEA2026_BIOLOGY_SCHEMA",json.dumps({
          "file":v["relative_path"],"uncompressed_bytes":v["original_uncompressed_bytes"],
          "header":v.get("column_header",[])[:70]},sort_keys=True),flush=True)

if __name__=="__main__":
    main()
