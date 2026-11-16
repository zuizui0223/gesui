#!/usr/bin/env python3
"""Official BMLUK 2020–2022 report archive: stage genuine 2022 FW1.1
annual river chemistry workbook, FW2 older 2019–2021 source ONLY.

Before seeing numeric results, record official site-year workbook provenance
and original archive member names, never claim a WWTP effect or 2022 benthos.
Temporary signed Actions artifact only, no Git original Excel rehosting.
"""
from __future__ import annotations
import argparse,hashlib,json,zipfile
from pathlib import Path
import requests

SOURCE="https://www.bmluk.gv.at/dam/jcr%3A24629723-d7da-465f-9c1f-048783b885a5/Wasserg%C3%BCte-Jahresbericht_2020-2022_Ergebnistabellen.zip"
MAX_ZIP_BYTES=4_000_000
SCHEMA="GESUI_20261116_BMLUK_ORIGINAL_2020_2022_FW1_2022_ANNUAL_SOURCE_STAGING_V1"

def sha(buf):return hashlib.sha256(buf).hexdigest()

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    r=requests.get(SOURCE,timeout=(15,80),headers={"User-Agent":"gesui-original-Austrian-government-WG2020-2022-source-inventory/20261116"})
    if r.status_code!=200 or not r.url.startswith("https://www.bmluk.gv.at/"):
        raise RuntimeError("BMLUK official 2020–22 report result workbook ZIP not available from original government")
    data=r.content
    if not (150000<len(data)<MAX_ZIP_BYTES and data[:2]==b"PK"):
        raise RuntimeError("Original national report ZIP response invalid size/magic")
    archive=a.output/"TEMP_BMLUK_OFFICIAL_SOURCE_2020_2022.zip"
    archive.write_bytes(data)
    source=[]
    staged=[]
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:
            raise RuntimeError("Original published Austrian report ZIP has an invalid CRC entry")
        for ent in z.infolist():
            if ent.is_dir():continue
            file=Path(ent.filename)
            ext=file.suffix.lower()
            record={"source_report_original_member":ent.filename,
              "filename":file.name,
              "size_uncompressed":ent.file_size,
              "source_extension":ext}
            if ext==".xlsx" and file.name.upper().startswith(("FW1","FW2","FW_1","FW_2")):
                if ent.file_size>3_000_000:
                    raise RuntimeError("Official annual river report source workbook too large")
                payload=z.read(ent.filename)
                if payload[:2]!=b"PK":
                    raise RuntimeError("Publisher native Excel workbook is not OpenXML")
                dest=a.output/file.name
                dest.write_bytes(payload)
                record["official_workbook_SHA256"]=sha(payload)
                record["temporarily_staged_source_file"]=file.name
                staged.append(record)
            source.append(record)
    if not any(z["filename"].upper().startswith(("FW1","FW_1")) for z in staged):
        raise RuntimeError("Missing official FW1 original 2022 annual chemistry report workbook")
    archive.unlink()
    manifest={"schema":SCHEMA,
      "original_BMLUK_source_page":"https://www.bmluk.gv.at/themen/wasser/wasserqualitaet/wassergute-jahresbericht-2020-2022.html",
      "original_public_ZIP_URL":SOURCE,
      "source_original_ZIP_SHA256":sha(data),
      "source_ZIP_bytes":len(data),
      "source_original_members":source,
      "river_FW1_FW2_original_workbooks":staged,
      "original_biology_FW2_if_present_is_2019_to2021_not_2022":True,
      "actual_2022_total_N_and_total_P_same_PKINT_not_yet_verified":True,
      "no_EDA_EEA_or_organism_effect_values_used":True}
    (a.output/"BMLUK_2020_2022_FW1_FW2_PUBLISHER_SOURCE_MANIFEST.json").write_text(
        json.dumps(manifest,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_BMLUK_2020_2022_OFFICIAL_FW1_SOURCE="+json.dumps(
      {k:v for k,v in manifest.items() if k not in ("source_original_members",)},ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True);main(p.parse_args())
