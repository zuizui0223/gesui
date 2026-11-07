#!/usr/bin/env python3
"""Source schema gate for 2023 European plant nutrient-affinity indicators.

Does not fit any source-plant association or change the 18 source-blind pairs.
Probes public Masaryk/JUICE European Ellenberg CSV and GBIF species taxonomy.
All sample species keys are fixed before looking at their trait indicator.
"""
from __future__ import annotations
from collections import Counter
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,io,json
import pandas as pd
import requests

URL="https://www.sci.muni.cz/botany/juice/NET/sources/EIV_Europe.csv"
GBIF="https://api.gbif.org/v1/species/"
YEAR_WINDOW="2022,2024"
MAX_BYTES=8_000_000
def get(url,cap):
    r=requests.get(url,timeout=(15,65),stream=True,
       headers={"User-Agent":"gesui-plantae-European-Ellenberg-EIV-source-gate/0.1",
                "Accept":"text/csv,application/json,*/*"})
    status=r.status_code
    r.raise_for_status()
    data=bytearray()
    for chunk in r.iter_content(131072):
        data.extend(chunk)
        if len(data)>cap:raise ValueError("Source was larger than fixed cap")
    r.close()
    return bytes(data),status,r.headers.get("Content-Type","")
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--plants",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    records=json.loads(args.plants.read_text())
    assert len(records["plants"])==36
    keys=sorted({int(k) for row in records["plants"] for k in row["specieskey_counts"]})
    if len(keys)!=1006:raise ValueError("Fixed pretrait GBIF species union changed")
    raw,status,mime=get(URL,MAX_BYTES)
    report={"source":URL,"source_status":status,"source_content_type":mime,
            "source_bytes":len(raw),"source_md5":hashlib.md5(raw).hexdigest(),
            "source_sha256":hashlib.sha256(raw).hexdigest(),
            "source_first_1000_characters":raw[:1000].decode("utf-8","replace")}
    decoded=None
    for enc in ("utf-8-sig","latin1"):
        try:
            text=raw.decode(enc)
            frame=pd.read_csv(io.StringIO(text),sep=None,engine="python",low_memory=False)
            decoded=enc
            break
        except Exception:continue
    if decoded is None:raise ValueError("EIV European source cannot be decoded")
    report.update(decoded_source_encoding=decoded,
                  parsed_rows=len(frame),
                  parsed_columns=list(map(str,frame.columns)),
                  first_three_species_rows=frame.head(3).astype(str).to_dict(orient="records"),
                  missing_by_column={str(k):int(frame[k].isna().sum()) for k in frame.columns})
    if len(frame)<1000 or len(frame.columns)<4:
        raise ValueError("Not the European indicator species table")
    sample=[keys[int(i*(len(keys)-1)/11)] for i in range(12)]
    taxa=[]
    for k in sample:
        try:
            b,st,m=get(GBIF+str(k),30_000)
            item=json.loads(b)
            taxa.append({"speciesKey":k,"status":st,
                         "scientificName":item.get("scientificName"),
                         "canonicalName":item.get("canonicalName"),
                         "rank":item.get("rank"),"status_taxon":item.get("taxonomicStatus"),
                         "kingdom":item.get("kingdom"),"phylum":item.get("phylum"),
                         "acceptedKey":item.get("acceptedKey")})
        except Exception as e:
            taxa.append({"speciesKey":k,"source_error":type(e).__name__+": "+str(e)[:180]})
    output={"schema":"gesui_2023_Ellenberg_Europe_plus_GBIF_taxon_source_gate_v01",
            "date_utc":datetime.now(timezone.utc).isoformat(),
            "indicator_source_reference":"Tichy et al 2023 JVS e13168 / Zenodo v2 DOI 10.5281/zenodo.7427088",
            "trait_source_licensed_original_not_mirrored":True,
            "original_new_GBIF_2022_2024_speciesKeys":len(keys),
            "EIV_Europe_csv_source":report,
            "frozen_12_speciesKey_sample_indices":sample,
            "GBIF_taxon_metadata_sample":taxa,
            "nutrient_affinity_is_ecological_optimum_not_site_P_concentration":True,
            "no_observation_effect_fitted":True}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(output,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_EIV_SOURCE_SCHEMA",json.dumps({k:v for k,v in report.items()
           if k not in ("first_three_species_rows","source_first_1000_characters","missing_by_column")},ensure_ascii=False),flush=True)
    print("GESUI_EIV_GBIF_TAXON_SAMPLE",json.dumps(taxa,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
