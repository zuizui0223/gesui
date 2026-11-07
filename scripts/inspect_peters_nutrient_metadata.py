#!/usr/bin/env python3
"""Inspect source-locked 2026 Peters TN/TP observational site/year metadata.

Download only the five small public CSV/text artifacts (no predicted fields);
verify file md5 as specified by 2026-10-29 figshare v5. Report schema/quality
before any response matching, with no source raw site catalogue redistribution.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import io
import json
from pathlib import Path
import re
import requests
import pandas as pd

ART="https://api.figshare.com/v2/articles/30546371"
EXPECTED_DOI="10.25400/lincolnuninz.30546371.v5"
NAMES=("!README.txt","TN.2_meta.csv","TP.2_meta.csv","TN.2_years.csv","TP.2_years.csv")
CAP=5_000_000
def get(url,max_bytes):
    r=requests.get(url,timeout=(15,90),stream=True,
                   headers={"User-Agent":"gesui-plants-nutrient-observation-source/0.1"})
    r.raise_for_status()
    data=bytearray()
    for b in r.iter_content(131072):
        data.extend(b)
        if len(data)>max_bytes:
            r.close()
            raise ValueError("Source exceeded fixed byte cap: "+url)
    r.close()
    return bytes(data)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    article=json.loads(get(ART,500_000))
    if article.get("doi")!=EXPECTED_DOI:
        raise ValueError("Figshare source version changed unexpectedly")
    files={x["name"]:x for x in article["files"]}
    if set(NAMES)-set(files):
        raise ValueError("Required measured-data meta/years files unavailable")
    info=[]
    for name in NAMES:
        f=files[name]
        blob=get(f["download_url"],CAP)
        if hashlib.md5(blob).hexdigest()!=f.get("supplied_md5"):
            raise ValueError("Figshare file md5 mismatch: "+name)
        rec={"name":name,"file_id":f["id"],"source_size":f["size"],
             "content_md5":hashlib.md5(blob).hexdigest(),
             "bytes_downloaded":len(blob),
             "data_source":"original_observational_site_metadata" if "meta" in name else "observation_year_inventory"}
        if name=="!README.txt":
            txt=blob.decode("utf-8","replace")
            terms=("TN.2_meta","TP.2_meta","TN.2_years","TP.2_years","Latitude","Longitude",
                   "HYRIV","lat","lon","data source")
            snippets=[]
            lines=txt.splitlines()
            for i,line in enumerate(lines):
                if any(q.lower() in line.lower() for q in terms):
                    snippets.append({"line":i+1,"excerpt":" ".join(lines[max(0,i-1):min(len(lines),i+3)])[:850]})
            rec["readme_source_schema_excerpt"]=snippets[:45]
            rec["readme_first_1700_chars"]=txt[:1700]
        else:
            df=pd.read_csv(io.BytesIO(blob),low_memory=False)
            rec.update(row_count=len(df),column_names=list(map(str,df.columns)),
                       column_dtypes={str(k):str(v) for k,v in df.dtypes.items()},
                       first_values={str(k):[str(t)[:90] for t in df[k].dropna().head(2).tolist()]
                           for k in df.columns if re.search(r"lat|lon|coord|time|year|region|source|unit|order|hyriv|date|id",str(k),re.I)},
                       null_counts={str(k):int(df[k].isna().sum()) for k in df.columns})
            if "meta" in name:
                rec["coordinate_candidates"]=[str(k) for k in df.columns if re.search(r"lat|lon|east|north|coord|x_|y_",str(k),re.I)]
                rec["river_network_candidates"]=[str(k) for k in df.columns if re.search(r"hyriv|reach|hydro|stream|river",str(k),re.I)]
        info.append(rec)
        print("GESUI_VERIFIED_PETERS2026_FILE",
              json.dumps({k:v for k,v in rec.items() if k not in ("first_values","null_counts",
                                                                 "readme_source_schema_excerpt","readme_first_1700_chars")},
                         ensure_ascii=False),flush=True)
    out={"schema":"gesui_peters2026_figshare_v5_observational_meta_schema_gate_v01",
         "timestamp_utc":datetime.now(timezone.utc).isoformat(),
         "source_doi":article["doi"],"article_id":article["id"],
         "article_published_date":article.get("published_date"),
         "article_license":article.get("license"),
         "file_count_article":len(article["files"]),
         "used_file_names":list(NAMES),
         "no_original_site_rows_redistributed":True,
         "no_nutrient_to_plant_effect_test":True,
         "file_schema":info}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_MEASURED_NUTRIENT_SCHEMA_GATE_PASS",json.dumps({
       "files":len(info),"article_id":article["id"],"version_doi":article["doi"]},sort_keys=True))
if __name__=="__main__":main()
