#!/usr/bin/env python3
"""2026 USGS public Farmington 5-stream PFAS source-only file-and-key audit.

DOI 10.5066/P139C2SN (Wesner et al., USGS CC0).
Target external five streams from ONE Farmington River watershed.
Published data contain water, sediment, biofilm, larval/adult insects,
tetragnathid spiders, stable isotopes, and wet mass. These are NOT
documented matched municipal wastewater effluent outfall up/down groups.

No PFAS numeric outcome effects opened. Do not transfer Grand River
microplastic effects into independent PFAS effect estimates.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,time,urllib.error,urllib.parse,urllib.request
from collections import Counter
from pathlib import Path

SOURCE_DOI="10.5066/P139C2SN"
SOURCE_STABLE_ITEM="69b80121b66b010d85582e8b"
API="https://www.sciencebase.gov/catalog/item/"
SOURCE="https://doi.org/"+SOURCE_DOI
TARGET_TITLE="Per- and polyfluoroalkyl substance concentrations in water"
UA="gesui-Farmington-five-stream-public-source-schema-audit/2026-11-18"
MAX_META_BYTES=6_000_000
MAX_TABULAR_BYTES=2_000_000
SCHEMA="GESUI_20261118_USGS_FARMINGTON_FIVE_STREAM_2020_23_PFAS_BIOTA_SOURCE_ONLY_V1"

def sha(b):return hashlib.sha256(b).hexdigest()
def get(url,limit=MAX_META_BYTES):
    if not (url.startswith("https://www.sciencebase.gov/") or
            url.startswith("https://sciencebase.gov/")):
        raise RuntimeError("Source fetch domain not USGS ScienceBase")
    err=None
    for i in range(3):
        try:
            rq=urllib.request.Request(url,headers={
              "User-Agent":UA,"Accept":"application/json,text/csv,text/tab-separated-values,*/*"})
            with urllib.request.urlopen(rq,timeout=60) as f:
                b=f.read(limit+1)
                if len(b)>limit:raise RuntimeError("Source exceeds byte ceiling, not ingested")
                return b
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
            err=e
            if isinstance(e,urllib.error.HTTPError) and e.code not in (429,500,502,503,504):
                break
            time.sleep(1+i*2)
    raise RuntimeError(f"USGS public source request failed: {str(err)[:300]}")

def site_csv_schema(b,name):
    if b[:2]==b"PK" or name.lower().endswith((".xlsx",".xls",".zip")):
        return {"status":"BINARY_TABLE_NOT_INSPECTED","original_file_name":name,
                "no_numeric_values_read":True}
    text=b.decode("utf-8-sig","replace")
    delim="\t" if name.lower().endswith((".tab",".tsv")) else ","
    try:
        rdr=csv.DictReader(io.StringIO(text),delimiter=delim)
        cols=list(rdr.fieldnames or [])
        if not cols:raise RuntimeError("Empty header")
        if any(len(k)>160 for k in cols):raise RuntimeError("Wrong delimiter or binary file")
        # This is only a *site-year-stage* SOURCE opportunity audit.
        relevant=[k for k in cols if re.search(
          r"(site|stream|reach|year|date|sample|stage|taxa|taxon|spider|insect|matrix|compound|analyte|depth|unit|latitude|longitude|wwtp|outfall)",k,re.I)]
        uniq={k:set() for k in relevant if re.search(
          r"(site|stream|reach|year|stage|taxa|taxon|matrix|unit|wwtp|outfall)",k,re.I)}
        n=0;blank=Counter()
        for row in rdr:
            if not row or not any(str(v or "").strip() for v in row.values()):continue
            n+=1
            for k in relevant:
                v=str(row.get(k) or "").strip()
                if not v:blank[k]+=1
                elif k in uniq and len(uniq[k])<5000:uniq[k].add(v)
        return {"status":"ORIGINAL_FULL_SMALL_TABULAR_SCHEMA_INSPECTED",
          "original_file_name":name,"original_rows":n,"original_n_columns":len(cols),
          "original_columns":cols[:120],"relevant_source_fields":relevant[:80],
          "source_field_nonempty_rows":{k:n-blank[k] for k in relevant},
          "source_unique_site_year_stage_counts":{k:len(v) for k,v in uniq.items()},
          "source_site_year_stage_levels":{k:sorted(v)[:100] for k,v in uniq.items()
            if len(v)<=75},
          "no_raw_numeric_PFAS_concentrations_fitted":True}
    except (ValueError,RuntimeError) as e:
        return {"status":"SOURCE_SCHEMA_PARSE_STOP","error":str(e)[:180],
                "file":name}

def run(args):
    args.out.mkdir(parents=True,exist_ok=True)
    root=API+SOURCE_STABLE_ITEM+"?format=json"
    report={"schema":SCHEMA,"dataset_DOI":SOURCE_DOI,
       "ScienceBase_catalog_item":SOURCE_STABLE_ITEM,
       "source_item_public_url":API+SOURCE_STABLE_ITEM,
       "original_publication_DOI":"10.1021/acs.est.5c16942",
       "source_five_streams_not_five_independent_watersheds":True,
       "source_is_PFAS_not_municipal_wastewater_P_or_N":True,
       "no_new_cross_basin_spider_transfer_effect_computed":True}
    try:
        b=get(root)
        meta=json.loads(b)
        report["source_metadata_SHA256"]=sha(b)
        report["source_metadata_id"]=meta.get("id")
        report["source_metadata_title"]=meta.get("title")
        if meta.get("id")!=SOURCE_STABLE_ITEM or TARGET_TITLE.lower() not in str(meta.get("title","")).lower():
            raise RuntimeError("USGS catalogue data item/title mismatch")
        files=meta.get("files",[])
        report["original_ScienceBase_public_file_count"]=len(files)
        report["published_file_descriptors"]=[]
        if not files:raise RuntimeError("USGS source metadata has no published public file descriptors")
        for f in files:
            rec={"name":f.get("name"),"size":f.get("size"),
                 "contentType":f.get("contentType"),
                 "downloadUri":f.get("downloadUri"),
                 "source_original_file_sha256":f.get("checksum"),
                 "original_source_metadata_fields":sorted(f)}
            uri=f.get("downloadUri")
            size=f.get("size")
            if (isinstance(uri,str) and uri.startswith("https://www.sciencebase.gov/") and
                isinstance(size,int) and 0<size<=MAX_TABULAR_BYTES and
                str(f.get("name","")).lower().endswith((".csv",".tsv",".tab"))):
                try:
                    raw=get(uri,limit=MAX_TABULAR_BYTES)
                    rec["actual_source_bytes_sha256"]=sha(raw)
                    rec["original_schema"]=site_csv_schema(raw,str(f.get("name")))
                except Exception as e:
                    rec["original_schema"]={"status":"REMOTE_FILE_RETRIEVAL_STOP",
                       "error_type":type(e).__name__,"error":str(e)[:200]}
            else:
                rec["original_schema"]={"status":"METADATA_ONLY_SIZE_FORMAT_OR_LOCATION_UNSUPPORTED"}
            report["published_file_descriptors"].append(rec)
        report["status"]="ORIGINAL_USGS_PUBLISHED_ITEM_AND_FILES_SOURCE_AUDITED"
    except Exception as e:
        report["status"]="SOURCE_ACCESS_OR_ID_STOP_NOT_BIOLOGICAL_NULL"
        report["error_type"]=type(e).__name__
        report["error_message"]=str(e)[:400]
    path=args.out/"FARMINGTON_USGS_2020_2023_FIVE_STREAM_ORIGINAL_PFAS_STAGE_SOURCE_QA.json"
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_FARMINGTON_PFAS_SOURCE="+json.dumps({
      "status":report["status"],
      "title":report.get("source_metadata_title"),
      "n_published_files":report.get("original_ScienceBase_public_file_count"),
      "errors":report.get("error_message")},ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())
