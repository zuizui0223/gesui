#!/usr/bin/env python3
"""USGS source-only catalog + real file schemas: industrial-wastewater PFAS.

2025 Hubbard/Kotalik original source doi:10.5066/P13BSPG4
USGS ScienceBase item 68af0ddcd4be0253bb2c7c6b.
Potential independent water→larval insect→emerged adult→riparian spider
matrix in one INDUSTRIAL effluent-receiving stream (not municipal WWTP).
Original 2026 Perrotta microbiome paper reuses same outfall and must not
be treated as independent river experiment.

Neither PFAS numeric concentrations nor ecology effect signs are opened.
"""
from __future__ import annotations
import argparse, collections,csv,hashlib,io,json,re,time,urllib.error,urllib.parse,urllib.request,zipfile
from pathlib import Path

BASE="https://www.sciencebase.gov"
ITEM="68af0ddcd4be0253bb2c7c6b"
URI=f"{BASE}/catalog/item/{ITEM}?format=json"
DOI="10.5066/P13BSPG4"
MAX_RESPONSE=15_000_000
MAX_INDIVIDUAL=6_000_000
SCHEMA="GESUI_20261118_KOTALIK_2025_INDUSTRIAL_WWTP_PFAS_2021_22_USGS_RAW_SOURCE_SCHEMA_V1"

def sha(b):return hashlib.sha256(b).hexdigest()
def fetch(uri,maxbytes=MAX_RESPONSE):
    if not (uri.startswith("https://www.sciencebase.gov/") or uri.startswith("https://sciencebase.gov/")):
        raise RuntimeError("Original source URL outside USGS ScienceBase")
    error=None
    for i in range(3):
        try:
            req=urllib.request.Request(uri,headers={
                "User-Agent":"gesui-independent-2025-USGS-industrial-PFAS-source-audit/2026-11-18",
                "Accept":"application/json,application/octet-stream,text/plain,*/*"})
            with urllib.request.urlopen(req,timeout=100) as res:
                out=res.read(maxbytes+1)
                if len(out)>maxbytes:raise RuntimeError("Original public file over source cap")
                return out
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
            error=e
            if isinstance(e,urllib.error.HTTPError) and e.code not in (429,500,502,503,504):break
            time.sleep(2+i*2)
    raise RuntimeError("Original USGS source unavailable "+type(error).__name__+":"+str(error)[:210])

def original_table(raw,filename):
    text=raw.decode("utf-8-sig","replace")
    if raw[:2]==b"PK" or filename.lower().endswith((".xls",".xlsx",".zip")):
        if raw[:2]!=b"PK":return {"kind":"unsupported_archive_binary"}
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            entries=[{"name":x.filename,"uncompressed_bytes":x.file_size} for x in z.infolist()]
        return {"kind":"original_archive_file_names_only","n":len(entries),"files":entries[:100],
                "no_spreadsheet_numeric_values_loaded":True}
    if not filename.lower().endswith((".csv",".tsv",".tab",".txt")):
        return {"kind":"other_original_text_or_binary","n_bytes":len(raw)}
    # USGS Soyoil original .txt data files are TAB-delimited, not comma-delimited.
    first=text.splitlines()[0] if text else ""
    delim="\t" if filename.lower().endswith((".tsv",".tab")) or
        (filename.lower().endswith(".txt") and "\t" in first) else ","
    reader=csv.DictReader(io.StringIO(text),delimiter=delim)
    fields=list(reader.fieldnames or [])
    if not fields or any(len(str(f))>250 for f in fields):
        return {"kind":"source_tabular_parse_failed","field_count":len(fields)}
    candidate=[f for f in fields if re.search(
        r"(site|river|stream|outfall|up|down|location|matrix|material|water|tissue|species|taxa|stage|sample|year|date|medium|type|unit|method|analyte|compound|pfas|wwtp|effluent)",f,re.I)]
    keys={f:collections.Counter() for f in candidate if re.search(
        r"(site|river|stream|location|up|down|matrix|material|tissue|species|taxa|stage|year|medium|type|unit|method|wwtp)",f,re.I)}
    seen=collections.Counter();n=0
    for row in reader:
        if not row or not any((str(v).strip() for v in row.values() if v is not None)):continue
        n+=1
        for f in candidate:
            val=str(row.get(f) or "").strip()
            if not val:continue
            seen[f]+=1
            if f in keys and len(keys[f])<1000:keys[f][val]+=1
    return {"kind":"true_original_csv_or_tsv",
       "original_data_rows_n":n,"original_col_count":len(fields),
       "original_columns":fields[:140],
       "source_site_year_stage_support_fields":candidate[:80],
       "source_nonblank_counts":dict(seen),
       "source_key_category_counts":{k:len(v) for k,v in keys.items()},
       "source_key_levels":{k:dict(v) for k,v in keys.items() if len(v)<=40},
       "PFAS_concentrations_not_fitted":True}

def run(args):
    args.out.mkdir(parents=True,exist_ok=True)
    r={"schema":SCHEMA,"source":"USGS_ScienceBase_public","source_DOI":DOI,
       "source_2025_publication_DOI":"10.1021/acs.est.5c04867",
       "source_2026_related_microbiome_DOI":"10.1021/acsesttox.6c00051",
       "source_different_publications_NOT_independent_river_replication":True,
       "industrial_food_processing_effluent_NOT_municipal_wwtp":True,
       "source_NOT_a_new_independent_ecological_effect_estimate":True}
    try:
        meta_b=fetch(URI)
        m=json.loads(meta_b)
        r["original_metadata_sha256"]=sha(meta_b)
        r["original_item_id"]=m.get("id")
        r["original_item_title"]=m.get("title")
        if m.get("id")!=ITEM or "Per- and polyfluoroalkyl" not in m.get("title",""):
            raise RuntimeError("USGS original DOI item/title changed")
        filemeta=m.get("files",[])
        if not filemeta:raise RuntimeError("Original USGS DOI contains zero files")
        r["original_source_file_count"]=len(filemeta)
        r["original_files"]=[]
        for f in filemeta:
            name=str(f.get("name") or "");size=f.get("size")
            rec={"source_name":name,"source_metadata_size":size,"source_content_type":f.get("contentType"),
                "source_checksum":f.get("checksum"),"source_uri_domain":str(urllib.parse.urlparse(str(f.get("downloadUri") or "")).netloc),
                "source_filename":name}
            uri=f.get("downloadUri")
            if (isinstance(uri,str) and uri.startswith("https://www.sciencebase.gov/") and
                isinstance(size,int) and 0<size<=MAX_INDIVIDUAL and
                name.lower().endswith((".csv",".tsv",".tab",".txt",".zip",".xlsx"))):
                try:
                    b=fetch(uri,maxbytes=MAX_INDIVIDUAL)
                    rec["download_status"]="ORIGINAL_FILE_ACTUALLY_RETRIEVED"
                    rec["original_file_sha256"]=sha(b)
                    rec["bytes"]=len(b)
                    rec["schema_or_archive_manifest"]=original_table(b,name)
                except Exception as e:
                    rec["download_status"]="STOP_ORIGINAL_FILE_FETCH_OR_SCHEMA"
                    rec["error_type"]=type(e).__name__
                    rec["error"]=str(e)[:230]
            else:rec["download_status"]="ORIGINAL_METADATA_ONLY_FILE_NOT_INGESTED"
            r["original_files"].append(rec)
        r["status"]="PASS_USGS_DOI_ITEM_FILE_SCHEMAS_SOURCE_ONLY_NO_NUMERIC_EFFECT"
    except Exception as e:
        r["status"]="STOP_ORIGINAL_USGS_ITEM_ACCESS_NOT_ECOLOGICAL_NULL"
        r["original_source_error_type"]=type(e).__name__
        r["original_source_error_summary"]=str(e)[:320]
    path=args.out/"INDUSTRIAL_PFAS_2021_2022_UPDOWN_INSECT_SPIDER_ORIGINAL_USGS_SOURCE_QA.json"
    path.write_text(json.dumps(r,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("GESUI_KOTALIK_2025_USGS_ORIGINAL_SOURCE="+json.dumps({
       "status":r["status"],"title":r.get("original_item_title"),
       "n_source_files":r.get("original_source_file_count"),
       "retrieved":[q["source_name"] for q in r.get("original_files",[])
          if q["download_status"]=="ORIGINAL_FILE_ACTUALLY_RETRIEVED"],
       "error":r.get("original_source_error_summary")},ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())
