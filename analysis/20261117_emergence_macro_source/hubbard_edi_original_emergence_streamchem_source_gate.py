#!/usr/bin/env python3
"""SOURCE ONLY: genuine EDI Hubbard Brook 2018-24 emergent aquatic insects
and long-term continuous stream chemistry source packet metadata.

Physically exported emergent insects (source counts) != riparian predator
uptake or wastewater-derived nutrients. Does not fit ecological outcomes.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,re,time,urllib.error,urllib.parse,urllib.request,xml.etree.ElementTree as ET
from pathlib import Path

BASE="https://pasta.lternet.edu"
GROUP="knb-lter-hbr"
SOURCES=[
 ("measured_emergence","455","1","495c55e97b678d14ae40c27296fd62cc"),
 ("continuous_stream_chemistry","208","14",None)]
UA="gesui-2026-11-17-EDI-hubbard-brook-verified-source-first"
META_LIMIT=2_000_000
FIRST_BYTES=60000

def request(url,limit=META_LIMIT,head=None):
    if not (url.startswith("https://pasta.lternet.edu/") or url.startswith("https://pasta-d.lternet.edu/")):
        raise RuntimeError("Only original public EDI archive endpoints allowed")
    req=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":"*/*",**(head or {})})
    last=None
    for a in range(3):
        try:
            with urllib.request.urlopen(req,timeout=65) as r:
                b=r.read(limit+1)
                if len(b)>limit:raise RuntimeError("Source material exceeds bounded request")
                return b,r.status,str(r.headers.get("Content-Type",""))
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
            last=e
            if isinstance(e,urllib.error.HTTPError) and e.code not in [429,500,502,503,504]:break
            time.sleep(2+2*a)
    raise RuntimeError("Original EDI source blocked: "+str(last))

def local(tag):return tag.split("}")[-1]
def meta_fields(root):
    # Keep original scientific attribute names, never infer units unreported.
    tables=[]
    for dt in root.iter():
        if local(dt.tag)!="dataTable":continue
        entity=next((str(e.text).strip() for e in dt.iter() if local(e.tag)=="entityName" and e.text),"unknown")
        columns=[str(e.text).strip() for e in dt.iter() if local(e.tag)=="attributeName" and e.text]
        physical=[str(e.text).strip() for e in dt.iter() if local(e.tag)=="objectName" and e.text]
        table={"original_entity_name":entity,"original_physical_file_names":physical,
               "attribute_names":columns,"n_attribute_fields":len(columns)}
        tables.append(table)
    return tables

def one(name,ident,rev,expected_doi):
    suffix=f"/{GROUP}/{ident}/{rev}"
    result={"source_key":name,"scope":GROUP,"identifier":ident,"revision":rev,
            "original_immutable_resource_map":BASE+"/package/eml"+suffix,
            "original_metadata":BASE+"/package/metadata/eml"+suffix,
            "expected DOI":"https://doi.org/10.6073/pasta/"+expected_doi if expected_doi else None,
            "source_only_no_WWTP_or_terrestrial_consumer":True}
    if expected_doi:
        raw,status,typ=request(BASE+"/package/doi/eml"+suffix,1000)
        assert status==200, "Original DOI did not return 200"
        doi=raw.decode("utf-8","replace").strip()
        result["original_DOI_response"]=doi
        if expected_doi not in doi:
            raise RuntimeError("Original EDI packet DOI and original publication disagree")
    meta,status,typ=request(BASE+"/package/metadata/eml"+suffix,META_LIMIT)
    root=ET.fromstring(meta)
    result["original_metadata_sha256"]=hashlib.sha256(meta).hexdigest()
    result["source_table_schema"]=meta_fields(root)
    mp,status,typ=request(BASE+"/package/eml"+suffix,500_000)
    text=mp.decode("utf-8","replace")
    # EDI serves resource-map line URLs for true immutable entities, not
    # just a paper citation. Retain those original links for reproducibility.
    links=sorted(set(re.findall(r'https://pasta(?:-d)?\.lternet\.edu/package/data/eml/[^\s<>"\']+',text)))
    result["original_data_entity_count"]=len(links)
    result["resource_map_sha256"]=hashlib.sha256(mp).hexdigest()
    if not links:
        raise RuntimeError("Genuine immutable EDI source did not advertise any data entities")
    scans=[]
    for link in links[:20]:
        rec={"source_entity_url":link}
        try:
            b,status,typ=request(link,FIRST_BYTES,{"Range":f"bytes=0-{FIRST_BYTES-1}"})
            rec["HTTP_status"]=status
            rec["original_content_type"]=typ
            rec["bounded_raw_bytes"]=len(b)
            rec["first_raw_bytes_SHA256"]=hashlib.sha256(b).hexdigest()
            header=b.decode("utf-8-sig","replace").splitlines()[:4]
            if header:
                sep="\t" if header[0].count("\t")>header[0].count(",") else ","
                # source only: don't release actual records, only data fields
                cols=next(csv.reader([header[0]],delimiter=sep),[])
                rec["original_header_fields"]=[k[:100] for k in cols[:90]]
                rec["first_line_n_columns"]=len(cols)
                rec["header_has_site_and_time"]=any(
                    any(t in c.lower() for t in ("site","stream","watershed","locat","wshd"))
                    for c in cols) and any(any(t in c.lower() for t in ("date","year","month","week","time")) for c in cols)
        except Exception as ex:
            rec["error_type"]=type(ex).__name__
            rec["error_summary"]=str(ex)[:200]
        scans.append(rec)
    result["entity_original_header_source_scans"]=scans
    result["status"]="PASS_ORIGINAL_EDI_SOURCE_METADATA_AND_ENTITY_IDS_NO_EFFECT" if len(scans)>0 else "STOP_EDI_SOURCE"
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--dataset",choices=("all","emergence","chemistry"),default="all")
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    results=[]
    for src in SOURCES:
        if a.dataset!="all" and (a.dataset=="emergence") != (src[0]=="measured_emergence"):continue
        try: r=one(*src)
        except Exception as ex:
            r={"source_key":src[0],"status":"STOP_ORIGINAL_EDI_SOURCE_ACCESS",
               "error_type":type(ex).__name__,"error_summary":str(ex)[:450]}
        results.append(r)
        print("EDI_SOURCE_GATE="+json.dumps({
            "dataset":src[0],"status":r["status"],
            "n_data_entity":r.get("original_data_entity_count"),
            "n_original_tables":len(r.get("source_table_schema",[]))},ensure_ascii=False),flush=True)
    result={"schema":"GESUI_20261117_HUBBARD_BROOK_MEASURED_AQUATIC_EMERGENCE_AND_STREAMCHEM_SOURCE_QA_V1",
      "datasets":results,"no_ecological_effects_or_causal_WWTP_transfer_fitted":True,
      "do_not_conflate_one_HBEF_watershed_with_independent_global_macro_sample":True}
    (a.out/f"HBEF_2018_24_EMERGENCE_CHEM_SOURCE_{a.dataset.upper()}.json").write_text(
      json.dumps(result,indent=2,ensure_ascii=False)+"\n")

if __name__=="__main__":main()
