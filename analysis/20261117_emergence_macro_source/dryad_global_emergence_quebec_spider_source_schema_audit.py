#!/usr/bin/env python3
"""Source-first audit of real *measured* aquatic emergence and riparian spider
public Dryad datasets; never mistake literature site counts or mere EPTO
occurrences for measured flux, or urbanization for WWTP discharge exposure.

No ecological effect testing, no wastewater join, no ecological conclusions.
Only published source records/columns/completeness and candidate matching keys.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,os,time,urllib.error,urllib.parse,urllib.request,zipfile
from pathlib import Path

SETS=[
 {"name":"nash_global_emergence_2023","doi":"10.5061/dryad.zw3r228cs",
  "citation":"https://doi.org/10.1111/geb.13700"},
 {"name":"gagnon_quebec_emergence_spiders_2026","doi":"10.5061/dryad.v6wwpzh8z",
  "citation":"https://doi.org/10.5061/dryad.v6wwpzh8z"}]
BASE="https://datadryad.org/api/v2"
UA="gesui-source-audit-public-aquatic-emergence-v1/2026-11-17"
MAX_DOWNLOAD=2_500_000

def sha(x):return hashlib.sha256(x).hexdigest()
def get(url,max_bytes=MAX_DOWNLOAD):
    if not url.startswith(BASE+"/"):raise RuntimeError("No unsupported/untrusted URL fetch: "+url[:180])
    last=None
    for k in range(3):
        try:
            req=urllib.request.Request(url,headers={"Accept":"application/json,application/octet-stream,*/*","User-Agent":UA})
            with urllib.request.urlopen(req,timeout=75) as res:
                raw=res.read(max_bytes+1)
                if len(raw)>max_bytes:raise RuntimeError("Unexpected source file larger than limit: "+str(max_bytes))
                return raw
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
            last=e
            if isinstance(e,urllib.error.HTTPError) and e.code not in [429,500,502,503,504]:
                break
            time.sleep(2+k*3)
    raise RuntimeError(f"Dryad API unavailable for {url}: {type(last).__name__}:{last}")
def js(url):return json.loads(get(url,500_000).decode("utf-8"))
def file_list(obj):
    # Dryad HAL encodes arrays in _embedded, often "stash:files" and "stash:versions".
    if isinstance(obj,list):return [v for v in obj if isinstance(v,dict)]
    if not isinstance(obj,dict):return []
    e=obj.get("_embedded",{})
    out=[]
    for k,v in e.items() if isinstance(e,dict) else []:
        if isinstance(v,list):out+=file_list(v)
        if isinstance(v,dict):out+=file_list(v)
    if not out:
        for k,v in obj.items():
            if k not in ("_embedded",) and isinstance(v,list) and any(isinstance(q,dict) for q in v):
                out+=file_list(v)
    return out

def redacted_columns(fieldnames):
    # Clean already-published scientific column names, but no contact fields/records.
    return [str(v)[:110] for v in fieldnames[:120] if str(v).strip()]
def csv_member(raw,filename):
    stream=io.StringIO(raw.decode("utf-8-sig","replace"))
    reader=csv.reader(stream)
    try:columns=redacted_columns(next(reader))
    except StopIteration:columns=[]
    rows=0
    counts=collections.Counter()
    key_indices={k:i for i,k in enumerate(columns) if
                 any(a in k.lower() for a in ("site","river","stream","location","latitude","longitude","imperv","emerg","spider","isotope","year","month","date","taxon","length","biomass"))}
    unique={x:set() for x in key_indices if any(k in x.lower() for k in ("site","river","stream","year","month"))}
    empties=collections.Counter()
    for row in reader:
        if not row or not any(z.strip() for z in row):continue
        rows+=1
        for k,i in key_indices.items():
            if i>=len(row) or not row[i].strip():empties[k]+=1
            if k in unique and i<len(row) and row[i].strip() and len(unique[k])<=5000:
                unique[k].add(row[i].strip())
    return {"kind":"csv","file":filename,"n_data_rows":rows,
      "n_columns":len(columns),"column_names":columns,
      "candidate_matching_keys":[k for k in key_indices],
      "unique_values_for_site_river_year_keys":{k:len(v) for k,v in unique.items()},
      "blank_rows_by_candidate_key":dict(empties)}

def inspect_one(name,contents):
    if name.lower().endswith(".zip") or contents[:2]==b"PK":
        with zipfile.ZipFile(io.BytesIO(contents)) as z:
            names=z.namelist()
            if any(".." in Path(n).parts or n.startswith("/") for n in names):
                raise RuntimeError("Potentially unsafe published ZIP path; inspecting only")
            if name.lower().endswith(".xlsx"):
                # Workbook metadata only: sheet names, no raw tables or results.
                import xml.etree.ElementTree as ET
                wb=ET.fromstring(z.read("xl/workbook.xml"))
                ns={"x":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
                sheetnames=[x.attrib.get("name") for x in wb.findall(".//x:sheet",ns)]
                return {"kind":"published_xlsx_metadata","file":name,"sheet_names":sheetnames,
                        "workbook_file_members_n":len(names),
                        "full_data_rows_NOT_read":True}
            out=[]
            for child in names:
                if child.lower().endswith(".csv") and z.getinfo(child).file_size<=4_000_000:
                    out.append(csv_member(z.read(child),f"{name}/{child}"))
                elif child.lower().endswith((".txt",".md")) and z.getinfo(child).file_size<70_000:
                    # do not distribute full README contents; metadata only
                    out.append({"file":f"{name}/{child}","kind":"readme_text","bytes":z.getinfo(child).file_size})
                else:
                    out.append({"file":f"{name}/{child}","kind":"not_inspected",
                                "uncompressed_bytes":z.getinfo(child).file_size})
            return {"kind":"zip_package","file":name,
                "n_members":len(names),"source_file_member_metadata":out}
    if name.lower().endswith(".csv"):
        return csv_member(contents,name)
    if name.lower().endswith((".md",".txt")):
        return {"kind":"readme_text","file":name,"bytes":len(contents),
                "header_lines":[x.strip()[:120] for x in contents.decode("utf-8-sig","replace").splitlines()
                  if x.strip().startswith("#")][:25]}
    return {"kind":"unknown_format","file":name,"bytes":len(contents)}

def one_set(m):
    enc=urllib.parse.quote("doi:"+m["doi"],safe="")
    api=BASE+"/datasets/"+enc
    dataset=js(api)
    versions=js(api+"/versions")
    # Published REST provider versions change shape; log public HAL metadata
    # before attempting a version/file download. No ecology magnitudes.
    print("DRYAD_PUBLIC_VERSION_RESPONSE="+json.dumps({
        "dataset":m["name"],
        "version_response":versions,
        "dataset_api_links":dataset.get("_links",{}),
        "dataset_keys":list(dataset.keys())},ensure_ascii=False)[:7000],flush=True)
    candidates=file_list(versions)
    ranked=[v for v in candidates if isinstance(v.get("id"),int)]
    if not ranked:
        # Dryad version metadata may be linked from the dataset or latestVersion
        for k in ("latestVersion","versionNumber","id"):
            if isinstance(dataset.get(k),int):
                ranked=[{"id":dataset[k],"versionNumber":dataset.get("versionNumber",0)}]
                break
    if not ranked:raise RuntimeError("No accessible published version IDs in Dryad API")
    version=max(ranked,key=lambda d:float(d.get("versionNumber",0) or 0))
    versionid=int(version["id"])
    filemeta=js(BASE+f"/versions/{versionid}/files")
    files=file_list(filemeta)
    if not files:raise RuntimeError("Published dataset has no file records via version/files")
    outputs=[]
    for f in files:
        fileid=f.get("id")
        name=str(f.get("path") or f.get("filename") or f.get("name") or "")
        size=f.get("size") or f.get("fileSize") or f.get("sizeInBytes") or 0
        if not isinstance(fileid,int):
            outputs.append({"source_file_metadata_keys":sorted(f.keys()),"error":"missing published file numeric id"})
            continue
        rec={"file_id":fileid,"published_name":name,"published_size_bytes":size}
        try:
            raw=get(BASE+f"/files/{fileid}/download",MAX_DOWNLOAD)
            rec["actual_download_sha256"]=sha(raw)
            rec["actual_bytes"]=len(raw)
            rec["inspection"]=inspect_one(name,raw)
        except Exception as ex:
            rec["error_type"]=type(ex).__name__
            rec["error_summary"]=str(ex)[:320]
        outputs.append(rec)
    return {
      "dataset":m["name"],"source_doi":m["doi"],"citation":m["citation"],
      "dryad_dataset_title":dataset.get("title"),
      "published_dataset_version_id":versionid,
      "published_file_count":len(files),
      "source_version_metadata_sha256":sha(json.dumps(versions,sort_keys=True).encode()),
      "published_files":outputs}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--dataset",choices=("all","nash","quebec"),default="all")
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    sets=[x for x in SETS if args.dataset=="all" or (args.dataset=="nash")==(x["name"]=="nash_global_emergence_2023")]
    data=[]
    for s in sets:
        try:
            outcome=one_set(s)
        except Exception as e:
            outcome={"dataset":s["name"],"source_doi":s["doi"],"status":"SOURCE_ACCESS_STOP",
                     "error_type":type(e).__name__,"error_summary":str(e)[:450]}
        else:
            outcome["status"]="SOURCE_FILE_SCHEMA_AUDIT_COMPLETED_NO_WWTP_OR_ECOLOGY_TEST"
        data.append(outcome)
        print("DRYAD_EMERGENCE_SOURCE="+json.dumps(
           {"dataset":s["name"],"status":outcome["status"],
            "n_files":outcome.get("published_file_count"),
            "errors":sum(bool(f.get("error_type")) for f in outcome.get("published_files",[]))},
           ensure_ascii=False),flush=True)
    result={"schema":"GESUI_20261117_PUBLISHED_AQUATIC_EMERGENCE_RIPARIAN_SPIDER_SOURCE_ONLY_AUDIT_V1",
      "dataset_records":data,"source_identity_not_wwtp_exposure":True,
      "source_not_a_terrestrial_subsidy_fit":True,
      "do_not_treat_global_occurrences_as_emergence_biomass":True,
      "no_causal_sensitivity_or_ecological_effect_computed":True}
    (args.out/f"DRYAD_AQUATIC_EMERGENCE_SOURCE_METADATA_{args.dataset.upper()}.json").write_text(
       json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")

if __name__=="__main__":main()
