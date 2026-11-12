#!/usr/bin/env python3
"""SOURCE-ONLY 2026-11-12 exact Wardlaw et al original Borealis Dataverse.

Identify full public DOI and original version/file entities; inspect native
headers, location/year/type fields without opening empirical outcome contrasts.
No invented wastewater sources, no mirrors, no credential bypass.
"""
import collections,csv,hashlib,io,json,os,pathlib,re,time,urllib.parse,urllib.request,zipfile
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_borealis_wardlaw_original_source"
ROOT.mkdir(parents=True,exist_ok=True)
DOI="doi:10.5683/SP4/U0MN9S"
BASE="https://borealisdata.ca"
def fetch(url,retries=3):
    last=None
    for i in range(retries):
        try:
            req=urllib.request.Request(url,headers={"Accept":"application/json,text/plain,text/csv,application/octet-stream,*/*",
                        "User-Agent":"AcademicBorealisFieldSourceAudit/1.0 (published publicly accessible research, no credentials)"})
            with urllib.request.urlopen(req,timeout=75) as r:
                b=r.read(60_000_000)
                if len(b)>=60_000_000:raise RuntimeError("Source unexpectedly large: audit required")
                return b,r.url,r.headers.get("Content-Type","")
        except Exception as e:
            last=e;time.sleep(min(2+2*i,8))
    raise RuntimeError("Official Borealis public data unavailable "+url+" "+str(last))
sha=lambda b:hashlib.sha256(b).hexdigest()
q=urllib.parse.urlencode({"persistentId":DOI})
url=BASE+"/api/datasets/:persistentId/?"+q
blob,final,mime=fetch(url)
j=json.loads(blob)
if j.get("status")!="OK":raise RuntimeError("Borealis original native dataset lookup not OK")
d=j.get("data",{})
meta={
 "source_DOI":DOI,
 "api_url":url,"source_lookup_response_sha256":sha(blob),
 "dataset_id":d.get("id"),"persistentId":d.get("persistentId"),
 "globalId":d.get("globalId"),"publisher":"Borealis Canadian Dataverse",
 "dataset_citation":d.get("citation"),
 "source_subject":None,"source_year_method":"observed 2022–2023 upstream/downstream of 8 Grand River WWTP from paper",
 "license":None,"files":[],
 "biological_outcome_coefficients_opened":False
}
if meta["persistentId"] and meta["persistentId"]!=DOI:raise RuntimeError("Dataset DOI mismatch")
latest=d.get("latestVersion",{})
meta["version"]={"number":latest.get("versionNumber"),"minor":latest.get("versionMinorNumber"),
       "state":latest.get("versionState"),"createTime":latest.get("createTime"),
       "releaseTime":latest.get("releaseTime")}
meta["license"]=latest.get("license",{}).get("name") if isinstance(latest.get("license"),dict) else latest.get("license")
meta["termsOfUse"]=latest.get("termsOfUse")
blocks=latest.get("metadataBlocks",{})
citation=blocks.get("citation",{}).get("fields",[])
fields={x.get("typeName"):x.get("value") for x in citation}
for k in ("title","author","dsDescription","subject","publication","dateOfDeposit"):
    v=fields.get(k)
    if v is not None:
        txt=json.dumps(v,ensure_ascii=False,default=str)
        meta[k]=txt[:4500]
files=latest.get("files") or []
print("Borealis original published DOI",DOI,"native dataset id",meta["dataset_id"],"raw files",len(files),flush=True)
for f in files:
    df=f.get("dataFile",{})
    fileid=df.get("id")
    filename=df.get("filename") or f.get("label")
    size=df.get("filesize")
    contentType=df.get("contentType")
    entry={"source_file_id":fileid,"filename":filename,"source_size":size,"content_type":contentType,
           "source_storage_md5":df.get("md5"),"version":df.get("version"),
           "restricted":f.get("restricted"),"description":f.get("description"),
           "download_status":"not_attempted"}
    if f.get("restricted"):
        entry["download_status"]="restricted_BY_publisher_STOP"
        meta["files"].append(entry);continue
    if not fileid:
        entry["download_status"]="no_source_id";meta["files"].append(entry);continue
    dl=BASE+"/api/access/datafile/"+str(fileid)
    try:
        raw,resolved,ct=fetch(dl)
        entry["download_status"]="public_original_bytes_accessible"
        entry["download_link"]=dl;entry["resolved_link"]=resolved
        entry["bytes"]=len(raw);entry["sha256"]=sha(raw)
        if df.get("md5"):entry["md5_verified"]=hashlib.md5(raw).hexdigest()==df["md5"]
        name=re.sub("[^A-Za-z0-9_.-]","_",filename)
        path=ROOT/("Borealis"+str(fileid)+"__"+name)
        path.write_bytes(raw)
        entry["stored_temp_original"]=path.name
        try:
            decoded=raw.decode("utf-8-sig")
            if filename.lower().endswith((".csv",".tsv",".txt")):
                separator="\t" if filename.lower().endswith(".tsv") else ","
                rdr=csv.DictReader(io.StringIO(decoded,newline=""),delimiter=separator)
                head=rdr.fieldnames or []
                candidate=[h for h in head if any(c in h.lower() for c in
                  ["site","station","treatment","wwtp","plant","stream","upstream","downstream","region","locat","year","date","sample","taxon","family","type","insect","spider"])]
                count=0
                distinct={k:set() for k in candidate}
                for row in rdr:
                    count+=1
                    for k in distinct:
                        if len(distinct[k])<2000:
                            distinct[k].add(str(row.get(k) or "").strip())
                entry["raw_table_header"]=head
                entry["raw_table_rows"]=count
                entry["native_key_field_cardinality"]={k:len([x for x in v if x]) for k,v in distinct.items()}
                entry["native_key_field_examples"]={k:sorted([x for x in v if x])[:10] for k,v in distinct.items()
                       if len(v)<=40}
        except UnicodeError:pass
        if zipfile.is_zipfile(io.BytesIO(raw)):
            entry["zip_members"]=[]
            with zipfile.ZipFile(io.BytesIO(raw)) as z:
                for x in z.infolist():
                    if x.is_dir():continue
                    entry["zip_members"].append({"name":x.filename,"size":x.file_size})
        # XLSX structural sheet/counts only, no column values besides first row.
        if filename.lower().endswith(".xlsx"):
            try:
                import openpyxl
                wb=openpyxl.load_workbook(path,read_only=True,data_only=True)
                entry["xlsx_sheets"]=[]
                for sheet in wb:
                    head=[str(v) if v is not None else "" for v in next(sheet.iter_rows(min_row=1,max_row=1,values_only=True))]
                    entry["xlsx_sheets"].append({"sheet":sheet.title,
                        "rowmax":sheet.max_row,"colmax":sheet.max_column,
                        "headers":head[:120]})
            except Exception as e:entry["xlsx_structure_error"]=str(e)[:250]
    except Exception as e:
        entry["download_status"]="no_public_binary_access_STOP"
        entry["error"]=type(e).__name__+": "+str(e)[:240]
    meta["files"].append(entry)
    summary={k:v for k,v in entry.items() if k not in ("raw_table_header",)}
    print("BOREALIS_ORIGINAL_FILE_SOURCE",json.dumps(summary,sort_keys=True,default=str),flush=True)
out=ROOT/"BOREALIS_WARDLAW_2026_8WWTP_SOURCE_ONLY_QA.json"
out.write_text(json.dumps(meta,indent=2,sort_keys=True,ensure_ascii=False))
print("BOREALIS_WARDLAW_SOURCE_GATE_START")
print(json.dumps({
 "doi":DOI,"dataset_id":meta["dataset_id"],"version":meta["version"],
 "license":meta["license"],"n_files":len(meta["files"]),
 "n_downloaded_original":sum(f["download_status"]=="public_original_bytes_accessible" for f in meta["files"]),
 "n_restricted":sum(f.get("restricted",False) for f in meta["files"]),
 "file_summaries":[{"filename":f["filename"],"bytes":f.get("bytes"),
   "sha256":f.get("sha256"),"md5_verified":f.get("md5_verified"),
   "original_rows":f.get("raw_table_rows"),
   "columns":f.get("raw_table_header"),
   "xlsx_sheets":f.get("xlsx_sheets"),"status":f.get("download_status")} for f in meta["files"]]
},sort_keys=True,default=str))
print("BOREALIS_WARDLAW_SOURCE_GATE_END")
