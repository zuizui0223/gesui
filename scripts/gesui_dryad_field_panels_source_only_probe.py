#!/usr/bin/env python3
"""SOURCE-ONLY intake of Dryad field panels.

No numerical ecological effects, model fit, significance, or result-informed
subset. Uses ONLY Dryad's original public API and the API-advertised file links.
All original file checksums and native schemas archived. Raw bytes in temporary
Actions artifact only (public CC0 data) to allow exact post-source replay.
"""
import csv,hashlib,io,json,os,pathlib,re,time,urllib.error,urllib.parse,urllib.request,zipfile,collections
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_dryad_field_panels_source"
ROOT.mkdir(exist_ok=True)
PANELS={"A_ON_29":"doi:10.5061/dryad.sbcc2frfk",
        "B_QC_20":"doi:10.5061/dryad.v6wwpzh8z"}
BASE="https://datadryad.org"
def http(url,tries=3):
    err=None
    for t in range(tries):
        try:
            req=urllib.request.Request(url,headers={
                "Accept":"application/json,text/plain,text/csv,application/octet-stream,*/*",
                "User-Agent":"AcademicFieldSourceProvenance/1.0 (read-only, public Dryad ecology dataset)"})
            with urllib.request.urlopen(req,timeout=75) as resp:
                return resp.read(),resp.geturl(),resp.headers.get("Content-Type","")
        except Exception as e:
            err=e;time.sleep(min(2+2*t,7))
    raise RuntimeError("Official Dryad source unavailable "+url+" "+repr(err))
def urlpath(obj, rel):
    links=obj.get("_links",{}) if isinstance(obj,dict) else {}
    for key,value in links.items():
        if rel.lower() in key.lower():
            if isinstance(value,dict) and value.get("href"):return urllib.parse.urljoin(BASE,value["href"])
            if isinstance(value,list) and value and isinstance(value[0],dict):return urllib.parse.urljoin(BASE,value[0]["href"])
    return None
def csv_outline(blob):
    try:
        s=blob.decode("utf-8-sig")
    except Exception:
        return {"status":"not_utf8_text"}
    try:
        reader=csv.DictReader(io.StringIO(s,newline=""))
        cols=reader.fieldnames or []
        if len(cols)<2:return {"status":"not_tabular_csv","head":s[:200]}
        candidate=[c for c in cols if re.search(r"(site|station|location|sampling|point|date|period|river|stream|visit|sinu|species|family|taxa|guild)",c,re.I)]
        present=collections.defaultdict(set)
        n=0
        for row in reader:
            n+=1
            for c in candidate:
                v=str(row.get(c,"")).strip()
                if len(present[c])<10000:present[c].add(v)
        return {"status":"table_read","n_data_rows":n,"columns":cols,
            "candidate_source_key_distinct_nonblank":{k:len({v for v in p if v}) for k,p in present.items()},
            "candidate_key_value_sample":{k:sorted({v for v in p if v})[:10] for k,p in present.items()
                    if len(p)<=60 and "date" not in k.lower()}}
    except Exception as e:return {"status":"parse_fail","error":repr(e)[:300]}
all_out={}
for label,doi in PANELS.items():
    encoded=urllib.parse.quote(doi,safe="")
    u=BASE+"/api/v2/datasets/"+encoded
    info={"dataset_DOI":doi,"dataset_api_url":u,"dataset_response":None,"files":[]}
    try:
        blob,actual,mime=http(u)
        dataset=json.loads(blob)
        info["dataset_response"]={"sha256":hashlib.sha256(blob).hexdigest(),
          "dataset_id":dataset.get("id"),"versionNumber":dataset.get("versionNumber"),
          "versionStatus":dataset.get("versionStatus"),"license":dataset.get("license"),
          "publicationDate":dataset.get("publicationDate"),"title":dataset.get("title"),
          "dataset_url":actual}
        vu=urlpath(dataset,"version")
        if not vu:
            vu=BASE+"/api/v2/datasets/"+encoded+"/versions"
        vb,va,vm=http(vu)
        v=json.loads(vb)
        # The linked version may be a single version or HAL collection.
        version=(v.get("_embedded",{}).get("stash:versions",[{}])[-1] if isinstance(v,dict) else {})
        if not version or not isinstance(version,dict):version=v
        info["version_request"]={"url":vu,"sha256":hashlib.sha256(vb).hexdigest(),
            "body_keys":list(v.keys()) if isinstance(v,dict) else [],
            "version_id":version.get("id"),"filelink":urlpath(version,"files")}
        fu=urlpath(version,"files")
        if not fu and version.get("id"):fu=BASE+"/api/v2/versions/"+str(version["id"])+"/files"
        if not fu:raise RuntimeError("Version file list URL absent; no guessing")
        fb,fa,fm=http(fu)
        fileobj=json.loads(fb)
        entries=fileobj.get("_embedded",{}).get("stash:files")
        if entries is None:
            entries=fileobj.get("_embedded",{}).get("files")
        if entries is None and isinstance(fileobj,list):entries=fileobj
        if entries is None:
            # Other HAL file relation names allowed, use only value-shaped arrays.
            for val in fileobj.get("_embedded",{}).values():
                if isinstance(val,list) and val and isinstance(val[0],dict) and "path" in val[0]:
                    entries=val;break
        if not isinstance(entries,list):
            raise RuntimeError("Dryad official file metadata unsupported, source keys "+str(list(fileobj.keys())))
        info["file_list"]={"url":fa,"sha256":hashlib.sha256(fb).hexdigest(),
            "n_listed":len(entries)}
        for ix,f in enumerate(entries):
            fname=f.get("path") or f.get("name") or f.get("filename") or "file_"+str(ix)
            link=urlpath(f,"download")
            if not link and f.get("id"):
                link=BASE+"/api/v2/files/"+str(f["id"])+"/download"
            if not link:raise RuntimeError("Missing official per-file download URL for "+str(fname))
            data,actual,mime=http(link)
            fn=(re.sub(r"[^A-Za-z0-9_.-]","_",fname))
            stored=ROOT/(label+"__"+fn)
            stored.write_bytes(data)
            desc={"name":fname,"file_id":f.get("id"),"file_metadata":{k:f.get(k) for k in ("path","size","digest","status","mimeType","description")},
                "official_link":link,"resolved_link":actual,
                "size":len(data),"sha256":hashlib.sha256(data).hexdigest(),
                "stored_as":stored.name,"mime":mime}
            if zipfile.is_zipfile(io.BytesIO(data)):
                with zipfile.ZipFile(io.BytesIO(data)) as z:
                    desc["zip_members"]=[]
                    for member in z.infolist():
                        if member.is_dir():continue
                        if member.file_size>150_000_000:
                            raise RuntimeError("Archive member unexpectedly large")
                        contents=z.read(member)
                        member_record={"name":member.filename,"bytes":len(contents),
                            "sha256":hashlib.sha256(contents).hexdigest()}
                        if member.filename.lower().endswith(".csv") and member.file_size<20_000_000:
                            member_record["source_schema"]=csv_outline(contents)
                        elif member.filename.lower().endswith((".md",".txt")) and member.file_size<100_000:
                            member_record["text_head"]=contents.decode("utf8","replace")[:900]
                        desc["zip_members"].append(member_record)
            elif fn.lower().endswith(".csv"):
                desc["source_schema"]=csv_outline(data)
            elif fn.lower().endswith((".md",".txt")):
                desc["text_head"]=data.decode("utf8","replace")[:2500]
            info["files"].append(desc)
    except Exception as e:
        info["status"]="SOURCE_API_OR_FORMAT_STOP"
        info["error"]=type(e).__name__+" "+str(e)[:1500]
    else:
        info["status"]="OFFICIAL_SOURCE_METADATA_AND_FILES_RETRIEVED"
    all_out[label]=info
    print("DRYAD_SOURCE_PANEL_RESULT_START "+label,flush=True)
    short={k:v for k,v in info.items() if k not in ("files",)}
    short["file_summaries"]=[{"name":f["name"],"bytes":f["size"],"sha256":f["sha256"],
       "members":[{"name":m["name"],"rows":m.get("source_schema",{}).get("n_data_rows"),
                   "cols":m.get("source_schema",{}).get("columns",[])[:28]} for m in f.get("zip_members",[])]}
                     for f in info["files"]]
    print(json.dumps(short,sort_keys=True,default=str),flush=True)
    print("DRYAD_SOURCE_PANEL_RESULT_END "+label,flush=True)
reportpath=ROOT/"DRYAD_RIPARIAN_FIELD_PANELS_SOURCE_ONLY_QA.json"
reportpath.write_text(json.dumps(all_out,indent=2,sort_keys=True))
print("Original Dryad field files are SOURCE only; no ecological effect estimated.",flush=True)
