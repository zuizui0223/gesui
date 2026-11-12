#!/usr/bin/env python3
"""Read ONLY official public source links. No ecology estimates.

Checks Dryad's documented public *full-dataset* GET route, Zenodo original
record(s) explicitly crosslinked by published Dryad landing page, and Borealis
Dataverse's public search API. No credentials, CAPTCHA, cookies, scraped
session, or download of unpublished/unauthorized material.
"""
import csv,hashlib,io,json,os,pathlib,re,time,urllib.parse,urllib.request,zipfile,tarfile
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_wwtp_consumer_public_source_gate"
ROOT.mkdir(parents=True,exist_ok=True)
orig={
 "A":{
  "doi":"doi:10.5061/dryad.sbcc2frfk","versionid":342267,
  "original_expected_sha256":{
   "sites_full.csv":"3b4b38303fc3e70dcd00411742a20da8d77147cccc1ac56239736d9249e8afef",
   "sites_bird.csv":"b4dc6834eb0b2d6a1906adf2a74dd1d8180765378e11f3e967302f271e4a5849",
   "emergence.csv":"e36990e83a8a5f79ad06ca097911e849e2e6f3414b030d20753cf4376ac80dd4",
   "emergence_seasonality.csv":"719552e0c0967f7e60db939fccb68a62d9cfff2d3f7daf2cff2",
   "bird.csv":"5486ebd00866a6295cf4041f7dbaa873406c9c8daa0cc03ca8fb882065e1679f",
   "environment.csv":"17c21017c5e385ecbd9e8fef89325854da533b1fbca2b754b1434f1b361d4bee",
   "benthic.csv":"e9c5bfbcd1078ca4fdde5a1c5859692a7e67f5a4631ad1126033ec769b81c74b",
   "benthic_classifications.csv":"e273087edbce5c8f6aca595899369346139dca2632924ff09b5745e2ab3c74c8",
   "README.md":"049eea49b0ded5ef613675d1c9daef925a558b5009a6df6d94f4b4a1295407d0"
  }
 },
 "B":{
  "doi":"doi:10.5061/dryad.v6wwpzh8z","versionid":430321,
  "original_expected_sha256":{
   "File_3.zip":"fafa089c757115a83aceb8d3f254f2b4898416e2b470f5a2de8b07a82c8880ef",
   "README.md":"b3d4eb6e5dc260b6a7fe82156f7c43cadb01e81bd9dd61ac35d5273d88f9572d"
  }
 }
}
sha=lambda b:hashlib.sha256(b).hexdigest()
def fetch(url,timeout=45):
    req=urllib.request.Request(url,headers={
        "User-Agent":"gesui-public-research-data-source-check/1.0",
        "Accept":"application/json,application/zip,text/plain,application/octet-stream,*/*"})
    with urllib.request.urlopen(req,timeout=timeout) as r:
        b=r.read(40_000_000)
        if len(b)>=40_000_000:raise RuntimeError("Unexpected large original binary")
        return b,r.url,r.headers.get("Content-Type","")
def limited_fetch(url):
    try:return {"data":fetch(url),"error":None}
    except Exception as e:return {"data":None,"error":type(e).__name__+" "+str(e)[:240]}
def inspect_contents(blob):
    """Only file names and source columns and row-count. Never response slopes."""
    files={}
    if zipfile.is_zipfile(io.BytesIO(blob)):
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            for m in z.infolist():
                if m.is_dir() or m.file_size>15_000_000:continue
                files[m.filename]=z.read(m)
    elif tarfile.is_tarfile(fileobj:=None) if False else False:
        pass
    else:
        # Catch tar.gz or tar through read-only BytesIO
        try:
            with tarfile.open(fileobj=io.BytesIO(blob),mode="r:*") as t:
                for m in t.getmembers():
                    if not m.isfile() or m.size>15_000_000:continue
                    files[m.name]=t.extractfile(m).read()
        except (tarfile.TarError,EOFError):
            pass
    return files
def schema(b,filename):
    if filename.lower().endswith(".csv"):
        try:
            text=b.decode("utf-8-sig")
            reader=csv.DictReader(io.StringIO(text,newline=""))
            head=reader.fieldnames or []
            if len(head)>1000:raise RuntimeError("Unexpected column count")
            n=sum(1 for _ in reader)
            return {"nrows":n,"columns":head}
        except Exception as e:return {"parse_error":repr(e)[:150]}
    if filename.lower().endswith((".md",".txt")):
        # Metadata only, no user-authored site outcomes.
        return {"first_600_readme_chars":b.decode("utf8","replace")[:600]}
    return {}
out={}
for label,rec in orig.items():
    doi=rec["doi"];enc=urllib.parse.quote(doi,safe="")
    url=f"https://datadryad.org/api/v2/datasets/{enc}/download"
    q=limited_fetch(url)
    o={"source":"Dryad official documented full-dataset public download API",
       "canonical_DOI":doi,"original_versionid_metadata":rec["versionid"],
       "original_public_route":url,
       "status":"not_verified","download_status":q["error"],"files":{}}
    if q["data"] is not None:
        b,actual,mime=q["data"]
        o["archive_sha256"]=sha(b);o["archive_bytes"]=len(b)
        o["resolved_url"]=actual;o["mime"]=mime
        files=inspect_contents(b)
        o["archive_members_n"]=len(files)
        for path,bb in files.items():
            fn=path.split("/")[-1]
            expected=rec["original_expected_sha256"].get(fn)
            verified=bool(expected and sha(bb)==expected)
            meta={"original_path":path,"bytes":len(bb),"sha256":sha(bb),
                  "expected_Dryad_original_sha256":expected,
                  "exact_same_original":verified,
                  "source_structure":schema(bb,fn)}
            o["files"][fn]=meta
            if verified:
                f=ROOT/(label+"__"+fn)
                f.write_bytes(bb)
        o["status"]="ORIGINAL_ALL_FILES_SHA256_VERIFIED" if len(o["files"])==len(rec["original_expected_sha256"]) and all(x["exact_same_original"] for x in o["files"].values()) else "PARTIAL_OR_DIFFERENT_VERSION_DO_NOT_ASSUME_MATCH"
    out["Dryad_"+label]=o
    print("DRYAD_CANONICAL_PUBLIC_SOURCE",label,json.dumps({k:v for k,v in o.items() if k!="files"}),flush=True)
    print("DRYAD_FILE_SUMMARY",label,json.dumps([{k:row.get(k) for k in ["original_path","bytes","exact_same_original","source_structure"]} for row in o["files"].values()]),flush=True)

# Author/record provenance via *only* public Zenodo entry specifically linked by Dryad.
for num in ["11186734","11186733"]:
    api="https://zenodo.org/api/records/"+num
    z=limited_fetch(api)
    meta={"public_record_url":api}
    if z["data"] is None:meta["status"]="not_accessible";meta["error"]=z["error"]
    else:
        raw,actual,mime=z["data"]
        j=json.loads(raw)
        fs=j.get("files",{})
        if isinstance(fs,list):entries=fs
        else:entries=fs.get("entries",{}) if isinstance(fs,dict) else {}
        if isinstance(entries,dict):entries=list(entries.values())
        meta.update({"status":"PUBLIC_METADATA","record_title":j.get("metadata",{}).get("title"),
           "record_doi":j.get("metadata",{}).get("doi"),
           "license":j.get("metadata",{}).get("license"),
           "files": [{"name":f.get("key"),"size":f.get("size"),"checksum":f.get("checksum"),
                      "public_link":f.get("links",{}).get("self"),
                      "is_published":f.get("is_head")} for f in entries]})
        # Only inspect original plausible data files, no invented mirror: link specifically
        # registered on the original Dryad landing page as work referencing it.
        available={}
        for f in entries:
            fname=f.get("key","").split("/")[-1]
            original_sha=orig["A"]["original_expected_sha256"].get(fname)
            if not original_sha:continue
            fu=f.get("links",{}).get("self")
            if not fu:continue
            got=limited_fetch(fu)
            if got["data"] is None:
                available[fname]={"error":got["error"]};continue
            original=got["data"][0]
            checked=sha(original)==original_sha
            available[fname]={"bytes":len(original),"original_sha256_match":checked,
                    "sha256":sha(original),"columns":schema(original,fname)}
            if checked:
                (ROOT/("Zenodo_verified_A__"+fname)).write_bytes(original)
        meta["comparison_to_original_Dryad"] = available
    out["Zenodo_"+num]=meta
    print("ZENODO_SOURCE",num,json.dumps(meta,sort_keys=True),flush=True)

# Public Borealis metadata search only (standard Dataverse open search API).
for term in ['Wardlaw', '"Cross-ecosystem transfer"','"microplastics"']:
    api="https://borealisdata.ca/api/search?"+urllib.parse.urlencode({"q":term,"type":"dataset","per_page":30})
    got=limited_fetch(api)
    k="BorealisSearch_"+term.replace(" ","_")
    if got["data"] is None:val={"error":got["error"]}
    else:
        try:
            x=json.loads(got["data"][0])
            data=x.get("data",{})
            vals=data.get("items",[]) if isinstance(data,dict) else []
            val={"total_count":data.get("total_count"),"items":[{
                 "title":y.get("name"),"datasetDOI":y.get("global_id"),
                 "url":y.get("url"),"id":y.get("entity_id")} for y in vals[:30]],
                 "response_sha256":sha(got["data"][0])}
        except Exception as e:val={"parse_error":str(e)[:150],"head":got["data"][0][:300].decode("utf8","replace")}
    out[k]=val
    print("BOREALIS_METADATA",k,json.dumps(val,sort_keys=True),flush=True)
dest=ROOT/"ORIGINAL_DRYAD_ZENODO_BOREALIS_SOURCE_CHECK_20261112.json"
dest.write_text(json.dumps(out,indent=2,sort_keys=True))
