#!/usr/bin/env python3
"""Source-only adult feeding trait acquisition; no environmental outcomes read."""
import argparse, collections, hashlib, io, json, pathlib, re, urllib.parse, zipfile
import requests
import openpyxl
URL="https://butterflytraits.github.io/European-Butterfly-Traits/European_%26_Maghreb_Butterfly_Trait_data_v1.2.zip"
DEF="https://butterflytraits.github.io/European-Butterfly-Traits/trait_variables_definitions_v1.2.pdf"
HOSTZIP="1cf6ed4ee41282a7f1be2e5cc7e1542a186b4fef5e1a46adf5861543c0ebbd94"
ALIAS_SHA="4370ff528290711d2083891e0a83c6ba4eb9bda0a9230a3829aae06779f530fa"
FAMILIES={"Hesperiidae","Papilionidae","Pieridae","Lycaenidae","Nymphalidae","Riodinidae"}
CODES={"Hf","Er","Sf","Hd","Sa","Dp","An","Mi"}
def sha(b):return hashlib.sha256(b).hexdigest()
def fetch(u,maxbytes=80*1024*1024):
    with requests.get(u,stream=True,timeout=(15,60),headers={"User-Agent":"Academic butterfly trait provenance audit"}) as r:
        r.raise_for_status()
        assert urllib.parse.urlparse(r.url).hostname=="butterflytraits.github.io",r.url
        b=bytearray()
        for c in r.iter_content(1<<16):
            b.extend(c);assert len(b)<=maxbytes
        return bytes(b),{"requested_url":u,"resolved_url":r.url,"bytes":len(b),"sha256":sha(b),"content_type":r.headers.get("Content-Type")}
def txt(v):return "" if v is None else str(v).strip()
def norm(v):return re.sub(r"[^a-z0-9]","",txt(v).lower())
def binomial(v):
    s=txt(v)
    if not re.fullmatch(r"[A-Za-z][a-zA-Z-]+\s+[a-z][a-z-]+",s):return None
    return " ".join(s.lower().split())
def source_roster(p):
    b=p.read_bytes();assert sha(b)==HOSTZIP
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        q=z.read("butterfly_usage_aliases.json");assert sha(q)==ALIAS_SHA
        r=json.loads(q)
    assert len(r)==303 and len({v["original_speciesKey"] for v in r})==303 and len({v["accepted_speciesKey"] for v in r})==303
    assert collections.Counter(v["status"] for v in r)=={"ACCEPTED":290,"DOUBTFUL":13}
    return r
def schema(wb):
    result=[]
    for s in wb:
        first=[]
        for i,row in enumerate(s.iter_rows(min_row=1,max_row=min(12,s.max_row),values_only=True),1):
            # Preserve header/schema only; not broad redistribution of the workbook.
            first.append({"row":i,"nonempty_cells":[[j+1,txt(v)[:160]] for j,v in enumerate(row) if v is not None][:150]})
        result.append({"sheet":s.title,"max_row":s.max_row,"max_column":s.max_column,"first_rows":first})
    return result
def main(host,out):
    out.mkdir(parents=True,exist_ok=True);r=source_roster(host)
    result={"status":"SOURCE_STOP","original_host_roster":len(r),"accepted_roster":290,"family_counts":dict(collections.Counter(v["family"] for v in r)),"no_environmental_or_butterfly_outcomes_read":True}
    try:
        pdf,meta=fetch(DEF);assert pdf.startswith(b"%PDF")
        (out/"official_trait_definitions.pdf").write_bytes(pdf);result["definitions"]=meta
        b,meta=fetch(URL);assert zipfile.is_zipfile(io.BytesIO(b))
        (out/"publisher_trait_v1.2.zip").write_bytes(b);result["publisher_zip"]=meta
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            names=z.namelist();assert len(names)==len(set(names))
            result["members"]=[dict(path=i.filename,bytes=i.file_size,sha256=sha(z.read(i))) for i in z.infolist() if not i.is_dir()]
            xlsx=[n for n in names if n.lower().endswith(".xlsx") and not pathlib.PurePosixPath(n).name.startswith("~$")]
            assert len(xlsx)>=1
            workbooks={}
            for n in xlsx:
                w=openpyxl.load_workbook(io.BytesIO(z.read(n)),read_only=True,data_only=True)
                workbooks[n]=schema(w);w.close()
            result["workbooks"]=workbooks
        result["status"]="PUBLISHER_BYTES_ACQUIRED_SCHEMA_REVIEW_REQUIRED"
    except Exception as e:result["error"]=type(e).__name__+": "+str(e)[:1500]
    (out/"source_receipt.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("ADULT TRAIT SOURCE "+json.dumps({k:v for k,v in result.items() if k!="workbooks"},ensure_ascii=False,allow_nan=False),flush=True)
    if result["status"]=="SOURCE_STOP":raise SystemExit(2)
if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--host",required=True,type=pathlib.Path);a.add_argument("--out",required=True,type=pathlib.Path);x=a.parse_args();main(x.host,x.out)
