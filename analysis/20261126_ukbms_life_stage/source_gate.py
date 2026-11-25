#!/usr/bin/env python3
"""UKBMS EIDC official source acquisition and schema census; no population-index analysis."""
import argparse, collections, csv, hashlib, io, json, pathlib, re, urllib.request, urllib.parse, zipfile
from html.parser import HTMLParser
SOURCES = {
    "site_indices_2024": "53a1c705-4437-48bb-9bcd-b68977f8edd9",
    "site_locations_2024": "d7256b49-f2e3-4ae7-907d-af729610c768",
}
HOSTS = {"catalogue.ceh.ac.uk", "data-package.ceh.ac.uk"}
LIMIT = 256 * 1024**2
UA = "gesui UKBMS academic source census/1.0"
class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.links=[]
    def handle_starttag(self, tag, attrs):
        if tag == "a":
            u=dict(attrs).get("href")
            if u: self.links.append(u)
def checked(u):
    p=urllib.parse.urlsplit(u)
    if p.scheme!="https" or p.hostname not in HOSTS:
        raise ValueError("Unexpected source host")
    return u
class Redirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        checked(newurl)
        return super().redirect_request(req,fp,code,msg,headers,newurl)
OPENER=urllib.request.build_opener(Redirect())
def retrieve(u, p, limit):
    checked(u); h=hashlib.sha256(); n=0
    try:
        with OPENER.open(urllib.request.Request(u,headers={"User-Agent":UA,"Accept-Encoding":"identity"}),timeout=90) as r:
            checked(r.geturl())
            with p.open("wb") as f:
                while True:
                    b=r.read(1024*1024)
                    if not b: break
                    n+=len(b)
                    if n>limit: raise RuntimeError("Source size limit")
                    h.update(b); f.write(b)
            return {"url":u,"final_url":r.geturl(),"http":r.status,"bytes":n,"sha256":h.hexdigest(),"content_type":r.headers.get("Content-Type")}
    except BaseException:
        p.unlink(missing_ok=True); raise
def inspect(p):
    rows=[]
    with zipfile.ZipFile(p) as z:
        assert len(z.namelist()) == len(set(z.namelist()))
        assert sum(x.file_size for x in z.infolist()) <= 2*1024**3
        for i in z.infolist():
            if i.is_dir(): continue
            h=hashlib.sha256()
            with z.open(i) as f:
                for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
            r={"member":i.filename,"bytes":i.file_size,"sha256":h.hexdigest()}
            if i.filename.lower().endswith((".csv",".tsv",".txt")):
                with z.open(i) as raw:
                    # Determine encoding from a bounded prefix; don't silently replace bytes.
                    prefix=raw.read(65536)
                encoding="utf-8-sig"
                try: prefix.decode(encoding)
                except UnicodeDecodeError: encoding="cp1252"
                with z.open(i) as raw, io.TextIOWrapper(raw,encoding=encoding,newline="") as f:
                    header=f.readline()
                    if not header: r["empty"]=True
                    else:
                        delim="\t" if header.count("\t")>header.count(",") else ","
                        if header.count(";")>max(header.count(","),header.count("\t")): delim=";"
                        cols=next(csv.reader([header],delimiter=delim))
                        if len(cols)>1:
                            reader=csv.reader(f,delimiter=delim)
                            ct=collections.Counter(); n=0; bad=0
                            for line in reader:
                                n+=1
                                if len(line)!=len(cols):bad+=1
                            r.update(columns=cols,encoding=encoding,delimiter=delim,data_rows=n,width_mismatch=bad)
                        else: r["tabular_schema_unconfirmed"]=True
            rows.append(r)
    return rows
def main(out):
    out.mkdir(parents=True,exist_ok=True); receipt={"status":"SOURCE_STOP","population_index_values_analyzed":False,"sources":{}}
    try:
        for name,uuid in SOURCES.items():
            item={"doi":"10.5285/"+uuid};receipt["sources"][name]=item
            try:
                u="https://catalogue.ceh.ac.uk/documents/"+uuid
                page=out/(name+"_metadata.html")
                item["metadata"]=retrieve(u,page,16*1024**2)
                parser=Links();parser.feed(page.read_text(encoding="utf-8",errors="strict"))
                links=[]
                for link in parser.links:
                    full=urllib.parse.urljoin(u,link);sp=urllib.parse.urlsplit(full)
                    if sp.scheme!="https" or sp.hostname not in HOSTS:continue
                    if uuid in full and ("/download/" in sp.path or "/data/" in sp.path):
                        links.append(full)
                item["official_data_links"]=list(dict.fromkeys(links))
                if not links:raise ValueError("No UUID-bound official download link")
                # The official landing-page link is controlling; only its own CEH redirect target is allowed.
                p=out/(name+".zip"); errs=[]
                for link in item["official_data_links"][:4]:
                    try:
                        item["download"]=retrieve(link,p,LIMIT)
                        if not zipfile.is_zipfile(p):raise ValueError("Download is not a ZIP archive")
                        item["members"]=inspect(p);break
                    except Exception as e:
                        errs.append(type(e).__name__+": "+str(e)[:300]);p.unlink(missing_ok=True)
                else:raise RuntimeError(str(errs))
                item["status"]="ACQUIRED_SCHEMA_REVIEW_REQUIRED"
            except Exception as e:
                item.update(status="SOURCE_STOP",error=type(e).__name__+": "+str(e)[:500])
        if all(x.get("status")=="ACQUIRED_SCHEMA_REVIEW_REQUIRED" for x in receipt["sources"].values()):
            receipt["status"]="ACQUIRED_SCHEMA_REVIEW_REQUIRED"
    finally:
        (out/"source_receipt.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n")
        print("UKBMS SOURCE "+json.dumps(receipt,ensure_ascii=False))
    if receipt["status"]=="SOURCE_STOP":raise SystemExit(1)
if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--out",type=pathlib.Path,required=True);main(a.parse_args().out)
