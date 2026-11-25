#!/usr/bin/env python3
"""Pinned public NutNet Figshare v4 acquisition; source schema only, no ecology model."""
import argparse,csv,hashlib,json,pathlib,urllib.parse,urllib.request,zipfile,io,gzip
ARTICLE=26412295
VERSION=4
DOI="10.6084/m9.figshare.26412295.v4"
LIMIT=300*1024**2
API=f"https://api.figshare.com/v2/articles/{ARTICLE}/versions/{VERSION}"
UA="gesui academic open data verification/1.0"
def check(u,initial=False):
 p=urllib.parse.urlsplit(u);h=p.hostname or ""
 fg=h=="figshare.com" or h.endswith(".figshare.com") or h=="figshareusercontent.com" or h.endswith(".figshareusercontent.com")
 cloud=h.endswith(".amazonaws.com") and (h.startswith("s3") or ".s3" in h)
 if p.scheme!="https" or p.username or p.password or p.port not in (None,443) or not (fg or (cloud and not initial)):
  raise ValueError("unexpected source host "+h)
 return u
class Redirect(urllib.request.HTTPRedirectHandler):
 def redirect_request(self,req,fp,code,msg,headers,newurl):
  check(newurl)
  return super().redirect_request(req,fp,code,msg,headers,newurl)
op=urllib.request.build_opener(Redirect())
def retrieve(u,p,maxbytes):
 h=hashlib.sha256();m=hashlib.md5();n=0
 try:
  with op.open(urllib.request.Request(check(u,initial=True),headers={"User-Agent":UA,"Accept-Encoding":"identity"}),timeout=90) as r,p.open("wb") as f:
   check(r.geturl())
   while True:
    b=r.read(1024**2)
    if not b:break
    n+=len(b)
    if n>maxbytes:raise ValueError("size limit")
    f.write(b);h.update(b);m.update(b)
   return {"bytes":n,"sha256":h.hexdigest(),"md5":m.hexdigest(),"host":urllib.parse.urlsplit(r.geturl()).hostname}
 except BaseException:
  p.unlink(missing_ok=True);raise
def inspect_csv(raw):
 reader=csv.reader(io.TextIOWrapper(raw,encoding="utf-8-sig",newline=""))
 h=next(reader,[])
 if len(h)==1 and "\t" in h[0]:
  return {"header_preview":h[0][:2000],"schema":"TSV_REVIEW_REQUIRED"}
 n=0;bad=0
 for x in reader:
  n+=1;bad+=len(x)!=len(h)
 return {"columns":h,"rows":n,"width_mismatch":bad}
def inspect(p):
 if zipfile.is_zipfile(p):
  with zipfile.ZipFile(p) as z:
   assert len(z.namelist())==len(set(z.namelist()))
   assert sum(x.file_size for x in z.infolist())<3*1024**3
   rows=[]
   for a in z.infolist():
    if a.is_dir():continue
    rec={"member":a.filename,"bytes":a.file_size}
    with z.open(a) as f:
     h=hashlib.sha256()
     for chunk in iter(lambda:f.read(1024**2),b""):h.update(chunk)
     rec["sha256"]=h.hexdigest()
    if a.filename.lower().endswith(".csv"):
     with z.open(a) as f:rec["csv"]=inspect_csv(f)
    rows.append(rec)
   return rows
 if p.suffix.lower()==".csv":
  with p.open("rb") as f:return inspect_csv(f)
 return {"schema_review_required":True,"extension":p.suffix}
def main(out):
 out.mkdir(parents=True,exist_ok=True);r={"status":"SOURCE_STOP","ecological_model_fitted":False,"article_id":ARTICLE,"version":VERSION,"files":[]}
 try:
  meta=out/"article.json";r["metadata"]=retrieve(API,meta,8*1024**2)
  j=json.loads(meta.read_text());assert int(j["id"])==ARTICLE and int(j["version"])==VERSION and j["doi"].lower()==DOI
  r["doi"]=j["doi"];r["title"]=j["title"];r["license"]=j.get("license");files=j["files"];assert files and sum(x["size"] for x in files)<=LIMIT
  for x in files:
   n=x["name"];assert pathlib.Path(n).name==n and n not in {"article.json","source_receipt.json"}
   p=out/n;receipt=retrieve(x["download_url"],p,LIMIT)
   assert receipt["bytes"]==x["size"]
   md=x.get("computed_md5") or x.get("supplied_md5");assert md and receipt["md5"]==md.lower()
   r["files"].append({"id":x["id"],"name":n,**receipt,"inventory":inspect(p)})
  r["status"]="ACQUIRED_SCHEMA_REVIEW_REQUIRED"
 except Exception as e:
  r["error"]=type(e).__name__+": "+str(e)[:300]
 finally:
  (out/"source_receipt.json").write_text(json.dumps(r,indent=2)+"\n");print(json.dumps(r))
 if r["status"]=="SOURCE_STOP":raise SystemExit(1)
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);main(p.parse_args().out)
