#!/usr/bin/env python3
"""Natural England LTMN public-source acquisition and workbook schema inventory.
No ecological outcomes fitted. Only official publication links, no auth bypass."""
import argparse, hashlib, json, pathlib, urllib.parse, urllib.request, re, io
from html.parser import HTMLParser
import openpyxl, xlrd
HOST="publications.naturalengland.org.uk"
PAGES={
 "availability":"4967983926214656","butterfly":"5692468885979136",
 "vegetation":"6331137473904640","soil":"5143675349237760"
}
# Official file hrefs confirmed in the publication pages, never guessed.
FILES={
 "butterfly_index":("butterfly","5074328274862080",".xls"),
 "butterfly_transects":("butterfly","6033934228652032",".xlsx"),
 "butterfly_visits":("butterfly","6020843017338880",".xls"),
 "plant_survey_plots":("vegetation","5611258838253568",".xlsx"),
 "plant_lumped_species":("vegetation","6565859343204352",".xlsx"),
 "plant_plot_locations":("vegetation","4506313758081024",".xlsx"),
 "soil_chemistry":("soil","5956833842561024",".xlsx"),
 "soil_plot_locations":("soil","4920861424091136",".xlsx")
}
UA="gesui academic ecological source audit/1.0"
MAX=64*1024**2
def checked(u, *, redirect=False):
 p=urllib.parse.urlsplit(u)
 allowed={HOST,"storage.googleapis.com"} if redirect else {HOST}
 if p.scheme!="https" or p.hostname not in allowed or p.username or p.password or p.port not in (None,443):
  raise ValueError("unapproved source host: "+str(p.hostname))
 return u
class Redirect(urllib.request.HTTPRedirectHandler):
 def redirect_request(self, req, fp, code, msg, headers, newurl):
  checked(newurl,redirect=True)
  old=urllib.parse.urlsplit(req.full_url); new=urllib.parse.urlsplit(newurl)
  if new.hostname=="storage.googleapis.com" and not (old.hostname==HOST and old.path.startswith("/file/") or old.hostname==new.hostname):
   raise ValueError("storage redirect is not from an official attachment")
  return super().redirect_request(req,fp,code,msg,headers,newurl)
OP=urllib.request.build_opener(Redirect())
class Links(HTMLParser):
 def __init__(self):super().__init__();self.links=[]
 def handle_starttag(self,tag,attrs):
  if tag=="a":
   u=dict(attrs).get("href")
   if u:self.links.append(u)
def download(u,p,limit):
 h=hashlib.sha256();n=0
 try:
  with OP.open(urllib.request.Request(checked(u),headers={"User-Agent":UA,"Accept-Encoding":"identity"}),timeout=90) as r,p.open("wb") as f:
   checked(r.geturl(),redirect=True)
   while True:
    b=r.read(1024**2)
    if not b:break
    n+=len(b)
    if n>limit:raise ValueError("file size limit")
    f.write(b);h.update(b)
   return {"url":u,"final_url":urllib.parse.urlunsplit(urllib.parse.urlsplit(r.geturl())._replace(query="",fragment="")),"bytes":n,"sha256":h.hexdigest(),"content_type":r.headers.get("Content-Type")}
 except BaseException:
  p.unlink(missing_ok=True);raise
def summarize(p):
 with p.open("rb") as f:magic=f.read(8)
 if p.suffix==".xlsx":
  if magic[:4]!=b"PK\x03\x04":raise ValueError("not OOXML workbook")
  w=openpyxl.load_workbook(p,read_only=True,data_only=True)
  try:return [{"sheet":s.title,"declared_rows":s.max_row,"declared_columns":s.max_column,
    "header_strings":[str(v)[:120] for v in next(s.iter_rows(min_row=1,max_row=1,values_only=True),()) if isinstance(v,str)]}
    for s in w.worksheets]
  finally:w.close()
 if magic!=bytes.fromhex("d0cf11e0a1b11ae1"):raise ValueError("not binary XLS workbook")
 w=xlrd.open_workbook(p,on_demand=True)
 try:return [{"sheet":s.name,"rows":s.nrows,"columns":s.ncols,
   "header_strings":[str(v)[:120] for v in s.row_values(0) if isinstance(v,str)] if s.nrows else []}
   for s in (w.sheet_by_index(i) for i in range(w.nsheets))]
 finally:w.release_resources()
def main(out):
 out.mkdir(parents=True,exist_ok=True);r={"status":"SOURCE_STOP","population_or_plant_values_analyzed":False,"pages":{},"files":{}}
 links={}
 try:
  for k,id in PAGES.items():
   u=f"https://{HOST}/publication/{id}";p=out/(k+"_publication.html")
   try:
    meta=download(u,p,8*1024**2);lp=Links();lp.feed(p.read_text(encoding="utf-8",errors="strict"))
    links[k]={urllib.parse.urljoin(u,v) for v in lp.links}
    r["pages"][k]=meta
   except Exception as e:r["pages"][k]={"status":"SOURCE_STOP","error":type(e).__name__+": "+str(e)[:400]}
  for k,(page,id,ext) in FILES.items():
   u=f"https://{HOST}/file/{id}";p=out/(k+ext)
   try:
    if u not in links.get(page,set()):raise ValueError("exact href absent from retrieved official publication")
    meta=download(u,p,MAX)
    r["files"][k]={"status":"SCHEMA_REVIEW_REQUIRED","download":meta,"sheets":summarize(p)}
   except Exception as e:
    p.unlink(missing_ok=True);r["files"][k]={"status":"SOURCE_STOP","error":type(e).__name__+": "+str(e)[:400]}
  if all(v["status"]=="SCHEMA_REVIEW_REQUIRED" for v in r["files"].values()):r["status"]="ACQUIRED_SCHEMA_REVIEW_REQUIRED"
  elif any(v["status"]=="SCHEMA_REVIEW_REQUIRED" for v in r["files"].values()):r["status"]="PARTIAL_ACQUISITION"
 finally:
  (out/"source_receipt.json").write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
  print(json.dumps(r,ensure_ascii=False))
 if r["status"]!="ACQUIRED_SCHEMA_REVIEW_REQUIRED":raise SystemExit(1)
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=pathlib.Path);main(p.parse_args().out)
