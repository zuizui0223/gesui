#!/usr/bin/env python3
"""Acquire exact publisher ant interaction supplement and expert-validated range source.
No butterfly occurrence or environmental response evaluated.
"""
import argparse,csv,hashlib,io,json,pathlib,re,urllib.parse,zipfile
import requests
import openpyxl
from bs4 import BeautifulSoup
FDOI="10.3897/nl.44.68993.suppl1"
RDOI="10.5061/dryad.wm37pvmq1"
HOSTS={"doi.org","api.crossref.org","api.datacite.org","nl.pensoft.net","pensoft.net","www.pensoft.net","public.pensoft.net","datadryad.org","www.datadryad.org","dryad-assets.s3.amazonaws.com","dryad-data.s3.amazonaws.com","api.gbif.org"}
def allowed(u):
 p=urllib.parse.urlparse(u);return p.scheme=="https" and (p.hostname in HOSTS or bool(p.hostname and p.hostname.endswith(".pensoft.net")) or bool(p.hostname and p.hostname.endswith(".datadryad.org")))
def get(u,limit=100*1024**2):
 assert allowed(u),u
 # Validate every redirect rather than fetching an unreviewed external host.
 original=u
 for _ in range(10):
  with requests.get(u,timeout=(15,90),stream=True,allow_redirects=False,headers={"User-Agent":"Academic butterfly-ant data provenance audit","Accept":"*/*"}) as r:
   if 300<=r.status_code<400:
    u=urllib.parse.urljoin(u,r.headers["Location"]);assert allowed(u),("unreviewed redirect",u);continue
   r.raise_for_status();b=bytearray()
   for c in r.iter_content(1<<16):
    b.extend(c);assert len(b)<=limit,"size limit"
   v=bytes(b)
   return v,{"requested":original,"resolved":u,"sha256":hashlib.sha256(v).hexdigest(),"bytes":len(v),"content_type":r.headers.get("Content-Type"),"content_disposition":r.headers.get("Content-Disposition")}
 raise RuntimeError("redirect loop")
def workbook(b):
 if not zipfile.is_zipfile(io.BytesIO(b)):return None
 with zipfile.ZipFile(io.BytesIO(b)) as z:
  if not {"[Content_Types].xml","xl/workbook.xml"}<=set(z.namelist()):return None
 w=openpyxl.load_workbook(io.BytesIO(b),read_only=True,data_only=True);result=[]
 for s in w:
  rows=list(s.iter_rows(min_row=1,max_row=min(s.max_row,8),values_only=True))
  result.append({"sheet":s.title,"rows":s.max_row,"columns":s.max_column,"first_rows":[[None if v is None else str(v)[:300] for v in a] for a in rows]})
 w.close();return result
def fiedler(out):
 errors=[];attempts=["https://doi.org/"+FDOI]
 try:
  b,info=get("https://api.crossref.org/works/"+FDOI,5*1024**2)
  data=json.loads(b)["message"];assert data["DOI"].lower()==FDOI.lower()
  (out/"fiedler_crossref.json").write_bytes(b)
  for u in [data.get("URL"),*[x.get("URL") for x in data.get("link",[])],data.get("resource",{}).get("primary",{}).get("URL")]:
   if u and allowed(u) and u not in attempts:attempts.append(u)
 except Exception as e:errors.append({"crossref":str(e)[:800]})
 for u in attempts:
  try:
   b,meta=get(u);w=workbook(b)
   if w is not None:
    (out/"fiedler_suppl1.xlsx").write_bytes(b);return {"status":"SOURCE_ACQUIRED_SCHEMA_REVIEW","file":meta,"workbook":w,"errors":errors}
   if b.startswith(b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"):
    (out/"fiedler_suppl1.xls").write_bytes(b);return {"status":"SOURCE_ACQUIRED_LEGACY_SCHEMA_REVIEW","file":meta,"errors":errors}
   # If DOI resolves to a landing page, use only an unambiguous link for the same supplement.
   soup=BeautifulSoup(b,"html.parser");candidates=[]
   for a in soup.find_all("a",href=True):
    h=urllib.parse.urljoin(meta["resolved"],a["href"]);label=" ".join(a.stripped_strings).lower()
    if not allowed(h):continue
    if ("suppl1" in h.lower() or "appendix s1" in label or h.lower().endswith((".xlsx",".xls"))) and ("68993" in h or "appendix s1" in label or "suppl1" in h):candidates.append(h)
   candidates=sorted(set(candidates));errors.append({"landing":meta,"links":candidates})
   if len(candidates)==1 and candidates[0] not in attempts:attempts.append(candidates[0])
  except Exception as e:errors.append({"url":u,"error":type(e).__name__+": "+str(e)[:800]})
 return {"status":"SOURCE_STOP","errors":errors}
def dryad(out):
 rec={"status":"SOURCE_STOP","errors":[]}
 root="https://datadryad.org/api/v2/datasets/"+urllib.parse.quote("doi:"+RDOI,safe="")
 b,meta=get(root,5*1024**2);d=json.loads(b);assert d.get("identifier","").lower()==("doi:"+RDOI).lower(),d.get("identifier")
 (out/"dryad_dataset.json").write_bytes(b);rec["dataset"]=meta;rec["title"]=d.get("title")
 links=d.get("_links",{});v=links.get("stash:version",{}).get("href");assert v,v
 u=urllib.parse.urljoin(root,v);b,meta=get(u,5*1024**2);ver=json.loads(b);(out/"dryad_version.json").write_bytes(b);rec["version"]=meta
 file_href=ver.get("_links",{}).get("stash:files",{}).get("href")
 assert file_href,"no exact version files link"
 u=urllib.parse.urljoin(u,file_href);b,meta=get(u,5*1024**2);files=json.loads(b);(out/"dryad_files.json").write_bytes(b)
 rows=files.get("_embedded",{}).get("stash:files",[]);assert rows
 expected={"1._Regional_species_lists_of_European_ants.csv","2._Binary_ranges_maps_of_European__ants.zip","README.txt"}
 rec["files_manifest"]=meta;rec["available_paths"]=[x.get("path") for x in rows];found={}
 for f in rows:
  path=f.get("path");name=pathlib.PurePosixPath(path).name
  if name not in expected:continue
  assert name not in found
  h=f.get("_links",{}).get("stash:download",{}).get("href")
  assert h,("no declared download link",name)
  b,inf=get(urllib.parse.urljoin(u,h),32*1024**2);assert len(b)==int(f["size"]),name
  if f.get("digest") and f.get("digestType"):
   alg=f["digestType"].lower().replace("-","");assert hashlib.new(alg,b).hexdigest().lower()==f["digest"].lower()
  if name.endswith(".zip"):
   assert zipfile.is_zipfile(io.BytesIO(b))
   with zipfile.ZipFile(io.BytesIO(b)) as z:
    items=z.infolist();inf["member_count"]=len(items);inf["first_members"]=[{"path":q.filename,"bytes":q.file_size} for q in items[:30]]
  elif name.endswith(".csv"):
   inf["first_rows"]=b.decode("utf-8-sig",errors="strict").splitlines()[:4]
  else:inf["readme"]=b.decode("utf-8-sig",errors="strict")[:12000]
  (out/name).write_bytes(b);found[name]=inf
 rec["received"]=found;assert set(found)==expected,(found.keys(),expected)
 rec["status"]="SOURCE_ACQUIRED_SCHEMA_REVIEW";return rec
def main(out):
 out.mkdir(parents=True,exist_ok=True);r={"no_outcomes_evaluated":True}
 try:r["Fiedler2021"]=fiedler(out)
 except Exception as e:r["Fiedler2021"]={"status":"SOURCE_STOP","error":type(e).__name__+": "+str(e)[:1000]}
 try:r["EuropeanAntRanges"]=dryad(out)
 except Exception as e:r["EuropeanAntRanges"]={"status":"SOURCE_STOP","error":type(e).__name__+": "+str(e)[:1000]}
 (out/"receipt.json").write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
 print("ANT SOURCE "+json.dumps(r,ensure_ascii=False),flush=True)
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);main(p.parse_args().out)
