#!/usr/bin/env python3
"""Resolve publisher's exact three declared workbook attachments, not narrative SI numbers."""
import sys,pathlib,io,json,hashlib,zipfile,urllib.parse,re
import requests,openpyxl
from bs4 import BeautifulSoup
a=pathlib.Path(sys.argv[1]);out=pathlib.Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
b=a.read_bytes();assert sha(b)=="509bded2f4568a838ec2127f6c87731bf9ad170cb6916528ed8b6e85a8edd337"
z=zipfile.ZipFile(io.BytesIO(b));html=z.read("landing_3819283111be.html");assert sha(html)=="3819283111be98ad3cfc1f65a46141cde7eecfff08889664cad3788e64ff387f"
links={}
for x in BeautifulSoup(html,"html.parser").find_all("a",href=True):
 u=x["href"];m=re.search(r"S0266467426100510sup(00[456])\.xlsx$",u)
 if m:
  assert urllib.parse.urlparse(u).scheme=="https" and urllib.parse.urlparse(u).hostname=="static.cambridge.org"
  links[m.group(1)]=u
assert set(links)=={"004","005","006"}
r={"status":"SOURCE_STOP","source_archive_sha256":sha(b),"files":{}}
try:
 for k,u in sorted(links.items()):
  with requests.get(u,timeout=(15,60),stream=True,allow_redirects=False) as resp:
   resp.raise_for_status();assert resp.status_code==200
   v=bytearray()
   for c in resp.iter_content(65536):
    v.extend(c);assert len(v)<10*1024**2
  v=bytes(v);p=out/(k+".xlsx");p.write_bytes(v)
  w=openpyxl.load_workbook(io.BytesIO(v),read_only=True,data_only=True)
  r["files"][k]={"url":u,"sha256":sha(v),"bytes":len(v),"sheets":[{"sheet":s.title,"rows":s.max_row,"cols":s.max_column,"first_rows":[list(map(lambda q:None if q is None else str(q),row)) for row in s.iter_rows(min_row=1,max_row=3,values_only=True)]} for s in w]}
  w.close()
 r["status"]="SOURCE_TABLES_ACQUIRED_SCHEMA_REVIEW"
finally:
 (out/"receipt.json").write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
 print("BAIT TABLES "+json.dumps(r,ensure_ascii=False),flush=True)
assert r["status"]=="SOURCE_TABLES_ACQUIRED_SCHEMA_REVIEW"
