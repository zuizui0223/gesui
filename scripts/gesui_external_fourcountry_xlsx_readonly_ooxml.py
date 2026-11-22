#!/usr/bin/env python3
"""Minimal stdlib OOXML row reader for original publisher XLSX.

Read-only: no modifications, no guessed schema, preserve original ANALYSIS_ID
as a normalized source string. Does not require external spreadsheet libraries.
"""
import math,posixpath,re,xml.etree.ElementTree as ET,zipfile
MAIN="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_DOC="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
REL_PKG="http://schemas.openxmlformats.org/package/2006/relationships"
def index(col):
 n=0
 for c in col:
  n=n*26+ord(c)-64
 return n-1
def shared(z):
 if "xl/sharedStrings.xml" not in z.namelist():return []
 root=ET.fromstring(z.read("xl/sharedStrings.xml"))
 return ["".join(x.text or "" for x in t.iter(f"{{{MAIN}}}t"))
         for t in root.findall(f"{{{MAIN}}}si")]
def opensheets(path):
 with zipfile.ZipFile(path) as z:
  root=ET.fromstring(z.read("xl/workbook.xml"))
  rels=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
  targets={e.attrib["Id"]:e.attrib["Target"] for e in rels.findall(f"{{{REL_PKG}}}Relationship")}
  sheets={}
  for e in root.findall(f".//{{{MAIN}}}sheets/{{{MAIN}}}sheet"):
   target=targets[e.attrib[f"{{{REL_DOC}}}id"]]
   target=target.lstrip("/") if target.startswith("/") else posixpath.normpath("xl/"+target)
   sheets[e.attrib["name"]]=target
  return sheets
def normalize_id(v):
 if v is None:return ""
 s=str(v).strip()
 # Excel numeric sample IDs may use integral XML values such as 1502.0.
 # Never cast through float (would destroy large identifiers).
 if re.fullmatch(r"\d+\.0+",s):s=s.split(".")[0]
 return s
def read_all(path,sheet):
 with zipfile.ZipFile(path) as z:
  worksheets=opensheets(path)
  if sheet not in worksheets:raise ValueError("Publisher source worksheet not present: "+sheet)
  ss=shared(z)
  records=[]
  with z.open(worksheets[sheet]) as stream:
   for ev,el in ET.iterparse(stream,events=("end",)):
    if el.tag!=f"{{{MAIN}}}row":continue
    values={}
    for cell in el.findall(f"{{{MAIN}}}c"):
     ref=cell.attrib.get("r","")
     m=re.match(r"([A-Z]+)",ref)
     if not m:continue
     ix=index(m.group(1))
     typ=cell.attrib.get("t","n")
     v=cell.find(f"{{{MAIN}}}v")
     if typ=="s" and v is not None and v.text is not None:
      try:out=ss[int(v.text)]
      except Exception:out=""
     elif typ=="inlineStr":
      out="".join(t.text or "" for t in cell.iter(f"{{{MAIN}}}t"))
     elif v is None:
      out=None
     else:
      out=v.text
     values[ix]=out
    records.append((int(el.attrib.get("r",len(records)+1)),values))
    el.clear()
 if not records:raise RuntimeError("Empty original worksheet")
 row,head=records.pop(0)
 if row!=1:raise RuntimeError("Original worksheet header row not row1")
 last=max(head)
 names=[str(head.get(j,"")).strip() for j in range(last+1)]
 if len(set(names))!=len(names) or not names or names[0]!="ANALYSIS_ID":
  raise RuntimeError("Publisher original analysis ID header invalid or repeated columns")
 out=[]
 for number,record in records:
  vals=[record.get(j) for j in range(len(names))]
  out.append(vals)
 return names,out
def as_float(v):
 try:
  if v is None or str(v).strip()=="":return math.nan
  return float(v)
 except (ValueError,TypeError):return math.nan
