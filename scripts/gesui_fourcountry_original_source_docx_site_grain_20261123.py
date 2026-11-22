#!/usr/bin/env python3
"""Original Loughborough four-country Word source-summaries sampling-grain audit.

EXACT source file MD5 and byte count, read-only OOXML. Export no full original
copyrighted source text. Only restrained metadata, headings, short method
snippets and TABLE STRUCTURE for site/replicate/temporal sampling audit.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,pathlib,re,zipfile,xml.etree.ElementTree as ET
SOURCE_MD5="08d81e4d97af34053b453a1c1e4f8c7d"
SOURCE_BYTES=279236
W="http://schemas.openxmlformats.org/wordprocessingml/2006/main"
R="http://schemas.openxmlformats.org/package/2006/relationships"
WORD="word/document.xml"
COUNTRIES={
 "Australia":[r"\baustralia\b",r"new south wales",r"\bNSW\b"],
 "Brazil":[r"\bbrazil\b",r"callisto"],
 "New Zealand":[r"new zealand",r"\bNZ\b",r"death",r"gres"],
 "UK":[r"united kingdom",r"\bUK\b",r"britain",r"environment agency"]}
TERMS=re.compile(
 r"(?:\bsite|\bsites|\brepeat|\bvisit|\byear|\bperiod|\bspatial|\btemporal|\bsample"
 r"|\bannual|\bmonthly|\btime\b|\bdate\b|\bkick\b|\bsurber\b|\bsource"
 r"|\bcode\b|\bidentifier|\bcatchment|\briver\b|\bstream|\bnumber\b|\bnumerical"
 r"|\bquantitative|\bregion|\bsubstrat|\bmonitor|\btributar)",
 re.I)
def sha(path,kind):
 h=hashlib.new(kind)
 with path.open("rb") as f:
  for b in iter(lambda:f.read(1048576),b""):h.update(b)
 return h.hexdigest()
def paragraph(p):
 return " ".join("".join(x.text or "" for x in p.iter(f"{{{W}}}t")).split())
def all_content(z):
 root=ET.fromstring(z.read(WORD))
 body=root.find(f"{{{W}}}body")
 paras=[];tables=[]
 for el in body:
  if el.tag==f"{{{W}}}p":
   t=paragraph(el)
   if not t:continue
   ppr=el.find(f"{{{W}}}pPr")
   style=None
   if ppr is not None:
    tag=ppr.find(f"{{{W}}}pStyle")
    if tag is not None:style=tag.attrib.get(f"{{{W}}}val")
   paras.append({"original_paragraph_index":len(paras)+1,
     "text":t,"style":style})
  elif el.tag==f"{{{W}}}tbl":
   vals=[]
   for tr in el.findall(f"{{{W}}}tr"):
    cells=[]
    for tc in tr.findall(f"{{{W}}}tc"):
     text=" ".join(paragraph(p) for p in tc.iter(f"{{{W}}}p") if paragraph(p))
     cells.append(text)
    vals.append(cells)
   tables.append(vals)
 return paras,tables
def method_snippets(paras):
 # Provide internal audit with limited original textual spans and paragraph IDs.
 # No long original source passages are exported or committed; for each country
 # capture at most 30 source-metadata excerpts capped at 160 chars each.
 output={}
 for country,patterns in COUNTRIES.items():
  candidates=[]
  for i,p in enumerate(paras):
   t=p["text"]
   if any(re.search(pattern,t,re.I) for pattern in patterns):
    around=paras[max(i-2,0):min(i+7,len(paras))]
    for candidate in around:
     if TERMS.search(candidate["text"]):
      candidates.append((candidate["original_paragraph_index"],candidate["text"]))
  out=[];seen=set()
  for index,text in candidates:
   if index in seen:continue
   seen.add(index)
   out.append({"paragraph_index":index,
       "method_excerpt_first_160_characters":text[:160],
       "original_characters":len(text)})
   if len(out)>=30:break
  output[country]=out
 return output
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--original-docx",type=pathlib.Path,required=True)
 ap.add_argument("--out",type=pathlib.Path,required=True)
 a=ap.parse_args()
 size=a.original_docx.stat().st_size
 if size!=SOURCE_BYTES or sha(a.original_docx,"md5")!=SOURCE_MD5:
  raise RuntimeError("Original DOI publisher DOCX exact bytes / MD5 FAIL")
 with zipfile.ZipFile(a.original_docx) as z:
  if WORD not in z.namelist():raise RuntimeError("Publisher DOCX missing OOXML document")
  paras,tables=all_content(z)
  links=[]
  if "word/_rels/document.xml.rels" in z.namelist():
   root=ET.fromstring(z.read("word/_rels/document.xml.rels"))
   links=[{"target":x.attrib.get("Target","")[:180],
           "type":x.attrib.get("Type","").split("/")[-1]}
          for x in root.findall(f"{{{R}}}Relationship")
          if x.attrib.get("TargetMode")=="External"]
 if not paras:raise RuntimeError("Publisher source DOCX has no paragraphs")
 counts={c:0 for c in COUNTRIES}
 for row in paras:
  for k,terms in COUNTRIES.items():
   if any(re.search(p,row["text"],re.I) for p in terms):counts[k]+=1
 data={"schema":"gesui_2023_figshare_original_fourcountry_document_sampling_source_QA_v1",
  "original_DOI":"10.17028/rd.lboro.24615570.v1",
  "original_docx_file_ID":43299384,
  "original_docx_basename":"MathersEtAl-2023-PairedVisualFineSediment-DataSourceSummaries.docx",
  "verified_publisher_original_MD5":SOURCE_MD5,
  "verified_original_source_bytes":SOURCE_BYTES,
  "original_source_SHA256":sha(a.original_docx,"sha256"),
  "original_document_structural_counts":{"paragraphs":len(paras),"tables":len(tables)},
  "original_source_paragraphs_mentioning_country":counts,
  "original_document_heading_labels":[
    {"paragraph_index":r["original_paragraph_index"],
     "style":r["style"],"heading_first_160_chars":r["text"][:160]}
    for r in paras if (r["style"] or "").lower().startswith(("heading","title"))][:80],
  "country_source_method_snippets":method_snippets(paras),
  "source_table_structures":[{"table_index":i+1,
    "rows":len(rows),"columns_first_row":len(rows[0]) if rows else 0,
    "top3_rows_first6_cells":[[v[:100] for v in row[:6]] for row in rows[:3]]}
    for i,rows in enumerate(tables)][:40],
  "original_external_docx_reference_links":links[:50],
  "full_original_publisher_text_exported":False,
  "source_site_grouping_verified":False,
  "source_summary_only_no_numerical_biological_model":True,
  "reason":"Manual source-based verification of site-grouping needed before group independent inference. These are traceable excerpts and original DOI metadata, not evidence of cross-country sample independence."}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(data,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
 print("FOURCOUNTRY_ORIGINAL_DOCX_SITE_GRAIN",json.dumps({
  "source_sha":data["original_source_SHA256"],
  "structure":data["original_document_structural_counts"],
  "headings":data["original_document_heading_labels"],
  "country_methods":data["country_source_method_snippets"],
  "tables":data["source_table_structures"],
  "references":links[:20]},ensure_ascii=False,sort_keys=True)[:25000],flush=True)
if __name__=="__main__":main()
