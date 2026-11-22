#!/usr/bin/env python3
"""Inspect ORIGINAL published Loughborough Figshare XLSX structure, source only.

Strict standard-library ZIP/XML metadata only: lists worksheet names, dimensions,
row/column counts and textual header cells; NEVER reads biological abundance,
sediment numerical values or performs any regression.
"""
import argparse,hashlib,json,os,pathlib,posixpath,re,xml.etree.ElementTree as ET,zipfile
ARTICLE=24615570
FILE_ID=43299387
EXPECTED_MD5="1e4d6e37b73c663406e5567dd79d9452"
EXPECTED_SIZE=4640542
MAIN="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_DOC="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL="http://schemas.openxmlformats.org/package/2006/relationships"
def as_hash(path,algo):
 h=hashlib.new(algo)
 with path.open("rb") as f:
  for buf in iter(lambda:f.read(1024*1024),b""):h.update(buf)
 return h.hexdigest()
def parse_shared(z):
 if "xl/sharedStrings.xml" not in z.namelist():return []
 root=ET.fromstring(z.read("xl/sharedStrings.xml"))
 out=[]
 for si in root.findall(f"{{{MAIN}}}si"):
  out.append("".join(t.text or "" for t in si.iter(f"{{{MAIN}}}t")))
 return out
def read_sheet_structure(z,path,strings):
 rows=0;colletters=set();headrows=[];dimension=None;numeric_cells=0;text_cells=0
 with z.open(path) as f:
  for e,node in ET.iterparse(f,events=("end",)):
   if node.tag==f"{{{MAIN}}}dimension":
    dimension=node.attrib.get("ref")
   if node.tag!=f"{{{MAIN}}}row":continue
   rows+=1
   if rows<=5:
    visible=[]
    for cell in node.findall(f"{{{MAIN}}}c"):
     ref=cell.attrib.get("r","")
     letter=re.match(r"[A-Z]+",ref)
     if letter:colletters.add(letter.group(0))
     t=cell.attrib.get("t","n")
     v=cell.find(f"{{{MAIN}}}v")
     if t=="s" and v is not None and v.text is not None:
      try:value=strings[int(v.text)]
      except (ValueError,IndexError):value="[unresolved text]"
      text_cells+=1
     elif t=="inlineStr":
      value="".join(x.text or "" for x in cell.iter(f"{{{MAIN}}}t"))
      text_cells+=1
     elif t in ("str","e") and v is not None:
      value=str(v.text or "") if t=="str" else "[cell error]"
      text_cells+=1
     else:
      value="[non-text value withheld]"
      numeric_cells+=1
     visible.append({"column":letter.group(0) if letter else ref,"text_or_placeholder":value[:180]})
    headrows.append({"original_excel_row":int(node.attrib.get("r",rows)),"textual_cells":visible[:170]})
   node.clear()
 return {"original_worksheet_internal_path":path,"worksheet_excel_dimension":dimension,
  "excel_physical_row_count":rows,
  "first_five_rows_text_only":headrows,"data_numeric_values_read":False}
def main():
 parser=argparse.ArgumentParser()
 parser.add_argument("--xlsx",required=True,type=pathlib.Path)
 parser.add_argument("--out",required=True,type=pathlib.Path)
 a=parser.parse_args()
 if not a.xlsx.exists():raise RuntimeError("Publisher original workbook unavailable")
 size=a.xlsx.stat().st_size
 if size!=EXPECTED_SIZE:raise RuntimeError(f"Publisher original byte size changed {size} != {EXPECTED_SIZE}")
 md5=as_hash(a.xlsx,"md5")
 if md5!=EXPECTED_MD5:raise RuntimeError("Original publisher MD5 mismatch")
 with zipfile.ZipFile(a.xlsx) as z:
  wb=ET.fromstring(z.read("xl/workbook.xml"))
  relationships=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
  targets={r.attrib["Id"]:r.attrib["Target"] for r in relationships.findall(f"{{{PKG_REL}}}Relationship")}
  ss=parse_shared(z)
  sheets=[]
  for sh in wb.findall(f".//{{{MAIN}}}sheets/{{{MAIN}}}sheet"):
   ident=sh.attrib[f"{{{REL_DOC}}}id"]
   target=targets[ident]
   target=target.lstrip("/") if target.startswith("/") else posixpath.normpath("xl/"+target)
   if target not in z.namelist():raise RuntimeError("Missing referenced original XLSX XML worksheet "+target)
   structure=read_sheet_structure(z,target,ss)
   structure["original_sheet_name"]=sh.attrib["name"]
   structure["original_sheet_state"]=sh.attrib.get("state","visible")
   sheets.append(structure)
 data={"schema":"gesui_McKenzie_2023_fourcountry_finesediment_XLSX_original_source_schema_v1",
  "publisher_data_DOI":"10.17028/rd.lboro.24615570.v1",
  "source_publication_DOI":"10.1111/gcb.17084",
  "figshare_article_id":ARTICLE,
  "file_id":FILE_ID,
  "original_source_file_name":"MathersEtAl-2023-PairedVisualFineSediment-Data.xlsx",
  "publisher_supplied_MD5":EXPECTED_MD5,"verified_actual_MD5":md5,
  "original_source_SHA256":as_hash(a.xlsx,"sha256"),
  "original_source_file_bytes":size,
  "original_source_worksheets":sheets,
  "worksheet_count":len(sheets),
  "original_quantitative_data_cells_accessed":False,
  "site_id_source_join_proven":False,
  "original_EPT_individual_fraction_proven":False,
  "holdout_model_GO":False,
  "next_gate":"Audit fields, taxon abundance and site sampling-unit identity using publisher source. XLSX source units are NOT EPA XEMBED."}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")
 print("EXTERNAL_4COUNTRY_XLSX_SOURCE_SCHEMA",json.dumps({
  "source_SHA256":data["original_source_SHA256"],
  "worksheets":[{"name":s["original_sheet_name"],"rows":s["excel_physical_row_count"],
                  "dimension":s["worksheet_excel_dimension"],
                  "first_text_row":s["first_five_rows_text_only"][0]["textual_cells"][:45] if s["first_five_rows_text_only"] else []}
                 for s in sheets]},sort_keys=True),flush=True)
if __name__=="__main__":main()
