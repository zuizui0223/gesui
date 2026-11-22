#!/usr/bin/env python3
"""Targeted original publisher sampling-design assertions with provenance receipts.

Study-year substrata, NZ study labels, UK spring/autumn repeat and sewage
screening are author-documented; sample-id prefix grouping is ONLY a proxy.
Source-only; never fits sediment vs EPT or infers independent physical sites.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,pathlib,re,statistics,zipfile
import gesui_fourcountry_original_source_docx_site_grain_20261123 as doc
import gesui_external_fourcountry_xlsx_readonly_ooxml as oo
from gesui_external_fourcountry_source_grain_20261122 import ORIGINAL_SHA,ORIGINAL_MD5,digest

def source_paragraphs(docx):
 with zipfile.ZipFile(docx) as z:
  paras,tables=doc.all_content(z)
 return {p["original_paragraph_index"]:p["text"] for p in paras},tables
def require(source,num,pattern,label):
 text=source.get(num,"")
 hit=re.search(pattern,text,re.I)
 if hit is None:raise RuntimeError("Original publisher source text does not support "+label)
 return {"source_original_paragraph":num,"original_concept":label,
   "matched_source_keyword_count":len(list(re.finditer(pattern,text,re.I))),
   "verified":True,
   "selected_very_short_source_context":text[max(0,hit.start()-20):min(len(text),hit.end()+45)][:100]}
def NZ_patterns(xlsx):
 header,rows=oo.read_all(xlsx,"NZ taxa")
 if len(rows)!=342:raise RuntimeError("Original NZ source n mismatch")
 codes=[oo.normalize_id(r[0]) for r in rows]
 source_labels=collections.defaultdict(list)
 for label in codes:
  parts=label.split(".")
  source_labels[parts[0]].append(parts)
 out={}
 for label,records in sorted(source_labels.items()):
  n=len(records)
  ndot=collections.Counter(len(x)-1 for x in records)
  firsttwo=collections.Counter(".".join(x[:2]) for x in records)
  duplicated={k:v for k,v in firsttwo.items() if v>=2}
  out[label]={"source_observation_rows":n,
    "source_id_dot_segments_distribution":dict(sorted(ndot.items())),
    "source_prefix_first_two_tokens_unique":len(firsttwo),
    "source_rows_sharing_repeated_two_token_prefix":sum(duplicated.values()),
    "number_repeated_two_token_groups":len(duplicated),
    "maximum_repeats_one_prefix":max(firsttwo.values()),
    "original_ID_examples_capped":[".".join(x)[:65] for x in records[:3]],
    "site_identity_proven":False}
 return out
def summary_count(table):
 for tbl in table:
  if not tbl or tbl[0][0].strip()!="Source reference":continue
  data=[]
  for row in tbl[1:]:
   s=row[0].strip()
   match=re.search(r"\bn\s*=\s*(\d+)",s,re.I)
   data.append({"source_name":s[:80],"publisher_original_summary_n":int(match.group(1)) if match else None,
     "original_sampling_calendar":row[1][:45] if len(row)>1 else ""})
  return data
 raise RuntimeError("No original NZ methods source summary table")
def main():
 ap=argparse.ArgumentParser()
 for k in ("original-docx","original-xlsx","out"):ap.add_argument("--"+k,type=pathlib.Path,required=True)
 a=ap.parse_args()
 if doc.sha(a.original_docx,"md5")!=doc.SOURCE_MD5 or a.original_docx.stat().st_size!=doc.SOURCE_BYTES:
  raise RuntimeError("Original four-country publisher docx checksum differs")
 if digest(a.original_xlsx,"sha256")!=ORIGINAL_SHA or digest(a.original_xlsx,"md5")!=ORIGINAL_MD5:
  raise RuntimeError("Original four-country publisher xlsx checksum differs")
 source,tables=source_paragraphs(a.original_docx)
 checks={
  "UK_environment_agency_sewage_discharge_downstream_exclusion":
   require(source,14,r"not downstream of discharge from a sewage treat", "EA source discharge downstream exclusion"),
  "UK_academic_source_both_spring_and_autumn_sampling":
   require(source,16,r"(?:all sites in the UK were sampled in spring and autumn)", "UK academic spring and autumn site sampling"),
  "UK_environment_agency_spring_or_autumn_screen":
   require(source,15,r"spring\s*\(March\s*-\s*May\)\s*or\s*autumn", "EA seasonal source filtering"),
  "Australia_reference_site_AUS_RIVAS_filter":
   require(source,6,r"Only sites considered to pass the AusRivAS assessment", "AU reference-condition sampling site filter"),
  "NZ_multiple_different_academic_studies":
   require(source,11,r"aggregate of several academic studies", "NZ multi-study amalgam source"),
  "NZ_experimental_before_only_filter":
   require(source,11,r"only the [‘']before[’'] data was retained", "NZ before-only original source selection")}
 prefixes=NZ_patterns(a.original_xlsx)
 source_table=summary_count(tables)
 if sum(v["source_observation_rows"] for v in prefixes.values())!=342:raise RuntimeError("Source NZ original ID counts drift")
 output={"schema":"gesui_fourcountry_original_source_design_sampling_EXACT_RELEVANT_SIGNALS_v1",
  "source_original_docx_SHA256":doc.sha(a.original_docx,"sha256"),
  "source_original_xlsx_SHA256":ORIGINAL_SHA,
  "author_documented_source_filters":checks,
  "new_zealand_publisher_method_table":source_table,
  "new_zealand_XLSX_first_token_ID_patterns":prefixes,
  "new_zealand_XLSX_sample_n":342,
  "new_zealand_publisher_original_table_n":sum(x["publisher_original_summary_n"] or 0 for x in source_table),
  "original_source_method_table_counts_vs_release_workbook_not_identical":True,
  "no_independent_physical_site_ID_in_original_publisher_XLSX":True,
  "no_case_specific_impact_of_sewage_can_be_claimed":True,
  "UK_environment_agency_source_not_downstream_sewage_screen":True,
  "UK_academic_programme_must_not_be_assumed_same_as_environment_agency":True,
  "NZ_prefix_site_key_proven":False,
  "site_independent_uncertainty_inference_GO":False,
  "no_taxon_outcomes_EPT_sediment_model_fitted":True}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(output,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
 print("FOURCOUNTRY_SOURCE_DESIGN_EXPLICIT_FLAGS",json.dumps({
  "documented_filters":checks,
  "NZ_sample_ID_patterns":prefixes,
  "NZ_doc_source_original_n":output["new_zealand_publisher_original_table_n"],
  "NZ_final_workbook_source_n":342,
  "source_site_independent_go":False},ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":main()
