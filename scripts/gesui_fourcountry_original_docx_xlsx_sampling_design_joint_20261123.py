#!/usr/bin/env python3
"""Original 4-country publisher METHODS vs source XLSX sample-code grain audit.

Read-only original source DOCX and XLSX, exact MD5/SHA source checks. Probes
physical site/repeat/stratification evidence without fitting any ecology.
Never asserts that NZ dotted code tokens are proven site IDs or dates.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,pathlib,re,statistics,zipfile
import gesui_fourcountry_original_source_docx_site_grain_20261123 as design
import gesui_external_fourcountry_xlsx_readonly_ooxml as xlsx
from gesui_external_fourcountry_source_grain_20261122 import SHEETS,ORIGINAL_SHA,ORIGINAL_MD5,digest

METHOD_PARA={
 "Australia":[5,6],
 "Brazil":[8,9],
 "New Zealand":[11,12],
 "UK":[14,15,16]}
SCHEMA_EXACT_ORIGINAL={"Australia":"Australia taxa","Brazil":"Brazil taxa",
                       "New Zealand":"NZ taxa","UK":"UK taxa"}
PATTERNS=[
 (r"(?:\brepeat\w*|\bmore than once\b|\bmultiple times\b|\bseveral years\b)","repeats"),
 (r"(?:\bsite[s]?\b|\blocation[s]?\b|\breach(?:es)?\b)","sites_reaches"),
 (r"(?:\bspring\b|\bautumn\b|\byear\w*\b|\bmonth\w*\b|\bsampling date\b)","dates_seasons"),
 (r"(?:\bfix\w* count\b|\bkick\b|\bsurber\b|\bmesh\b|\breplicat\w*\b)","sample_methods"),
 (r"(?:\bmonitor\w*\b|\bresearch project\w*\b|\bstud\w*\b|\bsource\b)","source_programme")]
def by_sentence(text):
 # Paragraphs in publisher appendix include abbreviations, so split
 # conservatively; excerpts <= 90 source chars per matched phrase.
 sents=re.split(r"(?<=[.;])\s+(?=[A-Z0-9])",text)
 return sents if len(sents)>1 else re.split(r"\s+,\s+",text)
def matched_excerpt(text,regex,span=95,maxhit=4):
 out=[]
 for match in re.finditer(regex,text,re.I):
  start=max(0,match.start()-35);end=min(len(text),match.end()+span-35)
  val=" ".join(text[start:end].split())
  if val not in out:out.append(val[:span])
  if len(out)>=maxhit:break
 return out
def paragraph_design(paras):
 lookup={p["original_paragraph_index"]:p["text"] for p in paras}
 out={}
 for country,ids in METHOD_PARA.items():
  snippets=[]
  for i in ids:
   txt=lookup.get(i,"")
   patterns={}
   for regex,label in PATTERNS:
    hits=matched_excerpt(txt,regex)
    if hits:patterns[label]=hits
   snippets.append({"source_paragraph_number":i,
      "source_paragraph_length":len(txt),"keyword_occurrence_by_method":{
        label:len(list(re.finditer(regex,txt,re.I))) for regex,label in PATTERNS},
      "short_original_method_context":patterns})
  out[country]=snippets
 return out
def id_structure(name,rows):
 ids=[xlsx.normalize_id(r[0]) for r in rows]
 if len(set(ids))!=len(ids):
  raise RuntimeError("Original publisher country ANALYSIS_ID not one-to-one "+name)
 numeric=sum(bool(re.fullmatch(r"\d+",k)) for k in ids)
 dotted=sum("." in k for k in ids)
 results={"original_n":len(ids),"unique_original_ANALYSIS_ID":len(set(ids)),
    "numeric_only_ANALYSIS_ID":numeric,"dotted_ANALYSIS_ID":dotted,
    "full_original_SITE_ID_column_present":False}
 if name=="New Zealand":
  grouped=collections.defaultdict(list)
  first=collections.Counter()
  tails_are_number=0
  for ident in ids:
   parts=ident.split(".")
   first[parts[0]]+=1
   if len(parts)>=3:
    prefix=".".join(parts[:2])
    grouped[prefix].append(parts[2:])
    tails_are_number+=int(bool(re.fullmatch(r"\d{4,6}",parts[-1])))
   else:
    grouped[parts[0]].append(parts[1:])
  gs=collections.Counter({key:len(group) for key,group in grouped.items()})
  results["NZ_source_prefix_first_token_counter"]=dict(sorted(first.items()))
  results["NZ_two_token_prefix_stats"]={
   "unique_prefix_groups":len(gs),
   "prefix_groups_with_2plus_rows":sum(c>=2 for c in gs.values()),
   "source_rows_in_repeated_prefix_groups":sum(c for c in gs.values() if c>=2),
   "median_rows_per_prefix_group":statistics.median(gs.values()),
   "max_rows_per_prefix_group":max(gs.values()),
   "sample_IDs_with_numeric_4_to6digit_suffix":tails_are_number,
   "unverified_full_site_identifier_status":"PROXY ONLY — no original docx explicit ANALYSIS_ID field decoding"}
  results["NZ_example_prefixes_capped"]=sorted(gs)[:8]
 else:
  results["independent_site_ID_status"]="STOP; original publisher source sample ID numeric sequential; no site metadata"
 return results
def tables_for_study_sources(tables):
 # Only source name/brief source sample-year heading; not reproduce table content
 if not tables:return {}
 selected=[]
 for tbl in tables:
  if not tbl:continue
  if "Source reference" not in tbl[0]:continue
  for row in tbl[1:]:
   ref=row[0].strip() if row else ""
   date=row[1].strip() if len(row)>1 else ""
   n=re.search(r"\bn\s*=\s*(\d+)",ref,re.I)
   selected.append({
    "original_reference_capped":ref[:95],
    "source_calendar_interval_capped":date[:58],
    "original_source_n":int(n.group(1)) if n else None,
    "taxon_count_protocol_differs_within_country":True})
 return {"NZ_original_study_subsources":selected,
  "sum_original_study_n":sum(x["original_source_n"] or 0 for x in selected)}
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--original-docx",type=pathlib.Path,required=True)
 p.add_argument("--original-xlsx",type=pathlib.Path,required=True)
 p.add_argument("--out",type=pathlib.Path,required=True)
 args=p.parse_args()
 if args.original_docx.stat().st_size!=design.SOURCE_BYTES or design.sha(args.original_docx,"md5")!=design.SOURCE_MD5:
  raise RuntimeError("Original Loughborough source report DOCX changed")
 if digest(args.original_xlsx,"sha256")!=ORIGINAL_SHA or digest(args.original_xlsx,"md5")!=ORIGINAL_MD5:
  raise RuntimeError("Original Loughborough original XLSX checksum changed")
 with zipfile.ZipFile(args.original_docx) as z:
  paras,tables=design.all_content(z)
 country_details=paragraph_design(paras)
 sampleids={}
 for name,(taxsheet,sed,expected) in SHEETS.items():
  header,rows=xlsx.read_all(args.original_xlsx,taxsheet)
  if len(rows)!=expected or header[0]!="ANALYSIS_ID":
   raise RuntimeError("Publisher original source taxon worksheet grain changed")
  sampleids[name]=id_structure(name,rows)
 nz_table=tables_for_study_sources(tables)
 output={"schema":"gesui_original_Loughborough_4country_document_xlsx_sample_unit_joint_QA_v1",
   "original_DOCX_SHA256":design.sha(args.original_docx,"sha256"),
   "original_XLSX_SHA256":ORIGINAL_SHA,
   "original_DOI":"10.17028/rd.lboro.24615570.v1",
   "country_source_method_provenance":country_details,
   "NZ_original_subsource_table_summary":nz_table,
   "original_ANALYSIS_ID_structure_by_country":sampleids,
   "zero_analyses_of_biological_sediment_associations":True,
   "independent_physical_SITE_ID_verified_any_country":False,
   "cohort_original_sample_n_not_physical_independent_site_n":6491,
   "NZ_dotted_prefix_grouping_for_information_not_confirmed_physical_site":True,
   "source_gate_site_independent_inference":"STOP",
   "reason":"Publisher workbook retains sample-level ANALYSIS_ID but not physical site ID; appendix mixes monitoring programmes, study years, taxon counts and samples. NZ dotted prefix may imply repeats but is not a validated sampling site key."}
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(output,sort_keys=True,indent=2)+"\n")
 print("FOURCOUNTRY_ORIGINAL_DOCX_XLSX_DESIGN_GRAIN",json.dumps({
  "source_docs":nz_table,"NZ_sample_id":sampleids["New Zealand"],
  "methods":country_details,
  "site_inference_STOP":True},ensure_ascii=False,sort_keys=True)[:20000],flush=True)
if __name__=="__main__":main()
