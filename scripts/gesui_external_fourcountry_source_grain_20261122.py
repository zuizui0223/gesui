#!/usr/bin/env python3
"""Original 6491 four-country source-grain, key coverage and measurement gate.

No individual-site model fit, no EPT classification from guessed names,
no correlations between taxa and sediment. Strict original publisher SHA+MD5.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,math,pathlib,re,statistics
import gesui_external_fourcountry_xlsx_readonly_ooxml as oo

ORIGINAL_SHA="a02573bef0e11ea19c9465a78ebecc130fb9c39d76d35a79863c2efd5debcd48"
ORIGINAL_MD5="1e4d6e37b73c663406e5567dd79d9452"
SHEETS={
 "Australia":("Australia taxa","Australia sediment",3099),
 "Brazil":("Brazil taxa","Brazil sediment",232),
 "New Zealand":("NZ taxa","Nz sediment",342),
 "UK":("UK taxa","UK Sediment",2818)}
def digest(path,kind):
 h=hashlib.new(kind)
 with open(path,"rb") as f:
  for b in iter(lambda:f.read(1048576),b""):h.update(b)
 return h.hexdigest()
def quantiles(x):
 x=sorted(v for v in x if math.isfinite(v))
 if not x:return {"valid":0}
 return {"valid":len(x),"min":x[0],"q10":x[int((len(x)-1)*.1)],
         "median":statistics.median(x),"q90":x[int((len(x)-1)*.9)],"max":x[-1]}
def idsource(keys):
 uniq=len(set(keys))
 bad=sum(not k for k in keys)
 types=collections.Counter()
 for k in keys:
  if re.fullmatch(r"\d+",k):types["numeric_integer_string"]+=1
  elif "." in k:types["dot_delimited_string"]+=1
  else:types["other_string"]+=1
 return {"rows":len(keys),"unique_ANALYSIS_ID":uniq,
         "duplicates":len(keys)-uniq,"blank":bad,
         "types":dict(types),"example_first_10":keys[:10],
         "no_independent_physical_site_ID_column":True}
def taxon_quality(rows):
 nrow=len(rows)
 if nrow==0:raise RuntimeError("No original taxon rows")
 nonfinite=0;nonzero=0;noninteger=0;negative=0;missing=0;sums=[];counts=[]
 for r in rows:
  v=[oo.as_float(x) for x in r[1:]]
  present=[x for x in v if math.isfinite(x)]
  missing+=len(v)-len(present)
  noninteger+=sum(abs(x-round(x))>1e-9 for x in present)
  negative+=sum(x<0 for x in present)
  nonzero+=sum(x>0 for x in present)
  nonfinite+=sum(not math.isfinite(x) for x in v)
  sums.append(sum(present))
  counts.append(sum(x>0 for x in present))
 return {"n_original_rows":nrow,"sum_all_numeric_taxon_columns":quantiles(sums),
   "number_nonzero_taxon_columns_per_sample":quantiles(counts),
   "invalid_numeric_or_missing_cells":nonfinite,
   "noninteger_numeric_taxon_cells":noninteger,
   "negative_numeric_taxon_cells":negative,
   "positive_numeric_taxon_cells":nonzero,
   "taxon_value_scale_observation_only_no_biological_model":True}
def inspect_one(path,country,spec):
 taxa_sheet,env_sheet,expected=spec
 names,rows=oo.read_all(path,taxa_sheet)
 enames,erows=oo.read_all(path,env_sheet)
 if len(rows)!=expected or len(erows)!=expected:
  raise RuntimeError("Original paper-reported paired sample counts changed "+country)
 if enames[0]!="ANALYSIS_ID" or len(enames)!=2:
  raise RuntimeError("Original sediment schema not ANALYSIS_ID + one visual fines indicator")
 ids=[oo.normalize_id(r[0]) for r in rows]
 envkeys=[oo.normalize_id(r[0]) for r in erows]
 status="GO"
 problems=[]
 if len(set(ids))!=len(ids) or len(set(envkeys))!=len(envkeys):
  status="STOP";problems.append("duplicate sample IDs: do not aggregate")
 if set(ids)!=set(envkeys):
  status="STOP";problems.append("original taxa and sediment ANALYSIS_ID do not form exact 1:1 match")
 fine={oo.normalize_id(r[0]):oo.as_float(r[1]) for r in erows}
 x=[fine[i] for i in ids]
 valid=[v for v in x if math.isfinite(v) and 0<=v<=100]
 hist=quantiles(valid)
 coverage=len(valid)
 distinct=len(set(round(v,8) for v in valid))
 if coverage<expected or distinct<25 or max(valid,default=0)-min(valid,default=0)<20:
  status="STOP";problems.append("sediment % original source gradient or coverage fails gate")
 reported=taxon_quality(rows)
 if reported["invalid_numeric_or_missing_cells"]>0:
  # Do NOT turn blank taxon cells into zero without publisher explanation
  problems.append("original taxon cells missing: biological denominator requires publisher meaning audit")
 if reported["negative_numeric_taxon_cells"]>0:
  status="STOP";problems.append("original negative taxon values")
 if status=="GO":status="SOURCE_SAMPLE_LINK_GO_SITE_UNKNOWN"
 return {
  "status":status,"source_problems":problems,
  "original_taxa_worksheet":taxa_sheet,"original_sediment_worksheet":env_sheet,
  "paper_reported_paired_n":expected,
  "original_taxon_fields":names[1:],
  "original_taxon_column_count":len(names)-1,
  "source_identifiers":{
    "taxon_analysis_id":idsource(ids),
    "sediment_analysis_id":idsource(envkeys),
    "exact_pairing_uid_set":set(ids)==set(envkeys),
    "same_rows_order":ids==envkeys,
    "unique_nonempty_pair_ID_count":len(set(ids)&set(envkeys))},
  "original_visual_sediment_field":enames[1],
  "visual_sediment_pct_range":hist,
  "distinct_visual_sediment_values":distinct,
  "taxon_value_quality":reported,
  "single_survey_independent_site_proven":False,
  "original_site_ID_not_exported_in_this_workbook":True,
  "source_next_taxonomy_crosswalk_needed":True}
def main():
 parser=argparse.ArgumentParser()
 parser.add_argument("--xlsx",type=pathlib.Path,required=True)
 parser.add_argument("--out",type=pathlib.Path,required=True)
 a=parser.parse_args()
 if digest(a.xlsx,"sha256")!=ORIGINAL_SHA or digest(a.xlsx,"md5")!=ORIGINAL_MD5:
  raise RuntimeError("Original publisher source 2023 Figshare XLSX bytes changed")
 results={k:inspect_one(a.xlsx,k,s) for k,s in SHEETS.items()}
 original_sums=sum(x["paper_reported_paired_n"] for x in results.values())
 if original_sums!=6491:raise RuntimeError("Published four-country sample count total changed")
 # Do NOT equate ANALYSIS_ID with independent site ID. No site is inferred
 # from dot-separated sample ID or country metadata.
 out={"schema":"gesui_external_McKenzie_fourcountry_6491_publisher_exact_sample_gate_v1",
  "original_source_SHA256":ORIGINAL_SHA,
  "publisher_original_MD5":ORIGINAL_MD5,
  "original_DOI":"10.17028/rd.lboro.24615570.v1",
  "gcb_publication_DOI":"10.1111/gcb.17084",
  "actual_country_sample_total":original_sums,
  "country_source_profiles":results,
  "actual_original_site_identifier_verified":False,
  "EPA_2023_field_XEMBED_equivalent_source_verified":False,
  "EPT_individual_fraction_crosswalk_not_yet_verified":True,
  "independent_site_grouping_possible":False,
  "inference_GO":False,
  "source_only_no_sediment_taxa_correlations_or_EPT_models":True}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
 print("FOURCOUNTRY_ORIGINAL_SOURCE_GRAIN_SUMMARY",json.dumps({
  k:{"status":v["status"],"n":v["paper_reported_paired_n"],
     "paired_ids":v["source_identifiers"]["unique_nonempty_pair_ID_count"],
     "analys_id":v["source_identifiers"]["taxon_analysis_id"],
     "original_visual_sediment_field":v["original_visual_sediment_field"],
     "visual_sediment_distribution":v["visual_sediment_pct_range"],
     "taxon_source_value_distribution":v["taxon_value_quality"]["sum_all_numeric_taxon_columns"],
     "taxon_noninteger_cells":v["taxon_value_quality"]["noninteger_numeric_taxon_cells"],
     "missing_taxon_cells":v["taxon_value_quality"]["invalid_numeric_or_missing_cells"],
     "source_problems":v["source_problems"]}
     for k,v in results.items()},sort_keys=True),flush=True)
if __name__=="__main__":main()
