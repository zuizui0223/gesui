#!/usr/bin/env python3
"""Independent four-country EPT ratio availability and UNKNOWN taxonomy mass QA.

Publisher original normalized *relative abundance* rows and visual fines. This
profiles measurement uncertainty only: DOES NOT fit sediment slope, correlation,
country transfer model, p-values, or claim sites are independent.
"""
import argparse,hashlib,json,math,pathlib,statistics
import gesui_external_fourcountry_xlsx_readonly_ooxml as src
from gesui_external_fourcountry_source_grain_20261122 import SHEETS,ORIGINAL_SHA,ORIGINAL_MD5,digest,quantiles
EPT={"Ephemeroptera","Plecoptera","Trichoptera"}
def percentile(x,p):
 return quantiles(x)
def run(a):
 raw=json.loads(a.source.read_text());tax=json.loads(a.taxonomy.read_text())
 if raw["schema"]!="gesui_external_McKenzie_fourcountry_6491_publisher_exact_sample_gate_v1":
  raise RuntimeError("Publisher sample-grain original source schema changed")
 if tax["schema"]!="gesui_original_McKenzie_4country_taxonomy_GBIF_EXACT_source_only_v1":
  raise RuntimeError("Original source-only GBIF taxonomy schema changed")
 if tax["Figshare_original_XLSX_SHA256"]!=ORIGINAL_SHA or raw["original_source_SHA256"]!=ORIGINAL_SHA:
  raise RuntimeError("Taxonomy and original EPA? No: external Figshare source SHA mismatch")
 if digest(a.xlsx,"sha256")!=ORIGINAL_SHA or digest(a.xlsx,"md5")!=ORIGINAL_MD5:
  raise RuntimeError("Source publisher-original XLSX bytes not equal SHA/MD5")
 mapping={x["original_label"]:x for x in tax["taxon_label_source_lookup"]}
 if len(mapping)!=302:raise RuntimeError("Original GBIF taxonomy label count changed")
 country={}
 for name,(sh,sed,n) in SHEETS.items():
  headers,rows=src.read_all(a.xlsx,sh)
  sh2,env=src.read_all(a.xlsx,sed)
  if len(headers)-1!=len(raw["country_source_profiles"][name]["original_taxon_fields"]):
   raise RuntimeError("Original publisher original taxon width changed "+name)
  if headers[1:]!=raw["country_source_profiles"][name]["original_taxon_fields"]:
   raise RuntimeError("Original source-only exact taxon column names changed")
  if len(rows)!=n or len(env)!=n:raise RuntimeError("Original country sample totals changed")
  if set(src.normalize_id(r[0]) for r in rows)!=set(src.normalize_id(r[0]) for r in env):
   raise RuntimeError("Original ANALYSIS_ID source sheet pairing changed")
  cols=headers[1:]
  classidx={k:[i for i,t in enumerate(cols) if mapping[t]["state"]==k] for k in ("EPT","NON_EPT","UNKNOWN")}
  EPTsum=[];unknown=[];nonEPT=[];fines=[];rows_gt_5=0;rows_zeropos=0
  fine={src.normalize_id(r[0]):src.as_float(r[1]) for r in env}
  for row in rows:
   vals=[src.as_float(v) for v in row[1:]]
   if len(vals)!=len(cols) or not all(math.isfinite(v) and v>=0 for v in vals):
    raise RuntimeError("Source abundance normalized cells invalid or missing "+name)
   total=sum(vals)
   if abs(total-1)>1e-9:raise RuntimeError("Original publisher taxon values are not normalized proportions")
   ept=sum(vals[i] for i in classidx["EPT"])
   unk=sum(vals[i] for i in classidx["UNKNOWN"])
   ne=sum(vals[i] for i in classidx["NON_EPT"])
   if abs(ept+unk+ne-1)>1e-9:raise RuntimeError("Source taxon exact-state mass not complete")
   f=fine[src.normalize_id(row[0])]
   if not math.isfinite(f) or not 0<=f<=100:raise RuntimeError("Source original fine cover missing")
   EPTsum.append(ept);unknown.append(unk);nonEPT.append(ne);fines.append(f)
   rows_gt_5+=int(unk>.05+1e-10)
   rows_zeropos+=int(ept==0)
  unknown_mean=float(statistics.mean(unknown))
  unknown_upper_q95=float(sorted(unknown)[round((len(unknown)-1)*.95)])
  gate=unknown_mean<=.02 and rows_gt_5/len(rows)<=.10
  country[name]={
   "source_original_sample_count":n,
   "publisher_SOURCE_ratio_sum_one_proven":True,
   "original_visual_fines_fraction_range":quantiles(fines),
   "taxon_header_classification":{
      "EPT_count":len(classidx["EPT"]),
      "NON_EPT_count":len(classidx["NON_EPT"]),
      "UNKNOWN_count":len(classidx["UNKNOWN"])},
   "unknown_taxon_fraction_of_source_relative_abundance":quantiles(unknown),
   "mean_unknown_taxon_fraction":unknown_mean,
   "fraction_samples_with_more_than_5percent_unknown":rows_gt_5/n,
   "source_EPT_known_lower_bound_relative_fraction":quantiles(EPTsum),
   "source_EPT_unknown_as_EPT_upper_bound_relative_fraction":quantiles([x+u for x,u in zip(EPTsum,unknown)]),
   "source_known_EPT_zero_fraction":rows_zeropos/n,
   "source_unknown_mass_gate_GO":gate,
   "physical_independent_SITE_ID_still_not_present":True,
   "source_no_sediment_EPT_slope_or_MSE_fitted":True}
 all_go=all(v["source_unknown_mass_gate_GO"] for v in country.values())
 result={
  "schema":"gesui_external_fourcountry_EPT_taxon_unknown_mass_source_gate_v1",
  "publisher_original_XLSX_SHA256":ORIGINAL_SHA,
  "taxonomy_original_source_json_sha256":hashlib.sha256(a.taxonomy.read_bytes()).hexdigest(),
  "country_profiles":country,
  "all_four_country_unknown_taxa_mass_GO":all_go,
  "country_heldout_descriptive_GO":all_go,
  "site_independent_inference_GO":False,
  "EPA_HUC8_history_specific_interaction_replication_GO":False,
  "numeric_benthic_sediment_associations_not_yet_fitted":True,
  "reason":"Known EPT relative mass is source-normalized taxon share. UNKNOWN mass bounded, not classified as non-EPT. No physical site ID or original individual count denominator in source workbook."}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print("FOURCOUNTRY_EPT_SOURCE_UNKNOWN_MASS_GATE",json.dumps({
  "all_country_GO":all_go,
  "country_summary":{k:{"n":v["source_original_sample_count"],
   "EPT_original_label_count":v["taxon_header_classification"]["EPT_count"],
   "unknown_original_taxon_label_count":v["taxon_header_classification"]["UNKNOWN_count"],
   "unknown_mass_mean":v["mean_unknown_taxon_fraction"],
   "more_than_5pct_unknown_samples":v["fraction_samples_with_more_than_5percent_unknown"],
   "known_EPT_relative_fraction":v["source_EPT_known_lower_bound_relative_fraction"],
   "country_GO":v["source_unknown_mass_gate_GO"]}
   for k,v in country.items()}},sort_keys=True),flush=True)
if __name__=="__main__":main()
