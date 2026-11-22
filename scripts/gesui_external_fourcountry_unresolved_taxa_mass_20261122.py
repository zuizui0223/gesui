#!/usr/bin/env python3
"""Profile unresolved original GBIF taxon *relative abundance mass* by name.

This is an objective source-label quality repair diagnostic, not biological
sediment association. Does not fit EPT outcome, test a slope or select countries.
"""
import argparse,collections,hashlib,json,pathlib,math
import gesui_external_fourcountry_xlsx_readonly_ooxml as x
from gesui_external_fourcountry_source_grain_20261122 import SHEETS,ORIGINAL_SHA,ORIGINAL_MD5,digest
def main():
 p=argparse.ArgumentParser()
 for name in ("xlsx","source","taxonomy","out"):
  p.add_argument("--"+name,required=True,type=pathlib.Path)
 a=p.parse_args()
 orig=json.loads(a.source.read_text());tax=json.loads(a.taxonomy.read_text())
 if orig["original_source_SHA256"]!=ORIGINAL_SHA or tax["Figshare_original_XLSX_SHA256"]!=ORIGINAL_SHA:
  raise ValueError("original publisher source and taxonomy SHA mismatch")
 if digest(a.xlsx,"sha256")!=ORIGINAL_SHA or digest(a.xlsx,"md5")!=ORIGINAL_MD5:
  raise ValueError("publisher original XLSX bytes not source-pinned")
 classification={z["original_label"]:z for z in tax["taxon_label_source_lookup"]}
 results={}
 for country,(t,s,expected) in SHEETS.items():
  h,rows=x.read_all(a.xlsx,t)
  if len(rows)!=expected:raise RuntimeError("paper original sample n changed")
  headers=h[1:]
  unresolved=[(i,z) for i,z in enumerate(headers) if classification[z]["state"]=="UNKNOWN"]
  summary=[]
  for i,label in unresolved:
   fractions=[x.as_float(row[i+1]) for row in rows]
   if not all(math.isfinite(v) and v>=0 for v in fractions):
    raise ValueError("unknown taxon numeric fraction invalid")
   entry=classification[label]
   summary.append({"original_taxon_header":label,"source_GBIF_match_type":entry.get("GBIF_match_type"),
      "source_GBIF_status":entry.get("GBIF_status"),
      "source_GBIF_taxonomic_rank":entry.get("source_taxonomic_rank"),
      "source_GBIF_taxonomic_order":entry.get("source_taxonomic_order"),
      "mean_relative_individual_fraction":sum(fractions)/len(rows),
      "fraction_samples_present":sum(v>0 for v in fractions)/len(rows),
      "share_of_all_unknown_mass_in_country":None})
  summary.sort(key=lambda z:-z["mean_relative_individual_fraction"])
  allunknown=sum(z["mean_relative_individual_fraction"] for z in summary)
  for item in summary:
   item["share_of_all_unknown_mass_in_country"]=item["mean_relative_individual_fraction"]/allunknown if allunknown else 0
  results[country]={"sample_count":expected,
    "total_unresolved_fraction_from_all_unknown_labels":allunknown,
    "all_unresolved_taxonomic_source_labels_ranked_by_relative_mass":summary,
    "no_within_country_site_sampling_ID_available":True}
 output={"schema":"gesui_original_McKenzie_fourcountry_unresolved_taxonomic_abundance_source_only_v1",
  "source_XLSX_SHA256":ORIGINAL_SHA,
  "source_GBIF_exact_QA_SHA256":hashlib.sha256(a.taxonomy.read_bytes()).hexdigest(),
  "countries":results,
  "no_sediment_EPT_covariate_associations_computed":True,
  "taxon_names_only_source_sanitization":True,
  "unresolved_names_not_silently_reclassified":True}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
 print("EXTERNAL_4COUNTRY_UNKNOWN_TAXON_SOURCE_MASS",json.dumps({
  k:{"mean_unknown_fraction":z["total_unresolved_fraction_from_all_unknown_labels"],
     "top12":[{"name":r["original_taxon_header"],
       "mean_fraction":r["mean_relative_individual_fraction"],
       "GBIF_status":r.get("source_GBIF_status")}
     for r in z["all_unresolved_taxonomic_source_labels_ranked_by_relative_mass"][:12]]}
     for k,z in results.items()},sort_keys=True),flush=True)
if __name__=="__main__":main()
