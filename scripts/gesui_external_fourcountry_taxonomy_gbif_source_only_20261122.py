#!/usr/bin/env python3
"""Taxonomic source-only validation: published Figshare 4-country taxa headers
against exact GBIF Backbone names and original Ephemeroptera/Plecoptera/Trichoptera ORDER.

Do not read biological abundance values, choose taxa from outcome signs, or infer
unmatched taxa are non-EPT. Unknown taxa remain UNKNOWN and count toward upper bound.
"""
from __future__ import annotations
import argparse,collections,concurrent.futures,datetime,hashlib,json,pathlib,time,urllib.parse,urllib.request,urllib.error

BACKBONE="https://api.gbif.org/v1/species/match"
EPT={"Ephemeroptera","Plecoptera","Trichoptera"}
def lookup(name):
 url=BACKBONE+"?"+urllib.parse.urlencode({"name":name,"strict":"true"})
 request=urllib.request.Request(url,headers={
  "User-Agent":"gesui-ecology-research-source-taxonomy-QA/20261122",
  "Accept":"application/json"})
 for i in range(4):
  try:
   with urllib.request.urlopen(request,timeout=35) as response:
    data=json.loads(response.read(100000))
   typ=data.get("matchType");rank=data.get("rank")
   canonical=data.get("canonicalName") or ""
   # GBIF scientificName may contain author; canonicalName is the
   # unambiguous original canonical match for exact family/order.
   strict=typ=="EXACT" and canonical.lower()==name.lower() and data.get("status") in ("ACCEPTED","SYNONYM")
   order=data.get("order")
   if strict and rank=="ORDER" and canonical in EPT:order=canonical
   state="EPT" if strict and order in EPT else ("NON_EPT" if strict and order else "UNKNOWN")
   return {"original_label":name,"state":state,"source_taxonomic_order":order if strict else None,
      "source_taxonomic_rank":rank if strict else None,
      "GBIF_match_type":typ,"GBIF_canonical_name":canonical,
      "GBIF_usageKey":data.get("usageKey") if strict else None,
      "GBIF_status":data.get("status") if strict else None,
      "http_status":200,"source_backbone":BACKBONE,"url":url}
  except urllib.error.HTTPError as e:
   err="HTTP_"+str(e.code)
   if e.code in (429,500,502,503,504) and i<3:
    time.sleep((i+1)*2);continue
   break
  except Exception as e:
   err=type(e).__name__+":"+str(e)[:160]
   if i<3:
    time.sleep(i+1);continue
   break
 return {"original_label":name,"state":"UNKNOWN","source_taxonomic_order":None,
   "GBIF_error":err,"source_backbone":BACKBONE,"url":url}
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--source",required=True,type=pathlib.Path)
 p.add_argument("--out",required=True,type=pathlib.Path)
 a=p.parse_args()
 src=json.loads(a.source.read_text())
 if src.get("schema")!="gesui_external_McKenzie_fourcountry_6491_publisher_exact_sample_gate_v1":
  raise RuntimeError("Unknown publisher source schema")
 if src["original_source_SHA256"]!="a02573bef0e11ea19c9465a78ebecc130fb9c39d76d35a79863c2efd5debcd48":
  raise RuntimeError("Original source XLSX hash changed")
 bycountry={k:set(v["original_taxon_fields"]) for k,v in src["country_source_profiles"].items()}
 all_names=sorted(set.union(*bycountry.values()))
 if len(all_names)<200:raise RuntimeError("Publisher original taxa columns unexpectedly incomplete")
 results=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
  for item in pool.map(lookup,all_names):
   results.append(item)
 valid={r["original_label"]:r for r in results}
 counts={k:collections.Counter(valid[x]["state"] for x in names)
         for k,names in bycountry.items()}
 overall=collections.Counter(x["state"] for x in results)
 output={
  "schema":"gesui_original_McKenzie_4country_taxonomy_GBIF_EXACT_source_only_v1",
  "Figshare_original_XLSX_SHA256":src["original_source_SHA256"],
  "source_original_Figshare_file_id":43299387,
  "GBIF_backbone_match_endpoint":BACKBONE,
  "lookup_date_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat(),
  "original_unique_taxon_labels":len(all_names),
  "overall_taxonomic_classification_counts":dict(overall),
  "by_country_header_taxonomy":{k:dict(v) for k,v in counts.items()},
  "taxon_label_source_lookup":results,
  "GBIF_EXACT_names_no_taxon_guessing":True,
  "DOUBTFUL_GBIF_matches_excluded":True,
  "source_original_OUTCOME_VALUES_NOT_READ":True,
  "EPT_fraction_numerical_analysis_GO":False,
  "next_gate":"Compute source numeric fraction occupied by UNKNOWN names per 4 countries before allowing EPT metrics."}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(output,sort_keys=True,indent=2)+"\n")
 print("EXTERNAL_4COUNTRY_GBIF_TAXONOMY_SOURCE",json.dumps({
   "unique_taxon_labels":len(all_names),"counts":dict(overall),
   "by_country":output["by_country_header_taxonomy"],
   "EPT_names":{k:[name for name in sorted(names) if valid[name]["state"]=="EPT"]
                 for k,names in bycountry.items()},
   "UNKNOWN_names":[r["original_label"] for r in results if r["state"]=="UNKNOWN"]
  },sort_keys=True),flush=True)
if __name__=="__main__":main()
