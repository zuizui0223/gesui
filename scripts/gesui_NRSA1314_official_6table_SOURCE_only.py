#!/usr/bin/env python3
"""SOURCE-only independent 2013–14 original US EPA NRSA cycle:
SiteInfo, WaterChem, benthic BETB, periphyton chlorophyll,
periphyton AFDM and measured PHab. Read field headers, original source SHA,
UID uniqueness and source frame overlap; ZERO biological outcome magnitudes.
"""
from __future__ import annotations
import csv,hashlib,json,os,pathlib,time,urllib.error,urllib.request
from collections import Counter,defaultdict
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui1314_independent_source"
ROOT.mkdir(exist_ok=True,parents=True)
HOST="https://www.epa.gov/sites/default/files/2019-04/"
FILES={
 "sites":"nrsa1314_siteinformation_wide_04292019.csv",
 "water":"nrsa1314_widechem_04232019.csv",
 "chla":"nrsa1314_widepchl_04232019.csv",
 "phab":"nrsa1314_phabmed_04232019.csv",
 "benthic":"nrsa1314_bentcnts_04232019.csv",
 "afdm":"nrsa1314_widepbio_04232019.csv"}
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def get(name,leaf):
 url=HOST+leaf
 target=ROOT/leaf
 for i in range(3):
  try:
   req=urllib.request.Request(url,headers={"User-Agent":"gesui-public-nrsa1314-independent-source/20261114"})
   with urllib.request.urlopen(req,timeout=120) as resp,target.open("wb") as w:
    typ=resp.headers.get("Content-Type")
    while b:=resp.read(1048576):w.write(b)
   if target.stat().st_size<100 or b"<html" in target.read_bytes()[:200].lower():
    raise RuntimeError("Not original EPA numerical table")
   return target,{"url":url,"sha256":sha(target),"nbytes":target.stat().st_size,"HTTP":200,"content_type":typ}
  except Exception as e:
   target.unlink(missing_ok=True)
   if isinstance(e,urllib.error.HTTPError) and e.code in (400,401,403,404,410):
    return None,{"url":url,"status":"UNAVAILABLE_NOT_ZERO","http":e.code,"error":str(e)[:160]}
   if i==2:return None,{"url":url,"status":"UNAVAILABLE_NOT_ZERO","error":type(e).__name__+":"+str(e)[:150]}
   time.sleep(2**i)
def inspect(p):
 with p.open(encoding="utf-8-sig",newline="") as f:
  r=csv.DictReader(f)
  cols=r.fieldnames
  if not cols or "UID" not in cols:raise RuntimeError("Original EPA 2013 UID header not found")
  uids=set();count=Counter();types=Counter();yrs=Counter();site_ids=set()
  co=Counter()
  for row in r:
   uid=row.get("UID","").strip()
   co["rows"]+=1
   if not uid:co["missing_UID_rows"]+=1;continue
   uids.add(uid);count[uid]+=1
   if "SAMPLE_TYPE" in row:types[row.get("SAMPLE_TYPE","").strip()]+=1
   if "YEAR" in row:yrs[row.get("YEAR","").strip()]+=1
   if "SITE_ID" in row:site_ids.add(row.get("SITE_ID","").strip())
  return {"rows":co["rows"],"UID_count":len(uids),"missing_UID_rows":co["missing_UID_rows"],
          "UIDs_with_repeated_rows":sum(v>1 for v in count.values()),
          "raw_header":cols,"source_sample_types":dict(types),
          "source_year_labels":dict(yrs),"distinct_SITE_ID_header":len(site_ids)},uids
def main():
 entries={};sets={}
 for k,v in FILES.items():
  p,record=get(k,v)
  if p is not None:
   attrs,s=inspect(p)
   record.update(attrs);sets[k]=s
   p.unlink()
   print("NRSA1314_SOURCE_AVAILABLE",k,"rows",attrs["rows"],"uniqueUID",attrs["UID_count"],
         "sampletypes",dict(attrs["source_sample_types"]),flush=True)
  else:print("NRSA1314_SOURCE_UNAVAILABLE",k,record.get("http"),record.get("error"),flush=True)
  entries[k]=record
 if len(sets)!=6:
  source={"available_table_names":sorted(sets),"unavailable_table_names":sorted(set(FILES)-set(sets)),
          "multi_table_source_gate":False}
 else:
  common=set.intersection(*sets.values())
  intersection_bio_water=set.intersection(*(sets[k] for k in ("sites","water","chla","phab","benthic")))
  source={"all_six_joinable_native_UIDs":len(common),
          "five_without_AFD​​M_joinable_native_UIDs":len(intersection_bio_water),
          "has_>=350_national_source_eligible":len(intersection_bio_water)>=350,
          "all_six_common_UID_SHA256":hashlib.sha256("\n".join(sorted(common)).encode()).hexdigest(),
          "joint_UIDs_source_only":sorted(common),
          "five_core_UIDs_source_only":sorted(intersection_bio_water),
          "multi_table_source_gate":len(intersection_bio_water)>=350}
 out={
  "schema":"gesui_NRSA1314_6originalEPA_tables_source-only_before_any_new_biology_v1",
  "UTC":"2026-11-14",
  "official_EPA_nars_source_page":"https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys",
  "original_table_receipts":entries,
  "original_source_only_join":source,
  "source_only_numerical_effects_not_read":True,
  "old_NR​​SA1819_interaction_known_a_priori":True,
  "no_source_selection_based_on_new_independent_effect_sign":True,
  "zero_rows_or_HTTP_404_NOT_biological_zero":True
 }
 p=ROOT/"NRSA1314_EPA_ORIGINAL_SIX_TABLE_SOURCE_GATE.json"
 p.write_text(json.dumps(out,indent=2,sort_keys=True))
 print("NRSA1314_SOURCE_GATE",json.dumps({k:v for k,v in source.items() if not isinstance(v,list)},sort_keys=True),flush=True)
if __name__=="__main__":main()
