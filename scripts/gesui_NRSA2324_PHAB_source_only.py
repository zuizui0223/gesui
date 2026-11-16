#!/usr/bin/env python3
"""SOURCE-only authentic EPA NRSA2324 Physical Habitat Metrics original CSV:
publisher SHA, native UID, one-to-one source mapping, metadata definitions,
and blinded physical covariate availability. NO nutrient/effect values fitted.
"""
import csv, hashlib, json, os, pathlib, time, urllib.request, urllib.error, collections

R=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_phab_source"
R.mkdir(parents=True,exist_ok=True)
PHAB="https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_physicalhabitat_metrics.csv"
META="https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_physicalhabitat_metrics_metadata.csv"
SITE="https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_siteinfo.csv"
SITE_SHA="f7d14cf1574d12828b3aa97bd7b08b95f22aec2115753d5b75514f8eed43d013"
# Frozen from the publisher's prior-cycle 2018-19 PHab interpretation table,
# not tuned after seeing 2023-24 nutrient or biomass responses.
CANDIDATE_GROUPS={
 "canopy_light_proxy":["XCDENMID","XPCAN"],
 "riparian_vegetation":["XGVEG","XFC_OHV"],
 "streambed_embedding":["XEMBED","XCEMBED"],
 "fine_sediment":["LITSB1_SF","LITSB1_FN"],
 "channel_width":["XWIDTH"],
}
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def download(url,name,expected=None,required=True):
 p=R/name;err=""
 for i in range(3):
  try:
   q=urllib.request.Request(url,headers={"User-Agent":"gesui-NRSA2324-physical-source-20261114","Accept":"text/csv,*/*"})
   with urllib.request.urlopen(q,timeout=105) as resp,p.open("wb") as f:
    while b:=resp.read(2**20):f.write(b)
   if p.stat().st_size<60 or b"<html" in p.read_bytes()[:200].lower():raise RuntimeError("empty or HTML not authentic CSV")
   dig=sha(p)
   if expected and dig!=expected:raise RuntimeError("Publisher original source SHA changed")
   return p,{"url":url,"source_sha256":dig,"source_bytes":p.stat().st_size,"http_status":200}
  except urllib.error.HTTPError as e:
   err=str(e)
   if e.code in (400,401,403,404,410):break
   if i<2:time.sleep(2**i)
  except Exception as e:
   err=type(e).__name__+":"+str(e)
   if i<2:time.sleep(2**i)
  p.unlink(missing_ok=True)
 if required:raise RuntimeError(name+" original publisher source unavailable "+err)
 return None,{"url":url,"source_error":err,"status":"UNAVAILABLE_NOT_ZERO"}
def rows(path):
 with path.open(encoding="utf-8-sig",newline="") as f:
  yield from csv.DictReader(f)
def strval(row,k):return str(row.get(k,"")).strip()
def main():
 frozen=json.loads(pathlib.Path(os.environ["FROZEN_2129_SOURCE"]).read_text())
 y=frozen["2023_24"]["uid_bridge_manifest"]
 valid=set(y["biomass"].values())&set(y["water"].values())&set(y["chla"].values())
 if len(valid)!=2129:raise RuntimeError("PR73 frozen original pre-numeric native UID cohort changed")
 sp,si=download(SITE,"nrsa2324_sites.csv",SITE_SHA)
 index={};bykey=collections.defaultdict(set)
 for r in rows(sp):
  uid=strval(r,"UID")
  if uid in index:raise RuntimeError("Original SiteInfo repeated unique UID")
  index[uid]=r
  k=tuple(strval(r,c) for c in ["SITE_ID","VISIT_NO","DATE_COL"])
  if all(k):bykey[k].add(uid)
 if not valid <= set(index):raise RuntimeError("Original 2129 native UID not fully supported by authentic SiteInfo")
 pp,pinfo=download(PHAB,"nrsa2324_phab.csv")
 matched=set();counts=collections.Counter();mapping={}
 originalcols=None
 metricfields=collections.Counter()
 for r in rows(pp):
  if originalcols is None:originalcols=list(r)
  counts["rows"]+=1
  uid=strval(r,"UID")
  k=tuple(strval(r,c) for c in ["SITE_ID","VISIT_NO","DATE_COL"])
  match=None
  if uid in valid:
   same=tuple(strval(index[uid],c) for c in ["SITE_ID","VISIT_NO","DATE_COL"])
   if all(k) and all(same) and k!=same:
    counts["UID_vs_site_date_conflict"]+=1
    continue
   match=uid;counts["exact_UID"]+=1
  elif all(k) and len(bykey[k]&valid)==1:
   match=next(iter(bykey[k]&valid));counts["source_site_visit_date_fallback"]+=1
  else:
   counts["unmatched"]+=1
  if match is None:continue
  matched.add(match)
  counts["matched_rows"]+=1
  mapping.setdefault(match,0);mapping[match]+=1
  for group,names in CANDIDATE_GROUPS.items():
   for name in names:
    if name in r and strval(r,name)!="":
     metricfields[name]+=1
 mi,meta=download(META,"nrsa2324_phab_metadata.csv",required=False)
 defs={}
 if mi:
  with mi.open(encoding="utf-8-sig",newline="") as f:
   rd=csv.reader(f)
   allrows=list(rd)
  for nm in {x for names in CANDIDATE_GROUPS.values() for x in names}:
   matched_rows=[row for row in allrows if nm.upper() in [str(v).strip().upper() for v in row[:3]]]
   if matched_rows:defs[nm]=matched_rows[:2]
  meta["header"]=allrows[0] if allrows else []
 result={
  "schema":"gesui_NR​​SA2324_ORIGINAL_PHAB_SOURCE_ONLY_20261114_v1",
  "original_2129_site_biomass_chla_water_source_UIDs":len(valid),
  "source_siteinfo_sha256":si["source_sha256"],
  "phab":pinfo,"phab_column_names":originalcols,
  "phab_rows_and_native_ID_QA":dict(counts),
  "native_site_visits_matched":len(matched),
  "matched_UIDs_with_multiple_PHAB_rows":sum(v>1 for v in mapping.values()),
  "HUC8_of_matched_visit_UIDs":len({strval(index[i],"HUC8") for i in matched if strval(index[i],"HUC8")}),
  "states_of_matched_visit_UIDs":len({strval(index[i],"PSTL_CODE") for i in matched if strval(index[i],"PSTL_CODE")}),
  "candidate_physical_predictor_predeclared_groups":CANDIDATE_GROUPS,
  "candidate_field_nonempty_count":dict(metricfields),
  "metadata":meta,"candidate_metadata_field_definitions":defs,
  "source_eligible_original_UIDs":sorted(matched),
  "no_measurement_numeric_analysis":True,
  "no_WWTP_causal_exposure":True,
  "light_interpretation":"Canopy cover is a field proxy for shade, not measured solar radiation.",
  "gate_phab_support_min400_100HUC8_20states":len(matched)>=400 and len({strval(index[i],"HUC8") for i in matched})>=100 and len({strval(index[i],"PSTL_CODE") for i in matched})>=20
 }
 out=R/"NRSA2324_PHAB_SOURCE_ONLY_QA_20261114.json"
 out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False))
 print("NRSA2324_PHAB_SOURCE_QA_START")
 print(json.dumps({k:v for k,v in result.items() if k not in ["source_eligible_original_UIDs","phab_column_names"]},ensure_ascii=False,sort_keys=True))
 print("NRSA2324_PHAB_COLUMNS",json.dumps(originalcols))
 print("NRSA2324_PHAB_SOURCE_QA_END")
if __name__=="__main__":main()
