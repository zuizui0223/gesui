#!/usr/bin/env python3
"""Blinded ID crosswalk + categorical lab QA for authentic 2018-19 and 2023-24
EPA NRSA biofilm and chlorophyll source tables; NO numeric effects."""
from __future__ import annotations
import collections,csv,hashlib,json,os,pathlib,time,urllib.request

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"nrsa_uid_biomass_crosswalk"
ROOT.mkdir(exist_ok=True,parents=True)
URLS={
 "2018_19":{
 "sites":"https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_SiteInfo.csv",
 "water":"https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_water_chemistry_chla_-_data.csv",
 "biomass":"https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_pbio_0.csv",
 "chla":"https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_PeriChla.csv"},
 "2023_24":{
 "sites":"https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_siteinfo.csv",
 "water":"https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_waterchemistry.csv",
 "biomass":"https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_periphytonbiomass.csv",
 "chla":"https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_periphytonchlorophyll-a.csv"}
}
EXPECTED={
 "2018_19_sites":"1b5b2119c881290191edeb47f1d80a7711f10440e9851574cfa86fd57be770e3",
 "2018_19_water":"cb9a3349ae3708834eab533e71423d74c1778bf76aba0737f0ae1bcc8f50baa1",
 "2018_19_biomass":"048b3b773c1e6a6513fb5f8a0bba7a04fdb8cdf4e9b2dce536740a1b38a0cfc",
 "2018_19_chla":"47a3663e9cc5e99e3b1a495bcb38775e2d5d610b7d00efbe43c778dd50b2e2a",
}
# 2018 biomass/chla source SHAs are checked against exact full hashes in previous
# SOURCE_ONLY manifest in the next action; this pass compares 2018 site/water pins.
EXPECTED={k:v for k,v in EXPECTED.items() if k in ("2018_19_sites","2018_19_water")}
def pull(url, label):
 p=ROOT/(label+".csv")
 req=urllib.request.Request(url,headers={"User-Agent":"gesui-id-audit/20261113"})
 for attempt in range(4):
  try:
   with urllib.request.urlopen(req,timeout=90) as f, p.open("wb") as g:
    while b:=f.read(2**20):g.write(b)
   sha=hashlib.sha256(p.read_bytes()).hexdigest()
   if label in EXPECTED and sha!=EXPECTED[label]:raise ValueError("2018 pinned source SHA mismatch")
   return p,sha
  except Exception:
   p.unlink(missing_ok=True)
   if attempt==3:raise
   time.sleep(2**attempt)
def rows(p):
 with p.open(encoding="utf-8-sig",newline="") as f:
  yield from csv.DictReader(f)
def key(r,fields):
 a=tuple(str(r.get(c,"")).strip() for c in fields)
 return a if all(a) else None
K={"site_visit_date":["SITE_ID","VISIT_NO","DATE_COL"],
   "site_visit":["SITE_ID","VISIT_NO"],
   "unique_id":["UNIQUE_ID"]}
report={}
for era,urls in URLS.items():
 D={};shas={}
 for label,url in urls.items():
  p,shas[label]=pull(url,era+"_"+label)
  D[label]=list(rows(p))
  print("Genuine",era,label,"rows",len(D[label]),"SHA",shas[label][:12],flush=True)
  p.unlink()
 sites=D["sites"]
 indexes={k:collections.defaultdict(set) for k in K}
 siteuids=collections.defaultdict(list)
 for r in sites:
  u=str(r.get("UID","")).strip()
  if not u:continue
  siteuids[u].append(r)
  for name,cols in K.items():
   v=key(r,cols)
   if v:indexes[name][v].add(u)
 summary={"source_sha256":shas,"site_unique_Uids":len(siteuids),
          "site_duplicated_Uids":sum(len(v)>1 for v in siteuids.values()),
          "tables":{},"uid_bridge_manifest":{}}
 for table in ["water","biomass","chla"]:
  dd=D[table]
  match=collections.Counter()
  unit=collections.Counter()
  unit_vol=collections.Counter()
  result_availability=collections.Counter()
  analyte=collections.Counter()
  sample=collections.Counter()
  flag=collections.Counter()
  mapped={}
  collisions=collections.Counter()
  for r in dd:
   uid=str(r.get("UID","")).strip()
   lookup={}
   for name,cols in K.items():
    k=key(r,cols)
    lookup[name]=indexes[name].get(k,set()) if k else set()
   direct=set([uid]) if uid in siteuids and len(siteuids[uid])==1 else set()
   matched=set(lookup["site_visit_date"])
   if len(direct)==1:
    canonical=uid;way="exact_UID"
   elif len(matched)==1:
    canonical=next(iter(matched));way="site_visit_date"
   elif len(lookup["site_visit"])==1:
    canonical=next(iter(lookup["site_visit"]));way="site_visit_ONLY"
   elif len(lookup["unique_id"])==1:
    canonical=next(iter(lookup["unique_id"]));way="unique_id_ONLY"
   else:
    canonical=None;way="AMBIGUOUS_OR_ABSENT"
   match[way]+=1
   if len(matched)==1 and len(direct)==1 and matched!=direct:
    collisions["UID_vs_site_visit_date_discordance"]+=1
   if len(lookup["site_visit_date"])>1:collisions["duplicate_site_visit_date_key"]+=1
   if canonical is not None and uid and len(mapped)<100000:
    if uid in mapped and mapped[uid]!=canonical:collisions["one_source_UID_multiple_dest"]+=1
    mapped[uid]=canonical
   unit[str(r.get("RESULT_UNITS","")).strip()]+=1
   unit_vol[str(r.get("RESULT_VOL_UNITS","")).strip()]+=1
   result_availability["RESULT_nonblank" if str(r.get("RESULT","")).strip() else "RESULT_blank"]+=1
   result_availability["RESULT_VOL_nonblank" if str(r.get("RESULT_VOL","")).strip() else "RESULT_VOL_blank"]+=1
   analyte[str(r.get("ANALYTE","")).strip()]+=1
   sample[str(r.get("SAMPLE_TYPE","")).strip()]+=1
   flag[str(r.get("NARS_FLAG","")).strip()]+=1
  summary["tables"][table]={
   "original_rows":len(dd),"unique_original_uid":len({str(z.get("UID","")).strip() for z in dd}),
   "join_path_rows":dict(match),"distinct_canonical_site_UID":len(set(mapped.values())),
   "key_collisions":dict(collisions),"RESULT_UNITS":dict(unit),
   "RESULT_VOL_UNITS":dict(unit_vol),"result_cell_availability":dict(result_availability),
   "ANALYTE":dict(analyte),"SAMPLE_TYPE":dict(sample),"NARS_FLAG":dict(flag),
   "direct_UID_site_intersection":len({str(z.get("UID","")).strip() for z in dd}&set(siteuids)),
  }
  summary["uid_bridge_manifest"][table]=mapped
 comb=set(summary["uid_bridge_manifest"]["biomass"].values()) &set(summary["uid_bridge_manifest"]["water"].values())
 both=comb &set(summary["uid_bridge_manifest"]["chla"].values())
 summary["source_only_canonical_UID_overlap"]={"biomass+water":len(comb),
        "biomass+water+chla":len(both)}
 report[era]=summary
# Store only source IDs + QA; no numeric biological responses and no effect sizes.
output=ROOT/"SOURCE_ONLY_NRSA_2ERA_UID_CONCORDANCE.json"
output.write_text(json.dumps(report,indent=2,sort_keys=True))
print("CROSSWALK_SUMMARY_START")
print(json.dumps({k:{n:v for n,v in val.items() if n!="uid_bridge_manifest"} for k,val in report.items()}))
print("CROSSWALK_SUMMARY_END")
