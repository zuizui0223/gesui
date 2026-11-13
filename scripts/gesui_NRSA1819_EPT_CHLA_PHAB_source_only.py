#!/usr/bin/env python3
"""SOURCE-only locked EPA NRSA 2018-19: EPT+river nutrients+areal chlorophyll
and physical habitat. All source identifiers from existing PR58, PR59 and
pre-effects PR73, no biological numeric outcome tested in this stage.
"""
from __future__ import annotations
import csv,hashlib,json,os,pathlib,time,urllib.request,collections
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui1819_triad_source"
ROOT.mkdir(parents=True,exist_ok=True)
LINKS={
 "sites":("https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_SiteInfo.csv","1b5b2119c881290191edeb47f1d80a7711f10440e9851574cfa86fd57be770e3"),
 "water":("https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_water_chemistry_chla_-_data.csv","cb9a3349ae3708834eab533e71423d74c1778bf76aba0737f0ae1bcc8f50baa1"),
 "benthos":("https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_benthic_macroinvertebrate_count_-_data.csv","4e9b3a0798b25f67c7a96d4800b94f86544bbd623abaf119358e067b158324a7"),
 "phab":("https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_physical_habitat_larger_set_of_metrics_-_data.csv","69b3b56962c04492423bd7eaf4d074c6ae8d253924d683ba12eb25fd857ce063"),
 "chla":("https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_PeriChla.csv","47a3663e9cc5e99ec54e581bda09c831604c0f267153f56c56ec6558221acf0a")
}
BENTHOS_FROZEN_SHA="fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91"
FROZEN_PRIOR_CHLA_DIRECT_CROSSWALK=2110
PHAB_FIELDS=("XPCAN","XGVEG","XFC_OHV","XWIDTH","XCDENMID","XEMBED","W1_HALL")
def get(url,key,expected):
 p=ROOT/(key+".csv")
 for k in range(3):
  try:
   req=urllib.request.Request(url,headers={"User-Agent":"gesui-2018-19-genuine-ecological-source-only/1.0"})
   with urllib.request.urlopen(req,timeout=120) as r,p.open("wb") as w:
    while chunk:=r.read(1<<20):w.write(chunk)
   sha=hashlib.sha256(p.read_bytes()).hexdigest()
   if sha!=expected:raise RuntimeError("Authentic EPA 2018/19 SHA mismatch "+key)
   return p,{"url":url,"sha256":sha,"bytes":p.stat().st_size}
  except Exception:
   p.unlink(missing_ok=True)
   if k==2:raise
   time.sleep(2**k)
def rows(p):
 with p.open(encoding="utf-8-sig",newline="") as f:yield from csv.DictReader(f)
def index(p):
 vals=set(); n=0;header=None;counts=collections.Counter()
 for r in rows(p):
  n+=1
  if header is None:header=list(r)
  uid=str(r.get("UID","")).strip()
  if uid:vals.add(uid);counts[uid]+=1
 return {"n_rows":n,"unique_UIDs":len(vals),"n_UID_with_repeated_rows":sum(x>1 for x in counts.values()),
         "columns":header},vals
def main():
 manifest=pathlib.Path(os.environ["NRSA_FROZEN_MANIFEST"])
 if hashlib.sha256(manifest.read_bytes()).hexdigest()!=BENTHOS_FROZEN_SHA:
  raise RuntimeError("Original NRSA pre-effect benthic UID list SHA changed")
 a=json.loads(manifest.read_text())
 M={r["uid"]:r for r in a["primary_berw_UID_subset"]}
 if len(M)!=1221:raise RuntimeError("Source-preselected native UID opportunity changed")
 cw=json.loads(pathlib.Path(os.environ["NRSA_BIOFILM_FROZEN_CROSSWALK"]).read_text())
 chlamap=cw["2018_19"]["uid_bridge_manifest"]["chla"]
 chla_canonical=set(chlamap.values())
 if len(chla_canonical)!=FROZEN_PRIOR_CHLA_DIRECT_CROSSWALK:raise RuntimeError("Older 2018 source-only Chla mapping changed")
 evidence={};sets={}
 for key,(url,sha) in LINKS.items():
  p,rec=get(url,key,sha)
  info,ids=index(p)
  rec.update(info)
  evidence[key]=rec;sets[key]=ids
  p.unlink()
 if len(sets["sites"]&set(M))!=1221 or len(sets["water"]&set(M))!=1221:
  raise RuntimeError("Original strict pre-effect benthic 1221 SiteInfo/water source UIDs not fully present")
 source={uid for uid in M if uid in sets["sites"] and uid in sets["water"] and uid in sets["phab"] and uid in chla_canonical and uid in sets["benthos"]}
 attrs={"native_survey_UID_source_opportunity":len(source),
        "distinct_HUC8":len({M[v]["huc8"] for v in source}),
        "distinct_states":len({M[v]["state"] for v in source}),
        "distinct_physical_site_ID":len({M[v]["site_id"] for v in source}),
        "source_phab_has_required_fields":{k:k in evidence["phab"]["columns"] for k in PHAB_FIELDS},
        "source_source_count_not_complete_numeric_quality_cases":True}
 if len(source)<600 or attrs["distinct_HUC8"]<300 or attrs["distinct_states"]<25:
  raise RuntimeError("New same-visit triad source supports insufficient")
 if not all(attrs["source_phab_has_required_fields"][k] for k in PHAB_FIELDS[:4]):
  raise RuntimeError("Source PHab broad canopy/riparian/width fields absent")
 out={"schema":"gesui_NRSA1819_EPT_Chla_PHab_true_native_source_gate_v1",
      "original_1221_source_benthic_SHA256":BENTHOS_FROZEN_SHA,
      "frozen_2018_Chla_unique_site_date_bridge_not_numeric_selected":True,
      "original_publisher_CSV_source_evidence":evidence,"pre_num_source_eligibility":attrs,
      "source_UIDs":sorted(source),
      "no_EPT_or_NP_or_chla_or_PHab_numeric_outcomes_read":True,
      "no_wastewater_source_discharge":True,
      "no_insect_emergence_or_pollination_service":True}
 path=ROOT/"NRSA1819_ORIGINAL_EPT_CHLA_PHAB_SOURCE_ONLY_QA.json"
 path.write_text(json.dumps(out,sort_keys=True,indent=2))
 print("NRSA1819_TRIPLE_ECOSYSTEM_SOURCE_GATE",json.dumps(attrs,sort_keys=True))
if __name__=="__main__":main()
