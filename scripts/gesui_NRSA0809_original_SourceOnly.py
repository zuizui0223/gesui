#!/usr/bin/env python3
"""SOURCE-ONLY authentic EPA 2008/09 NRSA five full original CSV files.
Publisher metadata, raw source SHA, original native UID and 4/5 table join,
identification of true survey methods. NO numerical ecological effects read.
"""
import csv,hashlib,json,os,pathlib,time,urllib.request,urllib.error
from collections import Counter
R=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui0809_source"
R.mkdir(parents=True,exist_ok=True)
SRC={
 "siteinfo":"https://www.epa.gov/sites/default/files/2015-09/siteinfo_0.csv",
 "water":"https://www.epa.gov/sites/default/files/2015-09/chem.csv",
 "phab":"https://www.epa.gov/sites/default/files/2015-09/phabmed.csv",
 "bentcount":"https://www.epa.gov/sites/default/files/2016-11/nrsa0809bentcts.csv",
 "bentmmi":"https://www.epa.gov/sites/default/files/2015-09/bentcond.csv",
}
META={
 "siteinfo":"https://www.epa.gov/sites/default/files/2015-09/siteinfo_0.txt",
 "water":"https://www.epa.gov/sites/default/files/2015-09/chem.txt",
 "phab":"https://www.epa.gov/sites/default/files/2015-09/phabmed.txt",
 "bentcount":"https://www.epa.gov/sites/default/files/2016-11/nrsa0809bentctsmet.txt",
 "bentmmi":"https://www.epa.gov/sites/default/files/2015-09/bentcond.txt",
}
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  while z:=f.read(1024*1024):h.update(z)
 return h.hexdigest()
def download(k,url):
 target=R/(k+".csv")
 for n in range(3):
  try:
   req=urllib.request.Request(url,headers={"User-Agent":"gesui-NRSA0809-source-original-byteQA/20261114"})
   with urllib.request.urlopen(req,timeout=120) as w,target.open("wb") as f:
    typ=w.headers.get("Content-Type","")
    while b:=w.read(1024*1024):f.write(b)
   if target.stat().st_size<100 or b"<html" in target.read_bytes()[:256].lower():
    raise RuntimeError("Source URL returned HTML / small, not authentic EPA numeric CSV")
   return target,{"url":url,"status":"ORIGINAL_SOURCE_OK","source_SHA256":sha(target),
                  "source_bytes":target.stat().st_size,"content_type":typ}
  except Exception as e:
   target.unlink(missing_ok=True)
   if isinstance(e,urllib.error.HTTPError) and e.code in (400,401,403,404,410):
    return None,{"url":url,"status":"UNAVAILABLE_NOT_ZERO","http_status":e.code,"error":str(e)[:200]}
   if n==2:return None,{"url":url,"status":"UNAVAILABLE_NOT_ZERO","error":str(e)[:200]}
   time.sleep(2**n)
def audit(p):
 with p.open(encoding="latin1",newline="") as f:
  rd=csv.DictReader(f)
  cols=rd.fieldnames
  if not cols:raise RuntimeError("Original 0809 CSV header empty")
  upcols={k.strip().upper():k for k in cols}
  if "UID" not in upcols:raise RuntimeError("Original 0809 original UID missing "+repr(cols[:20]))
  uids=set();bysite=Counter();samp=Counter();n=0;years=Counter();states=Counter();siteclasses=Counter()
  for row in rd:
   n+=1
   uid=(row.get(upcols["UID"]) or "").strip()
   if uid:uids.add(uid);bysite[uid]+=1
   if "SAMPLE_TYPE" in upcols:samp[(row.get(upcols["SAMPLE_TYPE"]) or "").strip()]+=1
   if "SAMPLE_CAT" in upcols:samp["CATEGORY_"+(row.get(upcols["SAMPLE_CAT"]) or "").strip()]+=1
   for k in ("YEAR","PSTL_CODE","SITE_CLASS","SITETYPE"):
    if k not in upcols:continue
    x=(row.get(upcols[k]) or "").strip()
    if k=="YEAR":years[x]+=1
    elif k=="PSTL_CODE":states[x]+=1
    else:siteclasses[k+":"+x]+=1
  return {"original_rows":n,"native_UIDs":len(uids),
          "native_UIDs_with_multirows":sum(v>1 for v in bysite.values()),
          "source_header":cols,"source_sample_types":dict(samp),
          "source_years":dict(years),"source_states":dict(states),
          "source_site_classes":dict(siteclasses)},uids
def main():
 results={};sets={}
 for k,u in SRC.items():
  path,info=download(k,u)
  if path:
   profile,ids=audit(path)
   info.update(profile);sets[k]=ids
   path.unlink()
   print("NRSA0809_SOURCE_ORIGINAL",k,"ROWS",profile["original_rows"],"UIDS",len(ids),
         "SAMPLE",info["source_sample_types"],"CLASS",info["source_site_classes"],flush=True)
  else:print("NRSA0809_SOURCE_FAIL",k,info.get("http_status"),info.get("error"),flush=True)
  results[k]=info
 metadata={}
 for k,u in META.items():
  try:
   req=urllib.request.Request(u,headers={"User-Agent":"gesui-original-NRSA0809-metadata"})
   with urllib.request.urlopen(req,timeout=50) as resp:b=resp.read(200000)
   txt=b.decode("latin1",errors="replace")
   lines=txt.splitlines()
   relevant=[line for line in lines if any(x in line.upper() for x in ("XEMBED","XCDENMID","NTL","PTL","TOTAL300","BERW","BETB","SAMP_TYPE","CHLA","AFDM","HUC8","SAMPLE_TYPE"))]
   metadata[k]={"url":u,"sha256":hashlib.sha256(b).hexdigest(),"bytes":len(b),
                "relevant_metadata_lines":relevant[:25]}
  except Exception as e:
   metadata[k]={"url":u,"status":"UNAVAILABLE","error":str(e)[:150]}
 if set(("siteinfo","water","phab","bentcount"))<=set(sets):
  raw_common=set.intersection(*(sets[k] for k in ("siteinfo","water","phab","bentcount")))
  mmi=set.intersection(*(sets[k] for k in ("siteinfo","water","phab","bentcount","bentmmi"))) if "bentmmi" in sets else set()
 else:raw_common=set();mmi=set()
 source={"native_UID_four_source_intersection":len(raw_common),
         "native_UID_five_with_MMI_intersection":len(mmi),
         "native_UID_four_source_set_sha256":hashlib.sha256("\n".join(sorted(raw_common)).encode()).hexdigest(),
         "native_UID_four_original_source_only":sorted(raw_common),
         "no_periphyton_Chla_or_AFDM_source_assumed":True}
 print("NRSA0809_FULL_ORIGINAL_SOURCE_GATE",json.dumps({k:v for k,v in source.items() if not isinstance(v,list)},sort_keys=True),flush=True)
 out={"schema":"gesui_NRSA0809_AUTHENTIC_EPA_WATER_BENTHOS_PHAB_MMI_SOURCE_ONLY_v1",
      "original_EPA_NARS_reference":"https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys",
      "original_data_sources":results,"original_metadata":metadata,
      "original_source_join":source,"not_numeric_biodiversity_effects":True,
      "not_periphyton_biomass_or_emergent_insect_mass":True,
      "no_2018_EPT_Chla_replication_if_2008_pigment_absent":True,
      "no_WWTP_N_P_effluent_or_ecosystem_service":True}
 p=R/"NRSA0809_EPA_ORIGINAL_SOURCE_ONLY_COUNT_AND_GRANULARITY.json"
 p.write_text(json.dumps(out,sort_keys=True,indent=2))
 print("NRSA0809_SOURCE_OUT",str(p),flush=True)
if __name__=="__main__":main()
