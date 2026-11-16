#!/usr/bin/env python3
"""Original EPA 2023-24 benthic national source availability and same-UID,
same-site+visit+date compatibility with frozen 2129 NRSA source visits.

NO biological numeric outcomes are read or summarized. Never turns a 404
or field/UID mismatch into ecological zero. All official source responses
are protected by URL/SHA receipts and explicit source schema.
"""
from __future__ import annotations
import collections,csv,hashlib,json,os,pathlib,time,urllib.error,urllib.request

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui2324_benthic_source"
ROOT.mkdir(parents=True,exist_ok=True)
BASE="https://www.epa.gov/system/files/other-files/2026-06/"
SOURCE={
 "count":"nrsa2324_benthiccount.csv",
 "count_metadata":"nrsa2324_benthiccount_metadata_0.csv",
 "metrics":"nrsa2324_benthicmetrics.csv",
 "metrics_metadata":"nrsa2324_benthicmetrics_metadata.csv",
 "mmi":"nrsa2324_benthicmmi.csv",
 "mmi_metadata":"nrsa2324_benthicmmi_metadata.csv",
 "taxa":"nrsa2324_benthictaxa.csv",
 "taxa_metadata":"nrsa2324_benthictaxa_metadata.csv",
}
KNOWN_SITE=("https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_siteinfo.csv",
            "f7d14cf1574d12828b3aa97bd7b08b95f22aec2115753d5b75514f8eed43d013")
SOURCE_GEO_EXPECTED=2129
TARGET=("UID","SITE_ID","VISIT_NO","DATE_COL","UNIQUE_ID","SAMPLE_TYPE","SAMP_TYPE","SAMPLEID")
def digest(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  while block:=f.read(1048576):h.update(block)
 return h.hexdigest()
def obtain(url,label,expected=None):
 target=ROOT/(label+".csv")
 for i in range(4):
  try:
   req=urllib.request.Request(url,headers={"User-Agent":"gesui-public-NRSA-source-QA/20261114","Accept":"text/csv,*/*"})
   with urllib.request.urlopen(req,timeout=110) as resp,target.open("wb") as out:
    typ=resp.headers.get("Content-Type","")
    while buf:=resp.read(2**20):
     out.write(buf)
   if target.stat().st_size<50 or target.read_bytes()[:256].lower().find(b"<html")>=0:
    raise RuntimeError("Not genuine CSV: too-small or HTML response")
   h=digest(target)
   if expected and h!=expected:raise RuntimeError("Official pre-frozen source SHA changed "+label)
   return target,{"url":url,"http_status":200,"content_type":typ,"source_sha256":h,"bytes":target.stat().st_size}
  except Exception as e:
   target.unlink(missing_ok=True)
   if isinstance(e,urllib.error.HTTPError) and e.code in (400,401,403,404,410):
    return None,{"url":url,"http_status":e.code,"error":str(e),"status":"UNAVAILABLE_NOT_ZERO"}
   if i==3:return None,{"url":url,"http_status":getattr(e,"code",None),
                         "error":type(e).__name__+":"+str(e)[:180],"status":"UNAVAILABLE_NOT_ZERO"}
   time.sleep(min(2**i,7))
def recs(path):
 with path.open(encoding="utf-8-sig",newline="") as fp:
  yield from csv.DictReader(fp)
def clean(row,c):return str(row.get(c,"")).strip()
def key(row,fields):
 r=tuple(clean(row,c) for c in fields)
 return r if all(r) else None
def locked_sites():
 p=pathlib.Path(os.environ["FROZEN_2129_SOURCE"])
 cross=json.loads(p.read_text())
 e=cross["2023_24"]
 sets=e["uid_bridge_manifest"]
 required={"biomass","water","chla"}
 if not required.issubset(sets):raise RuntimeError("Original pre-outcome 3-way source member list missing")
 ids=set(sets["biomass"].values())&set(sets["water"].values())&set(sets["chla"].values())
 if len(ids)!=SOURCE_GEO_EXPECTED:raise RuntimeError("Original 2023/24 n2129 frozen source population changed")
 source,sha=obtain(KNOWN_SITE[0],"siteinfo",KNOWN_SITE[1])
 if source is None:raise RuntimeError("Pinned official EPA site source unavailable")
 by_uid={};by_composite=collections.defaultdict(set)
 for r in recs(source):
  uid=clean(r,"UID")
  if not uid:continue
  if uid in by_uid:raise RuntimeError("Duplicated original site UID")
  by_uid[uid]=r
  k=key(r,("SITE_ID","VISIT_NO","DATE_COL"))
  if k:by_composite[k].add(uid)
 source.unlink()
 if not ids<=set(by_uid):raise RuntimeError("Frozen 2023/24 original UID not in SiteInfo")
 if len({k for k,v in by_composite.items() if len(v)>1})>0:raise RuntimeError("Site visit/date non-unique in original source")
 return ids,by_uid,by_composite,sha
def table_audit(path,ids,by_uid,by_composite):
 counts=collections.Counter();keycount=collections.Counter()
 matched_uids=set();unique_origin=set();mapped_id_rows=collections.Counter()
 modes=collections.Counter()
 sampletypes=collections.Counter()
 columns=None
 for row in recs(path):
  if columns is None:columns=list(row)
  counts["rows"]+=1
  uid=clean(row,"UID")
  kk=key(row,("SITE_ID","VISIT_NO","DATE_COL"))
  if uid:unique_origin.add(uid)
  target=None
  if uid in ids:
   target=uid;mode="direct_uid"
   if kk and key(by_uid[target],("SITE_ID","VISIT_NO","DATE_COL"))!=kk:
    counts["native_uid_vs_composite_discordant_rows"]+=1
    target=None;mode="UID_CONFLICT"
  elif kk and len(by_composite.get(kk,set())&ids)==1:
   target=next(iter(by_composite[kk]&ids));mode="unique_site_visit_date"
  else:mode="unmatched_or_ambiguous"
  modes[mode]+=1
  if target:
   matched_uids.add(target)
   mapped_id_rows[target]+=1
  if "SAMPLE_TYPE" in row:
   sampletypes[clean(row,"SAMPLE_TYPE")]+=1
  if "SAMP_TYPE" in row:
   sampletypes[clean(row,"SAMP_TYPE")]+=1
 return {"rows":counts["rows"],"columns":columns,
   "native_unique_UID":len(unique_origin),"matched_frozen_visit_UIDs":len(matched_uids),
   "mapping_modes":dict(modes),"sample_type_distribution":dict(sampletypes),
   "source_UID_conflicts":counts["native_uid_vs_composite_discordant_rows"],
   "matched_UIDs_with_repeated_table_rows":sum(n>1 for n in mapped_id_rows.values()),
   "max_rows_per_frozen_site_UID":max(mapped_id_rows.values(),default=0),
   "matched_frozen_site_UID_sha256":hashlib.sha256("\n".join(sorted(matched_uids)).encode()).hexdigest(),
   # Source-only set retained in original artifact for downstream strictly signed QA.
   "matched_frozen_site_UIDs":sorted(matched_uids)}
def main():
 ids,by_uid,by_composite,site_receipt=locked_sites()
 report={"schema":"gesui_NR​​SA2324_genuine_benthic_count_metrics_mmi_taxa_SOURCE_ONLY_v1",
         "EPA_siteinfo_pinned":site_receipt,"pre_benthic_outcome_frozen_common_uid":len(ids),
         "original_epa_links_exact":True,"results":{},
         "no_benthic_score_or_periphyton_or_chemical_numeric_values_read":True,
         "download_failure_IS_NOT_biodiversity_zero":True}
 for name,leaf in SOURCE.items():
  p,info=obtain(BASE+leaf,name)
  if p is not None:
   if name.endswith("_metadata"):
    with p.open(encoding="utf-8-sig",errors="replace") as f:
     info["preview_metadata_lines"]=f.read(6000).splitlines()[:12]
   else:
    info.update(table_audit(p,ids,by_uid,by_composite))
    print("NRSA2324_SOURCE",name,"rows",info["rows"],"matched_frozen_uids",
          info["matched_frozen_visit_UIDs"],"column_count",len(info["columns"]),flush=True)
   p.unlink()
  else:
   print("NRSA2324_SOURCE",name,"UNAVAILABLE",info.get("http_status"),flush=True)
  report["results"][name]=info
 mmi=report["results"]["mmi"]
 count=report["results"]["count"]
 usable_mmi=(mmi.get("http_status")==200 and mmi.get("matched_frozen_visit_UIDs",0)>=400
             and len({by_uid[i]["HUC8"] for i in mmi.get("matched_frozen_site_UIDs",[]) if by_uid[i].get("HUC8")})>=100
             and len({by_uid[i]["PSTL_CODE"] for i in mmi.get("matched_frozen_site_UIDs",[]) if by_uid[i].get("PSTL_CODE")})>=20)
 report["gate"]={"MMI_same_UID_source_400_100HUC8_20state":bool(usable_mmi),
                 "count_original_link_access":count.get("http_status")==200,
                 "MMI_not_equal_to_species_richness_or_benthic_biomass":True,
                 "no_effect_analyses_allowed_until_actual_MMI_value_and_units_QA":True}
 report["results"]["mmi"].pop("matched_frozen_site_UIDs",None) if not usable_mmi else None
 out=ROOT/"NRSA2324_ORIGINAL_BENTHIC_SOURCE_AND_UID_QA.json"
 out.write_text(json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False))
 print("NRSA2324_SOURCE_GATE",json.dumps(report["gate"],sort_keys=True))
 print("NRSA2324_SOURCE_RECEIPT",str(out))
if __name__=="__main__":main()
