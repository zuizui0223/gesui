#!/usr/bin/env python3
"""SOURCE ONLY: inspect 2013–14 true EPA six publisher CSVs (source SHA frozen),
original SiteInfo population PROB vs HAND, native UID original BERW vs BETB,
schema-verified TN/TP, Chla and AFDM physical units/categorical NARS flags.
No numeric biological outcome magnitudes or regression coefficients.
"""
from __future__ import annotations
import collections,hashlib,json,os,pathlib,urllib.request
import pandas as pd

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui1314_independent_semantics"
ROOT.mkdir(exist_ok=True,parents=True)
BASE="https://www.epa.gov/sites/default/files/2019-04/"
NAMES={
 "sites":"nrsa1314_siteinformation_wide_04292019.csv",
 "water":"nrsa1314_widechem_04232019.csv",
 "chla":"nrsa1314_widepchl_04232019.csv",
 "afdm":"nrsa1314_widepbio_04232019.csv",
 "phab":"nrsa1314_phabmed_04232019.csv",
 "benthic":"nrsa1314_bentcnts_04232019.csv",
}
manifest=json.loads(pathlib.Path(os.environ["NRSA1314_ORIGINAL_SIX_SOURCE_QA"]).read_text())
receipts=manifest["original_table_receipts"]
if not manifest["original_source_only_join"]["multi_table_source_gate"]:
 raise RuntimeError("Original 2013 independent 6-table source availability gate not passed")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fetch(name):
 p=ROOT/NAMES[name]
 req=urllib.request.Request(BASE+NAMES[name],headers={"User-Agent":"gesui-original-2013-independent-measurement-contract"})
 with urllib.request.urlopen(req,timeout=120) as r,p.open("wb") as f:
  while chunk:=r.read(2**20):f.write(chunk)
 if sha(p)!=receipts[name]["sha256"]:raise RuntimeError("Source original publisher 2013 CSV SHA differs: "+name)
 return p
def read(name):
 p=fetch(name)
 df=pd.read_csv(p,dtype=str,keep_default_na=False,encoding="latin1",low_memory=False)
 p.unlink()
 df.columns=df.columns.str.upper().str.strip()
 if df.UID.eq("").any():raise RuntimeError("Original 2013 blank native UID "+name)
 return df
data={name:read(name) for name in NAMES}
site=data["sites"];B=data["benthic"]
if site.UID.duplicated().any():raise RuntimeError("Original SiteInfo UID nonunique")
for name in ("water","phab","chla","afdm"):
 if data[name].UID.duplicated().any():raise RuntimeError(name+" source repeated UID; STOP to disambiguate before outcome")
CANDIDATES={
 "sites":["UID","SITE_ID","YEAR","DATE_COL","VISIT_NO","PSTL_CODE","AG_ECO9","HUC8","HUC_8","LAT_DD83","LON_DD83","STRAH_ORD","SITETYPE","WGT_TP_CORE","WGT_TP"],
 "water":["NTL_RESULT","PTL_RESULT","NTL_UNITS","PTL_UNITS","NTL_RESULT_UNITS","PTL_RESULT_UNITS","COND_RESULT_UNITS","NTL_NARS_FLAG","PTL_NARS_FLAG","COND_RESULT","COND_UNITS","COND_NARS_FLAG","SAMPLE_TYPE","DATE_COL","SITE_ID"],
 "chla":["SAMPLE_TYPE","RESULT","RESULT_VOL","UNITS","RESULT_UNITS","RESULT_UNIT","RESULT_VOL_UNITS","UNIT","PARAMETER","PARAMETER_CODE","NARS_FLAG","RESULT_FLAG","CHLA_RESULT","CHLA_RESULT_UNITS","CHLA_NARS_FLAG","CHLA_UNITS","DATE_COL","SITE_ID"],
 "afdm":["SAMPLE_TYPE","RESULT","RESULT_VOL","UNITS","RESULT_UNITS","RESULT_UNIT","RESULT_VOL_UNITS","UNIT","PARAMETER","PARAMETER_CODE","NARS_FLAG","RESULT_FLAG","AFDM_RESULT","AFDM_RESULT_UNITS","AFDM_NARS_FLAG","DATE_COL","SITE_ID"],
 "phab":["UID","SITETYPE","XEMBED","XCEMBED","XCDENMID","XPCAN","XGVEG","XFC_OHV","XWIDTH","W1_HALL","DATE_COL","SITE_ID"],
 "benthic":["UID","SITE_ID","SAMPLE_TYPE","ORDER","FAMILY","TAXA_ID","TOTAL300","TOTAL","IS_DISTINCT300","DATE_COL"]
}
profiles={}
for k,df in data.items():
 cols=list(df.columns)
 present=[x for x in CANDIDATES[k] if x in cols]
 profile={"row_n":len(df),"uids":df.UID.nunique(),
          "present_fields":present,
          "missing_fields":sorted(set(CANDIDATES[k])-set(cols)),
          "all_columns":cols}
 for nm in present:
  if nm in ("SAMPLE_TYPE","SITETYPE","YEAR","NTL_UNITS","PTL_UNITS","COND_UNITS",
            "UNITS","RESULT_UNITS","RESULT_UNIT","RESULT_VOL_UNITS","UNIT",
            "PARAMETER","PARAMETER_CODE","NARS_FLAG","RESULT_FLAG","CHLA_UNITS",
            "NTL_NARS_FLAG","PTL_NARS_FLAG","COND_NARS_FLAG",
            "NTL_RESULT_UNITS","PTL_RESULT_UNITS","COND_RESULT_UNITS",
            "CHLA_RESULT_UNITS","CHLA_NARS_FLAG","AFDM_RESULT_UNITS","AFDM_NARS_FLAG"):
   cnt=collections.Counter(df[nm].str.strip().str.upper())
   profile["categorical_"+nm]=dict(cnt.most_common(25))
 profiles[k]=profile
 print("NRSA1314_SCHEMA_FIELDS",k,json.dumps({z:v for z,v in profile.items() if z!="all_columns"},sort_keys=True),flush=True)
 print("NRSA1314_HEADER",k,json.dumps(cols),flush=True)
 if k in ("chla","afdm","sites","phab","water"):
  print("NRSA1314_MEASUREMENT_FIELD_HEADERS",k,json.dumps([x for x in cols if any(y in x for y in ("CHLA","AFDM","XEMBED","XCEMBED","PCHL","PBIO","RESULT","NARS_FLAG","WGT_TP","YEAR","DATE_COL","XGV"))]),flush=True)
b_wade=set(B.loc[B.SAMPLE_TYPE.str.upper().eq("BERW"),"UID"])
b_boat=set(B.loc[B.SAMPLE_TYPE.str.upper().eq("BETB"),"UID"])
frame=set.intersection(*(set(data[k].UID) for k in data))
wade=frame & b_wade
boat=frame & b_boat
prob=set(site.loc[site.SITETYPE.str.upper().eq("PROB"),"UID"]) if "SITETYPE" in site else set(frame)
eligible=wade & prob
if "WGT_TP_CORE" in site:
 weight=pd.to_numeric(site.set_index("UID").loc[sorted(eligible),"WGT_TP_CORE"],errors="coerce")
 eligible={u for u,v in zip(sorted(eligible),weight) if pd.notna(v) and v>0}
elif "WGT_TP" in site:
 weight=pd.to_numeric(site.set_index("UID").loc[sorted(eligible),"WGT_TP"],errors="coerce")
 eligible={u for u,v in zip(sorted(eligible),weight) if pd.notna(v) and v>0}
huc_col=next((x for x in ["HUC8","HUC_8"] if x in site),None)
grp=site.set_index("UID").loc[sorted(eligible)] if eligible else site.iloc[[]]
gate={
 "source_all_six_UID_intersection":len(frame),"BERW_all_six_native_visits":len(wade),
 "BETB_all_six_native_visits":len(boat),"BERW_and_BETB_same_uid_overlap":len(wade&boat),
 "original_source_site_PROB_uid_n":len(prob),"original_wade_PROB_positive_weight_opportunity":len(eligible),
 "survey_site_huc_field":huc_col,
 "survey_states":int(grp.PSTL_CODE.nunique()) if "PSTL_CODE" in grp else None,
 "native_HUC8_count":int(grp[huc_col].replace("",pd.NA).nunique()) if huc_col else None,
 "survey_year_groups":dict(collections.Counter(grp.YEAR)) if "YEAR" in grp else None,
 "original_siteID_count":int(grp.SITE_ID.nunique()) if "SITE_ID" in grp else None,
 "source_eligible_native_UIDs":sorted(eligible),
 "eligible_for_numeric_if_n350_HUC200_states20":len(eligible)>=350 and
   (grp[huc_col].replace("",pd.NA).nunique()>=200 if huc_col else False) and
   (grp.PSTL_CODE.nunique()>=20 if "PSTL_CODE" in grp else False)
}
# Retrieve exact 2019 companion official measurement metadata, not substitute external assumptions.
meta={}
for key,leaf in {"chla":"nrsa1314_widepchl_meta_04292019.txt",
                 "afdm":"nrsa1314_widepbio_meta_04292019.txt",
                 "water":"nrsa1314_widechem_meta_04292019.txt"}.items():
 try:
  req=urllib.request.Request(BASE+leaf,headers={"User-Agent":"gesui-independent-2013-nrsa-source-units"})
  with urllib.request.urlopen(req,timeout=35) as r:raw=r.read(150000)
  desc=raw.decode("latin1",errors="replace")
  lines=[ln for ln in desc.splitlines() if any(z in ln.upper() for z in ("CHLA_RESULT","AFDM_RESULT","NTL_RESULT","PTL_RESULT","UNIT","AREAL"))]
  meta[key]={"official_url":BASE+leaf,"sha256":hashlib.sha256(raw).hexdigest(),
             "original_description_matching_measurement_lines":lines[:24],
             "total_metadata_bytes":len(raw)}
  print("NRSA1314_OFFICIAL_METADATA",key,json.dumps(meta[key]),flush=True)
 except Exception as e:
  meta[key]={"official_url":BASE+leaf,"status":"UNAVAILABLE","error":type(e).__name__+":"+str(e)[:150]}
  print("NRSA1314_OFFICIAL_METADATA_UNAVAILABLE",key,str(e)[:100],flush=True)
out={"schema":"NRSA1314_ORIGINAL_EPA_2019_SIX_TABLE_PRENUMERIC_SEMANTIC_CONTRACT_v1",
     "publisher_measurement_metadata_source_only":meta,
     "original_six_source_receipts_sha256":{k:receipts[k]["sha256"] for k in NAMES},
     "original_source_only_manifest_sha256":sha(pathlib.Path(os.environ["NRSA1314_ORIGINAL_SIX_SOURCE_QA"])),
     "all_source_schemas":profiles,
     "original_source_gate":gate,
     "source_wade_BERW_vs_boat_BETB_not_pooled":True,
     "no_chla_afdm_ept_TNTP_measurement_values_examined":True,
     "no_wastewater_source_or_actual_services":True}
p=ROOT/"NRSA1314_SIX_ORIGINAL_MEASUREMENT_SEMANTICS_SOURCE_ONLY.json"
p.write_text(json.dumps(out,indent=2,sort_keys=True))
print("NRSA1314_INDEPENDENT_SOURCE_SEMANTIC_GATE",json.dumps({k:v for k,v in gate.items() if not isinstance(v,list)},sort_keys=True))
