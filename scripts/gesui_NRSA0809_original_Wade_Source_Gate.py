#!/usr/bin/env python3
"""EPA2008/09 SOURCE-blind exact wadeable BERWW×CHEMW+PHab+original SiteInfo.
Lock original native UID, selected source methods, ORIGINAL CHEM NTL/PTL
reported in ug/L, site HUC source and physical XEMBED source availability.
Do not read quantitative benthic or nutrient outcome magnitudes in this gate.
"""
from __future__ import annotations
import collections,hashlib,json,os,pathlib,urllib.request,time
import pandas as pd
R=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui0809_semantics"
R.mkdir(parents=True,exist_ok=True)
orig=json.loads(pathlib.Path(os.environ["NRSA0809_SOURCE_LOCK"]).read_text())
if orig["schema"]!="gesui_NRSA0809_AUTHENTIC_EPA_WATER_BENTHOS_PHAB_MMI_SOURCE_ONLY_v1":
 raise RuntimeError("Original blinded 2008 source schema changed")
SRC=orig["original_data_sources"]
def take(k):
 url=SRC[k]["url"];p=R/(k+".csv")
 for attempt in range(3):
  try:
   req=urllib.request.Request(url,headers={"User-Agent":"gesui-EPA0809-original-numeric-blind-sample-code-contract"})
   with urllib.request.urlopen(req,timeout=120) as w,p.open("wb") as f:
    while b:=w.read(1<<20):f.write(b)
   h=hashlib.sha256(p.read_bytes()).hexdigest()
   if h!=SRC[k]["source_SHA256"]:raise RuntimeError("Original source SHA changed "+k)
   break
  except Exception:
   p.unlink(missing_ok=True)
   if attempt==2:raise
   time.sleep(2**attempt)
 data=pd.read_csv(p,dtype=str,keep_default_na=False,encoding="latin1",low_memory=False)
 p.unlink();data.columns=data.columns.str.upper().str.strip()
 return data
D={k:take(k) for k in ("siteinfo","water","phab","bentcount","bentmmi")}
H={}
for name,x in D.items():
 cols=list(x)
 H[name]=cols
 print("NRSA0809_ORIGINAL_COLS",name,json.dumps(cols),flush=True)
for name,x in D.items():
 for k in ("SAMPLE_TYPE","SAMPLE_CAT","SAM_CODE","SITE_CLASS","USE_NRSA","SITETYPE",
           "YEAR","NARS_FLAG","STRATA","ECO9","FW_ECO9","FW_ECO3"):
  if k in x:
   v=x[k].astype(str).str.strip()
   print("NRSA0809_SOURCE_CATEGORY",name,k,json.dumps(v.value_counts().head(20).to_dict()),flush=True)
# Source method labels frozen before numeric response inspection. Do not pool BELGB/BERWW or CHEMB/CHEMW.
S=D["siteinfo"].copy();W=D["water"].copy();P=D["phab"].copy();B=D["bentcount"].copy();M=D["bentmmi"].copy()
for name,x in (("siteinfo",S),("phab",P),("bentmmi",M)):
 if x.UID.duplicated().any():raise RuntimeError("Original 2008 source repeated "+name+" native UID")
# Duplicate visit source chemistry rows require disambiguation by original
# CHEMW sample and source primary, not arbitrary dedup.
bw=B[B.SAMPLE_TYPE.str.upper().str.strip().eq("BERWW")].copy()
ww=W[W.SAMPLE_TYPE.str.upper().str.strip().eq("CHEMW")].copy()
# Source sample category P primary for wadeable, not D duplicate.
if "SAMPLE_CAT" in bw:bw=bw[bw.SAMPLE_CAT.str.strip().str.upper().eq("P")]
if "SAMPLE_CAT" in ww:ww=ww[ww.SAMPLE_CAT.str.strip().str.upper().eq("P")]
print("NRSA0809_SOURCE_METHOD_LOCK",
      json.dumps({"sample_uid_berww_primary":bw.UID.nunique(),"sample_uid_CHEMW_primary":ww.UID.nunique(),
        "CHEMW_UID_repeated_rows":int(ww.UID.duplicated().sum()),
        "BERWW_native_taxon_records":len(bw),
        "sample_no_BERWW":int(B.SAMPLE_TYPE.str.strip().ne("BERWW").sum())},sort_keys=True),flush=True)
if ww.UID.duplicated().any():raise RuntimeError("More than one original CHEMW primary row per UID; source grain ambiguous")
source=set.intersection(*(set(y.UID) for y in (S,ww,P,bw,M)))
# Source population PROB and USE_NRSA if explicitly a valid primary/source field.
prob=source
if "SITE_CLASS" in S:
 prob=prob&set(S.loc[S.SITE_CLASS.str.strip().str.upper().eq("PROB"),"UID"])
# USE_NRSA categorical semantics require original field 'Y/1/T' not yet assumed.
source_eligible=prob.copy()
site=S.set_index("UID").loc[sorted(source_eligible)] if source_eligible else S.iloc[[]]
huc_candidate=next((k for k in ("HUC8","HUC_8","HUC_CODE","HUC8_OLD") if k in S.columns),None)
state_candidate=next((k for k in ("PSTL_CODE","STATE","STATE_ID") if k in S.columns),None)
gate={"original_all_five_source_UID":len(source),
      "original_primary_wadeable_probable_UID":len(source_eligible),
      "source_site_huc8_field":huc_candidate,
      "source_state_field":state_candidate,
      "original_source_HUC8":int(site[huc_candidate].replace("",pd.NA).nunique()) if huc_candidate else None,
      "original_source_states":int(site[state_candidate].replace("",pd.NA).nunique()) if state_candidate else None,
      "original_UIDs":sorted(source_eligible),
      "original_source_ecological_numeric_eligible_by_ID_only":(
        len(source_eligible)>=300 and huc_candidate is not None and state_candidate is not None and
        site[huc_candidate].replace("",pd.NA).nunique()>=150 and
        site[state_candidate].replace("",pd.NA).nunique()>=20)}
print("NRSA0809_ORIGINAL_PRENUMERIC_COHORT",json.dumps({k:v for k,v in gate.items() if not isinstance(v,list)},sort_keys=True),flush=True)
p=R/"NRSA0809_ORIGINAL_WADEABLE_BENTHIC_CHEMW_PHAB_SEMANTICS.json"
result={
 "schema":"gesui_NRSA0809_original_wadeable_CHEMW_BERWW_PHAB_MMI_pre_numeric_v1",
 "original_source_manifest_sha256":hashlib.sha256(pathlib.Path(os.environ["NRSA0809_SOURCE_LOCK"]).read_bytes()).hexdigest(),
 "original_EPA_CSVs_SHA256":{k:SRC[k]["source_SHA256"] for k in D},
 "native_csv_headers":H,
 "original_wadeable_primary_source_gate":gate,
 "no_survey_period_2008_Chla_AFDM_original_table_listed":True,
 "original_CHEMW_NTL_PTL_both_ug_L_per_EPA_metadata":True,
 "no_source_numeric_NP_or_EPT_or_MMI_effects_read":True,
 "no_2018_wadeable_protocol_claim_without_source_verification":True,
 "not_wastewater_delivery_or_ecosystem_services":True}
p.write_text(json.dumps(result,indent=2,sort_keys=True))
