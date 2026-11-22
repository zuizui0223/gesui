#!/usr/bin/env python3
"""Original EPA NRSA2023/24 source-only wadeable EPT+NP+six PHab opportunity.

Reads ONLY site, water chemistry, habitat environmental values and original
sample availability (not benthic taxon counts or chlorophyll numerical outcome).
Two original target cohorts: all eligible BERW and original same-UID periphyton
observation / quality source subset. Never fit ecological outcomes here.
"""
from __future__ import annotations
import collections,hashlib,json,pathlib,argparse
import numpy as np,pandas as pd

SHA={
 "sites":"f7d14cf1574d12828b3aa97bd7b08b95f22aec2115753d5b75514f8eed43d013",
 "water":"b34d3f9c4b4bdc763cb5e53a856bbab329f1646d248f223f78e9404bc4560d71",
 "habitat":"8a6094bb87f8a57cf45a35a1272b137e6a9327cefd158f4bee88b025eaf9b2f4",
 "benthos":"5290d679f7eb444e97c8d45096beb4531b450186cef190580f41f4b9430b0805",
 "chla":"aec44b58806481f5ee53ce66d61c46b4b9fad344d2da2cc8312fceb3e1fb192e"
}
PHAB=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","XWIDTH"]
GEO=["LAT_DD83","LON_DD83","STRAH_ORD","HUC8","PSTL_CODE","AG_ECO9","YEAR","DATE_COL","SITE_ID","VISIT_NO"]
def sha(path):
 h=hashlib.sha256()
 with open(path,"rb") as f:
  for b in iter(lambda:f.read(1<<20),b""):h.update(b)
 return h.hexdigest()
def read_one(path,name,cols):
 src=pd.read_csv(path,usecols=lambda col:col in set(cols+["UID"]),dtype=str,keep_default_na=False,low_memory=False)
 if "UID" not in src:raise RuntimeError(name+" no UID")
 if src.UID.duplicated().any():
  raise RuntimeError(name+" duplicate original UID")
 return src.set_index("UID",drop=False)
def health(x):
 z=pd.to_numeric(x,errors="coerce")
 return z
def date(x):
 return pd.to_datetime(x,format="mixed",errors="coerce").dt.normalize()
def correct_units(df,key,kind):
 val=health(df[key+"_RESULT"])
 unit=df[key+"_UNITS"].astype(str).str.strip()
 if kind=="P":
  accepted=unit.isin(("mg/L","mg P/L","ug/L","ug P/L"))
  factor=np.where(unit.isin(("ug/L","ug P/L")),1/1000,1.)
 else:
  accepted=unit.isin(("mg/L","mg N/L"))
  factor=np.ones(len(df))
 output=val*factor
 return output.where(accepted&(output>0))
def build(args):
 for k,pin in SHA.items():
  path=getattr(args,k)
  if sha(path)!=pin:raise RuntimeError("Official EPA original "+k+" SHA changed")
 site=read_one(args.sites,"site",GEO)
 water=read_one(args.water,"water",["SITE_ID","VISIT_NO","DATE_COL","NTL_RESULT","NTL_UNITS","NTL_NARS_FLAG","PTL_RESULT","PTL_UNITS","PTL_NARS_FLAG","COND_RESULT"])
 hab=read_one(args.habitat,"phab",["SITE_ID","VISIT_NO","DATE_COL","PROTOCOL"]+PHAB)
 # Original BenthicCount has MULTIPLE taxon rows per UID, but source-only
 # survey sample method and native keys must be unambiguous.
 b=pd.read_csv(args.benthos,usecols=["UID","SITE_ID","VISIT_NO","DATE_COL","SAMPLE_TYPE"],
               dtype=str,keep_default_na=False,low_memory=False)
 if b.duplicated(["UID","SITE_ID","VISIT_NO","DATE_COL","SAMPLE_TYPE"]).any():pass
 method_conflict=b.groupby("UID").agg({
   "SITE_ID":"nunique","VISIT_NO":"nunique","DATE_COL":"nunique","SAMPLE_TYPE":"nunique"})
 if method_conflict.gt(1).any().any():raise RuntimeError("Conflicting taxon sample metadata per UID")
 bent=b.drop_duplicates("UID").set_index("UID")
 chla=read_one(args.chla,"chla",["SITE_ID","VISIT_NO","DATE_COL","NARS_FLAG"])
 common=set(site.index)&set(water.index)&set(hab.index)&set(bent.index)
 source_info={"original_EPA_sha256":SHA,
  "all_site_water_habitat_benthic_UID_overlap":len(common),
  "original_UID_site_count":len(site),
  "original_chlorophyll_UIDs":len(chla),
  "source_only_not_benthic_abundance_or_chlorophyll_values":True}
 common_ids=sorted(common)
 d=site.loc[common_ids].copy()
 w=water.loc[common_ids];p=hab.loc[common_ids];s=bent.loc[common_ids]
 miss={}
 for name,f in (("water",w),("habitat",p),("benthic",s)):
  bad=((d.SITE_ID.to_numpy()!=f.SITE_ID.to_numpy())|
       (d.VISIT_NO.to_numpy()!=f.VISIT_NO.to_numpy())|
       (date(d.DATE_COL).to_numpy()!=date(f.DATE_COL).to_numpy()))
  miss[name]=int(bad.sum())
  if bad.any():raise RuntimeError("Publisher key or survey date mismatch "+name)
 d["method_ok"]=s.SAMPLE_TYPE.eq("BERW").to_numpy() & p.PROTOCOL.eq("WADEABLE").to_numpy()
 d["TN_mgL"]=correct_units(w,"NTL","N").to_numpy()
 d["TP_mgL"]=correct_units(w,"PTL","P").to_numpy()
 d["COND"]=health(w.COND_RESULT).to_numpy()
 d["tn_flag_ok"]=w.NTL_NARS_FLAG.astype(str).str.strip().eq("").to_numpy()
 d["tp_flag_ok"]=w.PTL_NARS_FLAG.astype(str).str.strip().eq("").to_numpy()
 for var in PHAB:d[var]=health(p[var]).to_numpy()
 d["LAT_DD83"]=health(d.LAT_DD83);d["LON_DD83"]=health(d.LON_DD83)
 d["STRAH_ORD"]=health(d.STRAH_ORD)
 d["year_float"]=health(d.YEAR)
 t=date(d.DATE_COL)
 d["date_valid"]=t.notna()
 d["season_sin"]=np.sin(2*np.pi*t.dt.dayofyear/365.25)
 d["season_cos"]=np.cos(2*np.pi*t.dt.dayofyear/365.25)
 dh=d.HUC8.astype(str).str.strip()
 # National EPA SiteInfo HUC8 is a named string identifier (sometimes H-prefixed).
 # Match old frozen 2008/2013/2018 normalizer exactly; no geographic guessing.
 hprefix=dh.str.fullmatch(r"[Hh]\d{8}")
 normalized=dh.where(~hprefix,dh.str.slice(1))
 hvalid=normalized.str.fullmatch(r"\d{7,8}")
 d["HUC8_native"]=dh
 d["HUC8"]=normalized.str.zfill(8)
 d["huc_valid"]=hvalid
 if d.loc[hvalid,"HUC8"].nunique()!=dh.loc[hvalid].nunique():
  raise RuntimeError("HUC8 native labels collapsed under exact canonicalization")
 source_info["HUC8_native_format"]={
   "H_prefix_eight_digits":int(hprefix.sum()),
   "plain_seven_or_eight_digits":int((~hprefix & dh.str.fullmatch(r"\d{7,8}")).sum()),
   "other_unrecognized":int((~hvalid).sum()),
   "first_five_original_examples":sorted(dh.unique().tolist())[:5]}
 d["geo_valid"]=d.PSTL_CODE.astype(str).str.strip().ne("")&d.AG_ECO9.astype(str).str.strip().ne("")
 d["prequal"]=d.method_ok & d.tn_flag_ok & d.tp_flag_ok & d.date_valid & d.huc_valid & d.geo_valid
 required=["TN_mgL","TP_mgL","COND","LAT_DD83","LON_DD83","STRAH_ORD",
           "year_float","season_sin","season_cos"]+PHAB
 valid=d.loc[d.prequal,required].replace([np.inf,-np.inf],np.nan)
 d["numeric_ok"]=False
 d.loc[valid.index,"numeric_ok"]=valid.notna().all(axis=1)
 for col in ("COND","XWIDTH"):
  d.loc[valid.index,"numeric_ok"]&=valid[col].gt(0)
 d["eligible"]=d.prequal & d.numeric_ok
 # CHLA source quality presence only; no pigment/other ecological value read.
 d["periphyton_record"]=d.index.isin(chla.index)
 cflag=chla.NARS_FLAG.astype(str).str.strip()
 d["periphyton_quality_available"]=d.index.isin(set(cflag[cflag.eq("")].index))
 source_info["missing_original_key_or_date_matches"]=miss
 source_info["raw_method_count"]=s.SAMPLE_TYPE.value_counts().to_dict()
 source_info["sample_quality_projections"]={key:int(d[key].sum()) for key in
   ("method_ok","tn_flag_ok","tp_flag_ok","date_valid","huc_valid","geo_valid",
    "prequal","numeric_ok","eligible","periphyton_record","periphyton_quality_available")}
 original_quality=d[d.eligible].copy()
 both=original_quality[original_quality.periphyton_quality_available]
 def meta(sub):
  return {"UID_visits":int(len(sub)),"physical_SITE_ID":int(sub.SITE_ID.nunique()),
          "HUC8":int(sub.HUC8.nunique()),"states":int(sub.PSTL_CODE.nunique()),
          "eco9":int(sub.AG_ECO9.nunique()),
          "n_source_years":{str(k):int(v) for k,v in sub.YEAR.value_counts().items()},
          "original_huc8_sha256":hashlib.sha256("\n".join(sorted(set(sub.HUC8))).encode()).hexdigest(),
          "source_go":len(sub)>=500 and sub.HUC8.nunique()>=300 and sub.PSTL_CODE.nunique()>=30}
 source_info["2023_original_all_source_valid"]=meta(original_quality)
 source_info["2023_original_periphyton_eligible_source_valid"]=meta(both)
 source_info["original_sample_method_not_pooled"]=True
 source_info["unreported_TP_never_zero"]=True
 source_info["not_sewage_derived_TN_TP"]=True
 return source_info
def main():
 parser=argparse.ArgumentParser()
 for k in SHA:parser.add_argument("--"+k,type=pathlib.Path,required=True)
 parser.add_argument("--out",type=pathlib.Path,required=True)
 x=parser.parse_args();report=build(x)
 x.out.parent.mkdir(parents=True,exist_ok=True)
 x.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print("FOURCYCLE_2023_SOURCE_ONLY",json.dumps(report,sort_keys=True),flush=True)
if __name__=="__main__":main()
