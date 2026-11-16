#!/usr/bin/env python3
"""Source-first site-level EPA original EPT FFG/HABIT/PTV trait eligibility
across three original national NRSA survey cycles. Source identity audit
done in prior original EPA source-only workflow, no environmental effects
modeled here. Taxon traits are inherited autecology codes not actual feeding.
"""
from __future__ import annotations
import hashlib,json,os,pathlib,re
import numpy as np,pandas as pd

T=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))
R=T/"gesui_threecycle_traits"
R.mkdir(exist_ok=True,parents=True)
M=T/"gesui_threecycle_rarefied_taxa"/"AUDIT_ONLY_ORIGINAL_EPA_3CYCLE_TAXON_COUNTS.csv"
V=T/"gesui_threecycle_ept_transfer"/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv"
MAP=R/"AUDIT_ONLY_EPA_ORIGINAL_3CYCLE_TAXA_FFG_HABIT_PTV.csv"
TQA=json.loads(pathlib.Path(os.environ["FROZEN_TRUE_EPA_TRAITS_SOURCE_QA"]).read_text())
if TQA["schema"]!="gesui_NR​​SA_three_orig_EPA_autecology_trait_ID_name_ORDER_coverage_SOURCE_ONLY_v1" or not all(TQA["predeclared_trait_EPT_70pct_individual_coverage_gate"].values()):
 raise RuntimeError("Prior original EPA official three-era trait taxon source gate not passed")
count=pd.read_csv(M,dtype={"cycle":str,"UID":str,"key":str,"order":str},keep_default_na=False)
full=pd.read_csv(V,dtype={"source_cycle":str,"source_UID":str,"huc8":str,"state":str,"eco9":str})
trait=pd.read_csv(MAP,dtype=str,keep_default_na=False)
if len(full)!=2648 or full.duplicated(["source_cycle","source_UID"]).any():
 raise RuntimeError("Original three-era 2648 native source visits changed")
if trait.duplicated(["source_cycle","key"]).any():
 raise RuntimeError("Original official TAXA_ID×cycle has ambiguous FFG/HABIT/PTV")
if count.duplicated(["cycle","UID","key"]).any():raise RuntimeError("Original positive taxon key duplicate")
count["count"]=pd.to_numeric(count["count"],errors="raise")
if count["count"].le(0).any():raise RuntimeError("Original BERW nonpositive counted data")
count["known"]=count.known.astype(str).str.lower().isin(("true","yes","1"))
# Original m200 subset from ORIGINAL source N/known flags, not response sign.
stats=count.groupby(["cycle","UID"],sort=True).apply(
 lambda d:pd.Series({"N":float(d["count"].sum()),
                      "classified":float(d.loc[d.known,"count"].sum()/d["count"].sum())}),
 include_groups=False).reset_index().rename(columns={"cycle":"source_cycle","UID":"source_UID"})
valid=stats[(stats.N>=200)&(stats.classified>=.95-1e-10)]
if len(valid)!=2209:raise RuntimeError("Original m200 site source eligibility differs")
full=full.merge(valid,on=["source_cycle","source_UID"],how="inner",validate="one_to_one")
count=count.merge(valid[["source_cycle","source_UID"]].rename(
 columns={"source_cycle":"cycle","source_UID":"UID"}),on=["cycle","UID"],validate="many_to_one")
E=count[count.order.isin(("EPHEMEROPTERA","PLECOPTERA","TRICHOPTERA"))].copy()
joined=E.merge(trait[["source_cycle","key","order","FFG_norm","HABIT_norm","PTV_norm"]],
 left_on=["cycle","key","order"],right_on=["source_cycle","key","order"],how="left",validate="many_to_one")
if len(joined)!=len(E):raise RuntimeError("EPA original source trait join expanded or lost taxon rows")
def fractional(code,wanted,allowed):
 c=str(code).strip().upper()
 if not c:return None
 vals=[v.strip() for v in c.split(",") if v.strip()]
 if not vals:return None
 if any(v not in allowed for v in vals):return None
 # Duplicated mode labels are NOT extra trait roles, avoid doubled counts.
 unique=list(dict.fromkeys(vals))
 return float(wanted in unique)/len(unique)
FG={"CF","CG","PA","PI","PR","SC","SH"}
HAB={"BU","CB","CN","DV","PK","SK","SP","SW"}
joined["scraper_weight"]=joined.FFG_norm.map(lambda s:fractional(s,"SC",FG))
joined["clinger_weight"]=joined.HABIT_norm.map(lambda s:fractional(s,"CN",HAB))
joined["ptv_numeric"]=pd.to_numeric(joined.PTV_norm,errors="coerce")
if joined.ptv_numeric.dropna().lt(0).any() or joined.ptv_numeric.dropna().gt(10).any():
 raise RuntimeError("Original EPA PTV source outside 0..10, no field coercion")
answers=[]
for (cycle,uid),d in joined.groupby(["cycle","UID"],sort=True):
 n=d["count"].to_numpy(float)
 total=float(n.sum())
 if total<=0:continue
 z={"source_cycle":cycle,"source_UID":uid,
    "original_EPT_counted_individuals":total}
 for field,name in [("scraper_weight","EPT_scraper_FFG_fraction"),
                    ("clinger_weight","EPT_clinger_HABIT_fraction"),
                    ("ptv_numeric","EPT_mean_pollution_tolerance_PTV")]:
  v=d[field].to_numpy(float)
  good=np.isfinite(v)
  denom=float(n[good].sum())
  coverage=denom/total
  z[name+"_original_EPT_counted_coverage"]=coverage
  z[name]=float(np.dot(v[good],n[good])/denom) if denom>0 else np.nan
 answers.append(z)
source_site=pd.DataFrame(answers)
data=full.merge(source_site,on=["source_cycle","source_UID"],how="left",validate="one_to_one")
if len(data)!=2209:raise RuntimeError("Original m200 source sites duplicated in trait join")
outcome={
 "EPT_scraper_FFG_fraction",
 "EPT_clinger_HABIT_fraction",
 "EPT_mean_pollution_tolerance_PTV"}
for name in outcome:
 data.loc[data[name+"_original_EPT_counted_coverage"].lt(.8-1e-10),name]=np.nan
 if name!="EPT_mean_pollution_tolerance_PTV":
  if data[name].dropna().lt(0).any() or data[name].dropna().gt(1).any():
   raise RuntimeError("Fractional FFG/HABIT must fall within [0,1]")
 else:
  if data[name].dropna().lt(0).any() or data[name].dropna().gt(10).any():
   raise RuntimeError("PTV mean should be 0–10")
QA={}
for era in ("2008_09","2013_14","2018_19"):
 d=data[data.source_cycle.eq(era)]
 info={"source_m200_UID":len(d),
       "source_m200_HUC8":d.huc8.nunique(),"sites_with_any_EPT_count":int(d.original_EPT_counted_individuals.notna().sum()),
       "EPT_zero_count_sites_with_missing_not_zero_trait":int(d.original_EPT_counted_individuals.isna().sum()),
       "native_EPT_count_median_if_positive":float(d.original_EPT_counted_individuals.median()),
       "original_explicit_trait_80pct_or_more_site_eligibility":{}}
 for name in sorted(outcome):
  e=d[d[name].notna()]
  info["original_explicit_trait_80pct_or_more_site_eligibility"][name]={
    "n_native_UID":len(e),"n_HUC8":e.huc8.nunique(),"states":e.state.nunique(),
    "source_mean_response_only_not_effect":float(e[name].mean()) if len(e) else None,
    "source_min":float(e[name].min()) if len(e) else None,
    "source_max":float(e[name].max()) if len(e) else None,
    "original_UID_SET_SHA256":hashlib.sha256("\n".join(sorted(e.source_UID)).encode()).hexdigest(),
    "source_predefined_target_n250_HUC150_state20_gate":bool(len(e)>=250 and e.huc8.nunique()>=150 and e.state.nunique()>=20)}
 QA[era]=info
 print("EPA_ORIGINAL_3CYCLE_EPT_FUNCTIONAL_TRAITS_SITES_SOURCE_GATE",era,json.dumps(info,sort_keys=True),flush=True)
data.to_csv(R/"AUDIT_ONLY_3CYCLE_ORIGINAL_EPT_FFG_HABIT_PTV_SAMESITE.csv",index=False)
O={
 "schema":"gesui_EPA3cycle_SOURCEONLY_ORIGINAL_EPT_SC_CN_PTV_nativeUID_true_coverage_v1",
 "original_EPA_TAXA_ID_name_ORDER_FFG_HABIT_PTV_gate":TQA["predeclared_trait_EPT_70pct_individual_coverage_gate"],
 "official_EPA_original_trait_source_QA_SHA256":hashlib.sha256(pathlib.Path(os.environ["FROZEN_TRUE_EPA_TRAITS_SOURCE_QA"]).read_bytes()).hexdigest(),
 "source_original_2209_same_count_sites_as_PR82":True,
 "source_original_sites_by_cycle_and_trait":QA,
 "original_multi_FFG_HABIT_modes_split_uniformly_not_observed_diet":True,
 "zero_EPT_source_sites_no_trait_ratio_imputed":True,
 "per_site_real_EPT_trait_coverage_at_least80pct":True,
 "no_external_prediction_model_fitted":True,
 "source_instream_TNTP_not_WWTP":True}
(R/"NRSA_3CYCLE_ORIGINAL_EPT_FFG_HABIT_PTV_SOURCE_ONLY_SITES.json").write_text(json.dumps(O,indent=2,sort_keys=True))
