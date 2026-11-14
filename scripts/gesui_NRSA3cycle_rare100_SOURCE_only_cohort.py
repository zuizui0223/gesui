#!/usr/bin/env python3
"""Source-only m100 census from authentic EPA three-cycle TOTAL300 count
cohorts pre-audited in PR82. No new taxon-richness scores or heldout effects.
Same m200 IDs and expanded n>=100 IDs are explicitly separate populations.
"""
from __future__ import annotations
import collections,hashlib,json,os,pathlib
import numpy as np,pandas as pd

T=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))
R=T/"gesui_threecycle_rarefied_taxa"
old_counts=pd.read_csv(R/"AUDIT_ONLY_ORIGINAL_EPA_3CYCLE_TAXON_COUNTS.csv",
 dtype={"cycle":str,"UID":str,"key":str,"order":str},keep_default_na=False)
full2648=pd.read_csv(T/"gesui_threecycle_ept_transfer"/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv",
 dtype={"source_cycle":str,"source_UID":str,"huc8":str,"state":str,"eco9":str})
if len(full2648)!=2648:raise RuntimeError("Prior source pop changed")
old_counts["count"]=pd.to_numeric(old_counts["count"],errors="raise")
old_counts["known"]=old_counts.known.astype(str).str.lower().isin(("true","1","yes"))
g=old_counts.groupby(["cycle","UID"],sort=True)
q=g.apply(lambda x:pd.Series({
 "N":float(x["count"].sum()),
 "known_frac":float(x.loc[x.known,"count"].sum()/x["count"].sum())}),include_groups=False).reset_index()
if len(q)!=2648:raise RuntimeError("Source counts don't span original all 2648 original UIDs")
q=q.rename(columns={"cycle":"source_cycle","UID":"source_UID"})
x=full2648.merge(q,on=["source_cycle","source_UID"],validate="one_to_one")
# Reconstruct the unchanged original source-only m200 eligibility directly
# from signed original TOTAL300 counts; no new m200 numerical outcome used.
orig2209=x[(x.N>=200)&(x.known_frac>=.95-1e-10)].copy()
if len(orig2209)!=2209:raise RuntimeError("Original source-only n2209 m200 count gate changed")
for name,n in {"2008_09":476,"2013_14":957,"2018_19":776}.items():
 if orig2209.source_cycle.eq(name).sum()!=n:raise RuntimeError("Original source-m200 by cycle count changed "+name)
source_only={}
for era in ["2008_09","2013_14","2018_19"]:
 v=x[x.source_cycle.eq(era)].copy()
 M100=v[(v.N>=100)&(v.known_frac>=.95-1e-10)].copy()
 M200=orig2209[orig2209.source_cycle.eq(era)].copy()
 sets100=set(M100.source_UID);sets200=set(M200.source_UID)
 if not sets200<=sets100:raise RuntimeError("m200 original UID must always pass source m100 gate "+era)
 add=M100[~M100.source_UID.isin(sets200)].copy()
 source_only[era]={
  "original_full_n":len(v),"original_full_HUC8":v.huc8.nunique(),
  "original_m200_site_n":len(M200),"original_m200_HUC8":M200.huc8.nunique(),
  "same2209_m100_sensitivity_n":len(M200),
  "m100_expanded_source_n":len(M100),"m100_expanded_HUC8":M100.huc8.nunique(),
  "m100_expanded_states":M100.state.nunique(),
  "new_source_UID_100_not200_n":len(add),
  "new_source_UID_100_not200_HUC8":add.huc8.nunique(),
  "source_new_UID_N_quantiles":{str(z):float(a) for z,a in add.N.quantile([0,.25,.5,.75,1]).items()},
  "source_m100_fraction_of_full_visits":float(len(M100)/len(v)),
  "source_m200_subset_of_m100":True,
  "source_mean_N_m100":float(M100.N.mean()),
  "source_mean_N_same_m200":float(M200.merge(q,on=["source_cycle","source_UID"]).N.mean()),
  "m100_source_gate_target_min250_HUC150_states20":bool(len(M100)>=250 and M100.huc8.nunique()>=150 and M100.state.nunique()>=20),
  "m100_selected_original_native_UID_sha256":hashlib.sha256("\n".join(sorted(sets100)).encode()).hexdigest()}
 # Eligibility fixed regardless of model performance; publish IDs for join
 # only in runner ephemeral, report source dataset manifest contains no values.
 if len(M100)<250 or M100.huc8.nunique()<150 or M100.state.nunique()<20:
  raise RuntimeError("Pre-frozen original source m100 target eligibility failed "+era)
 print("ORIGINAL_EPA_M100_SOURCE_ONLY_CENSUS",era,json.dumps(source_only[era],sort_keys=True),flush=True)
O={
 "schema":"gesui_EPA_3cycle_NRSA_m100_vs_m200_sourceONLY_two_population_frames_v1",
 "original_source_native_UID_n":len(full2648),
 "original_m200_native_UID_n":len(orig2209),
 "original_source_counts_reproduced_same_genuine_BERW_BERWW":True,
 "m100_N_min100_taxonomy95_no_false_zero":True,
 "data_only_no_new_rarefied_richness_values_computed":True,
 "no_EPT_external_effect_models_fit":True,
 "source_frames":source_only,
 "no_WWTP_source_loads_or_benthic_biomass_services":True}
path=R/"NRSA_3ERA_RARE100_SOURCE_ONLY_COHORT_AND_ADDED_N_QA.json"
path.write_text(json.dumps(O,sort_keys=True,indent=2))
