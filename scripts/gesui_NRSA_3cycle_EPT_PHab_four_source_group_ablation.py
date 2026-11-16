#!/usr/bin/env python3
"""Posthoc explanatory PHab four source-defined factor-group ablations on
THREE real EPA NRSA cycles with zero HUC8 overlap train vs target.
No source-level feature or site selection after observing any group effect.
No fitted biological mechanism, WWTP source or ecosystem service causality.
"""
from __future__ import annotations
import json,os,pathlib
import numpy as np,pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.linear_model import Ridge
BASE=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_threecycle_ept_transfer"
OLD=json.loads(pathlib.Path(os.environ["ORIGINAL_THREECYCLE_TRANSPORT_RESULT"]).read_text())
if OLD["schema"]!="gesui_NRSA_2008_2013_2018_REAL_EPT_external_period_AND_zero_shared_HUC8_v1":
 raise RuntimeError("Original three-cycle full-PHab transport prior source not independently confirmed")
d=pd.read_csv(BASE/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv",
  dtype={"source_cycle":str,"source_UID":str,"huc8":str,"state":str,"eco9":str})
if len(d)!=2648:raise RuntimeError("Original source panel changed")
CYCLES=["2008_09","2013_14","2018_19"]
GEO=["lat","lon","stream_order","season_sin","season_cos","log_cond"]
NUTRIENT=["log_tn","log_tp"]
GROUPS={
 "substrate_embeddedness":["XEMBED"],
 "canopy_and_vegetation_cover":["XCDENMID","XPCAN","XFC_OHV"],
 "channel_width":["log_width"],
 "human_disturbance_index":["W1_HALL"]}
if sorted({c for ll in GROUPS.values() for c in ll})!=sorted(["XEMBED","XCDENMID","XPCAN","XFC_OHV","log_width","W1_HALL"]):
 raise RuntimeError("Original approved PHab factor group partition changed")
PHAB=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]
basecols=GEO+NUTRIENT
def fit(train,test,numeric):
 prep=ColumnTransformer(transformers=[
    ("num",StandardScaler(),numeric),
    ("cat",OneHotEncoder(handle_unknown="ignore",sparse_output=False,dtype=float),["state","eco9"])],
    remainder="drop",sparse_threshold=0)
 Xtrain=prep.fit_transform(train);Xtest=prep.transform(test)
 count=train.groupby("huc8").source_UID.transform("size").to_numpy(float)
 w=(len(train)/train.huc8.nunique())/count
 yp=Ridge(alpha=10,solver="svd",fit_intercept=True).fit(
   Xtrain,train.ept_fraction.to_numpy(float),sample_weight=w).predict(Xtest)
 if not np.isfinite(yp).all():raise RuntimeError("Ablation OOS source predictions not finite")
 return yp
def err(target,p):
 hu,iv=np.unique(target.huc8.to_numpy(str),return_inverse=True)
 y=target.ept_fraction.to_numpy(float)
 return np.bincount(iv,weights=(y-p)**2,minlength=len(hu))/np.bincount(iv,minlength=len(hu))
def compare(target,old,new):
 a=err(target,old);b=err(target,new)
 improvement=float(100*(a.mean()-b.mean())/a.mean())
 rng=np.random.default_rng(20261114)
 sample=rng.integers(0,len(a),size=(2000,len(a)))
 bs=100*(a[sample].mean(axis=1)-b[sample].mean(axis=1))/a[sample].mean(axis=1)
 return {"HUC8_equalMSE_gain_pct":improvement,
         "fixedOOF_HUC8_2000_bootstrap95":list(map(float,np.quantile(bs,[.025,.5,.975])))}
result={}
for cycle in CYCLES:
 target=d[d.source_cycle.eq(cycle)].reset_index(drop=True)
 base=d[~d.source_cycle.eq(cycle)].copy()
 train=base[~base.huc8.isin(set(target.huc8))].reset_index(drop=True)
 if set(train.huc8)&set(target.huc8):raise RuntimeError("HUC8 source leakage in four-group ablation")
 prior=OLD["leave_one_SURVEY_cycle_out_and_remove_ALL_shared_HUC8_from_train"][cycle]
 if prior["gate"]!="VALID_2CYCLE_TRAIN_ZERO_TARGET_HUC8" or len(train)!=prior["QA"]["n_training_original_UID_after_strict_HUC_exclusion"]:
  raise RuntimeError("Original source zero-HUC8 train cohort not preserved "+cycle)
 pbase=fit(train,target,basecols)
 pfull=fit(train,target,basecols+PHAB)
 expected=prior["transport_results"]["MSE_equal_HUC8_primary"]
 if abs(err(target,pbase).mean()-expected["B+NP"])>1e-9 or abs(err(target,pfull).mean()-expected["B+PHab+NP"])>1e-9:
  raise RuntimeError("Subgroup model changes original approved full-PHab and base MSE before studying new effects "+cycle)
 group={}
 for factor,cols in GROUPS.items():
  solo=fit(train,target,basecols+cols)
  rest=[c for c in PHAB if c not in cols]
  without=fit(train,target,basecols+rest)
  group[factor]={
   "source_predefined_original_fields":cols,
   "standalone_added_to_geo_and_NP":compare(target,pbase,solo),
   "conditional_added_back_after_other_three_groups":compare(target,without,pfull),
   "standalone_exact_HUC8_equalMSE":float(err(target,solo).mean()),
   "full_model_without_this_group_exact_HUC8_equalMSE":float(err(target,without).mean())}
 result[cycle]={
  "n_target_UIDs":len(target),"target_HUC8":target.huc8.nunique(),
  "n_train_UIDs":len(train),"train_HUC8":train.huc8.nunique(),
  "base_geo_NP_MSE":float(err(target,pbase).mean()),
  "full_geo_NP_PHab_MSE":float(err(target,pfull).mean()),
  "all_six_PHab_total_predictive_gain_vs_geoNP":compare(target,pbase,pfull),
  "source_defined_four_factor_groups":group}
out={
 "schema":"gesui_NR​​SA_3cycle_EPT_real_PHAB_fourFactorGroup_POSTHOC_sourceOriginal_v1",
 "original_regression_source_nUID":2648,
 "source_3cycle_external_validated_run":43566511111,
 "source_native_HUC8_canonicalization_used_without_new_outcome_selection":True,
 "all_4factor_groups_were_frozen_after_total_PHab_result_but_before_group_tests":True,
 "four_groups":GROUPS,
 "not_independent_experiment_or_causal_ablation":True,
 "not_WWTP_nutrients_or_primary_productivity":True,
 "nontarget_dummy_source_onehot_train_only":True,
 "within_and_across_original_3cycle_all_HUC8_disjoint":True,
 "results":result}
p=BASE/"NRSA_THREECYCLE_EPT_PHAB_FOUR_FACTOR_GROUP_ABLATION.json"
p.write_text(json.dumps(out,indent=2,sort_keys=True))
print("GESUI_THREECYCLE_PHAB_FOURGROUP_SOURCE_ABLATION_RESULT")
print(json.dumps({c:{"n_target":s["n_target_UIDs"],"HUC8":s["target_HUC8"],
                 "fullPHab":s["all_six_PHab_total_predictive_gain_vs_geoNP"],
                 "groups":{name:{"alone":r["standalone_added_to_geo_and_NP"],
                                  "conditional":r["conditional_added_back_after_other_three_groups"]}
                    for name,r in s["source_defined_four_factor_groups"].items()}}
 for c,s in result.items()},sort_keys=True))
