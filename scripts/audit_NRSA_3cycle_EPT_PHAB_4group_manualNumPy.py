#!/usr/bin/env python3
"""Independent manual NumPy four-group PHab vs original EPA ambient TN/TP
external era transport. All original EPA source cohorts reconstructed in
runner; verify all group OOS gains and original 2000 HUC bootstrap intervals.
Independent code not independent field replication, no mechanism attribution.
"""
from __future__ import annotations
import json,os,pathlib
import numpy as np,pandas as pd
BASE=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_threecycle_ept_transfer"
orig=json.loads(pathlib.Path(os.environ["ORIGINAL_FOUR_FACTOR_RESULT"]).read_text())
if orig["schema"]!="gesui_NR​​SA_3cycle_EPT_real_PHAB_fourFactorGroup_POSTHOC_sourceOriginal_v1":
 raise RuntimeError("Source original 4-group model schema changed")
d=pd.read_csv(BASE/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv",
 dtype={"source_cycle":str,"source_UID":str,"huc8":str,"state":str,"eco9":str})
if len(d)!=2648:raise RuntimeError("Original 3 national period source frame changed")
GEO=["lat","lon","stream_order","season_sin","season_cos","log_cond"]
NP=["log_tn","log_tp"]
PHAB=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]
GROUPS={
 "substrate_embeddedness":["XEMBED"],
 "canopy_and_vegetation_cover":["XCDENMID","XPCAN","XFC_OHV"],
 "channel_width":["log_width"],
 "human_disturbance_index":["W1_HALL"]}
if GROUPS!=orig["four_groups"]:raise RuntimeError("Original four-group definitions changed")
def manual(train,test,nums):
 Xtr=train[nums].to_numpy(float);Xte=test[nums].to_numpy(float)
 m=Xtr.mean(0);sd=Xtr.std(0);sd=np.where(sd>1e-12,sd,1.)
 catstate=np.array(sorted(train.state.unique()),dtype=str)
 cateco=np.array(sorted(train.eco9.unique()),dtype=str)
 def matrix(D,A):
  return np.column_stack([
   (A-m)/sd,
   (D.state.to_numpy(str)[:,None]==catstate[None,:]).astype(float),
   (D.eco9.to_numpy(str)[:,None]==cateco[None,:]).astype(float)])
 A=np.column_stack([np.ones(len(train)),matrix(train,Xtr)])
 B=np.column_stack([np.ones(len(test)),matrix(test,Xte)])
 counts=train.groupby("huc8").source_UID.transform("size").to_numpy(float)
 w=(len(train)/train.huc8.nunique())/counts
 R=np.eye(A.shape[1])*10;R[0,0]=0
 beta=np.linalg.solve(A.T@(w[:,None]*A)+R,A.T@(w*train.ept_fraction.to_numpy(float)))
 return B@beta
def loss(target,p):
 u,inv=np.unique(target.huc8.to_numpy(str),return_inverse=True)
 y=target.ept_fraction.to_numpy(float)
 return np.bincount(inv,weights=(y-p)**2,minlength=len(u))/np.bincount(inv,minlength=len(u))
def compare(target,a,b):
 L=loss(target,a);N=loss(target,b)
 imp=100*(L.mean()-N.mean())/L.mean()
 rng=np.random.default_rng(20261114)
 samples=rng.integers(0,len(L),size=(2000,len(L)))
 bas=L[samples].mean(1);alt=N[samples].mean(1)
 boot=list(map(float,np.quantile(100*(bas-alt)/bas,[.025,.5,.975])))
 return float(imp),boot
aud={}
for cycle in ["2008_09","2013_14","2018_19"]:
 test=d[d.source_cycle.eq(cycle)].reset_index(drop=True)
 train=d[~d.source_cycle.eq(cycle)]
 train=train[~train.huc8.isin(set(test.huc8))].reset_index(drop=True)
 if set(test.huc8)&set(train.huc8):raise RuntimeError("Independent source 3cycle train HUC leakage")
 expected=orig["results"][cycle]
 base=manual(train,test,GEO+NP)
 full=manual(train,test,GEO+NP+PHAB)
 if abs(loss(test,base).mean()-expected["base_geo_NP_MSE"])>2e-8 or abs(loss(test,full).mean()-expected["full_geo_NP_PHab_MSE"])>2e-8:
  raise RuntimeError("Independent four-group base/full source original 3cycle HUC mismatch "+cycle)
 res={}
 for label,grp in GROUPS.items():
  alone=manual(train,test,GEO+NP+grp)
  remaining=[x for x in PHAB if x not in grp]
  omit=manual(train,test,GEO+NP+remaining)
  a0,ci0=compare(test,base,alone)
  a1,ci1=compare(test,omit,full)
  previous=expected["source_defined_four_factor_groups"][label]
  for gain,boot,pref in [
   (a0,ci0,previous["standalone_added_to_geo_and_NP"]),
   (a1,ci1,previous["conditional_added_back_after_other_three_groups"])]:
   if abs(gain-pref["HUC8_equalMSE_gain_pct"])>2e-6:
    raise RuntimeError("Independent all four source group model gain mismatch "+cycle+"/"+label)
   if max(abs(x-y) for x,y in zip(boot,pref["fixedOOF_HUC8_2000_bootstrap95"]))>2e-6:
    raise RuntimeError("Independent all four source group HUC bootstrap mismatch "+cycle+"/"+label)
  res[label]={"standalone_gain":a0,"standalone_95":ci0,
              "conditional_gain":a1,"conditional_95":ci1}
 aud[cycle]={"n_target":len(test),"HUC8_target":test.huc8.nunique(),
             "n_train":len(train),"HUC8_train":train.huc8.nunique(),
             "all_4_groups_independent_matched":True,"groups":res}
out={
 "schema":"gesui_NR​​SA_three_era_PHAB_fourGroup_SOURCE_independentNumPyQA_v1",
 "all_eight_group_tests_each_of_three_periods_numerically_reproduced":True,
 "same_original_2648_EPA_source_native_UIDs":True,
 "zero_shared_HUC8_any_heldout_cycle":True,
 "independent_manual_normal_equation_not_original_sklearn_helper":True,
 "effect_sign_not_causal_substrate_or_sewage":True,
 "results":aud}
p=BASE/"NRSA_THREECYCLE_PHAB_FOURGROUP_INDEPENDENT_NUMPY_QA.json"
p.write_text(json.dumps(out,sort_keys=True,indent=2))
print("GESUI_THREECYCLE_PHAB_FOURGROUP_INDEPENDENT_NUMPY_SUCCESS")
print(json.dumps({
 c:{k:{"alone":v["standalone_gain"],"alone95":v["standalone_95"],
        "conditional":v["conditional_gain"],"conditional95":v["conditional_95"]}
        for k,v in data["groups"].items()} for c,data in aud.items()},sort_keys=True))
