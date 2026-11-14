#!/usr/bin/env python3
"""INDEPENDENT NUMERICAL AUDIT: manual numpy normal-equation Ridge10 fits
for genuine three EPA NRSA source cycles. No sklearn / original model
helpers are used. The original native per-UID source data are reconstituted
in the same runner, never exported. Catchment-equal train weights and
train-only category vocab; no shared HUC8 in leave-cycle-out train.
"""
from __future__ import annotations
import hashlib,json,os,pathlib
import numpy as np,pandas as pd
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_threecycle_ept_transfer"
original=json.loads(pathlib.Path(os.environ["ORIGINAL_THREECYCLE_TRANSPORT_RESULT"]).read_text())
if original["schema"]!="gesui_NRSA_2008_2013_2018_REAL_EPT_external_period_AND_zero_shared_HUC8_v1":
 raise RuntimeError("Original independent 3cycle transport output schema changed")
data=pd.read_csv(ROOT/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv",
  dtype={"source_cycle":str,"source_UID":str,"site_id":str,"huc8":str,"state":str,"eco9":str})
if len(data)!=2648 or data.duplicated(["source_cycle","source_UID"]).any():
 raise RuntimeError("Original three cycle/source UID panel invalid")
CYCLES=["2008_09","2013_14","2018_19"]
NUMERIC=["lat","lon","stream_order","season_sin","season_cos","log_cond"]
HAB=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]
NP=["log_tn","log_tp"]
D={
 "B":NUMERIC,"B+NP":NUMERIC+NP,
 "B+PHab":NUMERIC+HAB,"B+PHab+NP":NUMERIC+HAB+NP}
def original_manual_train_predict(train,target,fields):
 # All columns and categories determined by train ONLY. Identical numerical
 # contract to sklearn StandardScaler + OneHotEncoder, independent matrix
 # implementation, intercept unpenalized and weights included in normal eq.
 train_num=train[fields].to_numpy(float)
 test_num=target[fields].to_numpy(float)
 mu=train_num.mean(axis=0)
 sd=train_num.std(axis=0)
 sd=np.where(sd>1e-12,sd,1.)
 st=np.array(sorted(train.state.unique()),dtype=str)
 ec=np.array(sorted(train.eco9.unique()),dtype=str)
 def vec(x,nums):
  return np.column_stack([
   (nums-mu)/sd,
   (x.state.to_numpy(str)[:,None]==st[None,:]).astype(float),
   (x.eco9.to_numpy(str)[:,None]==ec[None,:]).astype(float)])
 zz=vec(train,train_num);te=vec(target,test_num)
 n=len(train)
 sizes=train.groupby("huc8").source_UID.transform("size").to_numpy(float)
 w=(n/train.huc8.nunique())/sizes
 if abs(np.mean(w)-1)>1e-8:raise RuntimeError("Independent catchment total training weights mismatch")
 A=np.column_stack([np.ones(n),zz])
 B=np.column_stack([np.ones(len(target)),te])
 reg=np.eye(A.shape[1])*10
 reg[0,0]=0.
 y=train.ept_fraction.to_numpy(float)
 # Symmetric weighted normal equation, with unpenalized intercept, matches
 # weighted sklearn Ridge10 with train-only unweighted feature transforms.
 beta=np.linalg.solve(A.T@(w[:,None]*A)+reg,A.T@(w*y))
 return B@beta
def group_squared_errors(t,p):
 y=t.ept_fraction.to_numpy(float)
 g=t.huc8.to_numpy(str)
 u,iv=np.unique(g,return_inverse=True)
 return np.bincount(iv,weights=(y-p)**2,minlength=len(u))/np.bincount(iv,minlength=len(u))
def manual_boot_huc(t,p0,p1):
 a=group_squared_errors(t,p0);b=group_squared_errors(t,p1)
 rng=np.random.default_rng(20261114)
 idx=rng.integers(0,len(a),size=(2000,len(a)))
 ag=a[idx].mean(axis=1);bg=b[idx].mean(axis=1)
 return list(map(float,np.quantile(100*(ag-bg)/ag,[.025,.5,.975])))
checks={}
for target_cycle in CYCLES:
 orig=original["leave_one_SURVEY_cycle_out_and_remove_ALL_shared_HUC8_from_train"][target_cycle]
 if orig["gate"]!="VALID_2CYCLE_TRAIN_ZERO_TARGET_HUC8":
  raise RuntimeError("Original fit lacks primary prespecified train/test gate for "+target_cycle)
 target=data[data.source_cycle.eq(target_cycle)].reset_index(drop=True)
 train=data[~data.source_cycle.eq(target_cycle)].copy()
 train=train[~train.huc8.isin(set(target.huc8))].reset_index(drop=True)
 assert set(target.huc8).isdisjoint(set(train.huc8))
 qa=orig["QA"]
 if (len(target)!=qa["n_heldout_original_UID"] or
     len(train)!=qa["n_training_original_UID_after_strict_HUC_exclusion"] or
     train.huc8.nunique()!=qa["n_training_distinct_HUC8_after_exclusion"]):
  raise RuntimeError("Original source native UID/HUC8 zero-overlap partition mismatch for "+target_cycle)
 y=train.ept_fraction.to_numpy(float)
 size=train.groupby("huc8").source_UID.transform("size").to_numpy(float)
 w=(len(train)/train.huc8.nunique())/size
 p={"intercept":np.full(len(target),float(np.average(y,weights=w)))}
 for k,feats in D.items():
  p[k]=original_manual_train_predict(train,target,feats)
 reference=orig["transport_results"]
 expected=reference["MSE_equal_HUC8_primary"]
 manual={k:float(group_squared_errors(target,values).mean()) for k,values in p.items()}
 for k,value in manual.items():
  if abs(value-expected[k])>1e-8:
   raise RuntimeError("Independent manual NumPy weighted Ridge HUC8 MSE differs: "+
      target_cycle+"/"+k+":"+str((value,expected[k])))
 n=100*(manual["B+PHab"]-manual["B+PHab+NP"])/manual["B+PHab"]
 h=100*(manual["B+NP"]-manual["B+PHab+NP"])/manual["B+NP"]
 if abs(n-reference["extra_NP_after_PHab_equalHUC8_gain_pct"])>1e-6 or abs(h-reference["extra_PHab_after_NP_equalHUC8_gain_pct"])>1e-6:
  raise RuntimeError("Independent original incremental gain source mismatch "+target_cycle)
 cn=manual_boot_huc(target,p["B+PHab"],p["B+PHab+NP"])
 ch=manual_boot_huc(target,p["B+NP"],p["B+PHab+NP"])
 for a,b in [(cn,reference["HUC8_bootstrap2000_fixedOOF95_NP_after_PHab"]),
             (ch,reference["HUC8_bootstrap2000_fixedOOF95_PHab_after_NP"])]:
  if max(abs(v-u) for v,u in zip(a,b))>2.e-6:
   raise RuntimeError("Independent original HUC bootstrap quantiles mismatch "+target_cycle)
 # Important niche: target-only groups NEVER seen in any of the other
 # two original survey periods, not just removed from the fitted train.
 others=data[~data.source_cycle.eq(target_cycle)]
 uniquehuc=set(target.huc8)-set(others.huc8)
 isnew=target.huc8.isin(uniquehuc).to_numpy()
 novel=target.loc[isnew].reset_index(drop=True)
 if len(novel)>=150 and novel.huc8.nunique()>=100:
  pn={k:val[isnew] for k,val in p.items()}
  novel_np=100*(np.mean(group_squared_errors(novel,pn["B+PHab"]))-
           np.mean(group_squared_errors(novel,pn["B+PHab+NP"])))/np.mean(group_squared_errors(novel,pn["B+PHab"]))
  novel_hab=100*(np.mean(group_squared_errors(novel,pn["B+NP"]))-
           np.mean(group_squared_errors(novel,pn["B+PHab+NP"])))/np.mean(group_squared_errors(novel,pn["B+NP"]))
  assert abs(novel_np-orig["novel_target_huc8_not_in_either_training_source_cycle"]["extra_NP_after_PHab_equalHUC8_gain_pct"])<2e-6
  assert abs(novel_hab-orig["novel_target_huc8_not_in_either_training_source_cycle"]["extra_PHab_after_NP_equalHUC8_gain_pct"])<2e-6
  novelresult={"n_UIDs":len(novel),"HUC8":novel.huc8.nunique(),
   "gain_NP_after_PHab_pct":novel_np,"gain_PHab_after_NP_pct":novel_hab,
   "HUC8_2000_bootstrap95_NP_gain":manual_boot_huc(novel,pn["B+PHab"],pn["B+PHab+NP"]),
   "HUC8_2000_bootstrap95_PHab_gain":manual_boot_huc(novel,pn["B+NP"],pn["B+PHab+NP"])}
 else:
  novelresult={"status":"PREDECLARED_N150_HUC100_FAIL","n_UIDs":len(novel),"HUC8":novel.huc8.nunique()}
 checks[target_cycle]={
  "n_train_no_testHUC":len(train),"n_trainHUC":train.huc8.nunique(),
  "n_target":len(target),"targetHUC":target.huc8.nunique(),
  "independent_finite_matrix_normal_equation_HUC8_MSE":manual,
  "gains":{"NP_after_PHab":n,"PHab_after_NP":h},
  "independent_fixedOOF_HUC8_bootstrap95_NP_gain":cn,
  "independent_fixedOOF_HUC8_bootstrap95_PHab_gain":ch,
  "unique_target_original_HUC8_sensitivity":novelresult,
  "match_sklearn_implementation":True}
out={
 "schema":"gesui_NRSA3cycle_HUC8_equal_weight_source_2008_2013_2018_INDEPENDENT_NUMPY_AUDIT_v1",
 "original_reconstructed_3cycle_UID_total":2648,
 "original_distinct_HUC8_union":original["original_HUC8_union"],
 "prior_direct_transport_source_run":43566511111,
 "all_three_left_cycles_refitted_with_zero_shared_HUC8":True,
 "independent_normal_equation_not_sklearn":True,
 "same_original_EPA_field_samples_not_independent_field_replication":True,
 "source_taxon_2008_95pct_source_quality_preserved":True,
 "no_WWTP_NP_delivered_loads_or_ecosystem_service_claim":True,
 "per_unknown_cycle":checks}
path=ROOT/"NRSA_THREECYCLE_EPT_HUC8_MANUAL_NUMPY_AND_EXCLUSIVE_NEW_BASINS_AUDIT.json"
path.write_text(json.dumps(out,indent=2,sort_keys=True))
print("GESUI_THREE_CYCLE_EXTERNAL_EPT_MANUAL_INDEPENDENT_NUMPY_QA_SUCCESS")
print(json.dumps({k:{"n_test":v["n_target"],"HUC8_test":v["targetHUC"],
 "n_train":v["n_train_no_testHUC"],"N_gain":v["gains"]["NP_after_PHab"],
 "PHab_gain":v["gains"]["PHab_after_NP"],
 "NP_CI":v["independent_fixedOOF_HUC8_bootstrap95_NP_gain"],
 "PHab_CI":v["independent_fixedOOF_HUC8_bootstrap95_PHab_gain"],
 "newHUC":v["unique_target_original_HUC8_sensitivity"]}
 for k,v in checks.items()},sort_keys=True))
