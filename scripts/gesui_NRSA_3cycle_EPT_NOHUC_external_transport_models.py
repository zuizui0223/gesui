#!/usr/bin/env python3
"""Three-cycle authentic EPA NRSA 2008/2013/2018 EPT share prediction
TRANSPORT across a completely held-out survey era AND zero HUC8 overlap.

HUC8-equal-weight training and evaluation is primary; within-survey spatial
CV, site-equal fit and original 2008 unknown-order upper share are sensitivity
checks. All source signed in previous EPA source-only and numerical audits.
NO newly chosen ecological outcomes, mass loads, causal effects or services.
"""
from __future__ import annotations
import collections,hashlib,json,os,pathlib
import numpy as np,pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.linear_model import Ridge

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_threecycle_ept_transfer"
source=json.loads((ROOT/"NRSA_THREECYCLE_ORIGINAL_SOURCES_AND_HUC8_OVERLAP_QA.json").read_text())
if source["schema"]!="gesui_NR​​SA_3source_cycles_REAL_BERW_EPT_ambientNP_six_PHAB_source_only_v1":
 raise RuntimeError("Original pre-prediction 3 cycle native source /HUC8 frame was not independently audited")
src=pd.read_csv(ROOT/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv",
 dtype={"source_UID":str,"site_id":str,"huc8":str,"state":str,"eco9":str,"source_cycle":str})
source_era=["2008_09","2013_14","2018_19"]
expect={"2008_09":(586,437,48),"2013_14":(1116,706,48),"2018_19":(946,603,47)}
for name,(a,h,s) in expect.items():
 d=src[src.source_cycle.eq(name)]
 if len(d)!=a or d.huc8.nunique()!=h or d.state.nunique()!=s:
  raise RuntimeError("Pre-source-frozen cohort changed "+name)
if len(src)!=2648 or src.duplicated(["source_cycle","source_UID"]).any():
 raise RuntimeError("Source native UIDs not unique within real EPA cycle")
cats=["state","eco9"]
base=["lat","lon","stream_order","season_sin","season_cos","log_cond"]
phys=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]
npv=["log_tn","log_tp"]
models={"B":base,"B+NP":base+npv,"B+PHab":base+phys,"B+PHab+NP":base+phys+npv}
if any(x not in src.columns for v in models.values() for x in v):
 raise RuntimeError("Original source common 3-cycle EPT PHab N/P fields changed")
if not np.isfinite(src[base+phys+npv+["ept_fraction"]].to_numpy(float)).all():
 raise RuntimeError("Numeric source quality missing in previously frozen 2648 survey UID")
# Source data categories have never been mapped by ecological outcome sign:
# train only learns any categorical vocabulary; unseen target gets all zeros.
def weights_equal_huc(train):
 cnt=train.huc8.value_counts()
 n=len(train);g=len(cnt)
 w=(n/g)/train.huc8.map(cnt)
 if not np.isfinite(w).all() or abs(float(w.mean())-1)>1e-9:
  raise RuntimeError("Catchment-equal source training weights invalid")
 return w.to_numpy(float)
def fit_predict(train,test,features,response_col="ept_fraction",catchment_weighted=True):
 prep=ColumnTransformer(
   transformers=[
    ("continuous",StandardScaler(),features),
    ("categorical",OneHotEncoder(handle_unknown="ignore",sparse_output=False,dtype=float),cats)],
   remainder="drop",sparse_threshold=0)
 trainX=prep.fit_transform(train)
 testX=prep.transform(test)
 if not np.isfinite(trainX).all() or not np.isfinite(testX).all():
  raise RuntimeError("Source-normalized train/test covariates invalid")
 y=train[response_col].to_numpy(float)
 w=weights_equal_huc(train) if catchment_weighted else None
 fit=Ridge(alpha=10,solver="svd",fit_intercept=True)
 fit.fit(trainX,y,sample_weight=w)
 pred=fit.predict(testX)
 if not np.isfinite(pred).all():
  raise RuntimeError("Nonfinite heldout external-era original EPT predictions")
 return pred
def preds(train,test,response_col="ept_fraction",weighted=True):
 out={}
 y=train[response_col].to_numpy(float)
 out["intercept"]=np.full(len(test),float(np.average(y,weights=weights_equal_huc(train))) if weighted else float(y.mean()))
 for key,features in models.items():
  out[key]=fit_predict(train,test,features,response_col,weighted)
 return out
def group_errors(D,ypred,y=None):
 Y=D.ept_fraction.to_numpy(float) if y is None else np.asarray(y,dtype=float)
 g=D.huc8.to_numpy(str)
 u,inv=np.unique(g,return_inverse=True)
 count=np.bincount(inv,minlength=len(u))
 err=np.bincount(inv,weights=(Y-ypred)**2,minlength=len(u))/count
 return {"huc_ids":u,"huc_squared_errors":err,"equal_huc_MSE":float(err.mean()),
         "equal_UID_MSE":float(np.mean((Y-ypred)**2))}
def gain(mses,b,a):
 return float(100*(mses[b]-mses[a])/mses[b])
def bootstrap(target,allpred,b,a,y=None):
 yy=target.ept_fraction.to_numpy(float) if y is None else np.asarray(y,dtype=float)
 eb=group_errors(target,allpred[b],yy)["huc_squared_errors"]
 ea=group_errors(target,allpred[a],yy)["huc_squared_errors"]
 if len(eb)<20:raise RuntimeError("Too few HUC8 blocks for bootstrap")
 rng=np.random.default_rng(20261114)
 reps=rng.integers(0,len(eb),size=(2000,len(eb)))
 mb=eb[reps].mean(axis=1);ma=ea[reps].mean(axis=1)
 vals=100*(mb-ma)/mb
 return list(map(float,np.quantile(vals,[.025,.5,.975])))
def score(target,P,compare=True,y=None):
 yy=target.ept_fraction.to_numpy(float) if y is None else np.asarray(y,dtype=float)
 cols=list(P)
 hucs={}
 visits={}
 for k,p in P.items():
  e=group_errors(target,p,yy)
  hucs[k]=e["equal_huc_MSE"]
  visits[k]=e["equal_UID_MSE"]
 Y=yy
 g=target.huc8.to_numpy(str)
 mean_by_huc=pd.DataFrame({"huc8":g,"y":Y}).groupby("huc8",sort=True).y.mean()
 population_mean=float(mean_by_huc.mean())
 tgtvar=group_errors(target,np.full(len(target),population_mean),Y)["equal_huc_MSE"]
 outcome={"n_UID_visits":len(target),"n_HUC8":target.huc8.nunique(),
  "source_mean_EPT_fraction_equal_HUC8":population_mean,
  "original_target_equalHUC8_variance_of_EPT_vs_own_mean":tgtvar,
  "MSE_equal_HUC8_primary":hucs,"MSE_equal_UID_secondary":visits,
  "relative_to_target_variance_descriptive_R2_equal_HUC8":
    {k:float(1-v/tgtvar) if tgtvar>0 else None for k,v in hucs.items()},
  "extra_NP_after_PHab_equalHUC8_gain_pct":gain(hucs,"B+PHab","B+PHab+NP"),
  "extra_PHab_after_NP_equalHUC8_gain_pct":gain(hucs,"B+NP","B+PHab+NP"),
  "extra_NP_after_PHab_equalUID_gain_pct":gain(visits,"B+PHab","B+PHab+NP"),
  "extra_PHab_after_NP_equalUID_gain_pct":gain(visits,"B+NP","B+PHab+NP")}
 if compare:
  outcome["HUC8_bootstrap2000_fixedOOF95_NP_after_PHab"]=bootstrap(target,P,"B+PHab","B+PHab+NP",Y)
  outcome["HUC8_bootstrap2000_fixedOOF95_PHab_after_NP"]=bootstrap(target,P,"B+NP","B+PHab+NP",Y)
 return outcome
def sufficient(n,g,minimum_n,minimum_g):
 return n>=minimum_n and g>=minimum_g
all_groups={name:set(src.loc[src.source_cycle.eq(name),"huc8"]) for name in source_era}
LOCO={}
for targetera in source_era:
 target=src[src.source_cycle.eq(targetera)].reset_index(drop=True).copy()
 targetHUC=all_groups[targetera]
 oth=src[~src.source_cycle.eq(targetera)].reset_index(drop=True).copy()
 excluded=oth.huc8.isin(targetHUC)
 train=oth.loc[~excluded].reset_index(drop=True).copy()
 overlap_exact=set(train.huc8)&set(target.huc8)
 if overlap_exact:raise RuntimeError("Training HUC8 leakage into unknown whole-cycle test "+targetera)
 qa={
  "leftout_cycle":targetera,"n_heldout_original_UID":len(target),
  "heldout_original_HUC8":len(targetHUC),
  "n_training_OTHER_2_period_UID_prior_HUC_exclusion":len(oth),
  "n_training_source_UID_removed_for_any_target_HUC8":int(excluded.sum()),
  "n_training_HUC8_duplicated_across_source_periods_before_exclusion":oth.huc8.nunique(),
  "n_training_original_UID_after_strict_HUC_exclusion":len(train),
  "n_training_distinct_HUC8_after_exclusion":train.huc8.nunique(),
  "training_and_heldout_HUC8_overlap_count":len(overlap_exact),
  "training_original_cycle_counts":train.source_cycle.value_counts().to_dict(),
  "state_levels_unseen_in_training":sorted(set(target.state)-set(train.state)),
  "eco9_levels_unseen_in_training":sorted(set(target.eco9)-set(train.eco9)),
  "no_heldout_cycle_year_entered_the_model":True,
  "train_only_category_fit_and_numeric_normalization":True}
 if not sufficient(len(train),train.huc8.nunique(),400,200) or not sufficient(len(target),len(targetHUC),300,150):
  LOCO[targetera]={"QA":qa,"gate":"SOURCE_COHORT_NOT_SUPPORTED_FOR_CROSS_CYCLE_TRANSPORT"}
  continue
 P=preds(train,target)
 primary=score(target,P)
 # Stringent target-only HUC8 not present in either other original cycle,
 # distinct from primary zero-shared-HUC training constraint.
 newHUC=targetHUC-set(oth.huc8)
 sub=target[target.huc8.isin(newHUC)].reset_index(drop=True).copy()
 if sufficient(len(sub),sub.huc8.nunique(),150,100):
  ids=target.huc8.isin(newHUC).to_numpy()
  psubset={k:pred[ids] for k,pred in P.items()}
  exclusive={"gate":"SUPPORTED_ORIGINAL_TARGET_ONLY_HUC8",
             **score(sub,psubset)}
 else:
  exclusive={"gate":"INSUFFICIENT_ORIGINAL_TARGET_ONLY_HUC8",
             "n_UID":len(sub),"n_HUC8":sub.huc8.nunique()}
 # Prespecified secondary fitting weights: UID-equal rather than HUC8-equal.
 visit_weight_fit=score(target,preds(train,target,weighted=False),compare=False)
 # Source original 2008 EPT shares have unknown taxon ORDER <5% counts.
 # Compare extrema without sourcing new insect identifiers.
 frame_upper=src.copy()
 is08=frame_upper.source_cycle.eq("2008_09")
 frame_upper.loc[is08,"ept_fraction"]=frame_upper.loc[is08,"ept_fraction_upper_bound"]
 train_upper=frame_upper[(~frame_upper.source_cycle.eq(targetera))&
                         (~frame_upper.huc8.isin(targetHUC))].reset_index(drop=True)
 test_upper=frame_upper[frame_upper.source_cycle.eq(targetera)].reset_index(drop=True)
 upper=score(test_upper,preds(train_upper,test_upper),compare=False)
 LOCO[targetera]={"QA":qa,"gate":"VALID_2CYCLE_TRAIN_ZERO_TARGET_HUC8",
  "transport_results":primary,
  "novel_target_huc8_not_in_either_training_source_cycle":exclusive,
  "visit_equal_train_HUC8_equal_target_sensitivity":{
   "gain_NP_after_PHab_pct":visit_weight_fit["extra_NP_after_PHab_equalHUC8_gain_pct"],
   "gain_PHab_after_NP_pct":visit_weight_fit["extra_PHab_after_NP_equalHUC8_gain_pct"]},
  "source_2008_taxon_UNKNOWN_allEPT_upper_bound_sensitivity":{
   "gain_NP_after_PHab_pct":upper["extra_NP_after_PHab_equalHUC8_gain_pct"],
   "gain_PHab_after_NP_pct":upper["extra_PHab_after_NP_equalHUC8_gain_pct"]}}
# Context only: site-level 5fold by HUC within each survey period; no
# leftout-cycle effect in training, no source year as a fit covariate.
INTERNAL={}
for era in source_era:
 target=src[src.source_cycle.eq(era)].reset_index(drop=True).copy()
 fold=np.array([int(hashlib.sha256(str(v).encode()).hexdigest()[:8],16)%5
                for v in target.huc8])
 if min(np.bincount(fold,minlength=5))<40:
  INTERNAL[era]={"gate":"INSUFFICIENT_SOURCE_HUC8_FOLDS","fold_n":np.bincount(fold,minlength=5).tolist()}
  continue
 pred={k:np.full(len(target),np.nan) for k in ["intercept"]+list(models)}
 foldinfo=[]
 for f in range(5):
  tr=target.loc[fold!=f].reset_index(drop=True)
  te=target.loc[fold==f].reset_index(drop=True)
  if set(tr.huc8)&set(te.huc8):raise RuntimeError("Same-cycle HUC8-fold training leakage")
  fitp=preds(tr,te)
  indices=(fold==f)
  if len(te)!=indices.sum():raise RuntimeError("Predictions not aligned to original HUC fold")
  for k in pred:pred[k][indices]=fitp[k]
  foldinfo.append({"heldout_n":len(te),"heldout_HUC8":te.huc8.nunique(),
                   "train_n":len(tr),"train_HUC8":tr.huc8.nunique()})
 if not all(np.isfinite(p).all() for p in pred.values()):
  raise RuntimeError("Missing within-cycle EPT out-of-HUC8 predictions")
 INTERNAL[era]={"gate":"WITHIN_CYCLE_FIVEFOLD_HUC8_COMPLETE",
  "folds":foldinfo,
  "same_source_environment_benchmark":score(target,pred)}
result={
 "schema":"gesui_NRSA_2008_2013_2018_REAL_EPT_external_period_AND_zero_shared_HUC8_v1",
 "source_preattached_protocol":"docs/NRSA_THREE_CYCLE_EPT_AMBIENT_NP_PHAB_HUC8_TRANSPORT_PROTOCOL_20261114.md",
 "original_3_period_native_UIDs":len(src),"original_HUC8_union":source["distinct_union_HUC8"],
 "source_sample_intersections_see_separate_audit":source["original_overlap_nHUC8"],
 "response":"original source BERWW/BERW TOTAL300 EPT counted share, 2008 conservative unknown-order lower bound",
 "strict_2008_source_order_validation_against_EPA_published_279_visits":True,
 "three_field_cycles_genuine_EPA_source":True,
 "training_only_standardization_and_onehot_categories":True,
 "no_year_as_predictor_across_unseen_period":True,
 "Ridge_fixed_alpha":10,"training_HUC8_weighted_primary":True,
 "catchment_equal_weight_target_MSE_primary":True,
 "source_2008_periphyton_Chla_AFDM_unavailable":True,
 "source_2013_2018_EPT_panels_intersect_PChla_sampling":True,
 "not_outcome_based_source_or_feature_selection":True,
 "bootstrap_is_fixed_heldout_error_resampling_not_fit_uncertainty":True,
 "not_WWTP_derived_NP_exposure_or_organic_mass_or_ecosystem_service":True,
 "leave_one_SURVEY_cycle_out_and_remove_ALL_shared_HUC8_from_train":LOCO,
 "within_each_cycle_sameHUC8_group5fold_context_ONLY":INTERNAL}
out=ROOT/"NRSA_THREECYCLE_EPT_NP_PHAB_NO_HUC_LEAK_EXTERNAL_TRANSPORT.json"
out.write_text(json.dumps(result,indent=2,sort_keys=True))
print("GESUI_THREECYCLE_EXTERNAL_EPT_TRANSPORT_RESULT")
print(json.dumps({
 "source_UID_n":len(src),"HUC8_union":source["distinct_union_HUC8"],
 "leaveout":{
  name:{"source":data["QA"],"gate":data["gate"],
   "primary":{k:v for k,v in data.get("transport_results",{}).items()
              if k in ("n_UID_visits","n_HUC8","MSE_equal_HUC8_primary",
                       "extra_NP_after_PHab_equalHUC8_gain_pct",
                       "extra_PHab_after_NP_equalHUC8_gain_pct",
                       "HUC8_bootstrap2000_fixedOOF95_NP_after_PHab",
                       "HUC8_bootstrap2000_fixedOOF95_PHab_after_NP")},
   "novelHUC":{k:v for k,v in data.get("novel_target_huc8_not_in_either_training_source_cycle",{}).items()
                   if k in ("gate","n_UID","n_HUC8","n_UID_visits",
                            "extra_NP_after_PHab_equalHUC8_gain_pct",
                            "extra_PHab_after_NP_equalHUC8_gain_pct")},
   "alt_training_weights":data.get("visit_equal_train_HUC8_equal_target_sensitivity"),
   "2008_taxon_upper":data.get("source_2008_taxon_UNKNOWN_allEPT_upper_bound_sensitivity")}
   for name,data in LOCO.items()},
 "within_cycle":{
  name:{"gate":data["gate"],"gains":{
    k:v for k,v in data.get("same_source_environment_benchmark",{}).items()
    if k in ("extra_NP_after_PHab_equalHUC8_gain_pct",
             "extra_PHab_after_NP_equalHUC8_gain_pct",
             "HUC8_bootstrap2000_fixedOOF95_NP_after_PHab",
             "HUC8_bootstrap2000_fixedOOF95_PHab_after_NP")}}
  for name,data in INTERNAL.items()}},sort_keys=True))
