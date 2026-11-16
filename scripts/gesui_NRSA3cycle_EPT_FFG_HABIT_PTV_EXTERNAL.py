#!/usr/bin/env python3
"""Original EPA three national survey periods, SOURCE matched EPT functional
feeding scraper (SC), clinger habit (CN), and source pollution tolerance
(PTV). Independent 2008/2013 source taxon autecology plus 2018 original
name+ID+order matched prior autecology; all real EPT individual counts.

Three environmental predictor groups: geography/state/ecoregion, real
in-stream measured TN+TP, six original PHab. For each target survey cycle,
remove all original target-cycle HUC8 before fitting other two cycles. Fixed
Ridge10, train-only preprocessing, equal-HUC8 fitting weights; HUC8 error
MSE and 2000 fixed-OOF group bootstrap. Observational not wastewater causal.
"""
from __future__ import annotations
import hashlib,json,os,pathlib
import numpy as np,pandas as pd
import statsmodels.api as sm
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.linear_model import Ridge

TEMP=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))
ROOT=TEMP/"gesui_threecycle_traits"
Q=json.loads(pathlib.Path(os.environ["ORIGINAL_EPT_TRAIT_SITE_SOURCE_QA"]).read_text())
if Q["schema"]!="gesui_EPA3cycle_SOURCEONLY_ORIGINAL_EPT_SC_CN_PTV_nativeUID_true_coverage_v1":
 raise RuntimeError("Original EPA EPT site trait coverage and source native UID not locked")
df=pd.read_csv(ROOT/"AUDIT_ONLY_3CYCLE_ORIGINAL_EPT_FFG_HABIT_PTV_SAMESITE.csv",
 dtype={"source_cycle":str,"source_UID":str,"huc8":str,"state":str,"eco9":str})
original=pd.read_csv(TEMP/"gesui_threecycle_ept_transfer"/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv",
 dtype={"source_cycle":str,"source_UID":str,"huc8":str})
if len(df)!=2209 or len(original)!=2648 or df.duplicated(["source_cycle","source_UID"]).any():
 raise RuntimeError("Original three-cycle native EPT site cohort changed")
source_huc={c:set(original.loc[original.source_cycle.eq(c),"huc8"]) for c in ("2008_09","2013_14","2018_19")}
FIELD={
 "scraper_FFG_within_EPT":"EPT_scraper_FFG_fraction",
 "clinger_habit_within_EPT":"EPT_clinger_HABIT_fraction",
 "source_pollution_tolerance_within_EPT":"EPT_mean_pollution_tolerance_PTV"}
GEO=["lat","lon","stream_order","season_sin","season_cos","log_cond"]
N=["log_tn","log_tp"]
H=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]
MAT={"B":GEO,"B+NP":GEO+N,"B+PHab":GEO+H,"B+PHab+NP":GEO+H+N}
def fit(train,test,outcome,columns):
 prep=ColumnTransformer([
  ("numeric",StandardScaler(),columns),
  ("strata",OneHotEncoder(handle_unknown="ignore",sparse_output=False,dtype=float),["state","eco9"])
 ],remainder="drop",sparse_threshold=0)
 x=prep.fit_transform(train);z=prep.transform(test)
 g=train.huc8.value_counts()
 w=(len(train)/len(g))/train.huc8.map(g)
 yh=train[outcome].to_numpy(float)
 p=Ridge(alpha=10,fit_intercept=True,solver="svd").fit(
   x,yh,sample_weight=w.to_numpy(float)).predict(z)
 if not np.isfinite(p).all():raise RuntimeError("Functional trait OOS prediction nonfinite")
 return p
def losses(test,y,p):
 ids,inv=np.unique(test.huc8.to_numpy(str),return_inverse=True)
 return np.bincount(inv,weights=(y-p)**2,minlength=len(ids))/np.bincount(inv,minlength=len(ids))
def compare(test,y,pred):
 e={k:losses(test,y,v) for k,v in pred.items()}
 mse={k:float(v.mean()) for k,v in e.items()}
 gain=lambda b,a:float(100*(mse[b]-mse[a])/mse[b])
 rng=np.random.default_rng(20261114)
 ix=rng.integers(0,len(e["B"]),size=(2000,len(e["B"])))
 def interval(b,a):
  l=e[b][ix].mean(axis=1);r=e[a][ix].mean(axis=1)
  return list(map(float,np.quantile(100*(l-r)/l,[.025,.5,.975])))
 return {
  "original_target_HUC8_equal_MSE":mse,
  "original_extra_PHab_AFTER_real_NP_percent":gain("B+NP","B+PHab+NP"),
  "original_extra_real_NP_AFTER_PHab_percent":gain("B+PHab","B+PHab+NP"),
  "HUC8_fixedOOF_bootstrap2000_PHab95":interval("B+NP","B+PHab+NP"),
  "HUC8_fixedOOF_bootstrap2000_NP95":interval("B+PHab","B+PHab+NP")}
def descriptive_ols(test,outcome):
 X=pd.concat([
  test[GEO+H+N].astype(float),
  pd.get_dummies(test[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
 for col in ("XEMBED","log_tn","log_tp"):
  X[col]=(X[col]-X[col].mean())/X[col].std(ddof=0)
 y=test[outcome].to_numpy(float)
 fit=sm.OLS(y,sm.add_constant(X.to_numpy(float),has_constant="add")).fit(
  cov_type="cluster",cov_kwds={"groups":test.huc8.to_numpy(str),"use_correction":True})
 return {
  col:{"source_trait_change_per_1SD_original_predictor":float(fit.params[X.columns.get_loc(col)+1]),
    "HUC8_cluster_descriptive95":list(map(float,fit.conf_int()[X.columns.get_loc(col)+1]))}
  for col in ("XEMBED","log_tn","log_tp")}
out={}
for cycle in ("2008_09","2013_14","2018_19"):
 results={}
 for name,col in FIELD.items():
  test=df[(df.source_cycle.eq(cycle))&df[col].notna()].reset_index(drop=True)
  train=df[(~df.source_cycle.eq(cycle))&
           (~df.huc8.isin(source_huc[cycle]))&df[col].notna()].reset_index(drop=True)
  overlap=set(train.huc8)&set(test.huc8)
  if overlap:raise RuntimeError("Original source HUC8 leak in trait external model")
  QA={"original_EPT_native_trait_response":name,
      "source_target_UID":len(test),"source_target_HUC8":test.huc8.nunique(),
      "source_target_states":test.state.nunique(),
      "source_other_two_cycle_train_UID_after_all_original_targetHUC8_excluded":len(train),
      "source_other_two_cycle_train_HUC8":train.huc8.nunique(),
      "source_train_states":train.state.nunique(),
      "train_cycle_counts":train.source_cycle.value_counts().to_dict(),
      "cross_era_train_test_HUC8_overlap":len(overlap)}
  prior=Q["source_original_sites_by_cycle_and_trait"][cycle]["original_explicit_trait_80pct_or_more_site_eligibility"][col]
  if [len(test),test.huc8.nunique(),test.state.nunique()] != [prior["n_native_UID"],prior["n_HUC8"],prior["states"]]:
   raise RuntimeError("Original source-locked EPT trait site cohort changed "+cycle+"/"+name)
  gate=(len(test)>=250 and test.huc8.nunique()>=150 and test.state.nunique()>=20 and
        len(train)>=300 and train.huc8.nunique()>=200 and train.state.nunique()>=15)
  QA["source_precommitted_fullHUC8_train_test_support"]=bool(gate)
  if not gate:
   results[name]={"QA":QA,"status":"SOURCE_GATE_FAIL_NO_TRAIT_REGRESSION"}
   continue
  p={k:fit(train,test,col,v) for k,v in MAT.items()}
  y=test[col].to_numpy(float)
  stats=compare(test,y,p)
  ols=descriptive_ols(test,col)
  results[name]={"QA":QA,"status":"OBSERVATIONAL_GENUINE_EPA_TRAIT_NUMERIC_FIT",
                 "true_original_trait_response_mean":float(y.mean()),
                 "true_original_trait_response_std":float(y.std(ddof=0)),
                 "external":stats,
                 "descriptive_same_cycle_nutrient_bed_standardized_partial_slopes":ols}
 out[cycle]=results
OUT={
 "schema":"gesui_EP​​A_3cycle_EPT_FFG_HABIT_PTV_function_original_externalHUC8_v1",
 "EPA_true_three_era_original_counted_EPT_autecology_traits":True,
 "source_site_QA_artifact":"NRSA_3CYCLE_ORIGINAL_EPT_FFG_HABIT_PTV_SOURCE_ONLY_SITES.json",
 "source_outcome_group_n2209_before_per_trait_coverage":True,
 "not_direct_feeding_behavior_or_evolved_pollution_tolerance":True,
 "not_sewage_derived_NP_load_or_emergent_insect_mass":True,
 "source_original_three_era_full_targetHUC8_excluded_from_other_era_training":True,
 "fixed_four_original_environment_predictor_models_and_trainonly_categories":True,
 "precommitted_3_trait_endpoints":FIELD,
 "all_source_gate_fails_negative_and_positive_results_retained":True,
 "per_untouched_external_cycle":out}
p=ROOT/"NRSA_3CYCLE_EPT_SCRAPER_CLINGER_PTV_TRUE_FUNCTION_EXTERNAL.json"
p.write_text(json.dumps(OUT,sort_keys=True,indent=2))
print("GESUI_THREE_ERA_EPA_ORIGINAL_TRAITS_FUNCTIONAL_EXTERNAL_RESULT")
print(json.dumps({era:{name:{"status":r["status"],
    "n":r["QA"]["source_target_UID"],"HUC8":r["QA"]["source_target_HUC8"],
    "train_n":r["QA"]["source_other_two_cycle_train_UID_after_all_original_targetHUC8_excluded"],
    "PHab_gain":r.get("external",{}).get("original_extra_PHab_AFTER_real_NP_percent"),
    "NP_gain":r.get("external",{}).get("original_extra_real_NP_AFTER_PHab_percent"),
    "PHab_boot95":r.get("external",{}).get("HUC8_fixedOOF_bootstrap2000_PHab95"),
    "NP_boot95":r.get("external",{}).get("HUC8_fixedOOF_bootstrap2000_NP95"),
    "withinEra_observed_partial_slopes":r.get("descriptive_same_cycle_nutrient_bed_standardized_partial_slopes")}
 for name,r in models.items()} for era,models in out.items()},sort_keys=True))
