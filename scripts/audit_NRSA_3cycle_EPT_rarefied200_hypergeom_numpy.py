#!/usr/bin/env python3
"""Independent numerical QA: scipy.stats.hypergeom mathematical rarefaction
and explicit weighted NumPy Ridge normal equations on authentic original
NRSA 2008/2013/2018 taxon counts (m=200). Compares ALL original 3x3x4
external period+zero original HUC8 predictions, and HUC8 cluster intervals.
Extends source-only original 2209 same UID cohort with source unknown-order
EPT upper-bound exploratory sensitivity; not independent field data.
"""
from __future__ import annotations
import hashlib,json,os,pathlib
import numpy as np,pandas as pd
from scipy.stats import hypergeom

TEMP=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))
R=TEMP/"gesui_threecycle_rarefied_taxa"
original=json.loads((R/"NRSA_THREECYCLE_RARE200_EPT_SHARE_RICHNESS_EXTERNAL_TRANSPORT.json").read_text())
if original["schema"]!="gesui_NR​​SA_3cycle_true_original_RARE200_EPT_taxa_vs_same_UID_share_cross-era_v1":
 raise RuntimeError("Original source-verified rare200 model not found")
dd=pd.read_csv(R/"AUDIT_ONLY_EPA_THREECYCLE_RARE200_SAME_UID.csv",
 dtype={"source_cycle":str,"source_UID":str,"huc8":str,"state":str,"eco9":str})
count=pd.read_csv(R/"AUDIT_ONLY_ORIGINAL_EPA_3CYCLE_TAXON_COUNTS.csv",
 dtype={"cycle":str,"UID":str,"key":str,"order":str},keep_default_na=False)
count["known"]=(count.known.astype(str).str.lower().isin(("true","1","yes")))
if len(dd)!=2209:raise RuntimeError("m200 original source cohort changed")
R200=200
check={}
observed=dd.set_index(["source_cycle","source_UID"],drop=False)
for (cycle,uid),grp in count.groupby(["cycle","UID"],sort=False):
 key=(cycle,uid)
 if key not in observed.index:continue
 N=int(round(grp["count"].sum()))
 if N<R200:raise RuntimeError("Source effort min200 contradicted its actual original individual counts")
 n=grp["count"].to_numpy(int)
 # Independent SciPy hypergeom probability of observing ZERO individuals of
 # one original native taxon ID in a draw of exactly 200 distinct individuals.
 # This is equivalent but separately implemented from original gammaln formula.
 probability_present=1.0-hypergeom.pmf(0,N,n,R200)
 if not np.isfinite(probability_present).all():
  raise RuntimeError("Independent hypergeom rarefied original taxon probability nonfinite")
 code=grp.order.isin(("EPHEMEROPTERA","PLECOPTERA","TRICHOPTERA")).to_numpy()
 known=grp.known.to_numpy(bool)
 valid=grp.key.ne("").to_numpy()
 values={
  "rarefied_200_EPT_taxonID":float(probability_present[code&valid].sum()),
  "rarefied_200_nonEPT_taxonID":float(probability_present[(~code)&known&valid].sum()),
  "rarefied_200_total_taxonID":float(probability_present[valid].sum()),
  "rarefied_200_EPT_upper_if_unknown_order_all_EPT":float(probability_present[code&valid].sum()+probability_present[(~known)&valid].sum())}
 for field,v in values.items():
  other=float(observed.loc[key,field])
  if abs(v-other)>5e-8:
   raise RuntimeError("Independent hypergeom probability differs from original source "+str((key,field,v,other)))
 check[cycle]=check.get(cycle,0)+1
if check!={"2008_09":476,"2013_14":957,"2018_19":776}:
 raise RuntimeError("Original 3-epoch hypergeom full cohort not independently checked "+str(check))
CYCLES=["2008_09","2013_14","2018_19"]
GEO=["lat","lon","stream_order","season_sin","season_cos","log_cond"]
NP=["log_tn","log_tp"]
PHAB=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]
XSET={"B":GEO,"B+NP":GEO+NP,"B+PHab":GEO+PHAB,"B+PHab+NP":GEO+PHAB+NP}
RESP={
 "EPT_rarefied_200_taxonID_log1p":"log1p_rare200_EPT_taxonID",
 "EPT_counted_fraction_same_sites":"ept_fraction",
 "nonEPT_rarefied_200_taxonID_log1p":"log1p_rare200_nonEPT_taxonID"}
full=pd.read_csv(TEMP/"gesui_threecycle_ept_transfer"/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv",
 dtype={"source_cycle":str,"source_UID":str,"huc8":str})
if len(full)!=2648:raise RuntimeError("Original full EPA source HUC8 exposure leaks or changed")
def manual(train,test,col,cols):
 a=train[cols].to_numpy(float)
 b=test[cols].to_numpy(float)
 mu=a.mean(0);sd=a.std(0);sd=np.where(sd>1e-12,sd,1)
 ec=np.array(sorted(train.eco9.unique()),dtype=str)
 st=np.array(sorted(train.state.unique()),dtype=str)
 def make(D,v):
  return np.column_stack([
   (v-mu)/sd,
   (D.state.to_numpy(str)[:,None]==st[None,:]).astype(float),
   (D.eco9.to_numpy(str)[:,None]==ec[None,:]).astype(float)])
 a0=make(train,a);b0=make(test,b)
 A=np.column_stack([np.ones(len(train)),a0])
 B=np.column_stack([np.ones(len(test)),b0])
 cnt=train.groupby("huc8").source_UID.transform("size").to_numpy(float)
 w=(len(train)/train.huc8.nunique())/cnt
 pen=np.eye(A.shape[1])*10;pen[0,0]=0
 y=train[col].to_numpy(float)
 coef=np.linalg.solve(A.T@(w[:,None]*A)+pen,A.T@(w*y))
 return B@coef
def loss(target,y,p):
 u,idx=np.unique(target.huc8.to_numpy(str),return_inverse=True)
 err=np.bincount(idx,weights=(y-p)**2,minlength=len(u))/np.bincount(idx,minlength=len(u))
 return err
def score(target,Y,p):
 e={k:loss(target,Y,pr) for k,pr in p.items()}
 mse={k:float(q.mean()) for k,q in e.items()}
 rng=np.random.default_rng(20261114)
 H=len(e["B"])
 j=rng.integers(0,H,size=(2000,H))
 def gain(b,a):
  return float(100*(mse[b]-mse[a])/mse[b])
 def boot(b,a):
  bas=e[b][j].mean(1);alt=e[a][j].mean(1)
  return list(map(float,np.quantile(100*(bas-alt)/bas,[.025,.5,.975])))
 return {
  "MSE":mse,
  "extra_PHab_pct":gain("B+NP","B+PHab+NP"),
  "extra_NP_pct":gain("B+PHab","B+PHab+NP"),
  "PHab_HUC8_boot95":boot("B+NP","B+PHab+NP"),
  "NP_HUC8_boot95":boot("B+PHab","B+PHab+NP")}
out={}
for era in CYCLES:
 target=dd[dd.source_cycle.eq(era)].reset_index(drop=True)
 train=dd[(~dd.source_cycle.eq(era)) &
          (~dd.huc8.isin(set(full.loc[full.source_cycle.eq(era),"huc8"])))].reset_index(drop=True)
 previous=original["external_whole_cycle_zero_original_HUC8"][era]
 if previous["gate"]!="TRUE_OUT_OF_SURVEY_AND_ALL_ORIGINAL_TARGET_HUC8":
  raise RuntimeError("Source original m200 three-cycle test didn't meet precommitted gate "+era)
 if len(train)!=previous["QA"]["train_original_after_excluding_ALL_target_period_native_HUC8"] or len(target)!=previous["QA"]["heldout_original_source_uid_count_m200"]:
  raise RuntimeError("Original source-only train/test n changed "+era)
 if set(train.huc8)&set(target.huc8):raise RuntimeError("HUC8 leak")
 results={}
 for nm,col in RESP.items():
  p={k:manual(train,target,col,f) for k,f in XSET.items()}
  q=score(target,target[col].to_numpy(float),p)
  old=previous["responses"][nm]
  for key,ref in [
    ("extra_PHab_pct",old["gain_PHab_added_AFTER_NP_pct"]),
    ("extra_NP_pct",old["gain_NP_added_AFTER_PHab_pct"])]:
   if abs(q[key]-ref)>3e-7:raise RuntimeError("Manual normal equation point gain different "+era+"/"+nm+"/"+key)
  for key,ref in [
    ("PHab_HUC8_boot95",old["HUC8_fixedOOF_bootstrap95_PHab_after_NP"]),
    ("NP_HUC8_boot95",old["HUC8_fixedOOF_bootstrap95_NP_after_PHab"])]:
   if max(abs(x-y) for x,y in zip(q[key],ref))>2e-6:
    raise RuntimeError("Manual normal equation HUC8 fixed-error bootstrap differs "+era+"/"+nm+"/"+key)
  for k,origm in old["MSE_equal_HUC8"].items():
   if abs(q["MSE"][k]-origm)>1e-9:
    raise RuntimeError("Independent 3-cycle m200 MSE mismatch "+era+"/"+nm+"/"+k)
  results[nm]={"source_original_UID_n":len(target),"source_test_HUC8":target.huc8.nunique(),
               "normal_equation_source_model_MSE":q["MSE"],
               "PHab_gain_pct":q["extra_PHab_pct"],
               "NP_gain_pct":q["extra_NP_pct"],
               "PHab_95":q["PHab_HUC8_boot95"],
               "NP_95":q["NP_HUC8_boot95"]}
 # Explicit source taxonomy uncertainty (source 2008 unknown ORDER ->
 # all such original native taxon IDs assumed EPT = UPPER bound).
 # Does not alter target site selection, historical stream nutrient data,
 # geographical split, baseline or numerical fit method.
 au=dd.copy()
 idx=au.source_cycle.eq("2008_09")
 au.loc[idx,"log1p_rare200_EPT_taxonID"]=au.loc[idx,"log1p_rare200_EPT_upper"]
 t=au[au.source_cycle.eq(era)].reset_index(drop=True)
 tr=au[(~au.source_cycle.eq(era))&
       (~au.huc8.isin(set(full.loc[full.source_cycle.eq(era),"huc8"])))].reset_index(drop=True)
 pred={k:manual(tr,t,"log1p_rare200_EPT_taxonID",f) for k,f in XSET.items()}
 upper=score(t,t.log1p_rare200_EPT_taxonID.to_numpy(float),pred)
 out[era]={"independently_recomputed_original_3response_models":results,
           "source_unknown_ORDER_all_such_original_taxa_as_EPT_UPPER_sensitivity":{
            "n_original_target_UID":len(t),"HUC8":t.huc8.nunique(),
            "extra_PHab_gain_pct":upper["extra_PHab_pct"],
            "extra_NP_gain_pct":upper["extra_NP_pct"],
            "PHab_HUC8_boot95":upper["PHab_HUC8_boot95"],
            "NP_HUC8_boot95":upper["NP_HUC8_boot95"]}}
AUDIT={
 "schema":"gesui_EP​​A_three_cycle_rarefied200_hypergeom_independentMatrix_and_unknownORDER_UPPER_QA_v1",
 "original_exact_2209_native_source_UID_same_as_source_only":True,
 "independent_Scipy_hypergeom_2209_original_source_visit_rarefaction_matches_gamma_version":True,
 "every_EPT_nonEPT_total_original_expected_taxonID_probability_checked":True,
 "independent_manual_Ridge10_HUC8_weighted_normal_equation_all_3x3x4_externals_verified":True,
 "all_full-period_target_HUC8_excluded_from_other-period_eligible_train":True,
 "same_original_EPA_field_observations_not_new_independent_field_sample":True,
 "source_taxon_unknown_order_UPPER_bound_is_hypothetical_not_reclassified_field_taxa":True,
 "no_WWTP_delivered_NP_source_or_ecosystem_service":True,
 "true_per_cycle":out}
path=R/"NRSA_THREECYCLE_RARE200_INDEPENDENT_HYPERGEOM_WEIGHTEDRIDGE_QA.json"
path.write_text(json.dumps(AUDIT,indent=2,sort_keys=True))
print("GESUI_THREECYCLE_RARE200_INDEPENDENT_HYPERGEOM_NUMPY_QA_SUCCESS")
print(json.dumps({k:{"n":v["independently_recomputed_original_3response_models"]["EPT_rarefied_200_taxonID_log1p"]["source_original_UID_n"],
 "rareEPT_taxa":{name:z for name,z in v["independently_recomputed_original_3response_models"]["EPT_rarefied_200_taxonID_log1p"].items()
                  if name in ("PHab_gain_pct","NP_gain_pct","PHab_95","NP_95")},
 "unknown_SOURCE_UPPER":v["source_unknown_ORDER_all_such_original_taxa_as_EPT_UPPER_sensitivity"]}
 for k,v in out.items()},sort_keys=True))
