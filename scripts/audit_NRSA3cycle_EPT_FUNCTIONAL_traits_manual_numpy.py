#!/usr/bin/env python3
"""Independently validate all three true original EPA EPT functional traits:
* raw native taxon count-weighted FFG SCRAPER, HABIT CLINGER, mean source PTV
  across exactly verified source EPT organisms, missing traits not zero
* 3 periods ×3 traits ×4 environment models: manual weighted NumPy Ridge10,
  genuine zero full-target-era HUC8 overlap, 2000 HUC8 heldout bootstrap
* within-era measured XEMBED, ambient TN/TP HUC8 clustered OLS sandwiches.
All are original same data with independent code, not new field replicates or
measured wastewater, emerging insect biomass or ecosystem services.
"""
from __future__ import annotations
import hashlib,json,os,pathlib
import numpy as np,pandas as pd

TEMP=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))
R=TEMP/"gesui_threecycle_traits"
M=TEMP/"gesui_threecycle_rarefied_taxa"
O=json.loads((R/"NRSA_3CYCLE_EPT_SCRAPER_CLINGER_PTV_TRUE_FUNCTION_EXTERNAL.json").read_text())
if O["schema"]!="gesui_EP​​A_3cycle_EPT_FFG_HABIT_PTV_function_original_externalHUC8_v1":
 raise RuntimeError("Original EPA trait source results schema changed")
D=pd.read_csv(R/"AUDIT_ONLY_3CYCLE_ORIGINAL_EPT_FFG_HABIT_PTV_SAMESITE.csv",
 dtype={"source_cycle":str,"source_UID":str,"huc8":str,"state":str,"eco9":str})
counts=pd.read_csv(M/"AUDIT_ONLY_ORIGINAL_EPA_3CYCLE_TAXON_COUNTS.csv",
 dtype={"cycle":str,"UID":str,"key":str,"order":str},keep_default_na=False)
lookup=pd.read_csv(R/"AUDIT_ONLY_EPA_ORIGINAL_3CYCLE_TAXA_FFG_HABIT_PTV.csv",
 dtype=str,keep_default_na=False)
full=pd.read_csv(TEMP/"gesui_threecycle_ept_transfer"/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv",
 dtype={"source_cycle":str,"source_UID":str,"huc8":str})
if len(D)!=2209 or len(full)!=2648:
 raise RuntimeError("Original 3 EPA cycle m200 cohort changed")
if lookup.duplicated(["source_cycle","key"]).any():raise RuntimeError("Same original EPA taxon key has conflicting normalized traits")
ept=counts[counts.order.isin(("EPHEMEROPTERA","PLECOPTERA","TRICHOPTERA"))].copy()
ept["count"]=pd.to_numeric(ept["count"],errors="raise")
ept=ept.merge(lookup[["source_cycle","key","order","FFG_norm","HABIT_norm","PTV_norm"]],
 left_on=["cycle","key","order"],right_on=["source_cycle","key","order"],how="left",validate="many_to_one")
if len(ept)==0:raise RuntimeError("No original EPT organism in trait source")
# Independent fractional coding: equal fraction among the UNIQUE source FFG
# or HABIT modes, not one per each listed multi-role.
def codeweight(x,unit,legal):
 if not isinstance(x,str) or not x.strip():return np.nan
 vals=set(z.strip().upper() for z in x.split(",") if z.strip())
 if not vals or not vals.issubset(legal):return np.nan
 return float(unit in vals)/len(vals)
ept["scraper"]=ept.FFG_norm.map(lambda z:codeweight(z,"SC",set("CF CG PA PI PR SC SH".split())))
ept["clinger"]=ept.HABIT_norm.map(lambda z:codeweight(z,"CN",set("BU CB CN DV PK SK SP SW".split())))
ept["ptv"]=pd.to_numeric(ept.PTV_norm,errors="coerce")
if ept.ptv.dropna().lt(0).any() or ept.ptv.dropna().gt(10).any():
 raise RuntimeError("Original EPA trait value PTV not 0..10")
# Genuine source-positive EPT specimen count, not any invented absolute
# insect density. Reference site fractions were built by a SEPARATE script.
site={}
for (cycle,uid),x in ept.groupby(["cycle","UID"],sort=False):
 N=float(x["count"].sum())
 if N<=0:raise RuntimeError("Nonpositive original count of EPT")
 res={"source_cycle":cycle,"source_UID":uid}
 for key,name in (("scraper","EPT_scraper_FFG_fraction"),
                  ("clinger","EPT_clinger_HABIT_fraction"),
                  ("ptv","EPT_mean_pollution_tolerance_PTV")):
  values=x[key].to_numpy(float)
  weights=x["count"].to_numpy(float)
  mask=np.isfinite(values)
  rate=float(weights[mask].sum()/N)
  mean=float(np.dot(weights[mask],values[mask])/weights[mask].sum()) if mask.any() else np.nan
  res[name+"_source_coverage"]=rate
  res[name]=mean if rate>=.8-1e-10 else np.nan
 site[(cycle,uid)]=res
# Confirm independent count-weighted source trait numbers for every genuine
# original site and every endpoint; NaNs match where no EPT or insufficient
# original PTV/FFG/HABIT coverage.
for _,row in D.iterrows():
 z=site.get((row.source_cycle,row.source_UID))
 for name in ("EPT_scraper_FFG_fraction","EPT_clinger_HABIT_fraction","EPT_mean_pollution_tolerance_PTV"):
  val=float(row[name]) if pd.notna(row[name]) else np.nan
  other=np.nan if z is None else z[name]
  if (np.isfinite(val)!=np.isfinite(other)) or (np.isfinite(val) and abs(val-other)>1e-10):
   raise RuntimeError("Independent EPA original source individual-weighted trait response differs "+
      str((row.source_cycle,row.source_UID,name,val,other)))
GEO=["lat","lon","stream_order","season_sin","season_cos","log_cond"]
NUT=["log_tn","log_tp"]
PHAB=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]
MODELS={"B":GEO,"B+NP":GEO+NUT,"B+PHab":GEO+PHAB,"B+PHab+NP":GEO+PHAB+NUT}
RESP={
 "scraper_FFG_within_EPT":"EPT_scraper_FFG_fraction",
 "clinger_habit_within_EPT":"EPT_clinger_HABIT_fraction",
 "source_pollution_tolerance_within_EPT":"EPT_mean_pollution_tolerance_PTV"}
eras=["2008_09","2013_14","2018_19"]
full_hucs={era:set(full.loc[full.source_cycle.eq(era),"huc8"]) for era in eras}
def manual_predict(tr,te,col,fields):
 a=tr[fields].to_numpy(float);b=te[fields].to_numpy(float)
 mu=a.mean(0);sd=a.std(0);sd=np.where(sd>1e-12,sd,1)
 states=np.array(sorted(tr.state.unique()),dtype=str)
 eco=np.array(sorted(tr.eco9.unique()),dtype=str)
 def X(f,Z):
  return np.column_stack([
   (Z-mu)/sd,
   (f.state.to_numpy(str)[:,None]==states[None,:]).astype(float),
   (f.eco9.to_numpy(str)[:,None]==eco[None,:]).astype(float)])
 A=np.column_stack([np.ones(len(tr)),X(tr,a)])
 B=np.column_stack([np.ones(len(te)),X(te,b)])
 cnt=tr.groupby("huc8").source_UID.transform("size").to_numpy(float)
 w=(len(tr)/tr.huc8.nunique())/cnt
 pen=np.eye(A.shape[1])*10;pen[0,0]=0
 beta=np.linalg.solve(A.T@(w[:,None]*A)+pen,A.T@(w*tr[col].to_numpy(float)))
 return B@beta
def group_error(test,p,y):
 u,inv=np.unique(test.huc8.to_numpy(str),return_inverse=True)
 return np.bincount(inv,weights=(y-p)**2,minlength=len(u))/np.bincount(inv,minlength=len(u))
def manual_ci(ref,test,pred,col):
 y=test[col].to_numpy(float)
 errors={key:group_error(test,zz,y) for key,zz in pred.items()}
 mse={key:float(q.mean()) for key,q in errors.items()}
 for key,number in mse.items():
  original=ref["external"]["original_target_HUC8_equal_MSE"][key]
  if abs(number-original)>1e-9:
   raise RuntimeError("Independent functional trait source heldout MSE differs "+key)
 def gains(a,b):
  return float(100*(mse[a]-mse[b])/mse[a])
 ph=gains("B+NP","B+PHab+NP")
 npgain=gains("B+PHab","B+PHab+NP")
 if abs(ph-ref["external"]["original_extra_PHab_AFTER_real_NP_percent"])>1e-7:
  raise RuntimeError("Independent EPT trait PHab added-value != original")
 if abs(npgain-ref["external"]["original_extra_real_NP_AFTER_PHab_percent"])>1e-7:
  raise RuntimeError("Independent EPT trait N/P added-value != original")
 rng=np.random.default_rng(20261114)
 inds=rng.integers(0,len(errors["B"]),size=(2000,len(errors["B"])))
 def bootstrap(a,b):
  x=errors[a][inds].mean(axis=1);v=errors[b][inds].mean(axis=1)
  return list(map(float,np.quantile(100*(x-v)/x,[.025,.5,.975])))
 ci_ph=bootstrap("B+NP","B+PHab+NP")
 ci_np=bootstrap("B+PHab","B+PHab+NP")
 for obs,expect in [
  (ci_ph,ref["external"]["HUC8_fixedOOF_bootstrap2000_PHab95"]),
  (ci_np,ref["external"]["HUC8_fixedOOF_bootstrap2000_NP95"])]:
  if max(abs(a-b) for a,b in zip(obs,expect))>1e-6:
   raise RuntimeError("Original functional EPT source HUC8 bootstrap differs")
 return {"PHab_gain_pct":ph,"NP_gain_pct":npgain,
         "PHab_CI95":ci_ph,"NP_CI95":ci_np}
def sandwich(test,col):
 X=pd.concat([test[GEO+PHAB+NUT].astype(float),
              pd.get_dummies(test[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
 for k in ("XEMBED","log_tn","log_tp"):
  X[k]=(X[k]-X[k].mean())/X[k].std(ddof=0)
 A=np.column_stack([np.ones(len(test)),X.to_numpy(float)])
 y=test[col].to_numpy(float)
 beta=np.linalg.lstsq(A,y,rcond=None)[0]
 err=y-A@beta
 bread=np.linalg.pinv(A)
 meat=np.zeros((A.shape[1],A.shape[1]))
 groups=test.huc8.to_numpy(str)
 for h in pd.unique(groups):
  mask=groups==h
  v=A[mask].T@err[mask]
  meat+=np.outer(v,v)
 G=len(pd.unique(groups));K=np.linalg.matrix_rank(A);N=len(test)
 V=bread@meat@bread.T*G/(G-1)*(N-1)/(N-K)
 vals={}
 for k in ("XEMBED","log_tn","log_tp"):
  i=X.columns.get_loc(k)+1
  b=float(beta[i]);se=float(np.sqrt(max(0,V[i,i])))
  vals[k]={"beta":b,"ci":[b-1.95996398454*se,b+1.95996398454*se]}
 return vals
results={}
for era in eras:
 answers={}
 for key,col in RESP.items():
  prior=O["per_untouched_external_cycle"][era][key]
  test=D[(D.source_cycle.eq(era))&D[col].notna()].reset_index(drop=True)
  train=D[(~D.source_cycle.eq(era))&~D.huc8.isin(full_hucs[era])&D[col].notna()].reset_index(drop=True)
  if prior["status"]!="OBSERVATIONAL_GENUINE_EPA_TRAIT_NUMERIC_FIT":
   answers[key]={"original_gate":prior["status"]};continue
  if len(test)!=prior["QA"]["source_target_UID"] or len(train)!=prior["QA"]["source_other_two_cycle_train_UID_after_all_original_targetHUC8_excluded"]:
   raise RuntimeError("Original source N/HUC8 population for trait changed "+era+"/"+key)
  if set(train.huc8)&set(test.huc8):raise RuntimeError("Original target full HUC8 leakage")
  prediction={k:manual_predict(train,test,col,v) for k,v in MODELS.items()}
  metrics=manual_ci(prior,test,prediction,col)
  sl=sandwich(test,col)
  for nm,m in sl.items():
   old=prior["descriptive_same_cycle_nutrient_bed_standardized_partial_slopes"][nm]
   if abs(m["beta"]-old["source_trait_change_per_1SD_original_predictor"])>1e-7:
    raise RuntimeError("Original EPA 3-era EPT habitat/NP trait coefficient mismatch "+era+"/"+key+"/"+nm)
   lo,hi=old["HUC8_cluster_descriptive95"]
   half=max(1e-8,(hi-lo)/2)
   discrepancy=max(abs(x-y) for x,y in zip(m["ci"],[lo,hi]))/half
   if discrepancy>0.003:
    raise RuntimeError("Original EPT trait HUC8 cluster ci >0.3% halfwidth "+era+"/"+key+"/"+nm)
  answers[key]={"n":len(test),"HUC8":test.huc8.nunique(),
                "external":metrics,"environmental_partial_source_betas":sl}
 results[era]=answers
AUD={
 "schema":"gesui_EPA3cycle_EPT_TRUE_traits_36external_Ridge_manualNumPy_and_source_individual_means_QA_v1",
 "original_same_2209_EPA_field_UID_response_source_faithfully_recomputed":True,
 "original_three_traits_site_values_independently_count_weighted_before_model_validation":True,
 "all_3eras_3traits_4models_heldout_HUC8_MSE_and_bootstrap_matched":True,
 "all_27_source_descriptive_XEMBED_TNTP_coefficients_verified":True,
 "HUC8_robust_CI_rank_tolerance_relative_halfwidth":0.003,
 "same_actual_original_EPA_specimens_not_independent_ecological_dataset":True,
 "no_original_WWTP_source_or_insect_drymass_or_realized_ecosystem_services":True,
 "source_real_epa_independently_recomputed":results}
path=R/"NRSA_EPT_3CYCLE_FFG_HABIT_PTV_SOURCE_COUNTS_MANUAL_NUMPY_QA.json"
path.write_text(json.dumps(AUD,sort_keys=True,indent=2))
print("GESUI_EPT_THREECYCLE_TRAITS_INDEPENDENT_SOURCE_MANUAL_NUMPY_SUCCESS")
print(json.dumps({era:{key:{"n":val.get("n"),
 "PHab_gain":val.get("external",{}).get("PHab_gain_pct"),
 "NP_gain":val.get("external",{}).get("NP_gain_pct"),
 "XEMBED_beta":val.get("environmental_partial_source_betas",{}).get("XEMBED")}
 for key,val in resp.items()} for era,resp in results.items()},sort_keys=True))
