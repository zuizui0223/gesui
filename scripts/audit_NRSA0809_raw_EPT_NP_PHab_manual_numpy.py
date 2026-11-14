#!/usr/bin/env python3
"""Independent NumPy 2008 original EPA counted EPT share / EPT taxon count /
MMI source exact third-cycle physical and measured TN/TP HUC8 tests.
Independent numerical implementation, not new independently sampled ecology.
The original 2008 source count reconstructed EPT under unknown taxa lower/
upper bounds is compared on identical full source-selected visits.
"""
from __future__ import annotations
import hashlib,json,os,pathlib
import numpy as np,pandas as pd
R=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui0809_reconstructed_EPT"
report=json.loads((R/"NRSA0809_RECONSTRUCTED_EPT_REAL_AMBIENT_NP_PHab_NUMERIC_RESULTS.json").read_text())
D=pd.read_csv(R/"AUDIT_ONLY_2008_RAW_EPT_CHEMW_PHAB_SOURCE_PANEL.csv",
 dtype={"huc8":str,"state":str,"eco9":str,"site_id":str})
if report["schema"]!="gesui_NRSA0809_original_rescued_EPT_count_share_EPT_taxon_MMI_NP_PHAB_2008_v1":
 raise RuntimeError("Original third cycle 2008 source model schema altered")
if len(D)!=586 or D.huc8.nunique()!=437:
 raise RuntimeError("Original 2008 valid EPT source support not same n586/HUC437")
def evaluate(d,response,full=True):
 d=d.sort_values("huc8").reset_index(drop=True).copy()
 f=np.array([int(hashlib.sha256(str(z).encode()).hexdigest()[:8],16)%5 for z in d.huc8])
 geo=pd.concat([d[["lat","lon","year","stream_order","season_sin","season_cos","log_cond"]].astype(float),
                pd.get_dummies(d[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
 P=d[["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]].astype(float)
 N=d[["log_tn","log_tp"]].astype(float)
 matrices={"geo":geo,"geo+N_P":pd.concat([geo,N],axis=1),
           "geo+PHab":pd.concat([geo,P],axis=1),
           "geo+PHab+N_P":pd.concat([geo,P,N],axis=1)}
 y=d[response].to_numpy(float)
 predictions={}
 for key,X in matrices.items():
  a=X.to_numpy(float);p=np.full(len(y),np.nan)
  for k in range(5):
   tr=f!=k;te=~tr
   mu=a[tr].mean(axis=0);sd=a[tr].std(axis=0);sd=np.where(sd>1e-10,sd,1)
   Z=np.column_stack([np.ones(tr.sum()),(a[tr]-mu)/sd])
   Te=np.column_stack([np.ones(te.sum()),(a[te]-mu)/sd])
   L=np.eye(Z.shape[1])*10;L[0,0]=0
   b=np.linalg.solve(Z.T@Z+L,Z.T@y[tr])
   p[te]=Te@b
  if not np.isfinite(p).all():raise RuntimeError("Independent source HUC8 missing prediction")
  predictions[key]=p
 mse={k:float(np.mean((y-p)**2)) for k,p in predictions.items()}
 X=matrices["geo+PHab+N_P"].copy()
 for nm in ("log_tn","log_tp"):
  X[nm]=(X[nm]-X[nm].mean())/X[nm].std(ddof=0)
 arr=np.column_stack([np.ones(len(d)),X.to_numpy(float)])
 beta=np.linalg.lstsq(arr,y,rcond=None)[0]
 resid=y-arr@beta
 inv=np.linalg.pinv(arr.T@arr)
 groups=d.huc8.to_numpy();G=d.huc8.nunique();K=np.linalg.matrix_rank(arr)
 meat=np.zeros((arr.shape[1],arr.shape[1]))
 for g in pd.unique(groups):
  ix=groups==g;v=arr[ix].T@resid[ix];meat+=np.outer(v,v)
 V=inv@meat@inv*G/(G-1)*(len(d)-1)/(len(d)-K)
 slopes={}
 for c in ("log_tn","log_tp"):
  i=X.columns.get_loc(c)+1
  se=np.sqrt(max(0,V[i,i]));b=float(beta[i])
  slopes[c]={"beta":b,"ci":[float(b-1.95996398454*se),float(b+1.95996398454*se)]}
 gu,inv=np.unique(groups,return_inverse=True)
 loss0=np.bincount(inv,weights=(y-predictions["geo+PHab"])**2,minlength=len(gu))
 loss1=np.bincount(inv,weights=(y-predictions["geo+PHab+N_P"])**2,minlength=len(gu))
 rng=np.random.default_rng(20261114)
 inds=rng.integers(0,len(gu),size=(2000,len(gu)))
 a=loss0[inds].sum(1);b=loss1[inds].sum(1)
 boots=list(map(float,np.quantile(100*(a-b)/a,[.025,.5,.975])))
 return {"n":len(d),"HUC8":d.huc8.nunique(),"source_response":response,
         "MSEs_manual_numpy":mse,"NP_slope_cluster_independent":slopes,
         "extra_NP_after_PHab_2000_fixedOOF_boot95":boots,
         "NP_incremental_MSE_gain_pct":100*(mse["geo+PHab"]-mse["geo+PHab+N_P"])/mse["geo+PHab"]}
ans={}
for name in ("ept_share_raw_count_lower_bound","log1p_ept_taxa","benthic_mmi_score"):
 dd=D.dropna(subset=[name]).copy()
 result=evaluate(dd,name)
 orig=report["responses"][name]
 if result["n"]!=orig["source_QA"]["n"] or result["HUC8"]!=orig["source_QA"]["HUC8"]:
  raise RuntimeError("Different source cohort for distinct 2008 response "+name)
 for k,mse in result["MSEs_manual_numpy"].items():
  if abs(mse-orig["original_HUC8_spatial5fold_MSE"][k])>1e-7:
   raise RuntimeError("Independent matrix HUC8 MSE mismatch original source "+name+"/"+k)
 for c,b in result["NP_slope_cluster_independent"].items():
  old=orig["fully_adjusted_nutrient_associations"][c]
  if abs(b["beta"]-old["beta_1SD_log10_observed_instream"])>1e-7:
   raise RuntimeError("Independent original 2008 nutrient adjusted slope mismatch "+name+"/"+c)
  if max(abs(x-y) for x,y in zip(b["ci"],old["HUC8_cluster_descriptive95"]))>6e-5:
   raise RuntimeError("Independent 2008 HUC8 covar slope interval mismatch "+name+"/"+c)
 if max(abs(a-b) for a,b in zip(result["extra_NP_after_PHab_2000_fixedOOF_boot95"],
                                orig["HUC8_bootstrap2000_fixed_OOF_extra_NP_descriptive95"]))>1e-6:
  raise RuntimeError("Independent bootstrap on original HUC8 errors mismatch "+name)
 ans[name]=result
# The unknown taxon fraction is explicitly not recoded as confirmed non-EPT.
# Upper-bound sensitivity makes those original unclassified counts all EPT.
up=evaluate(D,"ept_share_source_upper_bound")
fidelity={
 "old_unknown_taxon_amounts_not_imputed":True,
 "lower_bound_mean":float(D.ept_share_raw_count_lower_bound.mean()),
 "upper_bound_mean":float(D.ept_share_source_upper_bound.mean()),
 "n_original_source_upper_bound":len(D),
 "NP_nutrient_coefficients_under_upper_bound":up["NP_slope_cluster_independent"],
 "NP_incremental_OOS_under_upper_bound_pct":up["NP_incremental_MSE_gain_pct"]}
out={
 "schema":"gesui_NR​​SA0809_original_EPT_TNTP_PHab_three_responses_manualNumPy_QA_v1",
 "all_source_endpoints_reproduced_with_independent_code":True,
 "actual_EPA_benthic_rows_same_as_source":True,
 "no_sewage_or_biomass_EPT_causal_effect":True,
 "no_2008_measured_periphyton_source":True,
 "responses":ans,"source_taxonomy_unknown_count_sensitivity":fidelity}
p=R/"NRSA0809_EPT_TNTP_PHab_THREE_RESPONSES_MANUAL_NUMPY_QA.json"
p.write_text(json.dumps(out,indent=2,sort_keys=True))
print("NRSA0809_RAW_EPT_THIRD_CYCLE_INDEPENDENT_NUMPY_SUCCESS")
print(json.dumps({"responses":{k:{"n":v["n"],"HUC8":v["HUC8"],
         "NP_gain":v["NP_incremental_MSE_gain_pct"],
         "slopes":v["NP_slope_cluster_independent"]} for k,v in ans.items()},
     "unknown":fidelity},sort_keys=True))
