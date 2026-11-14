#!/usr/bin/env python3
"""Independent normal-equation and HUC8 covariance QA of genuine EPA2013/14
BERW EPT~PCHL×embeddedness, plus harmonized original EPA2018/19 BERW panel.

Source cycles are distinct repeated NATIONAL SURVEYS, not necessarily disjoint
river sites. No scientific claim of causal sewage, benthos mass production,
or ecological service. Preserve 2013 weak/null heldout predictive outcome.
"""
from __future__ import annotations
import hashlib,json,os,pathlib
import numpy as np,pandas as pd
from scipy.stats import norm

TEMP=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))
ROOT=TEMP/"gesui1314_independent_numeric"
source=json.loads((ROOT/"NRSA1314_INDEPENDENT_BERW_AREAL_PCHL_AFDM_EMBED_RESULT.json").read_text())
d=pd.read_csv(ROOT/"AUDIT_ONLY_NRSA1314_BERW_EPT_CHLA_SOURCE_PANEL.csv",
              dtype={"UID":str,"huc8":str,"state":str,"eco9":str,"site_id":str})
af=pd.read_csv(ROOT/"AUDIT_ONLY_NRSA1314_AREAL_AFDM_SOURCE_PANEL.csv",
               dtype={"UID":str,"huc8":str,"state":str,"eco9":str,"site_id":str})
if len(d)!=1116 or d.huc8.nunique()!=706 or len(af)!=1094 or af.huc8.nunique()!=693:
 raise RuntimeError("Original independent EPA 2013 source sample size changed")
old=pd.read_csv(TEMP/"gesui1819_triad_numeric"/"AUDIT_ONLY_2018_19_EPT_Chla_PHAB_SOURCE_PANEL.csv",
                dtype={"UID":str,"huc8":str,"state":str,"eco9":str,"site_id":str})
if len(old)!=1221:raise RuntimeError("Old original 2018 fixed source UID census changed")
mask=old.chla_flag.fillna("").eq("")&old.tn_mgL_rawflag.fillna("").eq("")&old.tp_mgL_rawflag.fillna("").eq("")
q=old[mask].rename(columns={"log_cond":"logcond","log_tn":"logTN","log_tp":"logTP",
                             "log_width":"logwidth","log_chla":"logchla"}).copy()
cov=["lat","lon","year","season_sin","season_cos","stream_order","logcond",
     "logTN","logTP","XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","logwidth","logchla"]
q=q.dropna(subset=cov+["ept_fraction"])
q=q[(q.huc8!="")&(q.state!="")&(q.eco9!="")].sort_values("UID").reset_index(drop=True)
if len(q)!=946 or q.huc8.nunique()!=603:
 raise RuntimeError("2018–19 matched BERW original Chla environment n946/603 changed")
# 2013 and 2018 use EXACT same continuous covariates; XGVEG is not shared
# because original 2013 source uses XGV and is not equivalent metric.
def design(D,exposure):
 base=pd.concat([D[cov[:-1]].astype(float),
                 pd.get_dummies(D[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
 base[exposure]=D[exposure].astype(float)
 return base
def manual_cluster(D,exposure):
 X=design(D,exposure)
 zz=(D[exposure].to_numpy(float)-D[exposure].mean())/D[exposure].std(ddof=0)
 emb=(D.XEMBED.to_numpy(float)-D.XEMBED.mean())/D.XEMBED.std(ddof=0)
 X["quad"]=zz*zz
 X["interaction"]=zz*emb
 YY=D.ept_fraction.to_numpy(float)
 A=np.column_stack([np.ones(len(YY)),X.to_numpy(float)])
 b=np.linalg.lstsq(A,YY,rcond=None)[0]
 err=YY-A@b
 inv=np.linalg.pinv(A.T@A)
 meat=np.zeros((A.shape[1],A.shape[1]))
 groups=D.huc8.to_numpy()
 for u in pd.unique(groups):
  idx=groups==u
  v=A[idx].T@err[idx]
  meat+=np.outer(v,v)
 G=D.huc8.nunique();N=len(YY);K=np.linalg.matrix_rank(A)
 V=inv@meat@inv*G/(G-1)*(N-1)/(N-K)
 beta=float(b[-1]);se=float(np.sqrt(max(0,V[-1,-1])))
 return {"beta":beta,"CI":[beta-1.959963984540054*se,beta+1.959963984540054*se],
         "p":float(2*norm.sf(abs(beta/se))) if se>0 else 1.}
def manual_cv(D,exposure,interaction):
 X=design(D,exposure)
 YY=D.ept_fraction.to_numpy(float);arr=X.to_numpy(float)
 ie=X.columns.get_loc(exposure);ii=X.columns.get_loc("XEMBED")
 group=D.huc8.to_numpy()
 folds=np.array([int(hashlib.sha256(str(h).encode()).hexdigest()[:8],16)%5 for h in group])
 pred=np.full(len(YY),np.nan)
 for k in range(5):
  tr=folds!=k;te=~tr
  mean=arr[tr].mean(0);std=arr[tr].std(0)
  Z=(arr-mean)/np.where(std>1e-10,std,1)
  extras=[Z[:,ie]**2]
  if interaction:extras.append(Z[:,ie]*Z[:,ii])
  B=np.column_stack([Z]+extras)
  mu=B[tr].mean(0);sd=B[tr].std(0)
  V=(B-mu)/np.where(sd>1e-10,sd,1)
  A=np.column_stack([np.ones(tr.sum()),V[tr]])
  test=np.column_stack([np.ones(te.sum()),V[te]])
  pen=np.eye(A.shape[1])*10;pen[0,0]=0
  coef=np.linalg.solve(A.T@A+pen,A.T@YY[tr])
  pred[te]=test@coef
 return float(np.mean((YY-pred)**2))
def evaluate(D,exp):
 d=D.sort_values("UID").reset_index(drop=True).copy()
 if len(d)<300 or d.huc8.nunique()<160 or d.state.nunique()<20:
  return {"status":"not_enough_within_original_2013_survey_subcohort","n":len(d),"HUC8":d.huc8.nunique()}
 b=manual_cluster(d,exp)
 m0=manual_cv(d,exp,False);m1=manual_cv(d,exp,True)
 return {"n":len(d),"states":d.state.nunique(),"HUC8":d.huc8.nunique(),
         "physical_site_ID_n":d.site_id.nunique(),
         "Chla_or_AFDM_z_units":"within_survey_source_1SD_of_log10_reported_value",
         "interaction":b,"manual_MSE_Q":m0,"manual_MSE_I":m1,
         "independent_HUC8_heldout_MSE_gain_pct":100*(m0-m1)/m0}
out={"schema":"gesui_NRSA1314_NR​​SA1819_HARMONIZED_EPT_CHLA_EMBED_MANUAL_NUMPY_VALIDATE_20261114_v1",
     "two_independent_survey_cycles_original_EPA":True,
     "both_2013_and_2018_BERW_standardized_300_sample":True,
     "same_six_PHAB_and_core_source_covariates":True,
     "2013_14_original_source_harmonized":evaluate(d,"logchla"),
     "2013_14_original_areal_AFDM_secondary":evaluate(af,"logafdm"),
     "2018_19_original_source_harmonized":evaluate(q,"logchla"),
     "2013_only_original_2013":evaluate(d[d.year==2013],"logchla"),
     "2013_only_original_2014":evaluate(d[d.year==2014],"logchla"),
     "2013_unique_physical_site_first_visit":evaluate(d.sort_values(["year","UID"]).drop_duplicates("site_id"),"logchla"),
     "2013_original_benthic_TOTAL300_ge250":evaluate(d[d.total_benthic_counted>=250],"logchla"),
     "2013_HUC8_exclusive_of_all_2018_valid_HUC8":evaluate(d[~d.huc8.isin(set(q.huc8))],"logchla"),
     "2018_HUC8_exclusive_of_all_2013_valid_HUC8":evaluate(q[~q.huc8.isin(set(d.huc8))],"logchla"),
     "no_independent_experimental_confirmation_or_WWTP_attribution":True,
     "outcome_based_site_rediscovery_or_drop_positive_sign_not_used":True}
for key,m,truth in [
 ("2013_14_original_source_harmonized",out["2013_14_original_source_harmonized"],source["2013_14_biofilmChla_EPT_embed_primary"]),
 ("2013_14_original_areal_AFDM_secondary",out["2013_14_original_areal_AFDM_secondary"],source["2013_14_arealAFDM_EPT_embed_secondary"])]:
 if not (m["n"]==truth["n"] and abs(m["interaction"]["beta"]-truth["interaction"]["beta"])<1e-7
    and abs(m["manual_MSE_Q"]-truth["MSE_Q"])<1e-7
    and abs(m["manual_MSE_I"]-truth["MSE_interaction"])<1e-7):
  raise RuntimeError(f"Independent numeric method disagrees with original real EPA 2013 cycle {key}")
 for a,b in zip(m["interaction"]["CI"],truth["interaction"]["ci"]):
  if abs(a-b)>5e-5:raise RuntimeError("Independent ecological field HUC8 cluster uncertainty mismatch")
# Across-cycle overlap in HUC8 does not make repeat surveys experimental or
# independently located. Disclose source catchment overlap explicitly.
out["shared_2013_2018_HUC8"]=len(set(d.huc8)&set(q.huc8))
out["separate_2013_HUC8"]=len(set(d.huc8)-set(q.huc8))
out["separate_2018_HUC8"]=len(set(q.huc8)-set(d.huc8))
p=ROOT/"NRSA1314_1819_HARMONIZED_EPT_CHLA_EMBED_INDEPENDENT_NUMPY.json"
p.write_text(json.dumps(out,indent=2,sort_keys=True))
print("NRSA1314_1819_HARMONIZED_INDEPENDENT_NUMPY_PASS",
      json.dumps({k:v for k,v in out.items() if k in (
       "2013_14_original_source_harmonized","2013_14_original_areal_AFDM_secondary",
       "2018_19_original_source_harmonized","2013_only_original_2013","2013_only_original_2014",
       "2013_unique_physical_site_first_visit","2013_original_benthic_TOTAL300_ge250",
       "2013_HUC8_exclusive_of_all_2018_valid_HUC8","2018_HUC8_exclusive_of_all_2013_valid_HUC8",
       "shared_2013_2018_HUC8","separate_2013_HUC8","separate_2018_HUC8")},sort_keys=True))
