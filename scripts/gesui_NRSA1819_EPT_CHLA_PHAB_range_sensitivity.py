#!/usr/bin/env python3
"""Source-range and fixed sensitivity of true 2018–19 NRSA EPT, Chla and
taxon-ID richness with genuine stream nutrient and PHab covariates.

Posthoc falsification only: first repeated site visit, >=250 counted larvae,
exclude original XGVEG>1, year-specific 2018/19. No p-hacking or adjustment
of original sample to obtain a desired significant sign. Source errors are
distinct from biological nulls. The independent main source script is rerun
to regenerate its unique native UID visit panel; no raw rows are published.
"""
from __future__ import annotations
import hashlib,json,os,pathlib
import numpy as np,pandas as pd
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui1819_triad_numeric"
p=ROOT/"AUDIT_ONLY_2018_19_EPT_Chla_PHAB_SOURCE_PANEL.csv"
orig=ROOT/"NRSA1819_SAME_UID_BIOTIC_TRIAD_PHAB_NP_REAL_RESULTS.json"
source=json.loads(orig.read_text())
df=pd.read_csv(p,dtype={"UID":str,"site_id":str,"huc8":str,"state":str,"eco9":str})
if len(df)!=1221:raise RuntimeError("Original source 1221 UID cohort was modified")
physical=["XPCAN","XGVEG","XFC_OHV","W1_HALL","XCDENMID","XEMBED","XWIDTH","log_width"]
qa={}
for name in physical:
 q=pd.to_numeric(df[name],errors="coerce")
 qa[name]={
  "n_numeric":int(np.isfinite(q).sum()),
  "negative_n":int((q<0).sum()),
  "zero_n":int((q==0).sum()),
  "gt100_n":int((q>100).sum()),
  "gt1_n":int((q>1).sum()),
  "min":float(q.min()) if np.isfinite(q).any() else None,
  "median":float(q.median()) if np.isfinite(q).any() else None,
  "max":float(q.max()) if np.isfinite(q).any() else None}
print("NRSA1819_PHAB_ORIGINAL_NUMERIC_QA",json.dumps(qa,sort_keys=True),flush=True)
if any(qa[n]["negative_n"]>0 for n in ["XPCAN","XGVEG","XFC_OHV","W1_HALL","XCDENMID","XEMBED","XWIDTH"]):
 raise RuntimeError("PHab has negative physical sentinel(s): STOP until original source metric definition independently resolved")

basecols=["lat","lon","year","stream_order","season_sin","season_cos","log_cond","log_tn","log_tp","log_chla",
          "ept_fraction","log_taxon_richness","log_total_subsample_count"]
strictmask=df.chla_flag.fillna("").eq("") & df.tn_mgL_rawflag.fillna("").eq("") & df.tp_mgL_rawflag.fillna("").eq("")
D=df[strictmask].dropna(subset=basecols+["XPCAN","XGVEG","XFC_OHV","W1_HALL","log_width","XCDENMID","XEMBED"]).copy()
D=D[(D.state!="")&(D.huc8!="")&(D.eco9!="")]
if len(D)!=946 or D.huc8.nunique()!=603:
 raise RuntimeError("Strict original n946 HUC8=603 primary physical+biota cohort changed")
D=D.sort_values(["year","UID"]).reset_index(drop=True)
subsets={
 "full_original_946":np.ones(len(D),dtype=bool),
 "first_native_SITE_ID_visit_only":~D.duplicated("site_id"),
 "benthic_counted_subsample_at_least250":D.total_benthic_subsample_count.ge(250).to_numpy(),
 "exclude_original_XGVEG_gt1":D.XGVEG.le(1).to_numpy(),
 "year2018_only":D.year.eq(2018).to_numpy(),
 "year2019_only":D.year.eq(2019).to_numpy()
}
cols_h=["XPCAN","XGVEG","XFC_OHV","W1_HALL","log_width","XCDENMID","XEMBED"]
def calc(D):
 fold=np.array([int(hashlib.sha256(str(x).encode()).hexdigest()[:8],16)%5 for x in D.huc8])
 common=pd.concat([
  D[["lat","lon","year","stream_order","season_sin","season_cos","log_cond"]].astype(float),
  pd.get_dummies(D[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
 physical=D[cols_h].astype(float)
 chem=D[["log_tn","log_tp"]].astype(float)
 rows={}
 for outcome,extra,effort in [("ept_fraction","log_chla",False),
                               ("log_taxon_richness","log_chla",True)]:
  base=common if not effort else common.assign(log_total_subsample_count=D.log_total_subsample_count.to_numpy(float))
  X=pd.concat([base,physical,chem],axis=1)
  XF=X.assign(**{extra:D[extra].to_numpy(float)})
  Y=D[outcome].to_numpy(float)
  def cv(x):
   a=x.to_numpy(float);out=np.full(len(Y),np.nan)
   for f in range(5):
    tr=fold!=f;te=~tr
    s=StandardScaler().fit(a[tr])
    z=s.transform(a)
    out[te]=Ridge(alpha=10,fit_intercept=True,solver="svd").fit(z[tr],Y[tr]).predict(z[te])
   return np.mean((Y-out)**2)
  a,b=cv(X),cv(XF)
  E=XF.copy()
  for c in ("log_tn","log_tp",extra):
   E[c]=(E[c]-E[c].mean())/E[c].std(ddof=0)
  m=sm.OLS(Y,sm.add_constant(E.to_numpy(float),has_constant="add")).fit(
    cov_type="cluster",cov_kwds={"groups":D.huc8.to_numpy(),"use_correction":True})
  def slope(k):
   i=E.columns.get_loc(k)+1
   return {"beta":float(m.params[i]),"cluster95":list(map(float,m.conf_int()[i]))}
  rows[outcome]={
   "N":len(D),"base_N_PHab_MSE":float(a),"Chla_added_MSE":float(b),
   "Chla_added_predictive_gain_pct":float(100*(a-b)/a),
   "signed_coefficients":{"log_Chla":slope("log_chla"),"log_TN":slope("log_tn"),"log_TP":slope("log_tp")}
  }
 return rows
ans={}
for name,mask in subsets.items():
 d=D.loc[mask].copy()
 info={"n":len(d),"physical_SITE_ID_n":d.site_id.nunique(),
       "HUC8":d.huc8.nunique(),"states":d.state.nunique()}
 if len(d)<300 or d.huc8.nunique()<150 or d.state.nunique()<20:
  info["status"]="SOURCE/SUPPORT_GATE_FAIL";ans[name]=info;continue
 info["status"]="SENSITIVITY_FIT_OBSERVATIONAL"
 info["responses"]=calc(d)
 ans[name]=info
R={
 "schema":"gesui_NRSA1819_EPT_TAXONRICHNESS_Chla_PHAB_SOURCE_ORIGINAL_RANGE_POSTHOC_FALSIFICATION_v1",
 "original_phab_QA":qa,
 "original_parent_source_946_strict":True,
 "original_observed_response_not_WWTP_or_real_service":True,
 "posthoc_sensitivity_not_confirmatory_independent_dataset":True,
 "source_zero_vs_missing_separate":True,
 "source_comparisons":ans
}
o=ROOT/"NRSA1819_EPT_Chla_PHAB_RANGE_AND_SENSITIVITY_QA.json"
o.write_text(json.dumps(R,indent=2,sort_keys=True))
print("NRSA1819_EPT_CHLA_PHAB_SENSITIVITY_SUCCESS")
print(json.dumps({"range_QA":qa,"sensitivities":ans},sort_keys=True))
