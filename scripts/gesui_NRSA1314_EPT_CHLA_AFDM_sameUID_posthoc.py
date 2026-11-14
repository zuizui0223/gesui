#!/usr/bin/env python3
"""2013/14 SOURCE-FROZEN true same NRSA visit: EPT counted fraction vs
original areal AFDM organic dry standing stock and periphyton Chla pigment,
genuine in-stream TN/TP and physical PHab. Posthoc nested comparisons after
independent cross-cycle XEMBED interaction already seen.
"""
from __future__ import annotations
import hashlib,json,os,pathlib
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui1314_independent_numeric"
original=json.loads((ROOT/"NRSA1314_INDEPENDENT_BERW_AREAL_PCHL_AFDM_EMBED_RESULT.json").read_text())
df=pd.read_csv(ROOT/"AUDIT_ONLY_NRSA1314_AREAL_AFDM_SOURCE_PANEL.csv",
 dtype={"UID":str,"huc8":str,"state":str,"eco9":str,"site_id":str}).sort_values("UID").reset_index(drop=True)
if len(df)!=1094 or df.huc8.nunique()!=693 or df.state.nunique()!=47:
 raise RuntimeError("Original pre-frozen 2013 true same-site areal AFDM–Chla–EPT sample differs")
cov_base=["lat","lon","year","season_sin","season_cos","stream_order","logcond",
          "XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","logwidth"]
cov_nutrient=["logTN","logTP"]
exposures=["logchla","logafdm"]
Y=df.ept_fraction.to_numpy(float)
groups=df.huc8.astype(str).to_numpy()
folds=np.asarray([int(hashlib.sha256(g.encode()).hexdigest()[:8],16)%5 for g in groups])
if min(np.bincount(folds,minlength=5))<100:
 raise RuntimeError("Original source HUC8 5fold validation collapsed")
geo=pd.concat([df[cov_base].astype(float),
               pd.get_dummies(df[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
water=df[cov_nutrient].astype(float)
resource=df[exposures].astype(float)
XB={
 "base":pd.concat([geo,water],axis=1),
 "plus_chla":pd.concat([geo,water,resource[["logchla"]]],axis=1),
 "plus_afdm":pd.concat([geo,water,resource[["logafdm"]]],axis=1),
 "plus_both":pd.concat([geo,water,resource],axis=1)
}
def CV(X,y):
 x=X.to_numpy(float);pred=np.full(len(y),np.nan)
 for k in range(5):
  tr=folds!=k;te=~tr
  scaler=StandardScaler().fit(x[tr])
  z=scaler.transform(x)
  pred[te]=Ridge(alpha=10,fit_intercept=True,solver="svd").fit(z[tr],y[tr]).predict(z[te])
 if not np.isfinite(pred).all():raise RuntimeError("Missing heldout authentic EPA response")
 return pred
def bootstrap(y,pa,pb):
 uniq,inv=np.unique(groups,return_inverse=True)
 a=np.bincount(inv,weights=(y-pa)**2,minlength=len(uniq))
 b=np.bincount(inv,weights=(y-pb)**2,minlength=len(uniq))
 rng=np.random.default_rng(20261114)
 samples=rng.integers(0,len(uniq),size=(2000,len(uniq)))
 den=a[samples].sum(1);alt=b[samples].sum(1)
 return list(map(float,np.quantile(100*(den-alt)/den,[.025,.5,.975])))
def cluster(X,y,fields):
 d=X.copy()
 for c in fields:
  d[c]=(d[c]-d[c].mean())/d[c].std(ddof=0)
 fit=sm.OLS(y,sm.add_constant(d.to_numpy(float),has_constant="add")).fit(
   cov_type="cluster",cov_kwds={"groups":groups,"use_correction":True})
 return {c:{"beta_per1SD_source_log10":float(fit.params[d.columns.get_loc(c)+1]),
            "HUC8_cluster_descriptive95":list(map(float,fit.conf_int()[d.columns.get_loc(c)+1]))} for c in fields}
pred={k:CV(x,Y) for k,x in XB.items()}
mse={k:float(np.mean((Y-v)**2)) for k,v in pred.items()}
comparisons={}
for b,a in [("base","plus_chla"),("base","plus_afdm"),
            ("plus_chla","plus_both"),("plus_afdm","plus_both"),
            ("base","plus_both")]:
 comparisons[b+"__"+a]={
   "improvement_MSE_pct":float(100*(mse[b]-mse[a])/mse[b]),
   "fixed_OOF_HUC8_bootstrap95":bootstrap(Y,pred[b],pred[a])}
slopes=cluster(XB["plus_both"],Y,["logTN","logTP","logchla","logafdm"])
outputs={}
for outcome in ("logafdm","logchla"):
 y=df[outcome].to_numpy(float)
 X0=geo
 XN=pd.concat([geo,water],axis=1)
 a=CV(X0,y);b=CV(XN,y)
 MSE0=float(np.mean((y-a)**2));MSE1=float(np.mean((y-b)**2))
 outputs[outcome]={
  "source_original_biomass_or_pigment":"areal standing stock, NOT measured production",
  "N_P_incremental_MSE_improvement_pct":float(100*(MSE0-MSE1)/MSE0),
  "fixedOOF_HUC8_bootstrap95_NP_increment":bootstrap(y,a,b),
  "N_P_adjusted_descriptive_slopes":cluster(XN,y,["logTN","logTP"])}
OUT={
 "schema":"gesui_NR​​SA1314_original_same_1094_EPT_Chla_AREAL_AFDM_NP_PHAB_nested_posthoc_v1",
 "source_verified_2013_14_AFD​​M_mg_cm2_Chla_ug_cm2":True,
 "original_source_1094_UIDs_HUC8_693_states47":True,
 "original_SOURCE_uid_census_SHA256":original["original_preselected_1203_UID_source_manifest_SHA256"],
 "all_models_identical_source_native_uid_cohort":True,
 "no_EPA_WWTP_NP_delivery_or_actual_ecosystem_service_measured":True,
 "field_contrast_posthoc_after_old_rho_and_moderation_known":True,
 "EPT_counted_fraction_not_total_benthic_biomass":True,
 "AFDM_not_annual_algal_production":True,
 "original_n":len(df),"HUC8":df.huc8.nunique(),"states":df.state.nunique(),
 "original_HUC8_5fold_MSE":mse,
 "EPT_resource_nested_heldout_predictivity":comparisons,
 "EPT_adjusted_both_source_1SD_slopes":slopes,
 "source_absolute_ranks":{
   "EPT_fraction_vs_logA​​FDM":float(spearmanr(df.ept_fraction,df.logafdm).statistic),
   "EPT_fraction_vs_logChla":float(spearmanr(df.ept_fraction,df.logchla).statistic),
   "logAFDM_vs_logChla":float(spearmanr(df.logafdm,df.logchla).statistic)},
 "resource_endpoints_vs_ambient_nutrients":outputs}
dst=ROOT/"NRSA1314_EPT_CHLA_AREAL_AFDM_NP_PHAB_NESTED_SOURCE.json"
dst.write_text(json.dumps(OUT,indent=2,sort_keys=True))
print("NRSA1314_EPT_AREAL_AFDM_CHLA_NP_SAME_UID_RESULT",
 json.dumps({"n":len(df),"HUC8":df.huc8.nunique(),
             "rho":OUT["source_absolute_ranks"],"MSE":mse,"comparisons":comparisons,
             "EPT_slopes":slopes,"resources":outputs},sort_keys=True))
