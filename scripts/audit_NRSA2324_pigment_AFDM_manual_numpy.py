#!/usr/bin/env python3
"""Independent manual NumPy source-data validation of NRSA2023/24 pigment-AFDM
allocation vs real stream nutrients. Source rows are reconstructed from pinned
EPA original files; numerical fits are independent of sklearn/statsmodels.
Confirmatory? NO: this is a posthoc contrast of originally known components.
"""
import hashlib,json,os,pathlib
import numpy as np,pandas as pd
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_biomass_numeric"
src=pd.read_csv(ROOT/"AUDIT_ONLY_2023_24_PANEL.csv",dtype={"UID":str,"huc8":str,"state":str,"eco9":str,"site_id":str})
saved=json.loads((ROOT/"NRSA2324_PIGMENT_PER_ORGANIC_MASS_EXPLORATORY.json").read_text())
keep=src.chla_flag.fillna("").eq("")&src.afdm_flag.fillna("").eq("")&src.tn_mgL_rawflag.fillna("").eq("")&src.tp_mgL_rawflag.fillna("").eq("")
target=["lat","lon","stream_order","year","season_sin","season_cos","log_cond","log_tn","log_tp","log_chla","log_afdm"]
D=src.loc[keep].dropna(subset=target).copy()
D=D[(D.huc8.notna())&(D.huc8!="")&(D.state.notna())&(D.state!="")&(D.eco9.notna())&(D.eco9!="")]
D=D.sort_values("UID").reset_index(drop=True)
if len(D)!=1828 or D.huc8.nunique()!=946 or len(src)!=2129:raise RuntimeError("Original 2023 frozen cohort changed")

D["log_ratio"]=D.log_chla-D.log_afdm
fixed=["lat","lon","stream_order","year","season_sin","season_cos","log_cond"]
dummy=pd.get_dummies(D[["state","eco9"]],drop_first=True,dtype=float)
X0=np.column_stack([D[fixed].to_numpy(dtype=float),dummy.to_numpy(dtype=float)])
XNP=np.column_stack([X0,D[["log_tn","log_tp"]].to_numpy(dtype=float)])
XNP_AFDM=np.column_stack([XNP,D.log_afdm.to_numpy(dtype=float)])
fold=np.asarray([int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5 for x in D.huc8.astype(str)])
def heldout(y,X):
 pred=np.full(len(y),np.nan)
 for k in range(5):
  tr=fold!=k;te=fold==k
  x=X.astype(float).copy()
  mu=x[tr].mean(0);sd=x[tr].std(0)
  x=(x-mu)/np.where(sd>1e-10,sd,1)
  Z=np.column_stack([np.ones(int(tr.sum())),x[tr]])
  T=np.column_stack([np.ones(int(te.sum())),x[te]])
  R=np.diag(np.r_[0.,np.repeat(10.,Z.shape[1]-1)])
  b=np.linalg.solve(Z.T@Z+R,Z.T@y[tr])
  pred[te]=T@b
 return float(np.mean((y-pred)**2))
def slopes(y, extra=None):
 x=XNP if extra is None else XNP_AFDM
 columns=len(x[0])
 a=np.column_stack([np.ones(len(y)),x.astype(float)])
 for pos in range(columns-len((["log_tn","log_tp","log_afdm"] if extra is not None else ["log_tn","log_tp"])),columns):
  sd=a[:,pos+1].std()
  a[:,pos+1]=(a[:,pos+1]-a[:,pos+1].mean())/sd
 coef=np.linalg.lstsq(a,y,rcond=None)[0]
 residual=y-a@coef
 G=D.huc8.nunique();N=len(y);K=np.linalg.matrix_rank(a)
 bread=np.linalg.pinv(a.T@a)
 meat=np.zeros((a.shape[1],a.shape[1]))
 gr=D.huc8.to_numpy()
 for g in pd.unique(gr):
  ix=(gr==g)
  v=a[ix].T@residual[ix]
  meat+=np.outer(v,v)
 V=bread@meat@bread*G/(G-1)*(N-1)/(N-K)
 names=["log_tn","log_tp","log_afdm"] if extra else ["log_tn","log_tp"]
 return {name:{"beta":float(coef[-len(names)+i]),
   "ci":[float(coef[-len(names)+i]-1.959963984540054*np.sqrt(max(0,V[-len(names)+i,-len(names)+i]))),
         float(coef[-len(names)+i]+1.959963984540054*np.sqrt(max(0,V[-len(names)+i,-len(names)+i])))]}
         for i,name in enumerate(names)}
out={"schema":"NRSA2023_24_PIGMENT_AFDM_NUMPY_MANUAL_INDEPENDENT_NUMERIC_QA_v1",
     "strict_source_native_UID_n":len(D),"HUC8":D.huc8.nunique(),
     "observational_not_WWTP_dose":True,"different_numeric_library_than_main":True,
     "per_response":{}}
for key,outname,ext in [("log_ratio","log10_pigment_per_AFDM",None),
                        ("log_chla","log10_Chla_original","log_afdm"),
                        ("log_afdm","log10_AFDM_original",None)]:
 Y=D[key].to_numpy(dtype=float)
 model=saved["responses"][outname]
 b=heldout(Y,X0);n=heldout(Y,XNP)
 eb=model["out_of_HUC8_fivefold_MSE"]
 if abs(b-eb["baseline"])>5e-9 or abs(n-eb["plus_NP"])>5e-9:
  raise RuntimeError("Independent manual NumPy Ridge heldout result differs: "+outname+str((b,n,eb)))
 original=model["posthoc_within_huc8_clustered_slopes"]
 slopes_new=slopes(Y,ext)
 for name,x in slopes_new.items():
  old=original[name]
  if abs(x["beta"]-old["beta_per_1SD_source_log10_predictor"])>5e-8:
   raise RuntimeError("Original adjusted OLS beta mismatch "+outname+"/"+name)
  if max(abs(a-b) for a,b in zip(x["ci"],old["HUC8_cluster_descriptive95"]))>5e-5:
   raise RuntimeError("Original cluster CI discrepancy "+outname+"/"+name)
 out["per_response"][outname]={"manual_ridge_mse_baseline":b,"manual_ridge_mse_NP":n,
                               "raw_data_manual_cluster_slopes":slopes_new}
# Identity: For identical cohort + full additive model, residualised ratio slope
# equals residualised pigment slope minus AFDM slope exactly (no independent replication).
s_chla=slopes(D.log_chla.to_numpy(float))
s_afdm=slopes(D.log_afdm.to_numpy(float))
s_ratio=slopes(D.log_ratio.to_numpy(float))
for name in ["log_tn","log_tp"]:
 err=s_ratio[name]["beta"]-(s_chla[name]["beta"]-s_afdm[name]["beta"])
 if abs(err)>1e-9:raise RuntimeError("Basic algebraic ratio identity failed")
out["algebraic_ratio_slope_checked_not_independent_ecology"]=True
dest=ROOT/"NRSA2324_PIGMENT_AFDM_MANUAL_NUMPY_QA.json"
dest.write_text(json.dumps(out,sort_keys=True,indent=2))
print("NRSA2324_PIGMENT_AFDM_INDEPENDENT_MANUAL_QA_SUCCESS",json.dumps({
 "n":len(D),"source_HUC8":D.huc8.nunique(),
 "ratio_MSE_base":out["per_response"]["log10_pigment_per_AFDM"]["manual_ridge_mse_baseline"],
 "ratio_MSE_NP":out["per_response"]["log10_pigment_per_AFDM"]["manual_ridge_mse_NP"],
 "ratio_TP_slope":s_ratio["log_tp"]["beta"],
 "algebraic_identity":True},sort_keys=True))
