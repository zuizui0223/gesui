#!/usr/bin/env python3
"""Independent source-SHA and manual NumPy heldout/cluster QA for NRSA2324
PHab, areal AFDM/Chla and observed TN/TP. Exactly same original sites and
models as prior PR76 primary, but no sklearn, statsmodels or their helpers.

Also yields descriptive HUC8-resampled intervals for fixed heldout contrasts,
not causal uncertainty or independent field replication.
"""
import hashlib,json,os,pathlib,urllib.request, time
import numpy as np,pandas as pd

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_biomass_numeric"
SRC_URL="https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_physicalhabitat_metrics.csv"
SRC_SHA="8a6094bb87f8a57cf45a35a1272b137e6a9327cefd158f4bee88b025eaf9b2f4"
original=json.loads((ROOT/"NRSA2324_PHAB_NUTRIENT_PIGMENT_AREAL_BIOMASS_POSTHOC.json").read_text())
raw=pd.read_csv(ROOT/"AUDIT_ONLY_2023_24_PANEL.csv",
                dtype={"UID":str,"state":str,"eco9":str,"huc8":str,"site_id":str})
p=ROOT/"independent_PHAB_source.csv"
for i in range(3):
 try:
  req=urllib.request.Request(SRC_URL,headers={"User-Agent":"gesui-independent-NRSA2324-PHab-NP-check"})
  with urllib.request.urlopen(req,timeout=110) as r,p.open("wb") as f:
   while block:=r.read(1048576):f.write(block)
  if hashlib.sha256(p.read_bytes()).hexdigest()!=SRC_SHA:raise RuntimeError("Official EPA PHab SHA changed")
  break
 except Exception:
  if i==2:raise
  time.sleep(2**i)
phab=pd.read_csv(p,dtype=str,keep_default_na=False).set_index("UID")
p.unlink()
if phab.index.duplicated().any() or len(raw)!=2129:
 raise RuntimeError("Source original UID duplication/length")
widecols=["XPCAN","XGVEG","XFC_OHV","log_width"]
densecols=widecols+["XCDENMID","XEMBED"]
for name in set(widecols+densecols)-{"log_width"}:
 phab[name]=pd.to_numeric(phab[name],errors="coerce")
phab["XWIDTH"]=pd.to_numeric(phab["XWIDTH"],errors="coerce")
phab["log_width"]=np.log10(phab.XWIDTH.where(phab.XWIDTH>0))
data=raw.set_index("UID").join(phab[densecols],how="left")
keep=data.chla_flag.fillna("").eq("")&data.afdm_flag.fillna("").eq("")&data.tn_mgL_rawflag.fillna("").eq("")&data.tp_mgL_rawflag.fillna("").eq("")
core=["log_tn","log_tp","log_cond","lat","lon","year","stream_order","season_sin","season_cos","log_chla","log_afdm"]
data=data[keep].dropna(subset=core)
data=data[(data.huc8!="")&(data.state!="")&(data.eco9!="")]
if len(data)!=1828:raise RuntimeError("Original physiological quality source cohort changed")
data["ratio"]=data.log_chla-data.log_afdm

def design(d,cols):
 common=pd.concat([d[["lat","lon","year","stream_order","season_sin","season_cos","log_cond"]].astype(float),
    pd.get_dummies(d[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
 phy=d[cols].astype(float);nut=d[["log_tn","log_tp"]].astype(float)
 return {"baseline":common,"nutrients":pd.concat([common,nut],axis=1),
   "physical":pd.concat([common,phy],axis=1),
   "physical+nutrients":pd.concat([common,phy,nut],axis=1)}
def fold(h):
 return int(hashlib.sha256(str(h).encode()).hexdigest()[:8],16)%5
def manual_cv(X,y,F):
 a=X.to_numpy(float);yp=np.full(len(y),np.nan)
 for k in range(5):
  t=F!=k;e=~t
  mu=a[t].mean(0);sd=a[t].std(0)
  sd=np.where(sd>1.e-10,sd,1)
  ztr=(a[t]-mu)/sd;zte=(a[e]-mu)/sd
  A=np.column_stack([np.ones(t.sum()),ztr])
  B=np.column_stack([np.ones(e.sum()),zte])
  pen=np.eye(A.shape[1])*10;pen[0,0]=0
  theta=np.linalg.solve(A.T@A+pen,A.T@y[t])
  yp[e]=B@theta
 if not np.isfinite(yp).all():raise RuntimeError("Missing independently heldout prediction")
 return yp
def manual_slope(X,y,groups):
 A=X.copy()
 for nm in ["log_tn","log_tp"]:
  A[nm]=(A[nm]-A[nm].mean())/A[nm].std(ddof=0)
 xx=np.column_stack([np.ones(len(y)),A.to_numpy(float)])
 b=np.linalg.lstsq(xx,y,rcond=None)[0]
 residual=y-xx@b
 inv=np.linalg.pinv(xx.T@xx)
 meat=np.zeros((xx.shape[1],xx.shape[1]))
 for g in pd.unique(groups):
  q=groups==g;v=xx[q].T@residual[q];meat+=np.outer(v,v)
 N=len(y);G=len(pd.unique(groups));K=np.linalg.matrix_rank(xx)
 V=inv@meat@inv*G/(G-1)*(N-1)/(N-K)
 out={}
 for nm in ["log_tn","log_tp"]:
  ix=A.columns.get_loc(nm)+1
  se=np.sqrt(max(V[ix,ix],0))
  out[nm]={"beta":float(b[ix]),"95":[float(b[ix]-1.95996398454*se),float(b[ix]+1.95996398454*se)]}
 return out
def bootstrap_gain(y,ph,phnp,blocks):
 ids,inv=np.unique(blocks,return_inverse=True)
 N=len(ids)
 cnt=np.bincount(inv,minlength=N)
 e0=np.bincount(inv,weights=(y-ph)**2,minlength=N)
 e1=np.bincount(inv,weights=(y-phnp)**2,minlength=N)
 rng=np.random.default_rng(20261114)
 vals=[]
 for i in range(2000):
  j=rng.integers(0,N,size=N)
  z=e0[j].sum()
  if z<=0:continue
  vals.append(100*(z-e1[j].sum())/z)
 return [float(x) for x in np.quantile(vals,[.025,.5,.975])]
result={"schema":"NRSA2324_PHab_NP_3responses_manualNumPy_independent_QA_and_descriptive_HUC8_gain_CI",
 "original_source_SHA":SRC_SHA,"original_same_native_UID_cohort":len(raw),
 "different_numeric_estimator":True,"not_new_field_data":True,
 "not_experimental_or_WWTP_attribution":True,"tiers":{}}
for tier,cols in (("wide",widecols),("dense_substrate",densecols)):
 d=data.dropna(subset=cols).sort_index().copy()
 d=d.reset_index()
 stored=original["models"][tier]
 if len(d)!=stored["complete_case_n"] or d.huc8.nunique()!=stored["n_source_HUC8"]:
  raise RuntimeError(f"Original {tier} fixed physical source sample changed")
 X=design(d,cols);groups=d.huc8.to_numpy()
 F=np.array([fold(x) for x in groups])
 results={}
 for resp,rawkey in [("log_pigment_per_AFDM","ratio"),("log_chla","log_chla"),("log_afdm","log_afdm")]:
  y=d[rawkey].to_numpy(float)
  predictions={m:manual_cv(x,y,F) for m,x in X.items()}
  got={m:float(np.mean((p-y)**2)) for m,p in predictions.items()}
  expected=stored["responses"][resp]["MSEs_HUC8_fivefold"]
  for k in got:
   if abs(got[k]-expected[k])>1e-6:
    raise RuntimeError(f"Independent manual Ridge disagreement: {tier}/{resp}/{k} {got[k]} vs {expected[k]}")
  slopes=manual_slope(X["physical+nutrients"],y,groups)
  old=stored["responses"][resp]["TP_TN_slope_PHab_control"]
  for nm in ("log_tp","log_tn"):
   if abs(slopes[nm]["beta"]-old[nm]["beta_per_1SD_log10_nutrient"])>1e-7:
    raise RuntimeError(f"Independent PHab OLS TP/TN slope differs: {tier}/{resp}/{nm}")
   for z,ref in zip(slopes[nm]["95"],old[nm]["HUC8_cluster_descriptive95"]):
    if abs(z-ref)>0.00005:raise RuntimeError("HUC8 sandwich uncertainty different")
  results[resp]={"manual_mse":got,"verified_slopes":slopes,
    "HUC8_bootstrap_2000_fixedOOF_NP_added_to_PHab_MSE_gain_95":
       bootstrap_gain(y,predictions["physical"],predictions["physical+nutrients"],groups),
    "HUC8_bootstrap_2000_fixedOOF_PHab_added_to_NP_MSE_gain_95":
       bootstrap_gain(y,predictions["nutrients"],predictions["physical+nutrients"],groups)}
 result["tiers"][tier]={"n":len(d),"HUC8":d.huc8.nunique(),"responses":results}
dest=ROOT/"NRSA2324_PHAB_NP_INDEPENDENT_MANUAL_QA.json"
dest.write_text(json.dumps(result,sort_keys=True,indent=2))
print("PHAB_INDEPENDENT_NUMPY_QA_SUCCESS")
print(json.dumps({tier:{"n":o["n"],"HUC8":o["HUC8"],"response":{m:{
 "NP_gain_HUC8_boot95":a["HUC8_bootstrap_2000_fixedOOF_NP_added_to_PHab_MSE_gain_95"],
 "PHab_gain_HUC8_boot95":a["HUC8_bootstrap_2000_fixedOOF_PHab_added_to_NP_MSE_gain_95"],
 "TP_beta":a["verified_slopes"]["log_tp"]["beta"]} for m,a in o["responses"].items()}}
 for tier,o in result["tiers"].items()},sort_keys=True))
