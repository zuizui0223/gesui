#!/usr/bin/env python3
"""Independent NumPy Ridge10, HUC8 cluster sandwich QA of actual EPA 2013/14
same-UID EPT vs areal AFDM and Chla, vs source measured ambient N/P.
This uses the SAME field data independently reconstructed in the runner, not
an independently sampled field replication or randomized wastewater effect.
"""
import hashlib,json,os,pathlib
import numpy as np,pandas as pd

R=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui1314_independent_numeric"
D=pd.read_csv(R/"AUDIT_ONLY_NRSA1314_AREAL_AFDM_SOURCE_PANEL.csv",
              dtype={"UID":str,"huc8":str,"state":str,"eco9":str,"site_id":str}).sort_values("UID").reset_index(drop=True)
original=json.loads((R/"NRSA1314_EPT_CHLA_AREAL_AFDM_NP_PHAB_NESTED_SOURCE.json").read_text())
if len(D)!=1094 or D.huc8.nunique()!=693 or original["original_n"]!=1094:
 raise RuntimeError("Genuine original EPA 2013 strict source n1094/693 modified")
a=["lat","lon","year","season_sin","season_cos","stream_order","logcond",
   "XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","logwidth"]
ch=["logTN","logTP"]
dum=pd.get_dummies(D[["state","eco9"]],drop_first=True,dtype=float)
geo=pd.concat([D[a].astype(float),dum],axis=1)
nu=D[ch].astype(float)
phys=D[["logchla","logafdm"]].astype(float)
models={
 "base":pd.concat([geo,nu],axis=1),
 "plus_chla":pd.concat([geo,nu,phys[["logchla"]]],axis=1),
 "plus_afdm":pd.concat([geo,nu,phys[["logafdm"]]],axis=1),
 "plus_both":pd.concat([geo,nu,phys],axis=1)}
fold=np.asarray([int(hashlib.sha256(str(x).encode()).hexdigest()[:8],16)%5 for x in D.huc8])
group=D.huc8.to_numpy()
def cv(X,y):
 A=X.to_numpy(float);P=np.full(len(y),np.nan)
 for k in range(5):
  t=fold!=k;e=~t
  m=A[t].mean(0);s=A[t].std(0);s=np.where(s>1e-12,s,1)
  u=(A[t]-m)/s;v=(A[e]-m)/s
  Q=np.column_stack([np.ones(t.sum()),u])
  Z=np.column_stack([np.ones(e.sum()),v])
  ridge=np.eye(Q.shape[1])*10;ridge[0,0]=0
  bet=np.linalg.solve(Q.T@Q+ridge,Q.T@y[t])
  P[e]=Z@bet
 if not np.isfinite(P).all():raise RuntimeError("Original source HUC8 OOF missing")
 return float(np.mean((y-P)**2))
def coeff(X,y,parts):
 C=X.copy()
 for v in parts:C[v]=(C[v]-C[v].mean())/C[v].std(ddof=0)
 mat=np.column_stack([np.ones(len(y)),C.to_numpy(float)])
 beta=np.linalg.lstsq(mat,y,rcond=None)[0];err=y-mat@beta
 bread=np.linalg.pinv(mat.T@mat)
 meat=np.zeros((mat.shape[1],mat.shape[1]))
 for h in pd.unique(group):
  t=group==h;v=mat[t].T@err[t];meat+=np.outer(v,v)
 G=len(pd.unique(group));N=len(y);K=np.linalg.matrix_rank(mat)
 V=bread@meat@bread*G/(G-1)*(N-1)/(N-K)
 out={}
 for v in parts:
  i=C.columns.get_loc(v)+1
  sd=np.sqrt(max(0,V[i,i]))
  out[v]={"beta":float(beta[i]),"95":[float(beta[i]-1.959963984540054*sd),float(beta[i]+1.959963984540054*sd)]}
 return out
Y=D.ept_fraction.to_numpy(float)
mse={k:cv(X,Y) for k,X in models.items()}
for k in mse:
 if abs(mse[k]-original["original_HUC8_5fold_MSE"][k])>1e-7:
  raise RuntimeError("Independent source-native original MSE mismatch "+k)
bet=coeff(models["plus_both"],Y,["logTN","logTP","logchla","logafdm"])
for name,obj in bet.items():
 orig=original["EPT_adjusted_both_source_1SD_slopes"][name]
 if abs(obj["beta"]-orig["beta_per1SD_source_log10"])>1e-7:
  raise RuntimeError("Independent original cross-layer coefficient mismatch "+name)
 if max(abs(x-y) for x,y in zip(obj["95"],orig["HUC8_cluster_descriptive95"]))>5e-5:
  raise RuntimeError("Independent HUC8 robust interval mismatch "+name)
resource={}
for name in ["logafdm","logchla"]:
 y=D[name].to_numpy(float)
 b=cv(geo,y);n=cv(pd.concat([geo,nu],axis=1),y)
 source=original["resource_endpoints_vs_ambient_nutrients"][name]
 gain=100*(b-n)/b
 if abs(gain-source["N_P_incremental_MSE_improvement_pct"])>1e-7:
  raise RuntimeError("Independent source resource-yield N/P CV disagreement "+name)
 resource[name]={"baseline_MSE":b,"plus_NP_MSE":n,"NP_OOS_gain_pct":gain}
O={"schema":"gesui_NRSA1314_original_same_1094_EPT_Chla_AFDM_ambientNP_independent_matrix_Numpy_QA",
   "original_1094_UID_HUC8_693_verified":True,
   "original_source_separate_numeric_implementation":True,
   "independent_ecological_dataset_claim":False,
   "not_WWTP_exposure_or_food_web_function":True,
   "per_resource_model_MSE":mse,
   "per_source_1SD_adjusted_EPT_slopes":bet,
   "same_field_resource_outcome_NP_prediction":resource}
out=R/"NRSA1314_TRUE_EPT_CHLA_AFDM_NP_MATRIX_QA.json"
out.write_text(json.dumps(O,indent=2,sort_keys=True))
print("NRSA1314_EPT_AFDM_CHLA_NP_INDEPENDENT_NUMPY_PASS",
      json.dumps({"MSE":mse,"EPT_partial":bet,"biofilm_y_NP":resource},sort_keys=True))
