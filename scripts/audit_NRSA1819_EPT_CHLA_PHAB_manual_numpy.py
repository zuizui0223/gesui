#!/usr/bin/env python3
"""Independent manual NumPy audit of original genuine NRSA 2018–19 same-UID
stream N/P, benthic EPT and taxon-ID richness, and periphyton areal Chla,
with real field PHab. Uses new code/matrix and ephemeral original source rows.

All 2 PHab tiers x 3 responses x 5 models heldout MSEs, clustered OLS
and fixed-OOF HUC8 bootstrap must agree. No independently sampled field data,
sewage emissions or demonstrated mediation/ecosystem services.
"""
from __future__ import annotations
import hashlib,json,os,pathlib
import numpy as np,pandas as pd
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui1819_triad_numeric"
PANEL=ROOT/"AUDIT_ONLY_2018_19_EPT_Chla_PHAB_SOURCE_PANEL.csv"
RES=ROOT/"NRSA1819_SAME_UID_BIOTIC_TRIAD_PHAB_NP_REAL_RESULTS.json"
orig=json.loads(RES.read_text())
if orig["schema"]!="gesui_NRSA1819_TRUE_SAME_UID_EPT_Chla_TAXONRICHNESS_PHAB_NP_OBSERVATIONAL_MACRO_v1":
 raise RuntimeError("Original ecological source-backed study schema changed")
src=pd.read_csv(PANEL,dtype={"UID":str,"huc8":str,"state":str,"eco9":str,"site_id":str})
if len(src)!=1221 or src.UID.duplicated().any():raise RuntimeError("Authentic EPA source-locked UID panel changed")
base_vars=["lat","lon","year","stream_order","season_sin","season_cos","log_cond","log_tn","log_tp","log_chla",
           "ept_fraction","log_taxon_richness","log_total_subsample_count"]
qual=src.chla_flag.fillna("").eq("")&src.tn_mgL_rawflag.fillna("").eq("")&src.tp_mgL_rawflag.fillna("").eq("")
src=src[qual].dropna(subset=base_vars)
src=src[(src.huc8!="")&(src.state!="")&(src.eco9!="")]
if len(src)!=947 or src.huc8.nunique()!=603:raise RuntimeError("Original same-visit strict 2018 physiological cohort changed")

PHABS=[
 ("broad",["XPCAN","XGVEG","XFC_OHV","W1_HALL"]),
 ("dense",["XPCAN","XGVEG","XFC_OHV","W1_HALL","log_width","XCDENMID","XEMBED"])]
RESP=[
 ("ept_fraction","log_chla",False),
 ("log_chla","ept_fraction",False),
 ("log_taxon_richness","log_chla",True)
]
def cv(X,Y,fold):
 arr=X.to_numpy(float)
 predictions=np.empty(len(Y))
 for k in range(5):
  tr=(fold!=k);te=~tr
  mu=arr[tr].mean(axis=0);sd=arr[tr].std(axis=0)
  sd=np.where(sd>1e-12,sd,1)
  A=np.column_stack([np.ones(tr.sum()),(arr[tr]-mu)/sd])
  C=np.column_stack([np.ones(te.sum()),(arr[te]-mu)/sd])
  R=np.eye(A.shape[1])*10;R[0,0]=0
  b=np.linalg.solve(A.T@A+R,A.T@Y[tr])
  predictions[te]=C@b
 return predictions
def manual_cluster(X,Y,groups,names):
 XX=X.copy()
 for name in names:
  XX[name]=(XX[name]-XX[name].mean())/XX[name].std(ddof=0)
 A=np.column_stack([np.ones(len(Y)),XX.to_numpy(float)])
 b=np.linalg.lstsq(A,Y,rcond=None)[0]
 resid=Y-A@b
 inv=np.linalg.pinv(A.T@A)
 meat=np.zeros((A.shape[1],A.shape[1]))
 for group in np.unique(groups):
  ix=(groups==group)
  v=A[ix].T@resid[ix]
  meat+=np.outer(v,v)
 N=len(Y);G=len(np.unique(groups));K=np.linalg.matrix_rank(A)
 V=inv@meat@inv*G/(G-1)*(N-1)/(N-K)
 ret={}
 for name in names:
  i=XX.columns.get_loc(name)+1
  se=np.sqrt(max(V[i,i],0))
  ret[name]={"beta":float(b[i]),"ci":[float(b[i]-1.95996398454*se),float(b[i]+1.95996398454*se)]}
 return ret
def boot95(Y,a,b,groups):
 ids,inv=np.unique(groups,return_inverse=True)
 e0=np.bincount(inv,weights=(Y-a)**2,minlength=len(ids))
 e1=np.bincount(inv,weights=(Y-b)**2,minlength=len(ids))
 rng=np.random.default_rng(20261114)
 samples=rng.integers(0,len(ids),size=(1600,len(ids)))
 p=e0[samples].sum(axis=1);q=e1[samples].sum(axis=1)
 return list(map(float,np.quantile(100*(p-q)/p,[.025,.5,.975])))
audited={"schema":"gesui_NRSA1819_TRUE_EPT_Chla_taxonRichness_PHAB_manualNumPy_OOS_cluster_qa_v1",
 "independent_formula_not_field_dataset":True,"no_EPA_WWTP_attribution":True,
 "observed_same_native_UIDs":947,"tiers":{}}
for tier,cols in PHABS:
 d=src.dropna(subset=cols).sort_values("UID").copy()
 expected=orig["results"][tier]
 if expected["status"]!="OBSERVATIONAL_SAME_UID_FITS_COMPLETE":
  raise RuntimeError("Original tier didn't fit, cannot audit")
 if len(d)!=expected["source_QA"]["n_joined_complete_PHAB"] or d.huc8.nunique()!=expected["source_QA"]["n_HUC8"]:
  raise RuntimeError("Source original PHab tier sample changed")
 fold=np.array([int(hashlib.sha256(str(i).encode()).hexdigest()[:8],16)%5 for i in d.huc8])
 groups=d.huc8.to_numpy()
 audit={}
 for outcome,extra,richen in RESP:
  bcols=["lat","lon","year","stream_order","season_sin","season_cos","log_cond"]
  if richen:bcols+=["log_total_subsample_count"]
  common=pd.concat([d[bcols].astype(float),
                    pd.get_dummies(d[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
  N=d[["log_tn","log_tp"]].astype(float)
  H=d[cols].astype(float)
  X={"base":common,"N":pd.concat([common,N],axis=1),"H":pd.concat([common,H],axis=1),
      "NH":pd.concat([common,H,N],axis=1)}
  X["NH+other_biota"]=X["NH"].assign(**{extra:d[extra].to_numpy(float)})
  Y=d[outcome].to_numpy(float)
  pred={name:cv(xx,Y,fold) for name,xx in X.items()}
  errors={name:float(np.mean((Y-values)**2)) for name,values in pred.items()}
  prev=expected["responses"][outcome]["mse_HUC8_group5fold"]
  for key in errors:
   if abs(errors[key]-prev[key])>1e-7:
    raise RuntimeError(f"Independent HUC8 OOF MSE failed {tier}/{outcome}/{key}: {errors[key]} old {prev[key]}")
  labels=["log_tn","log_tp",extra]
  slopes=manual_cluster(X["NH+other_biota"],Y,groups,labels)
  prevSlopes=expected["responses"][outcome]["clustered_descriptive_association"]
  for lab,val in slopes.items():
   old=prevSlopes[lab]
   if abs(val["beta"]-old["beta_1SD_original_predictor"])>1e-7:
    raise RuntimeError(f"Independent HUC8 slope failed {tier}/{outcome}/{lab}")
   if max(abs(a-b) for a,b in zip(val["ci"],old["HUC8_cluster_descriptive95"]))>4e-5:
    raise RuntimeError(f"Independent HUC8 cluster CI failed {tier}/{outcome}/{lab}")
  b0=boot95(Y,pred["H"],pred["NH"],groups)
  b1=boot95(Y,pred["NH"],pred["NH+other_biota"],groups)
  original95=expected["responses"][outcome]["HUC8_fixedheldout_bootstrap95"]
  for test,truth in [(b0,original95["N_after_PHAB_pct"]),(b1,original95["other_biota_after_N_PHAB_pct"])]:
   if max(abs(a-b) for a,b in zip(test,truth))>1e-6:
    raise RuntimeError(f"Independent source HUC8 fixed OOF bootstrap failed {tier}/{outcome}")
  audit[outcome]={"n":len(d),"manualMSE":errors,
                  "manual_Ols_coefficients_and_95":slopes,
                  "HUC8_bootstrap_N_after_PHab":b0,
                  "HUC8_bootstrap_other_biota_after_PHab_N":b1}
 audited["tiers"][tier]={"n":len(d),"HUC8":d.huc8.nunique(),"responses":audit}
out=ROOT/"NRSA1819_TRUE_EPT_CHLA_PHAB_MANUAL_NUMPY_QA.json"
out.write_text(json.dumps(audited,indent=2,sort_keys=True))
print("NRSA1819_TRUE_TRIAD_INDEPENDENT_MANUAL_QA_SUCCESS")
print(json.dumps({
 tier:{"n":v["n"],"HUC8":v["HUC8"],
 "responses":{k:{"slopes":a["manual_Ols_coefficients_and_95"],
                "N_after_H_boot95":a["HUC8_bootstrap_N_after_PHab"],
                "other_biota_after_NH_boot95":a["HUC8_bootstrap_other_biota_after_PHab_N"]}
             for k,a in v["responses"].items()}} for tier,v in audited["tiers"].items()},sort_keys=True))
