#!/usr/bin/env python3
"""Independent NumPy ridge/OLS check of five fixed PHab × Chla interactions
and predeclared habitat/range/year sensitivity on original NRSA1819 n946.

No helper functions from primary implementation used, same raw pinned EPA
survey measurements reconstructed separately by preceding GH workflow steps.
Observational, after known negative EPT–Chla main-effect, no causal wastewater
or ecosystem-service mechanisms, no independent field replication.
"""
from __future__ import annotations
import hashlib,json,os,pathlib
import numpy as np,pandas as pd
from scipy.stats import norm

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui1819_triad_numeric"
MAIN=ROOT/"NRSA1819_EPT_CHLA_PHAB_INTERACTION_SOURCE_SIGNED_EXPLORATORY.json"
SOURCE=ROOT/"AUDIT_ONLY_2018_19_EPT_Chla_PHAB_SOURCE_PANEL.csv"
main=json.loads(MAIN.read_text())
df=pd.read_csv(SOURCE,dtype={"UID":str,"huc8":str,"state":str,"eco9":str,"site_id":str})
if len(df)!=1221 or df.UID.duplicated().any():raise RuntimeError("Original 1221 visits changed")
num=["lat","lon","year","stream_order","season_sin","season_cos","log_cond",
     "log_tn","log_tp","XPCAN","XGVEG","XFC_OHV","W1_HALL","log_width",
     "XCDENMID","XEMBED","log_chla","ept_fraction","total_benthic_subsample_count"]
quality=(df.chla_flag.fillna("").eq("")&
         df.tn_mgL_rawflag.fillna("").eq("")&
         df.tp_mgL_rawflag.fillna("").eq(""))
d=df.loc[quality].dropna(subset=num).copy()
d=d[(d.state!="")&(d.eco9!="")&(d.huc8!="")].sort_values("UID").reset_index(drop=True)
if len(d)!=946 or d.huc8.nunique()!=603 or d.state.nunique()!=47:
 raise RuntimeError("Original n946/603 source after numeric QA mismatch")
if main["cohort"]!=946 or main["HUC8"]!=603:
 raise RuntimeError("Original first-stage interaction result invalid source")
mods=["XCDENMID","XEMBED","XPCAN","XGVEG","log_width"]
cov=["lat","lon","year","stream_order","season_sin","season_cos","log_cond",
     "log_tn","log_tp","XPCAN","XGVEG","XFC_OHV","W1_HALL","log_width",
     "XCDENMID","XEMBED","log_chla"]
dummy=pd.get_dummies(d[["state","eco9"]],drop_first=True,dtype=float)
controls=pd.concat([d[cov].astype(float),dummy],axis=1)
Y=d.ept_fraction.to_numpy(float)
group=d.huc8.to_numpy()
sc={k:(float(d[k].mean()),float(d[k].std(ddof=0))) for k in mods+["log_chla"]}
def model_ols(D,Xbase,field,fullscale=sc):
 A=Xbase.copy()
 ch=(D.log_chla.to_numpy(float)-fullscale["log_chla"][0])/fullscale["log_chla"][1]
 ec=(D[field].to_numpy(float)-fullscale[field][0])/fullscale[field][1]
 A["chla_z_squared"]=ch*ch
 A["chla_z_times_"+field]=ch*ec
 mat=np.column_stack([np.ones(len(D)),A.to_numpy(float)])
 y=D.ept_fraction.to_numpy(float)
 b=np.linalg.lstsq(mat,y,rcond=None)[0]
 residual=y-mat@b
 bread=np.linalg.pinv(mat.T@mat)
 meat=np.zeros((mat.shape[1],mat.shape[1]))
 groups=D.huc8.to_numpy()
 for g in np.unique(groups):
  loc=groups==g
  u=mat[loc].T@residual[loc]
  meat+=np.outer(u,u)
 G=len(np.unique(groups));N=len(y);K=np.linalg.matrix_rank(mat)
 V=bread@meat@bread*G/(G-1)*(N-1)/(N-K)
 k=A.columns.get_loc("chla_z_times_"+field)+1
 se=float(np.sqrt(max(V[k,k],0)))
 beta=float(b[k])
 return {"n":len(D),"HUC8":G,"beta_interaction_per1SD_x_1SD":beta,
         "HUC8_cluster_normal95":[beta-1.959963984540054*se,beta+1.959963984540054*se],
         "cluster_normal_p":float(2*norm.sf(abs(beta/se))) if se>0 else 1,
         "source":field}
def cv(D,mods):
 # Independent manual closed-form Ridge10, train-fold scaling and product
 # construction. No test outcomes used in training/standardization.
 fields=pd.concat([D[cov].astype(float),
                   pd.get_dummies(D[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
 # Repeated category columns matched by ordered values in full n946 model.
 x=fields.to_numpy(float)
 y=D.ept_fraction.to_numpy(float)
 gr=D.huc8.to_numpy()
 fold=np.array([int(hashlib.sha256(str(h).encode()).hexdigest()[:8],16)%5 for h in gr])
 out=np.full(len(y),np.nan)
 for k in range(5):
  tr=fold!=k;te=~tr
  if tr.sum()<200 or te.sum()<30:raise RuntimeError("Fold support insufficient in sensitivity")
  rawmu=x[tr].mean(0);rawsd=x[tr].std(0)
  z=(x-rawmu)/np.where(rawsd>1e-10,rawsd,1)
  ch=z[:,cov.index("log_chla")]
  interactions=[ch*ch]
  for m in mods:interactions.append(ch*z[:,cov.index(m)])
  mat=np.column_stack([z]+interactions)
  center=mat[tr].mean(0);scale=mat[tr].std(0)
  mat=(mat-center)/np.where(scale>1e-10,scale,1)
  train=np.column_stack([np.ones(tr.sum()),mat[tr]])
  test=np.column_stack([np.ones(te.sum()),mat[te]])
  pen=np.eye(train.shape[1])*10.;pen[0,0]=0
  beta=np.linalg.solve(train.T@train+pen,train.T@y[tr])
  out[te]=test@beta
 if not np.isfinite(out).all():raise RuntimeError("Missing heldout independent predictions")
 return float(np.mean((y-out)**2)),out

# Independent full cohort OLS for all 5 and Ridge for five single interaction
# and joint, verifying exact primary results without reusing any fit helper.
results={}
for mod in mods:
 value=model_ols(d,controls,mod)
 old=main["single_moderators"][mod]
 if abs(value["beta_interaction_per1SD_x_1SD"]-old["interaction_beta_per_1SD_Chla_and_1SD_environment"])>1e-7:
  raise RuntimeError("Independent manual PHab x Chla coefficient mismatch "+mod)
 if max(abs(a-b) for a,b in zip(value["HUC8_cluster_normal95"],old["HUC8_cluster_descriptive95"]))>5e-5:
  raise RuntimeError("Independent manual clustered interval mismatch "+mod)
 mse,_=cv(d,(mod,))
 if abs(mse-main["MSE_HUC8_disjoint_5fold"]["I_"+mod])>1e-7:
  raise RuntimeError("Independent manual whole-HUC8 out-of-fold MSE mismatch "+mod)
 results[mod]={"coefficient":value,"manual_HUC8_OOS_MSE":mse,
               "Holm5_original":old["Holm_5_family_adjusted_p_posthoc"],
               "original_95gain":old["HUC8_bootstrap_2000_fixedOOF95_extra_gain_pct"]}
ref,_=cv(d,())
if abs(ref-main["MSE_HUC8_disjoint_5fold"]["Q"])>1e-7:
 raise RuntimeError("Independent quadratic base MSE mismatch")
combined,_=cv(d,mods)
if abs(combined-main["MSE_HUC8_disjoint_5fold"]["Iall"])>1e-7:
 raise RuntimeError("Independent joint MSE mismatch")

# Sensitivity contrasts are posthoc relative to already measured main signal.
subs={
 "all_original_946":np.ones(len(d),dtype=bool),
 "first_native_SITE_ID_only":~d.sort_values(["year","UID"]).duplicated("site_id").reindex(d.index).to_numpy(),
 "benthic_counted_individuals_at_least250":d.total_benthic_subsample_count.ge(250).to_numpy(),
 "source_original_XCDENMID_at_most100":d.XCDENMID.le(100).to_numpy(),
 "source_original_XGVEG_at_most1":d.XGVEG.le(1).to_numpy(),
 "only_original_2018":d.year.eq(2018).to_numpy(),
 "only_original_2019":d.year.eq(2019).to_numpy()
}
robust={}
for label,mask in subs.items():
 subset=d.loc[mask].copy()
 if len(subset)<300 or subset.huc8.nunique()<180 or subset.state.nunique()<20:
  robust[label]={"n":len(subset),"HUC8":subset.huc8.nunique(),"status":"INSUFFICIENT_SUPPORT"}
  continue
 control=pd.concat([subset[cov].astype(float),
                    pd.get_dummies(subset[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
 est=model_ols(subset,control,"XEMBED",sc)
 q,pred0=cv(subset,())
 a,pred1=cv(subset,("XEMBED",))
 est["independent_HUC8_MSE_Q"]=q
 est["independent_HUC8_MSE_with_embed_interaction"]=a
 est["independent_heldout_gain_pct"]=float(100*(q-a)/q)
 est["n_distinct_physical_SITE_ID"]=subset.site_id.nunique()
 robust[label]=est
out={"schema":"NRSA1819_EPT_CHLA_5PHAB_CONTEXT_INDEPENDENT_NUMPY_PLUS_SEVEN_POSTHOC_QA_v1",
     "original_EPA_2018_source_n":1221,"strict_actual_946":True,"HUC8":603,
     "not_independent_field_replicate":True,"not_wastewater_causal":True,
     "original_preselected_habitat_interaction_family":mods,
     "all_five_independent_replicated":results,
     "independent_base_Q_MSE":ref,"independent_five_joint_MSE":combined,
     "original_parent_run_claim_matches":True,
     "embedding_seven_posthoc_sensitivities":robust}
path=ROOT/"NRSA1819_EPT_CHLA_PHAB_INTERACTION_INDEPENDENT_NUMPY_SENSITIVITY.json"
path.write_text(json.dumps(out,indent=2,sort_keys=True))
print("EPT_CHLA_CONTEXT_INDEPENDENT_NUMPY_SUCCESS")
print(json.dumps({"all_moderators":{k:{"beta":v["coefficient"]["beta_interaction_per1SD_x_1SD"],
                                     "p":v["coefficient"]["cluster_normal_p"],
                                     "Holm":v["Holm5_original"]} for k,v in results.items()},
                  "sensitivity":robust,
                  "MSE_Q_manual":ref,"MSE_Iall_manual":combined},sort_keys=True))
