#!/usr/bin/env python3
"""SOURCE-PINNED EPT–Chla embeddedness composition vs Chironomidae counted
share, EPT native taxon IDs and non-EPT native taxon IDs at same 946 EPA visits.
Posthoc after EPT moderation; prespecified 3 endpoints, all results retained.
"""
from __future__ import annotations
import collections,hashlib,json,math,os,pathlib,time,urllib.request
import numpy as np,pandas as pd
from scipy.stats import norm

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui1819_triad_numeric"
SRC=ROOT/"AUDIT_ONLY_2018_19_EPT_Chla_PHAB_SOURCE_PANEL.csv"
PARENT=ROOT/"NRSA1819_EPT_CHLA_PHAB_INTERACTION_SOURCE_SIGNED_EXPLORATORY.json"
BIO_URL="https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_benthic_macroinvertebrate_count_-_data.csv"
BIO_SHA="4e9b3a0798b25f67c7a96d4800b94f86544bbd623abaf119358e067b158324a7"
a=json.loads(PARENT.read_text())
if a["cohort"]!=946 or a["HUC8"]!=603 or "XEMBED" not in a["single_moderators"]:
 raise RuntimeError("Prior source-locked embeddedness moderation sample/result changed")
frame=pd.read_csv(SRC,dtype={"UID":str,"site_id":str,"huc8":str,"eco9":str,"state":str})
if len(frame)!=1221 or frame.UID.duplicated().any():raise RuntimeError("Original EPA source native UIDs changed")
p=ROOT/"EPA1819_BERW_original.csv"
for attempt in range(4):
 try:
  req=urllib.request.Request(BIO_URL,headers={"User-Agent":"gesui-EPT-CHIRONOMID-source-frozen/20261114"})
  with urllib.request.urlopen(req,timeout=120) as r,p.open("wb") as f:
   while b:=r.read(1<<20):f.write(b)
  if hashlib.sha256(p.read_bytes()).hexdigest()!=BIO_SHA:raise RuntimeError("Authentic EPA original benthic source SHA mismatch")
  break
 except Exception:
  p.unlink(missing_ok=True)
  if attempt==3:raise
  time.sleep(2**attempt)
B=pd.read_csv(p,dtype=str,keep_default_na=False,low_memory=False)
p.unlink()
B=B[B.UID.isin(set(frame.UID))&B.SAMPLE_TYPE.eq("BERW")].copy()
B["count"]=pd.to_numeric(B.TOTAL300,errors="coerce")
B=B[B["count"].gt(0)&np.isfinite(B["count"])].copy()
if B.TAXA_ID.eq("").any():raise RuntimeError("Original positive benthos record lacks taxon-ID")
B["is_EPT"]=B.ORDER.str.strip().str.lower().isin(("ephemeroptera","plecoptera","trichoptera"))
B["is_chironomid"]=B.FAMILY.str.strip().str.lower().eq("chironomidae")
if any(B.is_EPT&B.is_chironomid):raise RuntimeError("Source Chironomidae and EPT mutually exclusive taxonomy failed")
G=B.groupby("UID")
tot=G["count"].sum()
ept=(B["count"]*B.is_EPT.astype(float)).groupby(B.UID).sum()
chir=(B["count"]*B.is_chironomid.astype(float)).groupby(B.UID).sum()
rich_ept=B.loc[B.is_EPT].groupby("UID").TAXA_ID.nunique()
rich_nEPT=B.loc[~B.is_EPT].groupby("UID").TAXA_ID.nunique()
frame=frame.set_index("UID")
bio=pd.DataFrame(index=frame.index)
bio["original_total_count"]=tot.reindex(bio.index)
bio["source_EPT_share"]=ept.reindex(bio.index)/bio.original_total_count
bio["chironomid_fraction"]=(chir.reindex(bio.index,fill_value=0)/bio.original_total_count)
bio["EPT_taxonID_richness"]=rich_ept.reindex(bio.index,fill_value=0)
bio["nonEPT_taxonID_richness"]=rich_nEPT.reindex(bio.index,fill_value=0)
data=frame.join(bio)
if not np.allclose(data.ept_fraction.to_numpy(float),data.source_EPT_share.to_numpy(float),atol=1e-10,rtol=0):
 raise RuntimeError("Original independent source EPT share differs from PR77 original field data")
if not np.allclose(data.total_benthic_subsample_count.to_numpy(float),data.original_total_count.to_numpy(float),atol=1e-7,rtol=0):
 raise RuntimeError("Source original counted effort differs from PR77 original data")
flags=(data.chla_flag.fillna("").eq("")&data.tn_mgL_rawflag.fillna("").eq("")&data.tp_mgL_rawflag.fillna("").eq(""))
cols=["lat","lon","year","stream_order","season_sin","season_cos","log_cond",
      "log_tn","log_tp","XPCAN","XGVEG","XFC_OHV","W1_HALL","log_width",
      "XCDENMID","XEMBED","log_chla","ept_fraction","chironomid_fraction",
      "EPT_taxonID_richness","nonEPT_taxonID_richness"]
D=data.loc[flags].dropna(subset=cols).copy()
D=D[(D.huc8!="")&(D.state!="")&(D.eco9!="")].sort_index()
if len(D)!=946 or D.huc8.nunique()!=603 or D.state.nunique()!=47:
 raise RuntimeError("Source 946 original complete same native UID cohort changed")
D["log1p_EPT_taxonID"]=np.log1p(D.EPT_taxonID_richness)
D["log1p_nonEPT_taxonID"]=np.log1p(D.nonEPT_taxonID_richness)
D["log_count"]=np.log(D.original_total_count)
fixed=["lat","lon","year","stream_order","season_sin","season_cos","log_cond",
       "log_tn","log_tp","XPCAN","XGVEG","XFC_OHV","W1_HALL",
       "log_width","XCDENMID","XEMBED","log_chla"]
dummies=pd.get_dummies(D[["state","eco9"]],drop_first=True,dtype=float)
grp=D.huc8.astype(str).to_numpy()
fold=np.array([int(hashlib.sha256(str(z).encode()).hexdigest()[:8],16)%5 for z in grp])
if min(np.bincount(fold,minlength=5))<100:raise RuntimeError("5fold whole HUC8 eligible site counts collapse")
c=(D.log_chla.to_numpy(float)-D.log_chla.mean())/D.log_chla.std(ddof=0)
e=(D.XEMBED.to_numpy(float)-D.XEMBED.mean())/D.XEMBED.std(ddof=0)
c2=c*c;interaction=c*e

def original_model(outcome,richness=False):
 Y=D[outcome].to_numpy(float)
 base=D[fixed].astype(float).copy()
 if richness:base["log_count"]=D["log_count"].to_numpy(float)
 base=pd.concat([base,dummies],axis=1)
 origsize=base.shape[1]
 # All original nuisance covars and squaring of Chla controlled in both Q/I.
 X=base.copy()
 X["chla_z_squared"]=c2
 XI=X.assign(chla_z_x_XEMBED=interaction)
 def cv(Z,want_interaction):
  raw=Z.to_numpy(float)
  preds=np.empty(len(Y))
  for k in range(5):
   tr=fold!=k;te=~tr
   orig=raw[:,:origsize]
   mu=orig[tr].mean(0);sd=orig[tr].std(0)
   origz=(orig-mu)/np.where(sd>1e-10,sd,1)
   cc=origz[:,fixed.index("log_chla")]
   emb=origz[:,fixed.index("XEMBED")]
   more=[cc*cc]
   if want_interaction:more.append(cc*emb)
   X0=np.column_stack([origz]+more)
   zz=(X0-X0[tr].mean(0))/np.where(X0[tr].std(0)>1e-10,X0[tr].std(0),1)
   A=np.column_stack([np.ones(tr.sum()),zz[tr]])
   Bte=np.column_stack([np.ones(te.sum()),zz[te]])
   pen=np.eye(A.shape[1])*10;pen[0,0]=0
   beta=np.linalg.solve(A.T@A+pen,A.T@Y[tr])
   preds[te]=Bte@beta
  return preds
 # Cluster OLS descriptive association, using full sample source scaling.
 xx=np.column_stack([np.ones(len(Y)),XI.to_numpy(float)])
 b=np.linalg.lstsq(xx,Y,rcond=None)[0]
 residual=Y-xx@b
 meat=np.zeros((xx.shape[1],xx.shape[1]))
 for group in pd.unique(grp):
  ix=grp==group;u=xx[ix].T@residual[ix];meat+=np.outer(u,u)
 bread=np.linalg.pinv(xx.T@xx)
 N=len(Y);G=len(pd.unique(grp));K=np.linalg.matrix_rank(xx)
 var=bread@meat@bread*G/(G-1)*(N-1)/(N-K)
 ix=XI.columns.get_loc("chla_z_x_XEMBED")+1
 bet=float(b[ix]);se=float(np.sqrt(max(0,var[ix,ix])))
 p=float(2*norm.sf(abs(bet/se))) if se>0 else 1.
 pred_Q=cv(X,False);pred_I=cv(XI,True)
 mseQ=float(np.mean((Y-pred_Q)**2));mseI=float(np.mean((Y-pred_I)**2))
 Nblock,inv=np.unique(grp,return_inverse=True)
 err0=np.bincount(inv,weights=(Y-pred_Q)**2,minlength=len(Nblock))
 err1=np.bincount(inv,weights=(Y-pred_I)**2,minlength=len(Nblock))
 rng=np.random.default_rng(20261114)
 draws=rng.integers(0,len(Nblock),size=(2000,len(Nblock)))
 den=err0[draws].sum(1);num=err1[draws].sum(1)
 bs=100*(den-num)/den
 return {"response":outcome,"n":len(Y),
  "source_response_mean":float(Y.mean()),
  "beta_Chla_z_x_XEMBED":bet,
  "HUC8_cluster_descriptive95":[float(bet-1.959963984540054*se),float(bet+1.959963984540054*se)],
  "cluster_p_exploratory_raw":p,
  "geographic_HUC8_MSE_Q":mseQ,"geographic_HUC8_MSE_embed_interaction":mseI,
  "geographic_MSE_gain_pct":float(100*(mseQ-mseI)/mseQ),
  "HUC8_fixedOOF_2000_bootstrap95_gain_pct":list(map(float,np.quantile(bs,[.025,.5,.975])))}
RES=[
 original_model("chironomid_fraction",False),
 original_model("log1p_EPT_taxonID",True),
 original_model("log1p_nonEPT_taxonID",True)
]
names=sorted(range(len(RES)),key=lambda i:RES[i]["cluster_p_exploratory_raw"])
adjusted=0.
for i,j in enumerate(names):
 adjusted=max(adjusted,min(1.,(len(names)-i)*RES[j]["cluster_p_exploratory_raw"]))
 RES[j]["Holm_3_exploratory_adjusted_p"]=float(adjusted)
out={
 "schema":"gesui_NR​​SA1819_TRUE_EPT_Chla_embeddedness_original_Chironomid_taxonID_decomposition_v1",
 "source_benthos_original_SHA256":BIO_SHA,
 "original_946_EPT_fraction_independently_reconstructed":True,
 "original_n_946_HUC8_603_from_PR77_78":True,
 "posthoc_after_embed_interaction_seen":True,
 "no_sewage_NP_source_nor_biomass_drymass_nor_ecosystem_services":True,
 "Chironomidae_vs_EPT_is_compositional_relative_subsample":True,
 "taxonID_richness_not_species_richness":True,
 "original_taxon_source_QA":{
  "original_total_BERW_taxon_rows":len(B),"unique_taxon_UIDs":B.UID.nunique(),
  "native_positive_duplicate_UID_taxonID_rows":int(B.duplicated(["UID","TAXA_ID"]).sum()),
  "same_946_nonEPT_source_zero_count":int((D.nonEPT_taxonID_richness==0).sum()),
  "same_946_chironomid_source_zero_share":int((D.chironomid_fraction==0).sum())},
 "results":RES}
path=ROOT/"NRSA1819_EPT_CHLA_EMBEDDEDNESS_TAXON_DECOMPOSITION_ORIGINAL.json"
path.write_text(json.dumps(out,sort_keys=True,indent=2))
print("EPT_CHLA_EMBED_TAXON_DECOMPOSITION_RESULT")
print(json.dumps({"original_QA":out["original_taxon_source_QA"],"observed":RES},sort_keys=True))
