#!/usr/bin/env python3
"""POSTHOC robustness of original NRSA2018-19 riparian EPT effect modifications.
Original source-hash-pinned parent source model is run first in same runner;
this script uses its freshly verified four EPA original CSVs.
Checks EPA XGVEG >1, whether XGH+XGW can explain it, untrimmed/winsorized
versus primary prespecified 0..1 selection, with no new source year/species.
"""
import hashlib,json,math,os,pathlib
import numpy as np,pandas as pd
root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_nrsa_riparian_ept"
report=json.loads((root/"NRSA1819_EPT_TNTP_RIPARIAN_SAME_UID_HUC8_OUTCOME_REPORT_20261112.json").read_text())
m_path=pathlib.Path(os.environ["NRSA_FROZEN_MANIFEST"])
p_path=pathlib.Path(os.environ["NRSA_RIPARIAN_MANIFEST"])
assert hashlib.sha256(m_path.read_bytes()).hexdigest()=="fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91"
assert hashlib.sha256(p_path.read_bytes()).hexdigest()=="23cf9546b751d1eaedb9aeefcdb63daed1edf14d41bae9f551240252c90f2012"
m=pd.DataFrame(json.loads(m_path.read_text())["primary_berw_UID_subset"]).set_index("uid")
puids=set(json.loads(p_path.read_text())["eligible_UIDs"])
m=m.loc[sorted(puids)]
def read(n):
 d=pd.read_csv(root/(n+".csv"),dtype=str,keep_default_na=False,low_memory=False,encoding="utf-8-sig")
 d.columns=d.columns.str.upper()
 return d
B=read("benthos");S=read("sites").set_index("UID").reindex(m.index)
W=read("water").set_index("UID").reindex(m.index)
P=read("phab").set_index("UID").reindex(m.index)
B=B[B.UID.isin(puids)&B.SAMPLE_TYPE.eq("BERW")].copy()
B["n"]=pd.to_numeric(B.TOTAL300,errors="coerce").fillna(0)
B=B[B.n>0].copy()
is_ept=B.ORDER.str.strip().str.lower().isin(["ephemeroptera","plecoptera","trichoptera"])
frac=((B.n*is_ept).groupby(B.UID).sum()/B.groupby("UID").n.sum())
F=pd.DataFrame(index=m.index)
F["ept"]=frac.reindex(F.index)
F["count"]=B.groupby("UID").n.sum().reindex(F.index)
for n in ["state","eco9","huc8","site_id"]:
 F[n]=m[n]
for n in ["lat","lon","year"]:
 F[n]=pd.to_numeric(m[n],errors="coerce")
F["stream_order"]=pd.to_numeric(S.STRAH_ORD,errors="coerce")
F["logcond"]=np.log10(pd.to_numeric(W.COND_RESULT,errors="coerce"))
F["logTN"]=np.log10(pd.to_numeric(W.NTL_RESULT,errors="coerce"))
F["logTP"]=np.log10(pd.to_numeric(W.PTL_RESULT,errors="coerce")/1000)
for name,field in [("ground","XGVEG"),("canopy","XPCAN"),("overhang","XFC_OHV"),
                   ("human","W1_HALL"),("pipe","W1H_PIPE"),("xgh","XGH"),("xgw","XGW")]:
 F[name]=pd.to_numeric(P[field],errors="coerce")
F["flags_clear"]=(W.NTL_NARS_FLAG.str.strip()=="")&(W.PTL_NARS_FLAG.str.strip()=="")
XG=F.ground
over=F[(F.ground>1)&F.ground.notna()]
diff=(F.ground-F.xgh-F.xgw)
metricQA={
 "xgveg_exceeds_one_original":int((F.ground>1).sum()),
 "max_original_xgveg":float(F.ground.max()),
 "xgh_xgw_nonmissing":int((F.xgh.notna()&F.xgw.notna()).sum()),
 "xgveg_minus_xgh_plus_xgw_over1_abs_max":float(diff.loc[over.index].abs().max()) if len(over) else None,
 "xgveg_minus_xgh_plus_xgw_over1_median":float(diff.loc[over.index].median()) if len(over) else None,
 "xgh_xgw_sum_above1_n":int(((F.xgh+F.xgw)>1).sum()),
 "source_metadata_has_no_legal_range_for_xgveg":True
}
BASE=["lat","lon","year","stream_order","logcond","human"]
CATS=["state","eco9"]
REQ=BASE+["ept","logTN","logTP","ground","canopy","overhang"]
def reg(d,cols):
 x=pd.get_dummies(d[BASE+CATS],columns=CATS,drop_first=True,dtype=float)
 x=pd.concat([x,d[cols].astype(float)],axis=1)
 x.insert(0,"intercept",1.)
 return x
def fold(t):
 return int(hashlib.sha256(str(t).encode()).hexdigest()[:8],16)%5
def predict(d,cols):
 x=reg(d,cols).to_numpy(float);y=d.ept.to_numpy(float)
 folds=np.array([fold(x) for x in d.huc8])
 pp=np.zeros(len(y))
 for k in range(5):
    train=folds!=k;test=~train
    z=x.copy()
    mu=z[train,1:].mean(axis=0);sd=z[train,1:].std(axis=0)
    sd=np.where(sd<1e-12,1,sd)
    z[:,1:]=(z[:,1:]-mu)/sd
    zz=z[train];zt=z[test]
    pen=10*np.eye(zz.shape[1]);pen[0,0]=0
    b=np.linalg.solve(zz.T@zz+pen,zz.T@y[train])
    pp[test]=zt@b
 return float(np.mean((y-pp)**2))
def cluster_slope(d,cols):
 x=reg(d,cols)
 mat=x.to_numpy(float);y=d.ept.to_numpy(float)
 b=np.linalg.lstsq(mat,y,rcond=None)[0];err=y-mat@b
 bread=np.linalg.pinv(mat.T@mat);meat=np.zeros((len(b),len(b)))
 g=d.huc8.to_numpy(str)
 for cl in np.unique(g):
    s=mat[g==cl].T@err[g==cl]
    meat+=np.outer(s,s)
 G=len(np.unique(g));n=len(y);K=np.linalg.matrix_rank(mat)
 v=bread@meat@bread*G/(G-1)*(n-1)/(n-K)
 se=np.sqrt(np.maximum(v.diagonal(),0))
 out={}
 for c in ["logTP","logTN","ground","tp_ground_interaction"]:
  if c in x.columns:
   k=x.columns.get_loc(c)
   out[c]={"beta":float(b[k]),"descriptive95":[float(b[k]-1.96*se[k]),float(b[k]+1.96*se[k])]}
 return out
def evaluate(name,f):
 d=f.replace([np.inf,-np.inf],np.nan).dropna(subset=REQ).copy()
 d=d[d.canopy.between(0,1)&d.overhang.between(0,1)&(d.human>=0)]
 if len(d)<650:raise RuntimeError("Posthoc domain variant not supported")
 d["tp_ground_interaction"]=d.logTP*(d.ground-0.5)
 nutrients=["logTN","logTP"];habitat=["ground","canopy","overhang"]
 NH=nutrients+habitat;NHI=NH+["tp_ground_interaction"]
 errs={
  "N":predict(d,nutrients),
  "NH":predict(d,NH),
  "NHI":predict(d,NHI)}
 return {"n":len(d),"physical_sites":d.site_id.nunique(),
  "slope_NH":cluster_slope(d,NH),"slope_NHI":cluster_slope(d,NHI),
  "huc8_mse":errs,
  "habitat_additional_after_N_pct":100*(errs["N"]-errs["NH"])/errs["N"],
  "interaction_additional_after_NH_pct":100*(errs["NH"]-errs["NHI"])/errs["NH"]}
out={}
out["predefined_0_to1_XGVEG_drop"]=evaluate("drop",F[F.ground.between(0,1)])
out["all_original_XGVEG_including_gt1"]=evaluate("all",F)
C=F.copy();C["ground"]=C.ground.clip(upper=1)
out["XGVEG_winsorized_at1_posthoc"]=evaluate("cap",C)
q=F[F.flags_clear & F["count"].ge(250)]
out["clean_NARS_subsample250_all_raw_ground"]=evaluate("clean",q)
# Confirm main model source locked results, prevent changing it silently
s=report["HUC8_cluster_descriptive_OLS_fits"]["NHI"]["tp_ground_interaction"]["beta"]
if abs(s-out["predefined_0_to1_XGVEG_drop"]["slope_NHI"]["tp_ground_interaction"]["beta"])>1e-9:
 raise RuntimeError("Independent reimplementation failed to reproduce original interaction")
audit={"schema":"NRSA1819_EPT_TP_RIPARIAN_GROUND_SOURCE_VALUE_RANGE_POSTHOC_QA_v1",
       "source_hashes":report["source_sha256"],
       "original_NH_specification_reproduced":True,
       "numeric_ground_metric_QA":metricQA,
       "variants":out,
       "not_causal_modulation_or_pollination":True,
       "source_flora_and_butterfly_not_joined":True}
path=root/"NRSA1819_RIPARIAN_XGVEG_RANGE_INDEPENDENT_SENSITIVITY_20261112.json"
path.write_text(json.dumps(audit,indent=2,sort_keys=True))
print("NRSA_RIPARIAN_GROUND_RANGE_AUDIT_START")
print(json.dumps(audit,sort_keys=True))
print("NRSA_RIPARIAN_GROUND_RANGE_AUDIT_END")
