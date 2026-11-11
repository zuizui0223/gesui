#!/usr/bin/env python3
"""Frozen same-UID aquatic EPT × stream TN/TP × riparian vegetation bridge.

Exploratory and source-locked before reading the new numeric physical-habitat
responses. No butterflies, aquatic insect emergence, or wastewater attribution.
"""
import hashlib,json,math,os,pathlib,time,urllib.request
import numpy as np
import pandas as pd

PIN={
 "benthos":("https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_benthic_macroinvertebrate_count_-_data.csv",
   "4e9b3a0798b25f67c7a96d4800b94f86544bbd623abaf119358e067b158324a7"),
 "sites":("https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_SiteInfo.csv",
   "1b5b2119c881290191edeb47f1d80a7711f10440e9851574cfa86fd57be770e3"),
 "water":("https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_water_chemistry_chla_-_data.csv",
   "cb9a3349ae3708834eab533e71423d74c1778bf76aba0737f0ae1bcc8f50baa1"),
 "phab":("https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_physical_habitat_larger_set_of_metrics_-_data.csv",
   "69b3b56962c04492423bd7eaf4d074c6ae8d253924d683ba12eb25fd857ce063")
}
UID_SHA="fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91"
PHAB_UID_SHA="23cf9546b751d1eaedb9aeefcdb63daed1edf14d41bae9f551240252c90f2012"
root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_nrsa_riparian_ept"
root.mkdir(parents=True,exist_ok=True)
def sha(x):return hashlib.sha256(x).hexdigest()
uid_path=pathlib.Path(os.environ["NRSA_FROZEN_MANIFEST"])
phabuid_path=pathlib.Path(os.environ["NRSA_RIPARIAN_MANIFEST"])
if sha(uid_path.read_bytes())!=UID_SHA or sha(phabuid_path.read_bytes())!=PHAB_UID_SHA:
    raise RuntimeError("Original pre-outcome source cohort SHA mismatch")
M=pd.DataFrame(json.loads(uid_path.read_text())["primary_berw_UID_subset"]).set_index("uid")
phab_source=json.loads(phabuid_path.read_text())
pinned_uids=set(phab_source["eligible_UIDs"])
if len(M)!=1221 or not M.index.is_unique or len(pinned_uids)!=1220 or not pinned_uids.issubset(set(M.index)):
    raise RuntimeError("Original pre-numeric bridge samples changed")
for name,(url,expected) in PIN.items():
    path=root/(name+".csv")
    for attempt in range(4):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"gesui-source-pinned-riparian-20261112/1","Accept":"text/csv,*/*"})
            with urllib.request.urlopen(req,timeout=80) as r,path.open("wb") as f:
                while chunk:=r.read(1048576):f.write(chunk)
            if sha(path.read_bytes())!=expected:raise RuntimeError("EPA upstream source changed: "+name)
            break
        except Exception:
            path.unlink(missing_ok=True)
            if attempt==3:raise
            time.sleep(min(2**attempt+1,10))
    print("Verified original "+name+" SHA and byte length "+str(path.stat().st_size),flush=True)
def dataframe(name):
    t=pd.read_csv(root/(name+".csv"),dtype=str,keep_default_na=False,low_memory=False,encoding="utf-8-sig")
    t.columns=t.columns.str.upper()
    return t
B=dataframe("benthos");S=dataframe("sites");W=dataframe("water");P=dataframe("phab")
S=S[S.UID.isin(pinned_uids)].copy();W=W[W.UID.isin(pinned_uids)].copy();P=P[P.UID.isin(pinned_uids)].copy()
for name,t in [("sites",S),("water",W),("phab",P)]:
    if len(t)!=1220 or not t.UID.is_unique:raise RuntimeError(name+" not a unique exactly pinned UID join")
    t.set_index("UID",inplace=True)
M=M.reindex(sorted(pinned_uids))
S=S.reindex(M.index);W=W.reindex(M.index);P=P.reindex(M.index)
if not (S.SITE_ID.eq(M.site_id).all() and W.SITE_ID.eq(M.site_id).all() and P.SITE_ID.eq(M.site_id).all()):
    raise RuntimeError("Same-visit native site identity does not match across original files")
B=B[B.UID.isin(pinned_uids)&B.SAMPLE_TYPE.eq("BERW")].copy()
B["N"]=pd.to_numeric(B.TOTAL300,errors="coerce")
if (B.N.dropna()<0).any():raise RuntimeError("Negative original benthic count")
B=B[B.N>0].copy()
B["is_ept"]=B.ORDER.str.strip().str.lower().isin(["ephemeroptera","plecoptera","trichoptera"])
group=B.groupby("UID")
n=group.N.sum()
ept=(B.N*B.is_ept).groupby(B.UID).sum()
d=pd.DataFrame(index=M.index)
d["total300"]=n.reindex(d.index)
d["ept_fraction"]=(ept/n).reindex(d.index)
if d.total300.isna().any() or not d.ept_fraction.between(0,1).all():
    raise RuntimeError("Invalid source BERW EPT fraction")
d["state"]=M.state
d["eco9"]=M.eco9
d["huc8"]=M.huc8
d["site_id"]=M.site_id
for k in ["lat","lon","year"]:
    d[k]=pd.to_numeric(M[k],errors="coerce")
d["stream_order"]=pd.to_numeric(S.STRAH_ORD,errors="coerce")
d["logcond"]=np.log10(pd.to_numeric(W.COND_RESULT,errors="coerce").where(lambda x:x>0))
for k in ["NTL_RESULT","PTL_RESULT"]:
    if (pd.to_numeric(W[k],errors="coerce").dropna()<=0).any():
        raise RuntimeError("Non-positive measured source qualified EPA nutrient")
d["logTN"]=np.log10(pd.to_numeric(W.NTL_RESULT,errors="coerce"))
d["logTP"]=np.log10(pd.to_numeric(W.PTL_RESULT,errors="coerce")/1000)
if set(W.NTL_UNITS)!={"mg N/L"} or set(W.PTL_UNITS)!={"ug/L"}:
    raise RuntimeError("Source nutrient units changed")
for src,new in [("XGVEG","ground_cover"),("XPCAN","canopy_presence"),("XFC_OHV","overhang_cover"),
                 ("W1_HALL","human_disturbance"),("W1H_PIPE","pipe_index")]:
    d[new]=pd.to_numeric(P[src],errors="coerce")
d["tp_ground_interaction"]=d.logTP*(d.ground_cover-0.5)
d["clean_NP_flags"]=(W.NTL_NARS_FLAG.str.strip()=="")&(W.PTL_NARS_FLAG.str.strip()=="")
pre_n=len(d);quality={}
for v in ["ground_cover","canopy_presence","overhang_cover"]:
    z=d[v]
    quality[v]={"nonmissing":int(z.notna().sum()),"out_of_0_1":int((z.notna()&~z.between(0,1)).sum()),
                "min":float(z.min()),"median":float(z.median()),"max":float(z.max())}
    d.loc[~z.between(0,1),v]=np.nan
for v in ["human_disturbance","pipe_index"]:
    z=d[v]
    quality[v]={"nonmissing":int(z.notna().sum()),"negative":int((z.notna()&(z<0)).sum()),
                "min":float(z.min()),"median":float(z.median()),"max":float(z.max())}
    d.loc[z<0,v]=np.nan
required=["ept_fraction","lat","lon","year","stream_order","logcond","logTN","logTP",
          "ground_cover","canopy_presence","overhang_cover","human_disturbance"]
d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=required).copy()
d["tp_ground_interaction"]=d.logTP*(d.ground_cover-0.5)
if len(d)<650 or d.state.nunique()<25 or d.eco9.nunique()<6:
    raise RuntimeError("Frozen habitat-nutrient complete cases fail pre-outcome support")
def fold(huc):
    return int(hashlib.sha256(str(huc).encode()).hexdigest()[:8],16)%5
folds=np.array([fold(x) for x in d.huc8])
if min(np.bincount(folds,minlength=5))<50:
    raise RuntimeError("HUC8 source CV fold support insufficient")
base=["lat","lon","year","stream_order","logcond","human_disturbance"]
cat=pd.get_dummies(d[["state","eco9"]],drop_first=True,dtype=float)
X0=pd.concat([d[base].astype(float),cat],axis=1)
defs={
 "B":[],
 "N":["logTN","logTP"],
 "H":["ground_cover","canopy_presence","overhang_cover"],
 "NH":["logTN","logTP","ground_cover","canopy_presence","overhang_cover"],
 "NHI":["logTN","logTP","ground_cover","canopy_presence","overhang_cover","tp_ground_interaction"]
}
def matrix(cols):
    X=pd.concat([X0,d[cols].astype(float)],axis=1)
    X.insert(0,"intercept",1.)
    return X
X={name:matrix(cols) for name,cols in defs.items()}
Y=d.ept_fraction.to_numpy(float)
def cvpred(XX):
    features=XX.to_numpy(float)
    pred=np.empty(len(Y))
    for f in range(5):
        train=folds!=f;test=~train
        Z=features.copy()
        mu=Z[train,1:].mean(axis=0);sd=Z[train,1:].std(axis=0)
        sd=np.where(sd<1e-12,1.,sd)
        Z[:,1:]=(Z[:,1:]-mu)/sd
        A=Z[train];C=Z[test]
        penalty=10*np.eye(A.shape[1]);penalty[0,0]=0.
        beta=np.linalg.solve(A.T@A+penalty,A.T@Y[train])
        pred[test]=C@beta
    return pred
preds={name:cvpred(XX) for name,XX in X.items()}
e2={k:(Y-v)**2 for k,v in preds.items()}
mse={k:float(np.mean(v)) for k,v in e2.items()}
comparisons=[("B","N"),("B","H"),("N","NH"),("H","NH"),("NH","NHI")]
clusters=d.huc8.astype(str).to_numpy()
unq,inverse=np.unique(clusters,return_inverse=True)
m=len(unq)
rng=np.random.default_rng(20261112)
samples=rng.integers(0,m,size=(2000,m))
comparisons_results={}
for baseline,alt in comparisons:
    gain=100*(mse[baseline]-mse[alt])/mse[baseline]
    a=np.bincount(inverse,weights=e2[baseline],minlength=m)
    b=np.bincount(inverse,weights=e2[alt],minlength=m)
    bs_a=np.sum(a[samples],axis=1)
    bs_b=np.sum(b[samples],axis=1)
    dist=100*(bs_a-bs_b)/bs_a
    ci=np.percentile(dist,[2.5,97.5])
    comparisons_results[baseline+"→"+alt]={
      "site_equal_heldout_MSE_gain_pct":float(gain),
      "HUC8_cluster_bootstrap_descriptive95_gain_pct":ci.tolist(),
      "fraction_2000_resamples_positive":float(np.mean(dist>0))
    }
def ols_cluster(A,yy,cc):
    z=A.to_numpy(float)
    beta=np.linalg.lstsq(z,yy,rcond=None)[0]
    resid=yy-z@beta
    bread=np.linalg.pinv(z.T@z)
    meat=np.zeros((z.shape[1],z.shape[1]))
    for g in np.unique(cc):
        ind=cc==g;v=z[ind].T@resid[ind]
        meat+=np.outer(v,v)
    G=len(np.unique(cc));n=len(yy);rank=np.linalg.matrix_rank(z)
    factor=G/(G-1)*(n-1)/(n-rank)
    V=(bread@meat@bread)*factor
    se=np.sqrt(np.maximum(np.diag(V),0.))
    values={}
    for k in ["logTN","logTP","ground_cover","canopy_presence","overhang_cover","tp_ground_interaction"]:
        if k in A.columns:
            j=A.columns.get_loc(k);values[k]={"beta":float(beta[j]),
             "HUC8_cluster_descriptive95":[float(beta[j]-1.96*se[j]),float(beta[j]+1.96*se[j])]}
    return values
fits={k:ols_cluster(X[k],Y,clusters) for k in ["N","NH","NHI"]}
sens={}
for key,sel in {
    "clean_NARS":d.clean_NP_flags,
    "at_least_250_subsample":d.total300>=250,
    "first_original_site_visit_only":~d.sort_values(["site_id","year"]).duplicated("site_id").reindex(d.index),
    "clean_AND_250":d.clean_NP_flags&(d.total300>=250),
}.items():
    t=d.loc[sel].copy()
    if len(t)<450 or t.huc8.nunique()<250:
        sens[key]={"n":len(t),"gate":"insufficient"};continue
    dat=X["NH"].loc[t.index]
    sens[key]={"n":len(t),"sites":t.site_id.nunique(),
        "NH_slopes":ols_cluster(dat,t.ept_fraction.to_numpy(float),t.huc8.astype(str).to_numpy())}
out={
 "schema":"NRSA1819_SAME_UID_RIPARIAN_WATER_EPT_JOINT_OBSERVATIONAL_v1",
 "source_sha256":{n:expected for n,(_,expected) in PIN.items()},
 "source_manifest_SHA256":UID_SHA,"riparian_pre_numeric_UID_SHA256":PHAB_UID_SHA,
 "cohort":{"original_UIDs":1221,"source_fixed_riparian_UIDs":pre_n,
   "complete_case_visits":len(d),"distinct_physical_sites":d.site_id.nunique(),
   "states":d.state.nunique(),"eco9":d.eco9.nunique(),
   "distinct_HUC8":d.huc8.nunique(),"fold_counts":np.bincount(folds,minlength=5).tolist(),
   "clean_NP_NARS_flag_visits":int(d.clean_NP_flags.sum()),
   "at_least_250_benthic_individuals":int((d.total300>=250).sum())},
 "habitat_numeric_source_QA":quality,
 "mean_EPT_individual_fraction":float(Y.mean()),
 "covariance_qualifier":"same-time field co-organization, not water-to-land trophic flux, not WWTP attribution",
 "model_covariates":{"B":base+list(cat.columns),"nested_additions":defs},
 "HUC8_disjoint_5fold_ridge10_MSE":mse,
 "nested_out_of_region_predictive_contributions":comparisons_results,
 "HUC8_cluster_descriptive_OLS_fits":fits,
 "posthoc_fixed_cohort_sensitivities":sens,
 "limitations":[
   "Cross-sectional same-place associations cannot identify causal buffering or mediation.",
   "XGVEG/XPCAN/XFC_OHV are structural riparian covers, not species richness, butterfly hosts or flowering.",
   "Original W1H_PIPE is only proximity-weighted human disturbance index for pipes, not measured municipal wastewater.",
   "Europe EEA 2018/2022 WWTP-source EQR is different geography and outcome, not joined to these EPA river visits.",
   "No adult aquatic-insect emergence, predators' diets, butterfly visitation, pollination or ecosystem service measured."
 ]}
outpath=root/"NRSA1819_EPT_TNTP_RIPARIAN_SAME_UID_HUC8_OUTCOME_REPORT_20261112.json"
outpath.write_text(json.dumps(out,indent=2,sort_keys=True))
print("NRSA_RIPARIAN_EPT_BRIDGE_RESULT_START",flush=True)
print(json.dumps({
 "cohort":out["cohort"],"habitat_QA":quality,
 "original_MSE":mse,"comparisons":comparisons_results,
 "slopes":fits,"sensitivities":sens},sort_keys=True),flush=True)
print("NRSA_RIPARIAN_EPT_BRIDGE_RESULT_END",flush=True)
