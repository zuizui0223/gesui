#!/usr/bin/env python3
"""Exploratory NRSA 2018-19 FFG after 2026-11-11 original EPT results.

Fractionally allocate each count among comma-separated EPA FFGs to avoid
the invalid single-code-only metric. Same 1221 pre-outcome native UID cohort,
official input SHAs, original HUC8-disjoint 5fold and controls. Does NOT
measure emergence, nutrient subsidies or ecosystem services.
Must run after original source-locked benthic model, same GH runner temp.
"""
import collections,csv,hashlib,json,math,os,pathlib
import numpy as np,pandas as pd
root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"nrsa1819_functional"
manifest_path=pathlib.Path(os.environ["NRSA_FROZEN_MANIFEST"])
SHA="fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91"
if hashlib.sha256(manifest_path.read_bytes()).hexdigest()!=SHA:
    raise RuntimeError("Original native-UID frozen source hash mismatch")
source=json.loads(manifest_path.read_text())["primary_berw_UID_subset"]
m=pd.DataFrame(source).set_index("uid")
if len(m)!=1221 or not m.index.is_unique:raise RuntimeError("Original EPA fixed UID sample changed")
PRIMARY=["CG","CF","SH","SC","PR"]  # gatherer, filterer, shredder, scraper, predator
mass=collections.defaultdict(lambda:collections.Counter())
count=collections.Counter()
unknown_ffg=collections.Counter()
multiple_count=0
with (root/"benthos.csv").open(encoding="utf-8-sig",newline="") as f:
  for r in csv.DictReader(f):
    uid=r["UID"]
    if uid not in m.index or r["SAMPLE_TYPE"]!="BERW":continue
    try:n=float(r.get("TOTAL300") or "")
    except ValueError:continue
    if not math.isfinite(n) or n<=0:continue
    count[uid]+=n
    labels=list(dict.fromkeys(c.strip().upper() for c in (r.get("FFG") or "").split(",") if c.strip()))
    if len(labels)>1:multiple_count+=n
    known=[x for x in labels if x in PRIMARY]
    others=[x for x in labels if x not in PRIMARY]
    for t in others:unknown_ffg[t]+=n
    if known:
      # Do not give multi-function taxa more total mass than a single code taxon.
      for c in known:mass[uid][c]+=n*len(known)/len(labels)/len(known)
W=pd.read_csv(root/"water.csv",dtype=str,keep_default_na=False,low_memory=False).set_index("UID").reindex(m.index)
S=pd.read_csv(root/"sites.csv",dtype=str,keep_default_na=False,low_memory=False).set_index("UID").reindex(m.index)
f=pd.DataFrame(index=m.index)
f["n"]=pd.Series(count).reindex(m.index).fillna(0).astype(float)
if (f.n<=0).any():raise RuntimeError("Frozen UID lacks positive original sample count")
for c in PRIMARY:
    f["fraction_"+c]=pd.Series({u:mass[u][c] for u in count}).reindex(m.index).fillna(0)/f.n
f["assigned_fraction"]=f[["fraction_"+c for c in PRIMARY]].sum(axis=1)
if (f.assigned_fraction>1.000001).any():raise RuntimeError("Fractions exceed unity after composite FFG allocation")
ps=f[["fraction_"+c for c in PRIMARY]].div(f.assigned_fraction.replace(0,np.nan),axis=0)
f["ffg_evenness"]=-np.nansum(np.where(ps.to_numpy()>0,ps.to_numpy()*np.log(np.maximum(ps.to_numpy(),1e-15)),0),axis=1)/np.log(5)
f["year"]=pd.to_numeric(m.year,errors="coerce")
f["lat"]=pd.to_numeric(m.lat,errors="coerce");f["lon"]=pd.to_numeric(m.lon,errors="coerce")
f["eco9"]=m.eco9;f["state"]=m.state;f["huc8"]=m.huc8;f["site"]=m.site_id
f["stream_order"]=pd.to_numeric(S.STRAH_ORD,errors="coerce")
tn=pd.to_numeric(W.NTL_RESULT,errors="coerce")
tp=pd.to_numeric(W.PTL_RESULT,errors="coerce")
co=pd.to_numeric(W.COND_RESULT,errors="coerce")
f["logTN"]=np.log10(tn.where(tn>0))
f["logTP"]=np.log10((tp/1000).where(tp>0))
f["logcond"]=np.log10(co.where(co>0))
f["flag_clean"]=(W.NTL_NARS_FLAG.str.strip()=="")&(W.PTL_NARS_FLAG.str.strip()=="")
covars=["year","lat","lon","stream_order","logcond","logTN","logTP"]
d=f.replace([np.inf,-np.inf],np.nan).dropna(subset=covars+["ffg_evenness"]).copy()
if len(d)!=1220 or d.state.nunique()!=48 or d.eco9.nunique()!=9:
    raise RuntimeError("Same-site original model cohort unexpectedly changed")
cats=pd.get_dummies(d[["state","eco9"]],drop_first=True,dtype=float)
x0=np.column_stack([np.ones(len(d)),d[["lat","lon","year","stream_order","logcond"]].to_numpy(float),cats.to_numpy(float)])
x1=np.column_stack([x0,d[["logTN","logTP"]].to_numpy(float)])
labels=d.huc8.astype(str).to_numpy()
def fold(huc):
    return int(hashlib.sha256(huc.encode()).hexdigest()[:8],16)%5
splits=np.array([fold(t) for t in labels])
def fit_oos(y,x):
    pred=np.full(len(y),np.nan)
    for k in range(5):
        train=splits!=k;test=~train
        xx=x.copy()
        mu=xx[train,1:].mean(axis=0);sd=xx[train,1:].std(axis=0)
        sd[sd<1e-12]=1.
        xx[:,1:]=(xx[:,1:]-mu)/sd
        X=xx[train];Xt=xx[test]
        pen=np.eye(X.shape[1])*10.;pen[0,0]=0.
        beta=np.linalg.solve(X.T@X+pen,X.T@y[train])
        pred[test]=Xt@beta
    return float(np.mean((y-pred)**2))
def clustered_slopes(y,dat,groups):
    cx=pd.get_dummies(dat[["state","eco9"]],drop_first=True,dtype=float)
    X=np.column_stack([np.ones(len(dat)),dat[["lat","lon","year","stream_order","logcond","logTN","logTP"]].to_numpy(float),cx.to_numpy(float)])
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    u=y-X@beta;bread=np.linalg.pinv(X.T@X)
    meat=np.zeros((X.shape[1],X.shape[1]))
    for a in sorted(set(groups)):
        ind=groups==a
        v=X[ind].T@u[ind]
        meat+=np.outer(v,v)
    G=len(set(groups));N=len(y);K=np.linalg.matrix_rank(X)
    variance=bread@meat@bread*(G/(G-1)*(N-1)/(N-K))
    se=np.sqrt(np.maximum(variance.diagonal(),0))
    return {name:{"beta":float(beta[i]),"huc8_descriptive95":[float(beta[i]-1.96*se[i]),float(beta[i]+1.96*se[i])]} for i,name in [(6,"TN"),(7,"TP")]}
out={}
for key in ["fraction_"+c for c in PRIMARY]+["ffg_evenness"]:
    y=d[key].to_numpy(float)
    base=fit_oos(y,x0);full=fit_oos(y,x1)
    o={"sample_n":len(d),"mean":float(y.mean()),
       "5fold_huc8_MSE_baseline":base,"5fold_huc8_MSE_joint_nutrients":full,
       "5fold_huc8_MSE_gain_pct":100*(base-full)/base if base>0 else None,
       "joint_slope_per_10fold_nutrient":clustered_slopes(y,d,labels)}
    dd=d[(d.flag_clean)&(d.n>=250)]
    o["clean_flags_and_count250_posthoc"]={
        "n":len(dd),"slope":clustered_slopes(dd[key].to_numpy(float),dd,dd.huc8.astype(str).to_numpy())}
    out[key]=o
audit={
 "schema":"EPA_NRSA1819_FFG_MIXED_CODE_FRACTIONAL_ALLOCATION_POSTHOC_20261112_v1",
 "original_frozen_epa_UID_manifest_sha256":SHA,
 "previous_known_EPT_sign_not_used_for_source_selection":True,
 "multi_FFG_assignment_rule":"fractional 1/k among known FFG tokens, with unknown tokens receiving their allocated fraction of mass",
 "primary_codes":PRIMARY,
 "raw_positive_benthic_count_mass":sum(count.values()),
 "raw_multi_ffg_count_mass":multiple_count,
 "fraction_multiple_assignments":multiple_count/sum(count.values()),
 "unknown_label_mass_before_split":dict(unknown_ffg),
 "overall_mass_coverage_after_fractional_allocation":sum(sum(mass[u].values()) for u in count)/sum(count.values()),
 "n_valid_cases":len(d),"original_38_state_not_assumed":True,"states":d.state.nunique(),"eco9":d.eco9.nunique(),
 "fractional_ffg_total_by_site_max":float(f.assigned_fraction.max()),
 "no_direct_services_or_wastewater_exposure_measured":True,
 "guilds":out
}
path=root/"NRSA1819_FFG_FRACTIONAL_MULTILABEL_NP_POSTHOC_20261112.json"
path.write_text(json.dumps(audit,indent=2,sort_keys=True,default=float))
print("NRSA_FFG_FRACTIONAL_SOURCE_RESULT_START")
print(json.dumps({
 "n":len(d),"states":d.state.nunique(),
 "fraction_multi_ffg":audit["fraction_multiple_assignments"],
 "mass_coverage":audit["overall_mass_coverage_after_fractional_allocation"],
 "outputs":{k:{"mean":v["mean"],"heldout_gain_pct":v["5fold_huc8_MSE_gain_pct"],
       "slopes":v["joint_slope_per_10fold_nutrient"],
       "strict":v["clean_flags_and_count250_posthoc"]} for k,v in out.items()}
},sort_keys=True))
print("NRSA_FFG_FRACTIONAL_SOURCE_RESULT_END")
