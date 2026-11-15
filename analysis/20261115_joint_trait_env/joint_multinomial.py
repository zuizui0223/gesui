#!/usr/bin/env python3
"""Penalized conditional multinomial fourth-corner analysis on real EPA EPT counts.

Observed native EPA 2008/2013/2018 site×TAXA_ID count composition; official
EPA HABIT CN trait, within-order centered; environmental embeddedness×trait.
Baseline includes taxon-specific training prevalence and order×embedding,
trait×other measured chemistry/geography/PHab, and eco×trait/order controls.

This is a constrained, observable trait–environment joint-composition model,
NOT a full latent JSDM, a causal community assembly model, or sewage exposure.
"""
from __future__ import annotations
import argparse
import collections
import hashlib
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp

ROOT=Path(os.environ.get("RUNNER_TEMP","/tmp"))
R=ROOT/"gesui_threecycle_traits"
C=ROOT/"gesui_threecycle_rarefied_taxa"/"AUDIT_ONLY_ORIGINAL_EPA_3CYCLE_TAXON_COUNTS.csv"
T=R/"AUDIT_ONLY_EPA_ORIGINAL_3CYCLE_TAXA_FFG_HABIT_PTV.csv"
S=R/"AUDIT_ONLY_3CYCLE_ORIGINAL_EPT_FFG_HABIT_PTV_SAMESITE.csv"
ERAS=("2008_09","2013_14","2018_19")
FIELD="HABIT_norm"
HAB=set("BU CB CN DV PK SK SP SW".split())
EPT=("EPHEMEROPTERA","PLECOPTERA","TRICHOPTERA")
NUMERIC=["XEMBED","log_tn","log_tp","lat","lon","stream_order","log_cond",
         "season_sin","season_cos","XCDENMID","XPCAN","XFC_OHV",
         "W1_HALL","log_width","year"]
LAMBDA=.025
FOLDS=4

def score(value):
    c=str(value).strip().upper()
    if not c:return np.nan
    modes=list(dict.fromkeys(x.strip() for x in c.split(",") if x.strip()))
    if not modes or any(v not in HAB for v in modes):return np.nan
    return float("CN" in modes)/len(modes)

def source():
    sites=pd.read_csv(S,dtype={"source_cycle":str,"source_UID":str,"huc8":str,
                               "state":str,"eco9":str})
    if len(sites)!=2209 or sites.duplicated(["source_cycle","source_UID"]).any():
        raise ValueError("Original source n2209 rare200 sites drift")
    count=pd.read_csv(C,dtype={"cycle":str,"UID":str,"key":str,"order":str},
                      keep_default_na=False)
    count=count[count.order.isin(EPT)].copy()
    count["count"]=pd.to_numeric(count["count"],errors="raise")
    if count.duplicated(["cycle","UID","key"]).any() or not count["count"].gt(0).all():
        raise ValueError("Original source per-taxon counts invalid")
    trait=pd.read_csv(T,dtype=str,keep_default_na=False)
    if trait.duplicated(["source_cycle","key"]).any():
        raise ValueError("True EPA taxonomy/trait nonunique")
    count=count.merge(trait[["source_cycle","key","order",FIELD]],
        left_on=["cycle","key","order"],right_on=["source_cycle","key","order"],
        how="left",validate="many_to_one")
    count["clinger"]=count[FIELD].map(score)
    return sites,count

def matrices_for_era(sites,count,era):
    d=sites[sites.source_cycle==era].sort_values("source_UID").copy()
    if d.empty:raise ValueError("Unknown source survey era")
    c=count[(count.cycle==era)&count.UID.isin(set(d.source_UID))].copy()
    all_count=c.groupby("UID")["count"].sum()
    n_all=c.groupby("UID")["key"].nunique()
    known=c.loc[np.isfinite(c.clinger)].copy()
    known_count=known.groupby("UID")["count"].sum()
    n_known=known.groupby("UID")["key"].nunique()
    total=all_count.reindex(d.source_UID).fillna(0).to_numpy(float)
    matched=known_count.reindex(d.source_UID).fillna(0).to_numpy(float)
    n_k=n_known.reindex(d.source_UID).fillna(0).to_numpy(int)
    n_t=n_all.reindex(d.source_UID).fillna(0).to_numpy(int)
    # Post-result conservative quality sensitivity from prior source-locked PR98:
    # >=3 scored source taxa, >=80% of taxa and specimen count scored,
    # >=10 original counted EPT. This is NOT a fresh preregistration.
    good=(matched/np.maximum(total,1)>=.8-1e-10)&(n_k>=3)&\
         (n_k/np.maximum(n_t,1)>=.8-1e-10)&(total>=10)
    d=d.loc[good].copy().reset_index(drop=True)
    expect={"2008_09":(377,313,44),"2013_14":(767,526,47),
            "2018_19":(624,423,46)}[era]
    if (len(d),d.huc8.nunique(),d.state.nunique())!=expect:
        raise ValueError(f"Strict same-site source eligibility changed: {era}: "
                         f"{len(d),d.huc8.nunique(),d.state.nunique()} expected {expect}")
    c=known[known.UID.isin(set(d.source_UID))].copy()
    c["taxon"]=c["key"].astype(str)+"|"+c["order"].astype(str)
    tax=sorted(c.taxon.unique())
    index={k:j for j,k in enumerate(tax)}
    site_index={k:i for i,k in enumerate(d.source_UID)}
    Y=np.zeros((len(d),len(tax)),dtype=float)
    q=np.full(len(tax),np.nan)
    orders=np.empty(len(tax),dtype=object)
    for x in c.itertuples(index=False):
        i=site_index[x.UID];j=index[x.taxon]
        Y[i,j]+=float(x.count)
        if np.isfinite(q[j]) and abs(q[j]-float(x.clinger))>1e-12:
            raise ValueError("Trait assignment inconsistent for original taxon")
        q[j]=float(x.clinger);orders[j]=x.order
    if not np.isfinite(q).all() or (Y.sum(axis=1)<=0).any():
        raise ValueError("True EPA site×taxon map failed")
    # Genuine original count-based CWM must match source from PR83 on same UIDs.
    calc=(Y@q)/Y.sum(axis=1)
    if not np.allclose(calc,d.EPT_clinger_HABIT_fraction.to_numpy(float),atol=1e-10,rtol=0):
        raise ValueError("Old original EPA CWM count not reproduced")
    known_trait=Y>0
    if not np.all(np.sum(known_trait,axis=1)>=3):
        raise ValueError("Incidence source eligibility unexpectedly dropped")
    incidence=known_trait/known_trait.sum(axis=1,keepdims=True)
    abundance=Y/Y.sum(axis=1,keepdims=True)
    qq=np.empty_like(q)
    for o in EPT:
        ix=orders==o
        qq[ix]=q[ix]-q[ix].mean()
    qs=qq.std(ddof=0)
    if qs<1e-8:raise ValueError("No within-order functional trait variation")
    qq/=qs
    # Predict within-order trait alignment, controlling order response to E.
    Q=np.column_stack([qq,(orders==EPT[0]).astype(float),
                       (orders==EPT[1]).astype(float)])
    if np.linalg.matrix_rank(Q)<3:
        raise ValueError("Collapsed source EPT trait / order design")
    return d,Y,abundance,incidence,Q,dict(n_taxa=len(tax),
        source_scored_EPT_count_total=float(Y.sum()),original_strict_UID=len(d),
        HUC8=d.huc8.nunique(),states=d.state.nunique(),
        taxon_clinger_range=[float(q.min()),float(q.max())],
        within_order_trait_source_sd=float(qs))

def covariates(train,test):
    x=[];z=[]
    for name in NUMERIC:
        a=train[name].to_numpy(float);b=test[name].to_numpy(float)
        if not np.isfinite(a).all() or not np.isfinite(b).all():
            raise ValueError("Original source env missing "+name)
        sd=float(a.std(ddof=0))
        if sd<1e-10:sd=1.
        x.append((a-a.mean())/sd);z.append((b-a.mean())/sd)
    cats=train.eco9.astype(str).to_numpy()
    testcat=test.eco9.astype(str).to_numpy()
    for lv in sorted(set(cats))[1:]:
        x.append((cats==lv).astype(float));z.append((testcat==lv).astype(float))
    a=np.column_stack(x);b=np.column_stack(z)
    # Embeddedness is NUMERIC[0] for every train/test fold
    if not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError("Invalid covariate matrix")
    return a,b

def prevalence(F,groups):
    freq=collections.Counter(map(str,groups))
    w=np.array([1./freq[str(x)] for x in groups])
    w/=w.mean()
    prior=.0002
    p=(w[:,None]*F).sum(axis=0)+prior
    p/=p.sum()
    return np.log(p),w

def optimize(F,X,Q,base,w,include_embed,lam=LAMBDA):
    k,m=X.shape[1],Q.shape[1]
    eligible=np.ones((k,m),dtype=bool)
    if not include_embed:eligible[0,0]=False
    shape=(k,m);nparam=int(eligible.sum())
    w=w/w.sum()
    def unpack(z):
        v=np.zeros(shape);v[eligible]=z
        return v
    def calc(z):
        B=unpack(z)
        logits=base[None,:]+(X@B)@Q.T
        ls=logsumexp(logits,axis=1)
        p=np.exp(logits-ls[:,None])
        loss=float(np.dot(w,ls-(F*logits).sum(axis=1))+lam*np.sum(B**2))
        delta=(p-F)*w[:,None]
        grad=(X.T@(delta@Q)+2*lam*B)[eligible]
        return loss,grad
    result=minimize(calc,np.zeros(nparam),jac=True,method="L-BFGS-B",
                    options={"maxiter":240,"gtol":3e-6,"ftol":1e-11,"maxls":30})
    if not result.success and np.linalg.norm(result.jac)>2e-4:
        raise RuntimeError(f"Source multinomial fitting failed: {result.message} gradient {np.linalg.norm(result.jac)}")
    return unpack(result.x),dict(converged=bool(result.success),
        gradient_norm=float(np.linalg.norm(result.jac)),iterations=int(result.nit),
        penalty_strength=float(lam),train_objective=float(result.fun))

def predict(F,X,Q,logprev,B):
    logits=logprev[None,:]+(X@B)@Q.T
    return logsumexp(logits,axis=1)-(F*logits).sum(axis=1)

def evaluate(d,F,Q,typ):
    groups=d.huc8.astype(str).to_numpy()
    folds=np.array([int(hashlib.sha256(z.encode()).hexdigest()[:8],16)%FOLDS for z in groups])
    test_scores={"baseline":np.full(len(d),np.nan),"trait_embed":np.full(len(d),np.nan)}
    theta=[];fit_status=[]
    for f in range(FOLDS):
        tr=folds!=f;te=~tr
        if te.sum()<45 or tr.sum()<170 or len(set(groups[te]))<35:
            raise ValueError("Insufficient original HUC8 holdout support")
        Xt,Xv=covariates(d.loc[tr],d.loc[te])
        logprev,w=prevalence(F[tr],groups[tr])
        B0,s0=optimize(F[tr],Xt,Q,logprev,w,False)
        B1,s1=optimize(F[tr],Xt,Q,logprev,w,True)
        test_scores["baseline"][te]=predict(F[te],Xv,Q,logprev,B0)
        test_scores["trait_embed"][te]=predict(F[te],Xv,Q,logprev,B1)
        theta.append(float(B1[0,0]))
        fit_status.append(dict(fold=f,train_n=int(tr.sum()),test_n=int(te.sum()),
                               train_HUC8=len(set(groups[tr])),test_HUC8=len(set(groups[te])),
                               theta_trait_by_XEMBED=float(B1[0,0]),
                               null=s0,model=s1))
    if not all(np.isfinite(v).all() for v in test_scores.values()):
        raise ValueError("Missing original heldout taxon predictions")
    uniq,inv=np.unique(groups,return_inverse=True)
    n=np.bincount(inv)
    baseline=np.bincount(inv,weights=test_scores["baseline"])/n
    added=np.bincount(inv,weights=test_scores["trait_embed"])/n
    dloss=baseline-added
    obs=float(dloss.mean())
    bootstrap_ix=np.random.default_rng(20261115).integers(0,len(uniq),
                        size=(2000,len(uniq)))
    bb=dloss[bootstrap_ix].mean(axis=1)
    fullX,_=covariates(d,d)
    logprev,w=prevalence(F,groups)
    B,s=optimize(F,fullX,Q,logprev,w,True)
    return dict(type=typ,original_n=len(d),HUC8=len(uniq),
                test_folds=fit_status,theta_full_sample=float(B[0,0]),
                theta_fold_values=theta,theta_source_q_is_order_centered=True,
                CE_HUC8_baseline=float(baseline.mean()),
                CE_HUC8_trait_embed=float(added.mean()),
                CE_gain_nats_per_scored_taxon=float(obs),
                CE_gain_percent=float(100*obs/baseline.mean()),
                CE_gain_HUC8_fixed_error_bootstrap95=list(map(float,np.quantile(bb,[.025,.975]))),
                CE_percent_gain_HUC8_bootstrap95=list(map(float,
                    np.quantile(100*bb/baseline[bootstrap_ix].mean(axis=1),[.025,.975]))),
                posthoc_exploratory=True,no_causal_response=True,
                alpha_fixed_before_model=LAMBDA,model_spec="taxon train-prevalence + order-by-environment + within-order trait-by-other-environment; focal XEMBED-by-clinger added",
                full_fit_status=s)

def selftest():
    rng=np.random.default_rng(15)
    n,J=180,12
    e=rng.normal(size=n)
    q=np.r_[np.zeros(6),np.ones(6)]
    orders=np.r_[np.zeros(3),np.ones(3),np.zeros(3),np.ones(3)]
    Q=np.column_stack([q-q.mean(),orders,np.r_[np.ones(3),np.zeros(9)]])
    X=e[:,None]
    base=np.log(np.ones(J)/J)
    logits=-2*e[:,None]*(q-q.mean())[None,:]
    prob=np.exp(logits-logsumexp(logits,axis=1)[:,None])
    F=np.stack([rng.multinomial(200,p)/200 for p in prob])
    B,_=optimize(F,X,Q,base,np.ones(n),True,lam=.01)
    if not B[0,0]<-1:
        raise AssertionError("Synthetic trait-embeddedness sign test failed")
    b,_=optimize(F,X,Q,base,np.ones(n),False,lam=.01)
    if not predict(F,X,Q,base,B).mean()<predict(F,X,Q,base,b).mean():
        raise AssertionError("Synthetic heldout contrast direction failed")
    print("SYNTHETIC MULTINOMIAL TRAIT-INTERACTION TESTS PASS",B[0,0])

def run(out):
    sites,counts=source()
    results=dict(schema="GESUI_EPA3ERA_CONSTRAINED_MULTINOMIAL_TRAIT_EMBED_V1",
        not_full_JSDM=True,not_causal_or_WWTP_dose=True,
        source_taxon_count_not_field_individual_density=True,
        original_population_2209=True,posthoc_strict_trait_coverage=True,
        observation_scales=["counted_individual_conditional_composition","observed_scored_taxon_incidence_conditional_composition"],
        analysis_status="exploratory_posthoc",results={})
    for era in ERAS:
        d,Y,Fcount,Finc,Q,meta=matrices_for_era(sites,counts,era)
        results["results"][era]={"sample":meta,
           "count":evaluate(d,Fcount,Q,"counted_individual_composition"),
           "incidence":evaluate(d,Finc,Q,"observed_taxon_incidence_composition")}
        print("ERA_CONSTRAINED_MULTINOMIAL",era,json.dumps({
           kind:{k:v for k,v in results["results"][era][kind].items()
               if k in ("theta_full_sample","theta_fold_values","CE_gain_percent",
                        "CE_gain_HUC8_fixed_error_bootstrap95","CE_percent_gain_HUC8_bootstrap95")}
           for kind in ("count","incidence")},sort_keys=True),flush=True)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(results,sort_keys=True,indent=2,allow_nan=False))
    return results

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=ROOT/"gesui_ept_trait_joint"/"EPA3ERA_JOINT_TRAIT_EMBED_RESULTS.json")
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    if args.self_test:selftest()
    else:run(args.output)
