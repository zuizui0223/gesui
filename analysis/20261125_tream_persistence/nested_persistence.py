#!/usr/bin/env python3
"""TREAM four-window nested composition: source-locked incidence and environment."""
import argparse, collections, csv, hashlib, io, json, pathlib, sys, tempfile, zipfile
import numpy as np
import pandas as pd
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"20261125_tream_source"))
import acquire_tream as src
SEED=20261125; BOOT=2000
ORDERS={"Ephemeroptera","Plecoptera","Trichoptera","Diptera","Coleoptera","Odonata","Megaloptera","Lepidoptera","Hymenoptera","Neuroptera"}
SETS={"insect":ORDERS,"EPT":{"Ephemeroptera","Plecoptera","Trichoptera"}}
MISSING={"","na","nan","null","none"}
FIELDS=["ppt_mm_12moPrior","tmax_C_12moPrior","crop_perc_upstr","urban_perc_upstr"]
WS=[str(a)+"_"+str(b) for a,b in src.WINDOWS]
def digest(b):return hashlib.sha256(b).hexdigest()
def key(a):return (src.sitekey(a["site_id"]),src.year(a["year"]))
def win(y):return next((i for i,(a,b) in enumerate(src.WINDOWS) if a<=y<=b),None)
def native_rows(p):
    for _,r in src.rows(p):
        yield from r
def truth(v):return str(v).lower()=="true"
def finite(v):
    try:
        x=float(src.clean(v))
        return x if np.isfinite(x) else None
    except (ValueError,TypeError):return None
def signed(work):
    b=(work/"source.zip").read_bytes()
    assert digest(b)=="4abff94ca247f0b8c22e56c5527940e851d2925128440e4fbae2bcd4e2cd3b5e"
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        receipt=json.loads(z.read("source_receipt.json"))
        assert receipt["status"]=="SOURCE_STOP_KEY_RECONCILIATION"
        assert receipt["scoped_feasibility"]["scoped_key_acceptance"]
        assert receipt["scoped_feasibility"]["candidate_siteyears"]==6273
        return receipt,{n:z.read(n) for n in ("candidate_roster.csv","scoped_siteyear_label_support.csv","source_key_reconciliation.csv","scoped_source_coverage.csv")}
def make_sources(d,work,out):
    receipt,tables=signed(work)
    for n in src.NAMES:
        r=src.receive(n,receipt["table_metadata"][n],d/n)
        assert r["sha256"]==receipt["source"][n]["sha256"]
    # Replay original source-only gates, including the global failure, before outcomes.
    verify=out/"source_replay";verify.mkdir(exist_ok=True)
    src.out_global=verify
    source={"status":"SOURCE_STOP"}
    src.census(d,source)
    assert source["status"]==receipt["status"]
    assert source["key_anomaly_counts"]==receipt["key_anomaly_counts"]
    assert source["scoped_feasibility"]["scoped_key_acceptance"]
    assert source["calendar_opportunity_conservative"]==receipt["calendar_opportunity_conservative"]
    assert source["scoped_feasibility"]["nested_label_calendar_support"]==receipt["scoped_feasibility"]["nested_label_calendar_support"]
    for n,b in tables.items():
        assert (verify/n).read_bytes()==b,("Source-only table changed",n)
    roster=pd.read_csv(io.BytesIO(tables["candidate_roster.csv"]),dtype={"site_id":str,"project":str})
    roster=roster[roster.conservative_candidate.map(truth)].set_index("site_id")
    assert len(roster)==563 and "100000145" not in roster.index
    support=pd.read_csv(io.BytesIO(tables["scoped_siteyear_label_support.csv"]),dtype={"site_id":str,"project":str})
    assert len(support)==6273 and not support.duplicated(["site_id","year"]).any()
    K={(r.site_id,int(r.year)) for r in support.itertuples()}
    meta={};dates={}
    for a in native_rows(d/src.NAMES[1]):
        if src.clean(a.get("site_id")) not in roster.index:continue
        k=key(a)
        if k not in K:continue
        assert k not in meta and src.sitekey(a["site_id_wMissing"])==k[0] and src.year(a["year_wMissing"])==k[1]
        z=np.asarray([finite(a[f]) for f in FIELDS],float)
        assert np.isfinite(z).all()
        meta[k]=z;dates[k]=(src.clean(a["day"]),src.clean(a["month"]))
    assert set(meta)==K
    groups=collections.defaultdict(set)
    for a in native_rows(d/src.NAMES[2]):
        f=src.clean(a["Family"]);g=src.clean(a["Group"])
        if f.lower() not in MISSING and g.lower() not in MISSING:groups[f].add(g)
    conflict={f:sorted(gs) for f,gs in groups.items() if len(gs)>1}
    assert conflict==receipt["source_census"]["families_in_multiple_groups"]
    # No correction of family labels, and no use of numerical abundance above.
    family={kind:sorted(f for f,gs in groups.items() if len(gs)==1 and next(iter(gs)) in allowed) for kind,allowed in SETS.items()}
    fset={kind:set(names) for kind,names in family.items()}
    detected={kind:collections.defaultdict(set) for kind in SETS}
    invalid={kind:set() for kind in SETS}
    counts={kind:collections.defaultdict(collections.Counter) for kind in SETS}
    all_rows=collections.Counter()
    for a in native_rows(d/src.NAMES[2]):
        if src.clean(a.get("site_id")) not in roster.index:continue
        k=key(a)
        if k not in K:continue
        assert dates[k]==(src.clean(a["day"]),src.clean(a["month"]))
        all_rows[k]+=1
        g=src.clean(a["Group"]);f=src.clean(a["Family"])
        for kind,allowed in SETS.items():
            if g not in allowed:continue
            q=counts[kind][k];q["order_rows"]+=1
            if f not in fset[kind] or groups[f]!={g}:continue
            q["clean_rows"]+=1
            x=finite(a["abundance"])
            if x is None or x<0:
                invalid[kind].add(k);q["invalid_abundance_rows"]+=1
            elif x>0:
                detected[kind][k].add(f);q["positive_rows"]+=1
            else:q["zero_rows"]+=1
    for r in support.itertuples():
        k=(r.site_id,int(r.year));assert all_rows[k]==r.label_rows
        for kind,field1,field2 in [("insect","insect_label_rows","clean_insect_label_rows"),("EPT","EPT_label_rows","clean_EPT_label_rows")]:
            q=counts[kind][k]
            assert q["order_rows"]==getattr(r,field1) and q["clean_rows"]==getattr(r,field2),(kind,k,q)
    return receipt,roster,support,K,meta,groups,family,detected,invalid,counts
def build_panel(kind,receipt,roster,support,K,meta,groups,family,detected,invalid,counts,out):
    fam=family[kind];j={f:i for i,f in enumerate(fam)}
    eligible=collections.defaultdict(lambda:collections.defaultdict(list))
    label=collections.defaultdict(lambda:collections.defaultdict(int))
    year_rows=[]
    for k,y in sorted(K):
        site=k;key0=(k,y);t=win(y);assert t is not None
        q=counts[kind][key0]
        label[k][t]+=int(q["clean_rows"]>0)
        valid=bool(detected[kind][key0]) and key0 not in invalid[kind]
        if valid:eligible[k][t].append(y)
        year_rows.append({"universe":kind,"site_id":k,"year":y,"window":t,"valid_positive_year":valid,"detected_families":len(detected[kind][key0]),"invalid_year":key0 in invalid[kind],**{f:q[f] for f in ("order_rows","clean_rows","positive_rows","zero_rows","invalid_abundance_rows")}})
    lab=sorted(k for k in roster.index if min(label[k][t] for t in range(4))>=2)
    rec=receipt["scoped_feasibility"]["nested_label_calendar_support"][kind]
    assert len(lab)==rec["sites"] and digest(("\n".join(lab)+"\n").encode())==rec["roster_sha256"]
    ids=sorted(k for k in roster.index if min(len(eligible[k][t]) for t in range(4))>=2)
    assert set(ids)<=set(lab)
    X=np.zeros((len(ids),4,len(fam)));Z=np.zeros((len(ids),4,4))
    for i,k in enumerate(ids):
        for t in range(4):
            years=eligible[k][t]
            for y in years:
                f=detected[kind][(k,y)]
                X[i,t,[j[v] for v in f]]+=1/(len(f)*len(years))
            Z[i,t]=np.mean([meta[k,y] for y in years],axis=0)
    np.testing.assert_allclose(X.sum(axis=2),1,atol=1e-12,rtol=0)
    meta0=roster.loc[ids].copy()
    # Source labels fixed globally; one project belongs to exactly one country here.
    assert meta0.reset_index().groupby("project").country.nunique().max()==1
    cc=meta0.country.value_counts().to_dict()
    gate=len(ids)>=100 and len(cc)>=5 and meta0.project.nunique()>=10 and sum(v>=20 for v in cc.values())>=3
    cohorts=pd.DataFrame([dict(universe=kind,site_id=k,country=roster.loc[k,"country"],project=roster.loc[k,"project"],label_eligible=k in set(lab),positive_eligible=k in set(ids),**{"n_valid_"+WS[t]:len(eligible[k][t]) for t in range(4)}) for k in roster.index])
    cohorts.to_csv(out/(kind+"_roster.csv"),index=False)
    pd.DataFrame(year_rows).to_csv(out/(kind+"_year_support.csv"),index=False)
    used={(k,y) for k in ids for t in range(4) for y in eligible[k][t]}
    aggregate={name:sum(q[name] for k,q in counts[kind].items() if k in used) for name in ("order_rows","clean_rows","positive_rows","zero_rows","invalid_abundance_rows")}
    panel={"universe":kind,"source_label_sites":len(lab),"positive_sites":len(ids),"countries":cc,"projects":int(meta0.project.nunique()),"informative_projects":int((meta0.project.value_counts()>=2).sum()),"siteyears":len(used),"window_years":[sum(len(eligible[k][t]) for k in ids) for t in range(4)],"family_universe":len(fam),"family_universe_never_detected":int(np.sum(X.sum(axis=(0,1))==0)),"invalid_source_years":len(invalid[kind]),"no_positive_source_years":len([k for k in K if not detected[kind][k]]),"used_row_counts":aggregate,"roster_sha256":digest(("\n".join(ids)+"\n").encode()),"source_gate_pass":bool(gate)}
    # Save derived standardized family incidence, not the raw abundance rows.
    np.savez_compressed(out/(kind+"_panel.npz"),family=np.asarray(fam),order=np.asarray([next(iter(groups[f])) for f in fam]),site=np.asarray(ids),X=X,Z=Z,project=meta0.project.to_numpy(str),country=meta0.country.to_numpy(str))
    return X,Z,meta0,fam,panel
def project(X,fam,groups):
    B=np.zeros_like(X)
    for g in sorted({next(iter(groups[f])) for f in fam}):
        idx=[i for i,f in enumerate(fam) if groups[f]=={g}]
        B[:,:,idx]=X[:,:,idx].mean(axis=2,keepdims=True)
    W=X-B
    np.testing.assert_allclose(B+W,X,atol=1e-12,rtol=0)
    np.testing.assert_allclose(np.sum(B*W,axis=2),0,atol=1e-12,rtol=0)
    return {"total":X,"between_order":B,"within_order":W}
def stats(G):
    ev=(G[...,0,0]+2*G[...,0,1]+G[...,1,1])/4
    lv=(G[...,2,2]+2*G[...,2,3]+G[...,3,3])/4
    cov=G[...,:2,2:].sum(axis=(-2,-1))/4
    a=np.divide(2*cov,ev+lv,out=np.full_like(cov,np.nan),where=(ev+lv)>1e-20)
    denom=np.sqrt(np.maximum(ev,0)*np.maximum(lv,0))
    b=np.divide(cov,denom,out=np.full_like(cov,np.nan),where=denom>1e-20)
    v=np.trace(G,axis1=-2,axis2=-1)/4
    c=G[...,np.triu_indices(4,1)[0],np.triu_indices(4,1)[1]].mean(axis=-1)
    d=np.divide(c,v,out=np.full_like(v,np.nan),where=v>1e-20)
    return {"half_reproducibility":a,"half_cosine":b,"all_window_covariance_ratio":d,"mean_window_spatial_variance":v}
def scalar(v):return float(v) if np.isfinite(v) else None
def pers(Xs,meta,kind,out):
    ps=sorted(meta.project.unique());ii=[np.flatnonzero(meta.project.to_numpy()==p) for p in ps]
    n=len(ps);rng=np.random.default_rng(SEED);q=rng.exponential(size=(BOOT,n));q/=q.sum(axis=1,keepdims=True)
    equal=np.full((1,n),1/n)
    loo=np.ones((n,n))-np.eye(n);loo/=loo.sum(axis=1,keepdims=True)
    informative=np.asarray([len(x)>1 for x in ii])
    q0=np.concatenate((equal,q,loo))
    cache={};byproj=[];summ=[];raw={}
    for name,X in Xs.items():
        M=np.stack([X[idx].mean(axis=0) for idx in ii])
        G=np.stack([np.einsum("itf,isf->ts",X[idx],X[idx])/len(idx) for idx in ii])
        Gw=G-np.einsum("ptf,psf->pts",M,M)
        # Independent direct centering check for each project.
        for i,idx in enumerate(ii):
            Z=X[idx]-X[idx].mean(axis=0)
            np.testing.assert_allclose(Gw[i],np.einsum("itf,isf->ts",Z,Z)/len(idx),atol=1e-12,rtol=0)
            vals=stats(Gw[i])
            byproj.append({"universe":kind,"component":name,"project":ps[i],"sites":len(idx),**{k:scalar(v) for k,v in vals.items()}})
        for scope in ("all_space","within_project"):
            if scope=="within_project":GG=np.einsum("bp,pts->bts",q0,Gw)
            else:
                MM=np.einsum("bp,ptf->btf",q0,M)
                GG=np.einsum("bp,pts->bts",q0,G)-np.einsum("btf,bsf->bts",MM,MM)
            cache[scope,name]=GG;raw[scope,name]=stats(GG)
    for scope in ("all_space","within_project"):
        np.testing.assert_allclose(cache[scope,"total"],cache[scope,"between_order"]+cache[scope,"within_order"],atol=1e-12,rtol=0)
        for name in Xs:
            vals=raw[scope,name].copy()
            vals["variance_share"]=vals["mean_window_spatial_variance"]/raw[scope,"total"]["mean_window_spatial_variance"]
            for k,v in vals.items():
                boot=v[1:1+BOOT];valid=boot[np.isfinite(boot)];l=v[1+BOOT:][informative];l=l[np.isfinite(l)]
                lo,hi=np.quantile(valid,[.025,.975]) if len(valid) else (np.nan,np.nan)
                summ.append(dict(universe=kind,scope=scope,component=name,metric=k,estimate=scalar(v[0]),lo=scalar(lo),hi=scalar(hi),valid_draws=len(valid),loo_min=scalar(l.min()) if len(l) else None,loo_max=scalar(l.max()) if len(l) else None))
        for k in ("half_reproducibility","half_cosine","all_window_covariance_ratio"):
            v=raw[scope,"between_order"][k]-raw[scope,"within_order"][k];valid=v[1:1+BOOT];valid=valid[np.isfinite(valid)];l=v[1+BOOT:][informative];l=l[np.isfinite(l)]
            lo,hi=np.quantile(valid,[.025,.975]) if len(valid) else (np.nan,np.nan)
            summ.append(dict(universe=kind,scope=scope,component="between_minus_within",metric=k,estimate=scalar(v[0]),lo=scalar(lo),hi=scalar(hi),valid_draws=len(valid),loo_min=scalar(l.min()) if len(l) else None,loo_max=scalar(l.max()) if len(l) else None))
    pd.DataFrame(byproj).to_csv(out/(kind+"_project_persistence.csv"),index=False)
    return summ
def center(Y,idx):
    Y=Y.copy()
    for i in idx:Y[i]-=Y[i].mean(axis=0)
    return Y
def env_prediction(Xs,Z,meta,kind,out):
    projects=meta.project.to_numpy();countries=meta.country.to_numpy()
    ps=sorted(set(projects));idx=[np.flatnonzero(projects==p) for p in ps]
    informative=np.zeros(len(meta),bool)
    w=np.zeros(len(meta))
    for i in idx:
        if len(i)>=2:informative[i]=True
        w[i]=1/len(i)
    testc=sorted(c for c in set(countries) if sum((countries==c)&informative)>=10)
    if len(testc)<5:return [],{"status":"SOURCE_STOP_LT5_COUNTRIES","eligible_countries":testc}
    contrasts=[]
    def pair(A,mode,start):
        return (A[:,start]+A[:,start+1])/2 if mode=="state" else A[:,start+1]-A[:,start]
    for mode in ("state","change"):
        A=center(pair(Z,mode,0),idx);D=center(pair(Z,mode,2),idx)
        Ys={k:center(pair(v,mode,0),idx) for k,v in Xs.items()}
        Ts={k:center(pair(v,mode,2),idx) for k,v in Xs.items()}
        for c in testc:
            tr=countries!=c;te=countries==c
            ww=w[tr]
            # Site weights sum to number of training projects; fixed Ridge10.
            scale=np.sqrt(np.sum(A[tr]**2*ww[:,None],axis=0)/ww.sum())
            active=scale>1e-12;scale=np.where(active,scale,1.)
            a=A[tr]/scale;d=D[te]/scale
            a[:,~active]=0;d[:,~active]=0
            H=a.T@(ww[:,None]*a)+10*np.eye(a.shape[1]);pred={}
            for k in Xs:
                b=np.linalg.solve(H,a.T@(ww[:,None]*Ys[k][tr]))
                pred[k]=d@b
                target=Ts[k][te]
                e0=np.sum(target**2,axis=1);e1=np.sum((target-pred[k])**2,axis=1)
                wte=w[te]
                mse0=float(np.dot(wte,e0)/wte.sum());mse1=float(np.dot(wte,e1)/wte.sum())
                contrasts.append(dict(universe=kind,mode=mode,country=c,component=k,train_sites=int(tr.sum()),test_sites=int(te.sum()),test_projects=len(set(projects[te])),baseline_mse=mse0,added_mse=mse1,active_predictors=int(active.sum())))
            np.testing.assert_allclose(pred["total"],pred["between_order"]+pred["within_order"],atol=1e-12,rtol=0)
    data=pd.DataFrame(contrasts);data.to_csv(out/(kind+"_environment_country_losses.csv"),index=False)
    rng=np.random.default_rng(SEED);draw=rng.integers(0,len(testc),(BOOT,len(testc)));results=[]
    for k in Xs:
        gains={}
        for mode in ("state","change"):
            t=data[(data.component==k)&(data["mode"]==mode)].set_index("country").loc[testc]
            a=t.baseline_mse.to_numpy();b=t.added_mse.to_numpy()
            assert a.mean()>1e-20
            v=100*(a[draw].mean(axis=1)-b[draw].mean(axis=1))/a[draw].mean(axis=1)
            pt=100*(a.mean()-b.mean())/a.mean();gains[mode]=(pt,v)
            results.append(dict(universe=kind,component=k,mode=mode,gain=pt,lo=float(np.quantile(v,.025)),hi=float(np.quantile(v,.975)),improved_countries=int((a>b).sum()),countries=len(testc)))
        pt=gains["state"][0]-gains["change"][0];v=gains["state"][1]-gains["change"][1]
        results.append(dict(universe=kind,component=k,mode="state_minus_change",gain=pt,lo=float(np.quantile(v,.025)),hi=float(np.quantile(v,.975)),improved_countries=None,countries=len(testc)))
    return results,{"status":"COMPLETE","eligible_countries":testc,"ridge_alpha":10,"test_centering_not_absolute_prediction":True}
def main(work,out):
    out.mkdir(parents=True,exist_ok=True)
    result={"status":"STOP","source_scope":"isolated 563 site roster, global package STOP preserved","not_causal":True}
    with tempfile.TemporaryDirectory() as tmp:
        r,roster,support,K,meta,groups,family,detected,invalid,counts=make_sources(pathlib.Path(tmp),work,out)
        summaries=[];environment=[];panels={}
        for kind in SETS:
            X,Z,meta0,fam,panel=build_panel(kind,r,roster,support,K,meta,groups,family,detected,invalid,counts,out)
            panels[kind]=panel
            if not panel["source_gate_pass"]:
                panel["analysis_status"]="STOP_POSITIVE_SUPPORT"
                continue
            Xs=project(X,fam,groups)
            s=pers(Xs,meta0,kind,out);summaries+=s
            e,qa=env_prediction(Xs,Z,meta0,kind,out);environment+=e
            panel["environment_status"]=qa;panel["analysis_status"]="COMPLETE"
        pd.DataFrame(summaries).to_csv(out/"persistence_summary.csv",index=False)
        pd.DataFrame(environment).to_csv(out/"environment_summary.csv",index=False)
        result.update(status="COMPLETE" if panels["insect"]["source_gate_pass"] else "STOP_MAIN_POSITIVE_SUPPORT",panels=panels,persistence=summaries,environment=environment)
        (out/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
        print("TREAM ECOLOGY "+json.dumps({"status":result["status"],"panels":panels,"persistence":[s for s in summaries if s["scope"]=="within_project" and s["metric"]=="half_reproducibility"],"environment":environment},allow_nan=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--work",type=pathlib.Path,required=True);p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args();main(a.work,a.out)
