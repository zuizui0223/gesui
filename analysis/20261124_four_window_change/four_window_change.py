#!/usr/bin/env python3
"""Retrospective four-window lagged change test; fit without the 2023–25 response."""
import os
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
import sys, pathlib, argparse, hashlib, json, collections
import numpy as np, pandas as pd
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261124_consumer_forecast"))
import frozen_consumer_forecast as prior
A0,A1,A2,A3="2014-2016","2017-2019","2020-2022","2023-2025"
def read_new(work,stem,period):
    result={}
    for i in range(4):
        fs=list(work.rglob(f"{stem}_full_{i}.json"));assert len(fs)==1
        for r in json.loads(fs[0].read_text()):
            h=str(r["HYBAS_ID_L07"]);g=r["group"]
            assert r["period"]==period and (h,period,g) not in result
            result[h,period,g]=r
    assert len(result)==1332
    return result
def support(hs,fac,period):
    eligible=[];quality=[]
    for h in hs:
        rs=[fac[h,period,g] for g in ("Plantae",)+prior.FAMS]
        ok=all(r["source_status"]=="ok" for r in rs)
        p=rs[0]["source_facetted_occurrence_count"] if ok else 0
        b=sum(r["source_facetted_occurrence_count"] for r in rs[1:]) if ok else 0
        good=ok and p>=50 and b>=20
        eligible.append(good);quality.append(good and min(r["source_facet_coverage"] or 0 for r in rs)>=.8)
    return np.array(eligible),np.array(quality)
def resource(G,P0,P1):
    # Rows are butterfly species; weights sum to one over catalogue hosts.
    W=G/G.sum(1)[:,None];d=P1-P0
    D=d@W.T;V=np.abs(d)@W.T;Q=(d*d)@(W*W).T
    den=V*V-Q
    # Prior protocol support threshold is on unnormalized equal-weight host pairs.
    tol=1e-14/(G.sum(1)**2)
    C=np.divide(D*D-Q,den,out=np.full_like(D,np.nan),where=den>tol)
    C=np.clip(C,-1,1)
    return D,V,C,P1@W.T
def fit(x,y,xt,w):
    mu=x.mean(0);sd=x.std(0);sd[sd<1e-12]=1
    x=(x-mu)/sd;xt=(xt-mu)/sd
    m=prior.Ridge(alpha=10,fit_intercept=True).fit(x,y,sample_weight=w)
    p=m.predict(xt)
    X=np.column_stack((np.ones(len(x)),x));Xt=np.column_stack((np.ones(len(xt)),xt))
    pen=np.eye(X.shape[1])*10;pen[0,0]=0
    c=np.linalg.solve(X.T@(w[:,None]*X)+pen,X.T@(w*y))
    np.testing.assert_allclose(p,Xt@c,atol=1e-8,rtol=0)
    return p
def weights(h,groups):
    ct=collections.Counter(h);gh=collections.Counter(groups[list(ct)])
    w=np.array([1/(ct[i]*gh[groups[i]]) for i in h])
    return w/w.mean()
def unit():
    G=np.array([[1.,1.,1.,0.],[1.,0.,0.,0.]])
    p=np.array([[0.,1.,0.,0.],[0.,0.,0.,0.]])
    q=np.array([[1.,0.,.5,0.],[1.,1.,1.,0.]])
    D,V,C,L=resource(G,p,q)
    np.testing.assert_allclose(D[:,0],[1/6,1])
    np.testing.assert_allclose(V[:,0],[2.5/3,1])
    np.testing.assert_allclose(C[:,0],[-.5,1.])
    assert np.isnan(C[:,1]).all()
    np.testing.assert_allclose(weights(np.array([0,0,1,2]),np.array(["a","a","b"])),[.5,.5,1.,2.])
    rng=np.random.default_rng(1);x=rng.normal(size=(100,4));y=rng.normal(size=100)
    assert len(fit(x,y,x[:5],np.ones(100)))==5
    print("UNIT PASS: sign concordance, missing support, H6 weighting, independent ridge")
def train(work,out):
    hs,d1,d2,quality,fold,bk,pk,pf,bf,A,C,mapped,common,fac=prior.load_old(work)
    groups=d1.h6.to_numpy();earlier=read_new(work,"earlier",A0);fac.update(earlier)
    eg,eq=support(hs,fac,A0)
    ids=np.flatnonzero(eg);hs=[hs[i] for i in ids];d1=d1.iloc[ids];d2=d2.iloc[ids];quality=quality[ids]&eq[ids];fold=fold[ids];groups=groups[ids]
    if quality.sum()<40 or len(set(groups[quality]))<20:
        (out/"STOP.json").write_text(json.dumps({"status":"STOP_HISTORICAL_SUPPORT","quality_H7":int(quality.sum()),"quality_H6":len(set(groups[quality]))}))
        return
    z=[prior.arrays(hs,t,fac,bk,pk) for t in (A0,A1,A2)]
    B=[x[0] for x in z];P=[x[1] for x in z];PC=[x[2] for x in z]
    static=["W","area","riverflow","crop","urban","naturalout","EQR_original_site_equal_basin_period_mean"]
    panels=[("common_Lycaenidae",common,[("old",A),("Clarke",C)])]
    panels += [(fam,mapped&(bf==fam),[("old",A)]) for fam in prior.FAMS]
    rows=[];stops=[]
    for panel,mask,graphs in panels:
        kk=np.flatnonzero(mask);k=len(kk);keys=np.asarray(bk)[kk]
        fmask=bf==bf[kk[0]]
        good=[np.ones(len(hs),bool) for _ in range(3)]
        if panel=="common_Lycaenidae":
            for t in range(3):
                denom=B[t][:,mapped&fmask].sum(1)
                good[t]=(denom>0)&(B[t][:,mask].sum(1)>=.8*denom)
        rr=[{name:resource(G[mask],P[t],P[t+1]) for name,G in graphs} for t in (0,1)]
        identifiable=[np.logical_and.reduce([np.isfinite(v[2]) for v in r.values()]) for r in rr]
        def base(t,h,s):
            b0,b1=B[t],B[t+1];df=[d1,d2][t]
            # All features known at end of predictor interval, same transformations in both periods.
            g=np.column_stack((df[static].to_numpy(float),np.log1p(z[t+1][3]),np.log1p(z[t+1][4]),np.log1p(z[t+1][3])-np.log1p(z[t][3]),np.log1p(z[t+1][4])-np.log1p(z[t][4]),z[t+1][5],z[t+1][5]-z[t][5],P[t+1].sum(1),P[t+1].sum(1)-P[t].sum(1),b1[:,fmask].sum(1),b1[:,fmask].sum(1)-b0[:,fmask].sum(1)))
            return np.column_stack((b1[h,kk[s]],b1[h,kk[s]]-b0[h,kk[s]],g[h],np.eye(k)[s]))
        for f in range(5):
            trainmask=quality&(fold!=f)&good[0]&good[1]
            testmask=(fold==f)&good[0]&good[1]&good[2]
            for endpoint in ("direction","coherence"):
                # Both catalogues share identical rows in each interval. Train selection uses only its own predictor interval.
                m0=np.broadcast_to(trainmask[:,None],(len(hs),k)).copy()
                m1=np.broadcast_to(testmask[:,None],(len(hs),k)).copy()
                if endpoint=="coherence":
                    m0 &= identifiable[0]
                    m1 &= identifiable[0]&identifiable[1]
                h,s=np.where(m0);ht,st=np.where(m1)
                ntr=len(set(h));nh6=len(set(groups[h]));nsp=len(set(s))
                if ntr<30 or nh6<15 or nsp<(20 if panel=="common_Lycaenidae" else 5) or not len(ht):
                    stops.append(dict(panel=panel,fold=f,endpoint=endpoint,train_h7=ntr,train_h6=nh6,train_taxa=nsp,test_rows=len(ht),status="STOP_SUPPORT"));continue
                assert not set(groups[h])&set(groups[ht])
                x=base(0,h,s);xt=base(1,ht,st)
                y=B[2][h,kk[s]]-B[1][h,kk[s]]
                w=weights(h,groups);now=B[2][ht,kk[st]]
                for name,G0 in graphs:
                    G=G0[mask];D0,V0,C0,L0=rr[0][name];D1,V1,C1,L1=rr[1][name]
                    def pick(v,t):return v[h,s] if t==0 else v[ht,st]
                    if endpoint=="direction":
                        N=prior.nullgraph(G,PC[0],groups,pf,pk,np.unique(h))
                        dn0=resource(N,P[0],P[1])[0];dn1=resource(N,P[1],P[2])[0]
                        x0,x1=x,xt
                        ref=fit(np.column_stack((x0,pick(dn0,0))),y,np.column_stack((x1,pick(dn1,1))),w)
                        alt=fit(np.column_stack((x0,pick(D0,0))),y,np.column_stack((x1,pick(D1,1))),w)
                    else:
                        deg=G.sum(1)
                        x0=np.column_stack((x,pick(D0,0),pick(V0,0),pick(C0,0),pick(L0,0),np.log1p(deg[s])))
                        x1=np.column_stack((xt,pick(D1,1),pick(V1,1),pick(C1,1),pick(L1,1),np.log1p(deg[st])))
                        ref=fit(x0,y,x1,w)
                        alt=fit(np.column_stack((x0,pick(D0*C0,0))),y,np.column_stack((x1,pick(D1*C1,1))),w)
                    baseline=fit(x0,y,x1,w)
                    for i in range(len(ht)):
                        rows.append(dict(panel=panel,graph=name,endpoint=endpoint,fold=f,h=hs[ht[i]],h6=groups[ht[i]],key=keys[st[i]],historical_quality=bool(quality[ht[i]]),now=now[i],baseline=np.clip(baseline[i],-now[i],1-now[i]),reference=np.clip(ref[i],-now[i],1-now[i]),alternative=np.clip(alt[i],-now[i],1-now[i])))
    if not rows:
        (out/"STOP.json").write_text(json.dumps({"status":"STOP_MODEL_SUPPORT","details":stops},indent=2));return
    p=out/"change_predictions.csv";pd.DataFrame(rows).to_csv(p,index=False)
    (out/"train_stops.json").write_text(json.dumps(stops,indent=2))
    (out/"prediction_receipt.json").write_text(json.dumps({"sha256":prior.digest(p),"known_keys":bk,"rows":len(rows),"T3_read":False,"historical_quality_H7":int(quality.sum()),"historical_quality_H6":len(set(groups[quality]))},indent=2))
    print("PREDICTIONS SEALED",len(rows),prior.digest(p),"TRAIN STOPS",len(stops))
def score(work,out):
    receipt=json.loads((work/"prediction_receipt.json").read_text());p=work/"change_predictions.csv"
    assert receipt["T3_read"] is False and prior.digest(p)==receipt["sha256"]
    pred=pd.read_csv(p,dtype={"h":str,"h6":str,"key":str})
    fac=read_new(work,"third",A3);hs=sorted(set(pred.h));eg,eq=support(hs,fac,A3);quality=dict(zip(hs,eq));keys=receipt["known_keys"]
    y={};cov=[]
    for h,ok in zip(hs,eg):
        if not ok:continue
        bc=collections.Counter()
        for g in prior.FAMS:bc.update(fac[h,A3,g]["source_speciesKey_counts"])
        n=sum(bc.values());v=np.array([bc.get(k,0) for k in keys])
        for k,x in zip(keys,prior.incl(n,v,20)):y[h,k]=x
        cov.append({"h":h,"known_key_share":v.sum()/n})
    pred["future"]=[y.get((h,k),np.nan) for h,k in zip(pred.h,pred.key)]
    pred["quality"]=pred.historical_quality&pred.h.map(quality)
    pred=pred[pred.future.notna()].copy()
    for name in ("baseline","reference","alternative"):pred["e_"+name]=(pred.future-pred.now-pred[name])**2
    pred["e_no_change"]=(pred.future-pred.now)**2
    rows=[];losses=[]
    for sample,d in [("quality",pred[pred.quality]),("all",pred)]:
        for (panel,endpoint),g in d.groupby(["panel","endpoint"]):
            if panel=="common_Lycaenidae":
                a=g[g.graph=="old"][["h","key","fold"]].sort_values(["h","key"]).reset_index(drop=True)
                b=g[g.graph=="Clarke"][["h","key","fold"]].sort_values(["h","key"]).reset_index(drop=True)
                assert a.equals(b)
            for graph,q in g.groupby("graph"):
                nH=q.h.nunique();n=q.h6.nunique();nsp=q.key.nunique()
                if nH<40 or n<20 or nsp<(20 if panel=="common_Lycaenidae" else 5) or q.fold.nunique()<5:
                    rows.append(dict(sample=sample,panel=panel,endpoint=endpoint,graph=graph,H7=nH,H6=n,taxa=nsp,status="STOP_SUPPORT"));continue
                v=q.groupby(["h6","h"])[["e_baseline","e_reference","e_alternative","e_no_change"]].mean().groupby("h6").mean()
                m=v.to_numpy();mu=m.mean(0);assert min(mu)>0
                ii=np.random.default_rng(20261124).integers(0,n,(5000,n));bb=m[ii].mean(1)
                refs=[("structured_null",1),("baseline",0)] if endpoint=="direction" else [("no_interaction",1)]
                for label,j in refs+[("no_change",3)]:
                    u=100*(bb[:,j]-bb[:,2])/bb[:,j]
                    rows.append(dict(sample=sample,panel=panel,endpoint=endpoint,graph=graph,H7=nH,H6=n,taxa=nsp,status="FIT",reference=label,gain=100*(mu[j]-mu[2])/mu[j],lo=np.quantile(u,.025),hi=np.quantile(u,.975)))
                for h,r in v.iterrows():losses.append(dict(sample=sample,panel=panel,endpoint=endpoint,graph=graph,h6=h,**r.to_dict()))
    pd.DataFrame(rows).to_csv(out/"results.csv",index=False);pd.DataFrame(losses).to_csv(out/"losses.csv",index=False);pd.DataFrame(cov).to_csv(out/"coverage.csv",index=False)
    print(pd.DataFrame(rows).round(4).to_string(index=False))
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("stage",choices=["unit","train","score"]);p.add_argument("--work",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);a=p.parse_args()
    if a.stage=="unit":unit()
    else:
        a.out.mkdir(parents=True,exist_ok=True);{"train":train,"score":score}[a.stage](a.work,a.out)
