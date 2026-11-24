#!/usr/bin/env python3
"""Frozen historical forecast: plant resource -> later butterfly recording composition.
Stage train has no third-period species counts. Stage score cannot refit.
"""
import os
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
import argparse, hashlib, json, zipfile, pathlib, collections
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.linear_model import Ridge
ROOT=pathlib.Path(__file__).resolve().parents[2]
OLD=ROOT/"analysis/20261124_lycaenidae_resource_network"
FAMS=("Lycaenidae","Nymphalidae","Pieridae","Papilionidae","Hesperiidae")
T0,T1,T2="2017-2019","2020-2022","2023-2025"
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def jzip(p,m):
    with zipfile.ZipFile(p) as z:return json.loads(z.read(m))
def incl(n,c,m):
    n=np.asarray(n,float);c=np.asarray(c,float)
    assert np.all(n>=m) and np.all(c>=0) and np.all(c<=n)
    miss=np.ones_like(c,float)
    for i in range(m):miss*=np.maximum(0,1-c/(n-i))
    return 1-miss
def hsh(s):
    h=2166136261
    for c in s:h=((h^ord(c))*16777619)&0xffffffff
    return h
def load_old(work):
    m=json.loads((OLD/"SOURCE_MANIFEST.json").read_text())
    for x in m["sources"]:
        if x["name"] in ("gesui_host_source.zip","gesui_gbif_facets.zip"):
            assert digest(work/x["name"])==x["sha256"]
    hp=work/"gesui_host_source.zip"
    alias=jzip(hp,"butterfly_usage_aliases.json")
    edges=jzip(hp,"verified_candidate_host_edges.json")
    cross=jzip(hp,"plant_name_crosswalk.json")
    meta={str(x["original_speciesKey"]):x for x in alias}
    bk=sorted(meta);pk0=sorted({str(x["plantKey"]) for x in edges})
    assert len(bk)==303 and len(pk0)==990 and len(edges)==1565
    pf={str(x["accepted_speciesKey"]):x["match"].get("family","UNKNOWN") for x in cross if x["accepted_speciesKey"] is not None}
    rec=json.loads((OLD/"clarke_expanded_mapping_receipt.json").read_text())
    for f,k in [("clarke_edges.csv","original_edge_sha256"),("clarke_expanded_edges.csv","expanded_edge_sha256"),("clarke_expanded_plant_keys.csv","plant_keys_sha256")]:
        assert digest(OLD/f)==rec[k]
    extra=pd.read_csv(OLD/"clarke_expanded_plant_keys.csv",dtype={"plantKey":str})
    pk=pk0+extra.plantKey.tolist();assert len(pk)==len(set(pk))==1176
    pf.update(dict(zip(extra.plantKey,extra.family)))
    bidx={k:i for i,k in enumerate(bk)};pidx={k:i for i,k in enumerate(pk)}
    A=np.zeros((303,1176));C=A.copy();R=A.copy()
    for e in edges:A[bidx[str(e["butterflyKey"])],pidx[str(e["plantKey"])]]=1
    for f,g in [("clarke_edges.csv",R),("clarke_expanded_edges.csv",C)]:
        for e in pd.read_csv(OLD/f,dtype={"butterflyKey":str,"plantKey":str}).itertuples():
            if e.status<=3:g[bidx[e.butterflyKey],pidx[e.plantKey]]=1
    assert A.sum()==1565 and (C>=R).all()
    bf=np.array([meta[k]["family"] for k in bk])
    accepted=np.array([meta[k]["status"]=="ACCEPTED" for k in bk])
    oldmapped=accepted&(A.sum(1)>0)
    common=oldmapped&(bf=="Lycaenidae")&(R.sum(1)>0)
    assert common.sum()==46
    raw=jzip(work/"gesui_gbif_facets.zip","GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json")
    assert len(raw)==2808 and all(x["source_status"]=="ok" for x in raw)
    facets={(str(x["HYBAS_ID_L07"]),x["period"],x["group"]):x for x in raw}
    assert len(facets)==2808
    panel=OLD/"gesui_rebuilt_network.csv"
    assert digest(panel)=="544823b8b06b2d89a0e637a258a0b0cd578fe56c3b0347fb26115a4fc19a3870"
    df=pd.read_csv(panel,dtype={"h":str,"h6":str},keep_default_na=False)
    d0=df[df.period==T0].set_index("h");d1=df[df.period==T1].set_index("h")
    hs=sorted(set(d0.index)&set(d1.index));d0=d0.loc[hs];d1=d1.loc[hs]
    assert len(hs)==222 and (d0.h6==d1.h6).all()
    quality=(d0.minimum_species_coverage>=.8)&(d1.minimum_species_coverage>=.8)
    assert quality.sum()==173 and d0[quality].h6.nunique()==77
    hg=sorted(set(d0[quality].h6),key=lambda s:(hsh(s),s))
    hg+=sorted(set(d0.h6)-set(hg),key=lambda s:(hsh(s),s))
    folds=d0.h6.map({g:i%5 for i,g in enumerate(hg)}).to_numpy()
    return hs,d0,d1,quality.to_numpy(),folds,bk,pk,pf,bf,A,C,oldmapped,common,facets
def arrays(hs,period,fac,bk,pk):
    B=[];P=[];pc=[];nB=[];nP=[];ps=[];match=[]
    for h in hs:
        rs=[fac[h,period,g] for g in FAMS]
        r=fac[h,period,"Plantae"]
        assert all(x["source_status"]=="ok" for x in rs+[r])
        bc=collections.Counter()
        for x in rs:bc.update(x["source_speciesKey_counts"])
        p=r["source_speciesKey_counts"];nb=sum(bc.values());np_=sum(p.values())
        b=np.array([bc.get(k,0) for k in bk]);c=np.array([p.get(k,0) for k in pk])
        B.append(incl(nb,b,20));P.append(incl(np_,c,50));pc.append(c);nB.append(nb);nP.append(np_)
        ps.append(incl(np_,list(p.values()),50).sum());match.append(b.sum()/nb)
    return tuple(np.asarray(v) for v in (B,P,pc,nB,nP,ps,match))
def nullgraph(A,pc,groups,pf,pk,tr):
    labs=sorted(set(groups[tr]))
    z=np.array([(pc[tr[groups[tr]==h]]>0).any(0) for h in labs]).T.astype(float)
    norm=np.linalg.norm(z,axis=1);good=norm>0
    cl=np.full(len(pk),-1);cl[good]=KMeans(n_clusters=6,random_state=20261114,n_init=10).fit_predict(z[good]/norm[good,None])
    bins=np.digitize((pc[tr]>0).sum(0),[0,1,5,20,100,468],right=True)
    ss=np.array([pf[k]+"|"+str(b)+"|"+str(c) for k,b,c in zip(pk,bins,cl)])
    N=np.zeros_like(A)
    for s in sorted(set(ss)):
        ix=np.flatnonzero(ss==s);N[:,ix]=A[:,ix].mean(1)[:,None]
    np.testing.assert_allclose(A.sum(1),N.sum(1),atol=1e-10,rtol=0)
    return N
def model(x,y,xt,w):
    mu=x.mean(0);sd=x.std(0);sd[sd<1e-12]=1
    x=(x-mu)/sd;xt=(xt-mu)/sd
    model=Ridge(alpha=10,fit_intercept=True).fit(x,y,sample_weight=w)
    p=model.predict(xt)
    X=np.column_stack((np.ones(len(x)),x));Xt=np.column_stack((np.ones(len(xt)),xt))
    pen=np.eye(X.shape[1])*10;pen[0,0]=0
    co=np.linalg.solve(X.T@(w[:,None]*X)+pen,X.T@(w*y))
    np.testing.assert_allclose(p,Xt@co,atol=1e-8,rtol=0)
    return np.clip(p,0,1),float(co[-1]/sd[-1])
def train(work,out):
    hs,d0,d1,quality,fold,bk,pk,pf,bf,A,C,oldmapped,common,fac=load_old(work)
    B0,P0,pc0,nB0,nP0,ps0,_=arrays(hs,T0,fac,bk,pk)
    B1,P1,pc1,nB1,nP1,ps1,_=arrays(hs,T1,fac,bk,pk)
    groups=d0.h6.to_numpy();geo=["W","area","riverflow","crop","urban","naturalout","EQR_original_site_equal_basin_period_mean"]
    def base(B,P,nb,np_,ps,d,ix,mask):
        k=int(mask.sum());fam=B[:,bf==bf[mask][0]].sum(1)
        v=np.column_stack((d[geo].to_numpy(float),np.log1p(nb),np.log1p(np_),ps,P.sum(1),fam))
        return np.column_stack((B[ix][:,mask].reshape(-1),np.repeat(v[ix],k,axis=0),np.tile(np.eye(k),(len(ix),1))))
    rows=[];fits=[]
    panels=[("common_Lycaenidae",common,[("old",A),("Clarke",C)])]
    panels += [(f,oldmapped&(bf==f),[("old",A)]) for f in FAMS]
    for tag,mask,graphs in panels:
        k=int(mask.sum())
        assert k>=5
        for foldid in range(5):
            tr=np.flatnonzero(quality&(fold!=foldid));te=np.flatnonzero(fold==foldid)
            assert len(set(groups[tr]))>=15 and not set(groups[tr])&set(groups[te])
            if tag=="common_Lycaenidae":
                den=B0[:,oldmapped&(bf=="Lycaenidae")].sum(1);den1=B1[:,oldmapped&(bf=="Lycaenidae")].sum(1)
                good=(B0[:,mask].sum(1)>=.8*den)&(B1[:,mask].sum(1)>=.8*den1)&(den>0)&(den1>0)
                tr=tr[good[tr]];te=te[good[te]]
            assert len(set(groups[tr]))>=15
            x=base(B0,P0,nB0,nP0,ps0,d0,tr,mask)
            xt=base(B1,P1,nB1,nP1,ps1,d1,te,mask)
            y=B1[tr][:,mask].ravel()
            cnt=collections.Counter(groups[tr]);w=np.repeat([1/cnt[g]/k for g in groups[tr]],k);w=w/w.mean()
            pred0,_=model(x,y,xt,w)
            for name,G0 in graphs:
                G=G0[mask];N=nullgraph(G,pc0,groups,pf,pk,tr);deg=G.sum(1)
                assert np.all(deg>0)
                v0=(P0@G.T)/deg;v1=(P1@G.T)/deg
                u0=(P0@N.T)/deg;u1=(P1@N.T)/deg
                pn,cn=model(np.column_stack((x,u0[tr].ravel())),y,np.column_stack((xt,u1[te].ravel())),w)
                pg,cg=model(np.column_stack((x,v0[tr].ravel())),y,np.column_stack((xt,v1[te].ravel())),w)
                fits.append(dict(panel=tag,graph=name,fold=foldid,train_h7=len(tr),train_h6=len(set(groups[tr])),test_h7=len(te),species=k,host_coef=cg,null_coef=cn))
                for i,h in enumerate(te):
                    for q,key in enumerate(np.asarray(bk)[mask]):
                        j=i*k+q
                        rows.append(dict(panel=tag,graph=name,fold=foldid,h=hs[h],h6=groups[h],key=key,prior_quality=bool(quality[h]),p0=pred0[j],pnull=pn[j],phost=pg[j],persistence=B1[h,np.flatnonzero(mask)[q]]))
    p=out/"frozen_predictions.csv";pd.DataFrame(rows).to_csv(p,index=False)
    pd.DataFrame(fits).to_csv(out/"historical_fits.csv",index=False)
    (out/"prediction_receipt.json").write_text(json.dumps({"sha256":digest(p),"rows":len(rows),"historical_only":True,"test_period_species_counts_read":False,"known_butterfly_keys":bk},indent=2)+"\n")
    print("FROZEN BEFORE T2 SCORE",len(rows),digest(p))
def score(work,out):
    receipt=json.loads((work/"prediction_receipt.json").read_text());p=work/"frozen_predictions.csv"
    assert receipt["historical_only"] and not receipt["test_period_species_counts_read"] and digest(p)==receipt["sha256"]
    pred=pd.read_csv(p,dtype={"h":str,"h6":str,"key":str})
    fac={}
    for i in range(4):
        ff=list(work.rglob(f"third_full_{i}.json"));assert len(ff)==1
        for x in json.loads(ff[0].read_text()):
            k=(str(x["HYBAS_ID_L07"]),x["period"],x["group"]);assert k not in fac;fac[k]=x
    assert len(fac)==1332
    gatefile=list(work.rglob("third_period_support.json"));assert len(gatefile)==1
    gate={x["h"]:x for x in json.loads(gatefile[0].read_text())}
    y={};cover=[]
    for h,z in gate.items():
        if not z["rarefaction_support"]:continue
        bc=collections.Counter()
        for g in FAMS:bc.update(fac[h,T2,g]["source_speciesKey_counts"])
        n=sum(bc.values());assert n==z["butterfly_facetted_n"] and n>=20
        keys=receipt["known_butterfly_keys"]
        for k,v in zip(keys,incl(n,[bc.get(k,0) for k in keys],20)):y[h,k]=v
        cover.append(dict(h=h,h6=z["h6"],prior_key_record_share=sum(bc.get(k,0) for k in keys)/n))
    pred["y"]=[y.get((h,k),np.nan) for h,k in zip(pred.h,pred.key)]
    pred["quality"]=[gate[h]["three_period_quality"] for h in pred.h]
    pred=pred[pred.y.notna()].copy()
    for k in ("p0","pnull","phost","persistence"):pred["e_"+k]=(pred.y-pred[k])**2
    rows=[];losses=[]
    for sample,s in [("quality",pred[pred.quality]),("all",pred)]:
        for (panel,graph),g in s.groupby(["panel","graph"]):
            v=g.groupby(["h6","h"])[["e_p0","e_pnull","e_phost","e_persistence"]].mean().groupby("h6").mean()
            n=len(v);nH=g.h.nunique()
            if n<20 or nH<40:
                rows.append(dict(sample=sample,panel=panel,graph=graph,H7=nH,H6=n,status="STOP_SUPPORT"));continue
            m=v.to_numpy();mu=m.mean(0)
            assert np.all(mu>0)
            ids=np.random.default_rng(20261124).integers(0,n,(5000,n));b=m[ids].mean(1)
            for ref,ix in [("baseline",0),("structured_null",1),("persistence",3)]:
                q=100*(b[:,ix]-b[:,2])/b[:,ix]
                rows.append(dict(sample=sample,panel=panel,graph=graph,H7=nH,H6=n,status="FIT",reference=ref,gain=100*(mu[ix]-mu[2])/mu[ix],lo=np.quantile(q,.025),hi=np.quantile(q,.975),host_mse=mu[2],reference_mse=mu[ix]))
            for h,r in v.iterrows():losses.append(dict(sample=sample,panel=panel,graph=graph,h6=h,**r.to_dict()))
    pd.DataFrame(rows).to_csv(out/"forecast_results.csv",index=False)
    pd.DataFrame(losses).to_csv(out/"forecast_losses.csv",index=False)
    pd.DataFrame(cover).to_csv(out/"prior_key_coverage.csv",index=False)
    print(pd.DataFrame(rows).round(4).to_string(index=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("stage",choices=["train","score"]);p.add_argument("--work",type=pathlib.Path,required=True);p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    {"train":train,"score":score}[a.stage](a.work,a.out)
