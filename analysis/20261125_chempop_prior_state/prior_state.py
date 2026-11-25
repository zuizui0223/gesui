#!/usr/bin/env python3
"""Four equal windows of original matched English monitoring; source gate before new outcomes."""
import argparse,collections,hashlib,json,pathlib,sys,tempfile,zipfile
import numpy as np
import pandas as pd
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261125_chempop_resolution"))
import chemistry_order as prev
wq,old,f,src=prev.wq,prev.old,prev.f,prev.src
WINDOWS=((2011,2012),(2013,2014),(2015,2016),(2017,2018))
BACKGROUND=[c for c in f.COVARIATES if c!="year"]
NUTRIENTS=["chem_"+k for k in wq.NUTRIENTS]
MODELS=["P","PS","PB","PBE","PBC","PBCE","PSB","PSBE","PSBC","PSBCE"]
RESPONSES=["total","between_order","within_order"]
CONTRASTS=[("PSBCE","PSBC"),("PSBCE","PSBE"),("PSBCE","PBCE"),("PS","P"),("PB","P"),("PSB","PS"),("PBE","PB"),("PBC","PB"),("PBCE","PBC"),("PBCE","PBE"),("PSBE","PSB"),("PSBC","PSB"),("PSB","PB"),("PSBE","PBE"),("PSBC","PBC")]
SEED=20261125;BOOT=2000
def digest(b):return hashlib.sha256(b).hexdigest()
def roster_sha(ids):return digest(("\n".join(sorted(ids))+"\n").encode())
def win(y):return next((i for i,(a,b) in enumerate(WINDOWS) if a<=y<=b),None)
def finite(a):return bool(np.isfinite(np.asarray(a,float)).all())
def center(A,region):
    Z=A.copy()
    for p in sorted(set(region)):
        ix=region==p;Z[ix]-=Z[ix].mean(axis=0)
    return Z
def donor(Y,region):
    M=np.empty_like(Y)
    for p in sorted(set(region)):
        i=np.flatnonzero(region==p);assert len(i)>=2
        M[i]=(Y[i].sum(axis=0)-Y[i])/(len(i)-1)
    return M
def source_inputs(work,out):
    orig=prev.receipt(work/"chemistry.zip","18c764318c2143a40565b8d58be350d91c7b8e3b371c9705038401db8730ff43")
    prior=prev.receipt(work/"resolution.zip","9c04699d4fab6ba27b3e5ba99e6c17670b0c69bb912584dbe9a697391185edfc")
    with tempfile.TemporaryDirectory() as tmp:
        d=pathlib.Path(tmp)
        ret=src.retrieve(src.EIDC_URL,d/"original.zip",src.SOURCE_DOWNLOAD_LIMIT)
        scores,si=src.score_reference(d)
        assert si["SHA256"]==prior["score_sha"]==wq.EXPECTED_SCORE_SHA and len(scores)==92
        with zipfile.ZipFile(d/"original.zip") as z:
            names=z.namelist();assert len(names)==len(set(names))==326
            expected=prior["source_35_members"];assert len(expected)==35
            for k,v in expected.items():assert digest(z.read(k))==v,k
            raw,weights,orders,qa=f.construct_eligible(z,names,scores)
            assert len(raw)==7836 and qa["source_gate"]["eligible_site_n"]==707
            chem,_=wq.chemical_summaries(z,names,{(r["region"],r["site"],r["sample_date"]) for r in raw})
        joint,gate=wq.source_join_raw(raw,weights,chem)
        assert gate==orig["matched_source_gate"]==prior["cohort"]
        df,joined=wq.site_year_matched(joint,scores)
        assert joined==orig["matched_siteyear_source"]
        df0,A,_,families,s,_=f.site_year_response([v[0] for v in joint],[v[1] for v in joint],scores)
        assert df[["site_region_id","year"]].equals(df0[["site_region_id","year"]])
        assert (len(df),df.site_region_id.nunique(),df.region.nunique())==(2624,591,7)
        assert len(families)==85
        ref,oq=old.order_ref(families,s,orders);assert oq==prior["order_source"]
    df=df.reset_index(drop=True)
    df["site_region_id"]=df.site_region_id.astype(str);df["site"]=df.site.astype(str);df["region"]=df.region.astype(str)
    df["window"]=[win(int(y)) for y in df.year]
    assert not df.duplicated(["site_region_id","year"]).any()
    # Original sample method counts only, not numeric abundance values.
    meth=collections.defaultdict(collections.Counter)
    counts=collections.Counter()
    for r,_,_ in joint:
        k=(str(r["site_region_id"]),int(r["year"]))
        meth[k][str(r["method"])]+=1;counts[k]+=1
    assert len(meth)==len(df) and sum(counts.values())==5000
    ambiguous=set(df.groupby("site").region.nunique().loc[lambda x:x>1].index)
    rows=[];keep={}
    for sid,g in df.groupby("site_region_id",sort=True):
        assert g.region.nunique()==g.site.nunique()==1
        sample=g[g.window.notna()]
        yr=[sorted(int(x) for x in sample.loc[sample.window==t,"year"]) for t in range(4)]
        m=collections.Counter()
        for v in sample.itertuples():m.update(meth[sid,int(v.year)])
        eligible=all(yr) and str(g.site.iloc[0]) not in ambiguous and g.edf_pct.max()-g.edf_pct.min()<1e-10 and finite(sample[BACKGROUND+NUTRIENTS+["edf_log1p"]].to_numpy())
        rows.append(dict(site_region_id=sid,native_site=str(g.site.iloc[0]),region=str(g.region.iloc[0]),source_eligible=bool(eligible),native_id_ambiguous=str(g.site.iloc[0]) in ambiguous,edf_static=bool(g.edf_pct.max()-g.edf_pct.min()<1e-10),stable_abundance_method=len(m)==1,methods=";".join(sorted(m)),source_matched_samples=sum(m.values()),**{"years_"+str(t):";".join(map(str,yr[t])) for t in range(4)}))
    ro=pd.DataFrame(rows);ro.to_csv(out/"source_roster.csv",index=False)
    gates={}
    for name in ("main","stable_recording_method"):
        k=ro[ro.source_eligible].copy()
        if name!="main":k=k[k.stable_abundance_method]
        before=k.region.value_counts().to_dict()
        k=k[k.region.isin([p for p,n in before.items() if n>=10])]
        ids=sorted(k.site_region_id)
        idx=df.site_region_id.isin(ids)&df.window.notna()
        sites=k.region.value_counts().to_dict()
        g=dict(candidate_sites_before_region_filter=int(sum(before.values())),candidate_sites_by_region=before,sites=len(ids),regions=len(sites),site_counts=sites,siteyears=int(idx.sum()),window_siteyears=[int((idx&(df.window==t)).sum()) for t in range(4)],original_biological_analyses=int(sum(counts[(r.site_region_id,int(r.year))] for r in df[idx].itertuples())),roster_sha256=roster_sha(ids))
        g["pass"]=bool(g["sites"]>=100 and g["regions"]>=5)
        gates[name]=g;keep[name]=ids
        k.to_csv(out/(name+"_source_roster.csv"),index=False)
    receipt={"status":"SOURCE_GATE_COMPLETE","original_cohort":gate,"source_35_members":expected,"score_SHA256":si["SHA256"],"wrapper_SHA256":ret["SHA256"],"order_source":oq,"windows":WINDOWS,"cross_region_ambiguous_native_sites":sorted(ambiguous),"cohorts":gates,"new_outcomes_evaluated":False}
    (out/"source_receipt.json").write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    print("CHEMPOP REPEAT SOURCE "+json.dumps(receipt,allow_nan=False),flush=True)
    return df,A,s,ref,meth,keep,receipt
def build(df,A,s,ref,meth,ids):
    # Called only after source gate. Equal years within each calendar window.
    all_methods=sorted(set(x for m in meth.values() for x in m))
    scores=np.column_stack((s,ref,s-ref));Y=np.empty((len(ids),4,3));X=np.empty((len(ids),4,len(BACKGROUND)+1+len(NUTRIENTS)))
    methods=np.zeros((len(ids),4,len(all_methods)));regions=[]
    for i,sid in enumerate(ids):
        g=df[df.site_region_id==sid];regions.append(str(g.region.iloc[0]))
        for t in range(4):
            h=g[g.window==t];assert len(h)
            Y[i,t]=A[h.index.to_numpy()].mean(axis=0)@scores
            X[i,t]=h[BACKGROUND+["edf_log1p"]+NUTRIENTS].to_numpy(float).mean(axis=0)
            for r in h.itertuples():
                m=meth[sid,int(r.year)]
                for j,k in enumerate(all_methods):methods[i,t,j]+=m[k]/(sum(m.values())*len(h))
    np.testing.assert_allclose(Y[:,:,0],Y[:,:,1]+Y[:,:,2],atol=1e-12,rtol=0)
    assert finite(X) and finite(Y)
    np.testing.assert_allclose(methods.sum(axis=2),1,atol=1e-12,rtol=0)
    return Y,X,methods,np.array(regions,str),all_methods
def fit_ridge(A,D,R,w):
    rms=np.sqrt(np.sum(w[:,None]*A*A,axis=0)/w.sum());active=rms>1e-12
    scale=np.where(active,rms,1.)
    a=A/scale;d=D/scale;a[:,~active]=0;d[:,~active]=0
    H=a.T@(w[:,None]*a)+10*np.eye(a.shape[1])
    b=np.linalg.solve(H,a.T@(w[:,None]*R))
    a0=np.vstack([a*np.sqrt(w)[:,None],np.sqrt(10)*np.eye(a.shape[1])]);r0=np.vstack([R*np.sqrt(w)[:,None],np.zeros((a.shape[1],3))])
    b0=np.linalg.lstsq(a0,r0,rcond=None)[0]
    err=float(np.max(np.abs(d@b-d@b0)));assert err<1e-10
    return d@b,int(active.sum()),err
def evaluate(cohort,Y,X,methods,region,all_methods,ids,out):
    regions=sorted(set(region));P0=donor(Y[:,0],region);P2=donor(Y[:,2],region);S0=Y[:,0];S2=Y[:,2]
    Z0=center(X[:,0],region);Z2=center(X[:,2],region);M0=center(methods[:,0],region);M2=center(methods[:,2],region)
    pred={m:np.empty_like(Y[:,3]) for m in MODELS};fold=[]
    nb=len(BACKGROUND);E=[nb];C=list(range(nb+1,nb+1+len(NUTRIENTS)))
    for p in regions:
        tr=region!=p;te=~tr
        w=np.array([1/sum(region==q) for q in region[tr]])
        delta=S0[tr,0]-P0[tr,0];v=Y[tr,1,0]-P0[tr,0];den=np.dot(w,delta*delta)
        own=float(np.clip(np.dot(w,delta*v)/den,0,1)) if den>1e-20 else 0.
        # Only predictor-time training vocabulary. All retained columns centered within region.
        vocab=np.any(methods[tr,0]>0,axis=0)
        base0={"P":P0[tr],"PS":(1-own)*P0[tr]+own*S0[tr]}
        base2={"P":P2[te],"PS":(1-own)*P2[te]+own*S2[te]}
        for label in ("P","PS"):
            pred[label][te]=base2[label]
            R=Y[tr,1]-base0[label]
            for suffix,extra in [("B",[]),("BE",E),("BC",C),("BCE",E+C)]:
                cols=list(range(nb))+extra
                a=np.column_stack((Z0[tr][:,cols],M0[tr][:,vocab]));d=np.column_stack((Z2[te][:,cols],M2[te][:,vocab]))
                correction,nactive,err=fit_ridge(a,d,R,w)
                pred[label+suffix][te]=base2[label]+correction
                fold.append(dict(cohort=cohort,heldout_region=p,train_sites=int(tr.sum()),test_sites=int(te.sum()),train_regions=len(regions)-1,model=label+suffix,own_weight=own,active_columns=nactive,training_method_vocabulary=";".join(m for m,b in zip(all_methods,vocab) if b),solver_max_error=err))
    losses={};site=[]
    for m,z in pred.items():
        assert finite(z);np.testing.assert_allclose(z[:,0],z[:,1]+z[:,2],atol=1e-10,rtol=0)
        error=(Y[:,3]-z)**2
        for j,name in enumerate(RESPONSES):
            losses[name,m]=np.array([error[region==p,j].mean() for p in regions])
            site.extend(dict(cohort=cohort,site_region_id=k,region=p,response=name,model=m,squared_loss=float(e)) for k,p,e in zip(ids,region,error[:,j]))
    pd.DataFrame(site).to_csv(out/(cohort+"_site_losses.csv"),index=False)
    pd.DataFrame(fold).to_csv(out/(cohort+"_folds.csv"),index=False)
    rng=np.random.default_rng(SEED);ix=rng.integers(0,len(regions),(BOOT,len(regions)));loo=np.ones((len(regions),len(regions)),bool);np.fill_diagonal(loo,False)
    summaries=[];regional=[]
    for resp in RESPONSES:
        for i,p in enumerate(regions):
            regional.append(dict(cohort=cohort,region=p,response=resp,sites=int(sum(region==p)),**{m:float(losses[resp,m][i]) for m in MODELS}))
        for new,base in CONTRASTS:
            a=losses[resp,base];b=losses[resp,new];assert a.mean()>1e-20
            den=a[ix].mean(axis=1);valid=den>1e-20
            v=100*(a[ix].mean(axis=1)[valid]-b[ix].mean(axis=1)[valid])/den[valid]
            lv=np.array([100*(a[i].mean()-b[i].mean())/a[i].mean() for i in loo if a[i].mean()>1e-20])
            summaries.append(dict(cohort=cohort,response=resp,contrast=new+"_vs_"+base,gain=float(100*(a.mean()-b.mean())/a.mean()),lo=float(np.quantile(v,.025)),hi=float(np.quantile(v,.975)),valid_bootstrap_draws=int(valid.sum()),regions=len(regions),improved_regions=int(sum(a>b)),loo_min=float(lv.min()),loo_max=float(lv.max()),baseline_mse=float(a.mean()),added_mse=float(b.mean())))
    pd.DataFrame(regional).to_csv(out/(cohort+"_region_losses.csv"),index=False)
    return summaries
def main(work,out):
    out.mkdir(parents=True,exist_ok=True)
    df,A,s,ref,meth,keep,receipt=source_inputs(work,out)
    summaries=[]
    for cohort,ids in keep.items():
        if not receipt["cohorts"][cohort]["pass"]:continue
        Y,X,M,region,methods=build(df,A,s,ref,meth,ids)
        summaries+=evaluate(cohort,Y,X,M,region,methods,ids,out)
    pd.DataFrame(summaries).to_csv(out/"summary.csv",index=False)
    primary=next((v for v in summaries if v["cohort"]=="main" and v["response"]=="within_order" and v["contrast"]=="PSBCE_vs_PSBC"),None)
    result=dict(status="COMPLETE" if receipt["cohorts"]["main"]["pass"] else "SOURCE_STOP_MAIN",source=receipt["cohorts"],summary=summaries,primary=primary,primary_pass=bool(primary and primary["lo"]>0),posthoc=True,not_causal=True)
    (out/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("CHEMPOP REPEAT ECOLOGY "+json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--work",type=pathlib.Path,required=True);p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args();main(a.work,a.out)
