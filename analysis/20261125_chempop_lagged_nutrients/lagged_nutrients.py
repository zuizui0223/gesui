#!/usr/bin/env python3
"""Lagged nutrient-summary increment on exact PR170 repeated-site support."""
import argparse, json, pathlib, sys, zipfile
import numpy as np
import pandas as pd
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261125_chempop_prior_state"))
import prior_state as p
SEED=20261125
BOOT=2000
CONTRASTS=[("PSBCEH","PSBCE"),("PBCEH","PBCE"),("PSBCED","PSBCE"),("PSBCEDH","PSBCED"),("PSBCE","PBCE")]
def replay(work,out):
    out.mkdir(parents=True,exist_ok=True)
    zpath=work/"repeat.zip"
    assert p.digest(zpath.read_bytes())=="350182684a24a9ee095fcf5020a2ffa755da4d8349c7ea635ed9468db21f265a"
    with zipfile.ZipFile(zpath) as z:
        saved=json.loads(z.read("results.json"))
        assert saved["status"]=="COMPLETE"
        assert saved["source"]["main"]["roster_sha256"]=="5dd31f4839292160872bc40f3aaa55fdcd496924463deac803d3514261ae0fd3"
        assert not saved["source"]["stable_recording_method"]["pass"]
        df,A,s,ref,meth,keep,receipt=p.source_inputs(work,out)
        assert receipt["cohorts"]==saved["source"]
        ids=keep["main"]
        assert len(ids)==131
        Y,X,M,region,methods=p.build(df,A,s,ref,meth,ids)
        rep=out/"previous_replay";rep.mkdir(exist_ok=True)
        summ=p.evaluate("main",Y,X,M,region,methods,ids,rep)
        keys=["cohort","site_region_id","region","response","model"]
        import io
        a=pd.read_csv(io.BytesIO(z.read("main_site_losses.csv")),dtype={"site_region_id":str}).sort_values(keys).reset_index(drop=True)
        b=pd.read_csv(rep/"main_site_losses.csv",dtype={"site_region_id":str}).sort_values(keys).reset_index(drop=True)
        assert len(a)==len(b)==3930 and a[keys].equals(b[keys])
        err=float(np.max(np.abs(a.squared_loss-b.squared_loss)));assert err<1e-12
        keys2=["cohort","response","contrast"]
        a=pd.DataFrame(saved["summary"]).sort_values(keys2).reset_index(drop=True)
        b=pd.DataFrame(summ).sort_values(keys2).reset_index(drop=True)
        assert len(a)==len(b)==45 and a[keys2].equals(b[keys2])
        fields=[k for k in a if k not in keys2]
        np.testing.assert_allclose(a[fields].to_numpy(float),b[fields].to_numpy(float),atol=1e-9,rtol=0)
    return Y,X,M,region,methods,ids,dict(site_loss_replay_max_error=err,site_loss_rows=3930,contrast_rows=45,source=receipt["cohorts"])
def rms(A,w):
    return np.sqrt(np.sum(w[:,None]*A*A,axis=0)/w.sum())
def project(A,D,H,J,w):
    """Residualize additional block in training, using no outcome response."""
    scale=rms(A,w);active=scale>1e-12
    a=A[:,active]/scale[active];d=D[:,active]/scale[active]
    sw=np.sqrt(w)
    beta,_,rank,_=np.linalg.lstsq(a*sw[:,None],H*sw[:,None],rcond=1e-10)
    h=H-a@beta;j=J-d@beta
    err=float(np.max(np.abs(a.T@(w[:,None]*h)))) if a.shape[1] else 0.
    # rank-aware least squares, check retained orthogonality numerically
    assert err<1e-7,(rank,err)
    before=rms(H,w);after=rms(h,w)
    frac=np.divide(after**2,before**2,out=np.zeros_like(after),where=before>1e-12)
    return h,j,dict(rank=int(rank),columns=int(a.shape[1]),raw_rms=before.tolist(),residual_rms=after.tolist(),residual_fraction=frac.tolist(),max_orthogonality_error=err)
def ridge(A,D,R,w):
    # Evaluate identical fixed fit on training and heldout predictors.
    z,n,e=p.fit_ridge(A,np.vstack([A,D]),R,w)
    return z[:len(A)],z[len(A):],dict(active_columns=n,solver_max_error=e)
def evaluate(Y,X,M,region,methods,ids,out):
    names=["PBCE","PSBCE","PBCEH","PSBCEH","PSBCED","PSBCEDH"]
    pred={k:np.empty_like(Y[:,3]) for k in names}
    regions=sorted(set(region));diagnostics=[]
    donors=[p.donor(Y[:,t],region) for t in (1,2)]
    recent=[p.center(X[:,t],region) for t in (1,2)]
    meth=[p.center(M[:,t],region) for t in (1,2)]
    nchem=len(p.NUTRIENTS)
    H=[p.center(X[:,t,-nchem:],region) for t in (0,1)]
    older=[p.center(Y[:,t,1:],region) for t in (0,1)]
    own=[p.center(Y[:,t,1:],region) for t in (1,2)]
    pp=[p.center(d[:,1:],region) for d in donors]
    for reg in regions:
        tr=region!=reg;te=~tr
        w=np.array([1/sum(region==r) for r in region[tr]])
        delta=Y[tr,1,0]-donors[0][tr,0];v=Y[tr,2,0]-donors[0][tr,0]
        den=np.dot(w,delta*delta);alpha=float(np.clip(np.dot(w,delta*v)/den,0,1)) if den>1e-20 else 0.
        vocab=np.any(M[tr,1]>0,axis=0)
        a=np.column_stack((recent[0][tr],meth[0][tr][:,vocab]))
        d=np.column_stack((recent[1][te],meth[1][te][:,vocab]))
        for label in ("P","PS"):
            b0=donors[0][tr].copy();b1=donors[1][te].copy()
            q0=np.column_stack((a,pp[0][tr]));q1=np.column_stack((d,pp[1][te]))
            if label=="PS":
                b0=(1-alpha)*b0+alpha*Y[tr,1]
                b1=(1-alpha)*b1+alpha*Y[te,2]
                q0=np.column_stack((q0,own[0][tr]));q1=np.column_stack((q1,own[1][te]))
            c0,c1,diag=ridge(a,d,Y[tr,2]-b0,w)
            base=label+"BCE"
            train=b0+c0;test=b1+c1
            pred[base][te]=test
            diagnostics.append(dict(region=reg,model=base,own_weight=alpha if label=="PS" else 0,train_sites=int(tr.sum()),test_sites=int(te.sum()),**diag))
            h,j,diag0=project(q0,q1,H[0][tr],H[1][te],w)
            g0,g1,diag1=ridge(h,j,Y[tr,2]-train,w)
            pred[base+"H"][te]=test+g1
            diagnostics.append(dict(region=reg,model=base+"H",**diag0,**diag1))
            if label=="PS":
                old0,old1,diag0=project(q0,q1,older[0][tr],older[1][te],w)
                g0,g1,diag1=ridge(old0,old1,Y[tr,2]-train,w)
                dt=train+g0;dv=test+g1
                pred[base+"D"][te]=dv
                diagnostics.append(dict(region=reg,model=base+"D",**diag0,**diag1))
                h,j,diag0=project(np.column_stack((q0,older[0][tr])),np.column_stack((q1,older[1][te])),H[0][tr],H[1][te],w)
                _,g1,diag1=ridge(h,j,Y[tr,2]-dt,w)
                pred[base+"DH"][te]=dv+g1
                diagnostics.append(dict(region=reg,model=base+"DH",**diag0,**diag1))
    site=[];losses={};region_rows=[]
    for name,z in pred.items():
        assert np.isfinite(z).all()
        np.testing.assert_allclose(z[:,0],z[:,1]+z[:,2],atol=1e-10,rtol=0)
        loss=(Y[:,3]-z)**2
        for i,r in enumerate(p.RESPONSES):
            losses[r,name]=np.asarray([loss[region==reg,i].mean() for reg in regions])
            site.extend(dict(site_region_id=k,region=reg,response=r,model=name,squared_loss=float(e)) for k,reg,e in zip(ids,region,loss[:,i]))
    for i,r in enumerate(p.RESPONSES):
        for j,reg in enumerate(regions):
            region_rows.append(dict(response=r,region=reg,sites=int(sum(region==reg)),**{k:float(losses[r,k][j]) for k in names}))
    rng=np.random.default_rng(SEED);ix=rng.integers(0,len(regions),(BOOT,len(regions)));loo=~np.eye(len(regions),dtype=bool)
    summary=[]
    for r in p.RESPONSES:
        for new,base in CONTRASTS:
            a=losses[r,base];b=losses[r,new]
            assert a.mean()>1e-20
            vals=100*(a[ix].mean(axis=1)-b[ix].mean(axis=1))/a[ix].mean(axis=1)
            lv=[100*(a[q].mean()-b[q].mean())/a[q].mean() for q in loo]
            summary.append(dict(response=r,contrast=new+"_vs_"+base,gain=float(100*(a.mean()-b.mean())/a.mean()),lo=float(np.quantile(vals,.025)),hi=float(np.quantile(vals,.975)),improved_regions=int(sum(a>b)),regions=len(regions),loo_min=float(min(lv)),loo_max=float(max(lv)),base_mse=float(a.mean()),new_mse=float(b.mean())))
    pd.DataFrame(site).to_csv(out/"site_losses.csv",index=False)
    pd.DataFrame(region_rows).to_csv(out/"region_losses.csv",index=False)
    pd.DataFrame(summary).to_csv(out/"summary.csv",index=False)
    (out/"fold_diagnostics.json").write_text(json.dumps(diagnostics,indent=2,allow_nan=False)+"\n")
    return summary,diagnostics
def main(work,out):
    Y,X,M,region,methods,ids,qa=replay(work,out)
    summary,diagnostics=evaluate(Y,X,M,region,methods,ids,out)
    primary=next(x for x in summary if x["response"]=="total" and x["contrast"]=="PSBCEH_vs_PSBCE")
    result=dict(status="COMPLETE",qa=qa,primary=primary,primary_pass=primary["lo"]>0,summary=summary,fold_diagnostics=diagnostics,not_causal=True,posthoc=True)
    (out/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("CHEMPOP LAGGED NUTRIENTS "+json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--work",type=pathlib.Path,required=True);a.add_argument("--out",type=pathlib.Path,required=True);v=a.parse_args();main(v.work,v.out)
