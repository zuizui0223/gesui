#!/usr/bin/env python3
"""Posthoc strict fold-specific early reference; family mass versus identity. No actual resource biomass, local feeding, or sewage causality."""
from pathlib import Path
import sys
P=Path(sys.argv[1])
prior=P/'gesui_strict_reference_falsifier.py'
exec(compile(prior.read_text().split('rows=[];foldrows=[];')[0],str(prior),'exec'))
fam=np.array([pf[k] for k in pk]);families=sorted(set(fam)); assert len(families)==71

def pieces(C,groups):
    mass=np.zeros(len(q));ident=mass.copy()
    for g in groups:
        x,y,c=PA[:,g],PC[:,g],C[:,g];m0,m1=x.sum(1),y.sum(1);t0,t1=(x*c).sum(1),(y*c).sum(1)
        u0=np.divide(t0,m0,out=np.zeros_like(t0),where=m0>1e-14);u1=np.divide(t1,m1,out=np.zeros_like(t1),where=m1>1e-14)
        u0=np.where(m0>1e-14,u0,u1);u1=np.where(m1>1e-14,u1,u0)
        mass+=(m1-m0)*.5*(u0+u1);ident+=(u1-u0)*.5*(m0+m1)
    return mass,ident

def endpoints(M,mode):
    for f in families: np.testing.assert_allclose(M[:,fam==f].sum(1),0,atol=1e-10,rtol=0)
    weights=BA if mode=='early_fixed' else .5*(BA+BB)
    C=weights@M;dp=PC-PA;allp=(C*dp).sum(1)
    mass,identity=pieces(C,[fam==f for f in families]); fm,fi=pieces(C,[fab]);nm,ni=pieces(C,[fam==f for f in families if f!='Fabaceae'])
    qt,ic=pieces(C,[np.ones(len(pk),bool)])
    fabpart=(C[:,fab]*dp[:,fab]).sum(1);nonpart=(C[:,~fab]*dp[:,~fab]).sum(1)
    obs=(weights@A[L]*dp).sum(1);null=(weights@(A[L]-M)*dp).sum(1)
    for v in [mass+identity,fm+fi+nm+ni,qt+ic,fabpart+nonpart,obs-null]:np.testing.assert_allclose(v,allp,atol=1e-12,rtol=0)
    return {'plant_delta':allp,'all_mass':qt,'all_identity':ic,'family_mass':mass,'within_family_identity':identity,'fab_delta':fabpart,'nonfab_delta':nonpart,'fab_mass':fm,'fab_identity':fi,'nonfab_mass':nm,'nonfab_identity':ni,'observed_graph_delta':obs,'null_graph_delta':null}
rows=[];contrasts=[];lossrows=[]
for tag,mask in [('quality',np.asarray(quality)),('all',np.ones(len(q),bool))]:
 ids=np.flatnonzero(mask)
 for mode in ['symmetric','early_fixed']:
    losses=defaultdict(list);coefs=defaultdict(list);err=0
    for f in sorted(set(fold[ids])):
        tr=ids[fold[ids]!=f];te=ids[fold[ids]==f];assert not set(H6[tr])&set(H6[te]);M,_=reference(tr);ys=endpoints(M,mode)
        early=(BA@M*PA).sum(1);X=np.column_stack((X0,early));X1=np.column_stack((X,E0));cat=np.array([f'{c}::{m}' for c,m in zip(country,method)])
        for name,y in ys.items():
            p0,e0,_=fit(tr,te,y,X,cat);p1,e1,c=fit(tr,te,y,X1,cat);err=max(err,e0,e1);coefs[name].append(c)
            for g in sorted(set(H6[te])):
                z=H6[te]==g;l0=float(np.mean((y[te][z]-p0[z])**2));l1=float(np.mean((y[te][z]-p1[z])**2));losses[name].append((g,l0,l1));lossrows.append(dict(sample=tag,consumer=mode,endpoint=name,h6=g,mse0=l0,mse1=l1))
    boot={};gains={}
    for name,ls in losses.items():
        df=pd.DataFrame(ls,columns=['h6','mse0','mse1']).sort_values('h6');assert not df.h6.duplicated().any();m=df[['mse0','mse1']].to_numpy();l0,l1=m.mean(0)
        idx=np.random.default_rng(20261124).integers(0,len(m),(5000,len(m)));v=m[idx].mean(1);bs=100*(v[:,0]-v[:,1])/v[:,0];boot[name]=bs;lo,hi=np.quantile(bs,[.025,.975]);gain=100*(l0-l1)/l0;gains[name]=gain
        rows.append(dict(sample=tag,consumer=mode,endpoint=name,H7=len(ids),H6=len(m),gain=gain,lo=lo,hi=hi,coef_min=min(coefs[name]),coef_max=max(coefs[name]),independent_prediction_max_diff=err))
    for x,y in [('all_identity','all_mass'),('within_family_identity','family_mass'),('fab_delta','nonfab_delta'),('plant_delta','null_graph_delta')]:
        bs=boot[x]-boot[y];lo,hi=np.quantile(bs,[.025,.975]);b0,b1=np.quantile(bs,[.025/4,1-.025/4]);contrasts.append(dict(sample=tag,consumer=mode,A=x,B=y,diff=gains[x]-gains[y],lo=lo,hi=hi,bonf4_lo=b0,bonf4_hi=b1))
for name,v in [('results',rows),('contrasts',contrasts),('losses',lossrows)]:pd.DataFrame(v).to_csv(P/f'gesui_family_strict_{name}.csv',index=False)
print(pd.DataFrame(rows).round(4).to_string(index=False));print(pd.DataFrame(contrasts).round(4).to_string(index=False))
