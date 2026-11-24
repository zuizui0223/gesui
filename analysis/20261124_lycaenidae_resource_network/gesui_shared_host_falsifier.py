#!/usr/bin/env python3
"""Exploratory catalogue shared-host topology; strict training-only early reference, not observed competition or sewage effect."""
from pathlib import Path
import sys
P=Path(sys.argv[1] if len(sys.argv)>1 else '/mnt/data/gesui_sync_replay')
prior=P/'gesui_strict_reference_falsifier.py'
exec(compile(prior.read_text().split('rows=[];foldrows=[];')[0],str(prior),'exec'))
AL=A[L]; mapped=AL.sum(1)>0; AL=AL[mapped]; Lidx=np.flatnonzero(L)[mapped]; assert len(Lidx)==52
Dg=AL.sum(0).astype(int); ii,jj=np.triu_indices(len(AL),1); pairgraph=AL[ii]*AL[jj]
assert len(ii)==1326 and np.array_equal(pairgraph.sum(0),Dg*(Dg-1)/2)
def pairprob(ix):
 counts=B[ix][:,Lidx]; n=np.asarray(NB)[ix,None]; p=PB[ix][:,Lidx]
 z=p[:,ii]+p[:,jj]-incl(n,counts[:,ii]+counts[:,jj],20)
 assert z.min()>-1e-12 and z.max()<=1+1e-12
 z=np.maximum(z,0)
 k=np.arange(0,len(ix),31); j=np.arange(0,len(ii),89)
 c1=counts[k][:,ii[j]];c2=counts[k][:,jj[j]]
 expected=1-hypergeom.pmf(0,n[k],c1,20)-hypergeom.pmf(0,n[k],c2,20)+hypergeom.pmf(0,n[k],c1+c2,20)
 np.testing.assert_allclose(z[np.ix_(k,j)],expected,atol=1e-10,rtol=0)
 return z
Q0,Q1=pairprob(a),pairprob(b)
S0,S1=BA[:,mapped]@AL,BB[:,mapped]@AL
J0,J1=Q0@pairgraph,Q1@pairgraph
static=pd.DataFrame(dict(plantKey=pk,degree=Dg,family=[pf[k] for k in pk],name=[pn[k] for k in pk],potential_pairs=pairgraph.sum(0).astype(int)))
static.to_csv(P/'gesui_shared_host_degree.csv',index=False)
print('degree counts',static.degree.value_counts().sort_index().to_dict(),flush=True)
def groups(train):
 ix=a[train]; hs=sorted(set(d.h6.iloc[ix]));assert len(hs)>=15
 freq=(Pcount[ix]>0).sum(0); bins=np.digitize(freq,[0,1,5,20,100,468],right=True)
 z=np.array([(Pcount[ix[d.h6.iloc[ix].to_numpy()==h]]>0).any(0) for h in hs]).T.astype(float); norm=np.linalg.norm(z,axis=1);lab=np.full(len(pk),-1);good=norm>0
 lab[good]=KMeans(n_clusters=6,random_state=20261114,n_init=10).fit_predict(z[good]/norm[good,None])
 st=np.array([pf[k]+'|'+str(f)+'|'+str(g) for k,f,g in zip(pk,bins,lab)])
 deg=np.array([s+'|degree'+str(g) for s,g in zip(st,Dg)])
 g0=[np.flatnonzero(st==s) for s in sorted(set(st))];g1=[np.flatnonzero(deg==s) for s in sorted(set(deg))]
 assert all(len(set(Dg[j]))==1 for j in g1)
 return g0,g1

def avg(C,gs):
 N=np.empty_like(C)
 for j in gs:N[:,j]=C[:,j].mean(1)[:,None]
 np.testing.assert_allclose(N.sum(1),C.sum(1),atol=1e-9,rtol=0)
 return N

def decomp(C,g0,g1):
 n0,n1=avg(C,g0),avg(C,g1); total=C-n0; degree=n1-n0; identity=C-n1
 np.testing.assert_allclose(total,degree+identity,atol=1e-12,rtol=0)
 return total,degree,identity
rows=[];lossrows=[];contrasts=[];foldmeta=[]
for tag,mask in [('quality',np.asarray(quality)),('all',np.ones(len(q),bool))]:
 ids=np.flatnonzero(mask);cat=np.array([f'{c}::{m}' for c,m in zip(country,method)])
 for mode in ['early_fixed','symmetric']:
  losses=defaultdict(list);coef=defaultdict(list);err=0.
  for f in sorted(set(fold[ids])):
   tr=ids[fold[ids]!=f];te=ids[fold[ids]==f];assert not set(H6[tr])&set(H6[te]);g0,g1=groups(tr)
   M,_=reference(tr); np.testing.assert_allclose(S0-avg(S0,g0),BA@M,atol=1e-12,rtol=0)
   sa=S0 if mode=='early_fixed' else .5*(S0+S1);ja=J0 if mode=='early_fixed' else .5*(J0+J1)
   early_s=(PA*(S0-avg(S0,g0))).sum(1);early_j=(PA*(J0-avg(J0,g0))).sum(1)
   sd,dd,idn=decomp(sa,g0,g1);pd0,pdd,pid=decomp(ja,g0,g1);dp=PC-PA;d2=np.broadcast_to((Dg*(Dg-1)/2).astype(float)[None,:],PA.shape).copy()
   ys={ 'single_total':(dp*sd).sum(1),'single_degree':(dp*dd).sum(1),'single_identity':(dp*idn).sum(1), 'pair_total':(dp*pd0).sum(1),'pair_degree':(dp*pdd).sum(1),'pair_identity':(dp*pid).sum(1), 'pair_uniform_degree':(dp*(d2-avg(d2,g0))).sum(1), 'pair_observed':(dp*ja).sum(1),'pair_null':(dp*avg(ja,g0)).sum(1)}
   for stem in ['single','pair']:np.testing.assert_allclose(ys[stem+'_total'],ys[stem+'_degree']+ys[stem+'_identity'],atol=1e-12,rtol=0)
   np.testing.assert_allclose(ys['pair_total'],ys['pair_observed']-ys['pair_null'],atol=1e-12,rtol=0)
   X=np.column_stack((X0,early_s,early_j));X1=np.column_stack((X,E0))
   foldmeta.append(dict(sample=tag,consumer=mode,fold=f,train_H7=len(tr),test_H7=len(te),train_H6=len(set(H6[tr])),test_H6=len(set(H6[te])),strata=len(g0),degree_strata=len(g1),movable_shared_keys=sum(len(j) for j in g1 if len(j)>1 and Dg[j[0]]>=2)))
   for name,y in ys.items():
    p0,e0,_=fit(tr,te,y,X,cat);p1,e1,c=fit(tr,te,y,X1,cat);err=max(err,e0,e1);coef[name].append(c)
    for g in sorted(set(H6[te])):
     z=H6[te]==g;l0=float(np.mean((y[te][z]-p0[z])**2));l1=float(np.mean((y[te][z]-p1[z])**2));losses[name].append((g,l0,l1));lossrows.append(dict(sample=tag,consumer=mode,fold=f,endpoint=name,h6=g,mse0=l0,mse1=l1))
  boots={};gains={}
  for name,ls in losses.items():
   df=pd.DataFrame(ls,columns=['h6','mse0','mse1']).sort_values('h6');assert not df.h6.duplicated().any();m=df[['mse0','mse1']].to_numpy();idx=np.random.default_rng(20261124).integers(0,len(m),(5000,len(m)));v=m[idx].mean(1);bs=100*(v[:,0]-v[:,1])/v[:,0];lo,hi=np.quantile(bs,[.025,.975]);boots[name]=bs
   l0,l1=m.mean(0);gain=100*(l0-l1)/l0;gains[name]=gain;rows.append(dict(sample=tag,consumer=mode,endpoint=name,H7=len(ids),H6=len(m),gain=gain,lo=lo,hi=hi,coef_min=min(coef[name]),coef_max=max(coef[name]),independent_prediction_max_diff=err,mse0=l0,mse1=l1))
  for x,y in [('pair_total','single_total'),('pair_degree','pair_identity'),('single_degree','single_identity'),('pair_total','pair_uniform_degree')]:
   bs=boots[x]-boots[y];lo,hi=np.quantile(bs,[.025,.975]);b0,b1=np.quantile(bs,[.025/4,1-.025/4]);contrasts.append(dict(sample=tag,consumer=mode,A=x,B=y,diff=gains[x]-gains[y],lo=lo,hi=hi,bonf4_lo=b0,bonf4_hi=b1))
for name,v in [('results',rows),('contrasts',contrasts),('losses',lossrows),('folds',foldmeta)]:pd.DataFrame(v).to_csv(P/f'gesui_shared_host_{name}.csv',index=False)
print(pd.DataFrame(rows)[['sample','consumer','endpoint','gain','lo','hi','coef_min','coef_max']].round(4).to_string(index=False));print(pd.DataFrame(contrasts).round(4).to_string(index=False))
