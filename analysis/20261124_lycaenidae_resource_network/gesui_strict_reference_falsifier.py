#!/usr/bin/env python3
"""Posthoc training-only early geographic null and early-consumer-fixed plant term.
No new field data. Fold-specific reference defines train/test target; no causal claims."""
from pathlib import Path
import sys
P=Path(sys.argv[1] if len(sys.argv)>1 else '/mnt/data/gesui_sync_replay')
prior=P/'gesui_fab_resource_consumer.py'
exec(compile(prior.read_text().split('base=q[K].to_numpy(float);ee=q[E].to_numpy(float)')[0],str(prior),'exec'))
from collections import defaultdict
BA,BB=PB[a][:,L],PB[b][:,L];PA,PC=PP[a],PP[b];E0=q[E].to_numpy(float);X0=q[K].to_numpy(float);H6=q.h6.to_numpy();country=q.country_holdout.to_numpy();method=q.method.replace('','__UNKNOWN__').to_numpy()
def reference(train):
 ix=a[train]; hs=sorted(set(d.h6.iloc[ix]));assert len(hs)>=15
 freq=(Pcount[ix]>0).sum(0);bins=np.digitize(freq,[0,1,5,20,100,468],right=True)
 z=np.array([(Pcount[ix[d.h6.iloc[ix].to_numpy()==h]]>0).any(0) for h in hs]).T.astype(float);norm=np.linalg.norm(z,axis=1);lab=np.full(len(pk),-1);good=norm>0
 lab[good]=KMeans(n_clusters=6,random_state=20261114,n_init=10).fit_predict(z[good]/norm[good,None])
 strata=np.array([pf[k]+'|'+str(f)+'|'+str(g) for k,f,g in zip(pk,bins,lab)])
 # Exact mean of the within-stratum permutations, not 499-draw approximation.
 N0=np.zeros_like(A[L]);AL=A[L]
 for s in sorted(set(strata)):
  j=np.flatnonzero(strata==s);N0[:,j]=AL[:,j].mean(1)[:,None]
 np.testing.assert_allclose(N0.sum(1),AL.sum(1),atol=1e-10,rtol=0)
 return AL-N0,dict(ref_H7=len(train),ref_H6=len(hs),nonzero_keys=int(good.sum()),strata=len(set(strata)))
def components(M):
 u=BA@M;v=BB@M;dp=PC-PA;dc=v-u;early=np.sum(u*PA,1);late=np.sum(v*PC,1)
 p=np.sum(u*dp,1);c=np.sum(dc*PA,1);i=np.sum(dc*dp,1)
 np.testing.assert_allclose(late-early,p+c+i,atol=1e-12,rtol=0)
 return early,{'plant_symmetric':p+.5*i,'consumer_symmetric':c+.5*i,'plant_early_consumer_fixed':p,'consumer_early_plant_fixed':c,'interaction':i,'total_delta':late-early}
def fit(tr,te,y,x,cat):
 cats=sorted(set(cat[tr]));cc=np.column_stack([cat==k for k in cats]).astype(float);ctr=collections.Counter(H6[tr]);w=np.array([1/ctr[g] for g in H6[tr]]);w/=w.mean()
 mu=x[tr].mean(0);sd=x[tr].std(0);sd[sd<1e-12]=1;X=np.column_stack((np.ones(len(x)),(x-mu)/sd,cc));pen=np.eye(X.shape[1])*10;pen[0,0]=0
 c=np.linalg.solve(X[tr].T@(w[:,None]*X[tr])+pen,X[tr].T@(w*y[tr]));p=X[te]@c;other=Ridge(alpha=10,fit_intercept=True).fit(X[tr,1:],y[tr],sample_weight=w).predict(X[te,1:]);err=float(np.max(np.abs(other-p)));assert err<1e-8
 return p,err,float(c[x.shape[1]]/sd[-1])
rows=[];foldrows=[];lossrows=[];contrastrows=[];M,gm=reference(np.arange(len(q)));ge,gy=components(M)
for tag,mask in [('quality',np.asarray(quality)),('all',np.ones(len(q),bool))]:
 ids=np.flatnonzero(mask)
 for mode in ['H6','country']:
  splits=[]
  if mode=='H6':
   for k in sorted(set(fold[ids])):splits.append((str(k),ids[fold[ids]!=k],ids[fold[ids]==k]))
   cat=np.array([f'{c}::{m}' for c,m in zip(country,method)])
  else:
   for k in sorted(set(country[ids])):
    te=ids[country[ids]==k];tr=ids[(country[ids]!=k)&~np.isin(H6[ids],H6[te])]
    if len(te)>=10 and len(set(H6[te]))>=5 and len(tr)>=30 and len(set(H6[tr]))>=15:splits.append((k,tr,te))
   cat=method
  for refname in ['early_all_paired_exact','training_early_exact']:
   losses=defaultdict(list);maxerr=defaultdict(float);coef=defaultdict(list)
   for k,tr,te in splits:
    assert not set(H6[tr])&set(H6[te])
    if refname=='training_early_exact':M,meta=reference(tr);y0,ys=components(M)
    else:y0,ys=ge,gy;meta=gm
    X=np.column_stack((X0,y0));X1=np.column_stack((X,E0));foldrows.append(dict(sample=tag,mode=mode,reference=refname,fold=k,train_H7=len(tr),test_H7=len(te),train_H6=len(set(H6[tr])),test_H6=len(set(H6[te])),**meta))
    for name,y in ys.items():
     p0,e0,_=fit(tr,te,y,X,cat);p1,e1,c=fit(tr,te,y,X1,cat);maxerr[name]=max(maxerr[name],e0,e1);coef[name].append(c)
     for g in sorted(set(H6[te])):
      z=H6[te]==g;v0=float(np.mean((y[te][z]-p0[z])**2));v1=float(np.mean((y[te][z]-p1[z])**2));losses[name].append((k,g,v0,v1));lossrows.append(dict(sample=tag,mode=mode,reference=refname,fold=k,h6=g,endpoint=name,mse0=v0,mse1=v1))
   boots={}
   for name,ls in losses.items():
    df=pd.DataFrame(ls,columns=['fold','h6','mse0','mse1'])
    if mode=='H6':
     assert not df.h6.duplicated().any();m=df[['mse0','mse1']].to_numpy();idx=np.random.default_rng(20261124).integers(0,len(m),(5000,len(m)));v=m[idx].mean(1);bs=100*(v[:,0]-v[:,1])/v[:,0];lo,hi=np.quantile(bs,[.025,.975]);boots[name]=bs
    else:
     m=df.groupby('fold')[['mse0','mse1']].mean().to_numpy();lo=hi=np.nan
     for k,g in df.groupby('fold'):
      l0,l1=g[['mse0','mse1']].mean();rows.append(dict(sample=tag,mode=mode,reference=refname,endpoint=name,fold=k,H6=g.h6.nunique(),gain=100*(l0-l1)/l0,mse0=l0,mse1=l1,lo=np.nan,hi=np.nan))
    l0,l1=m.mean(0);rows.append(dict(sample=tag,mode=mode,reference=refname,endpoint=name,fold='ALL',H6=df.h6.nunique(),gain=100*(l0-l1)/l0,mse0=l0,mse1=l1,lo=lo,hi=hi,coef_min=min(coef[name]),coef_max=max(coef[name]),independent_prediction_max_diff=maxerr[name]))
   if mode=='H6':
    for x,y in [('plant_symmetric','consumer_symmetric'),('plant_early_consumer_fixed','consumer_early_plant_fixed'),('plant_early_consumer_fixed','interaction')]:
     diff=boots[x]-boots[y];lo,hi=np.quantile(diff,[.025,.975]);lo3,hi3=np.quantile(diff,[.025/3,1-.025/3]);contrastrows.append(dict(sample=tag,reference=refname,A=x,B=y,lo=lo,hi=hi,bonf3_lo=lo3,bonf3_hi=hi3))
for name,records in [('results',rows),('folds',foldrows),('losses',lossrows),('contrasts',contrastrows)]:pd.DataFrame(records).to_csv(P/f'gesui_strict_reference_{name}.csv',index=False)
r=pd.DataFrame(rows);print(r[(r.reference=='training_early_exact') & ((r.fold=='ALL') | ((r.mode=='country') & (r.endpoint=='plant_early_consumer_fixed')))][['sample','mode','endpoint','fold','H6','gain','lo','hi','coef_min','coef_max']].round(3).to_string(index=False));print(pd.DataFrame(contrastrows).round(3).to_string(index=False));
