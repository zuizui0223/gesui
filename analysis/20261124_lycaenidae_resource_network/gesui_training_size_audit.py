#!/usr/bin/env python3
"""Fixed-target training restriction diagnostic; not catalogue replication.
Hold identical early-only full-training reference, retained test H7, fixed common40 graph.
Random training subsets match the retained training H7 count AND exact H6 size multiset.
Only regression training samples change. Draw distributions are descriptive, not null p-values.
"""
from pathlib import Path
import sys
import numpy as np,pandas as pd
P=Path(sys.argv[1]);prior=P/'gesui_catalogue_attrition_audit.py'
exec(compile(prior.read_text().split('rows=[];losses=[];cat=')[0],str(prior),'exec'))
cat=np.array([f'{c}::{m}' for c,m in zip(country,method)])
R=200;rng=np.random.default_rng(20261124)
rows=[];details=[];design=[];out=[]
def random_training(pool,sizes):
 groups={g:pool[H6[pool]==g] for g in sorted(set(H6[pool]))}
 res=[]
 for n in sorted(sizes,reverse=True):
  elig=[g for g,v in groups.items() if len(v)>=n]
  assert elig, 'cluster-size match impossible'
  g=rng.choice(elig);res.extend(rng.choice(groups.pop(g),n,replace=False).tolist())
 res=np.array(sorted(res));assert len(res)==sum(sizes) and sorted(collections.Counter(H6[res]).values())==sorted(sizes)
 return res
for sample,sm in [('quality',np.asarray(quality)),('all',np.ones(len(q),bool))]:
 ids=np.flatnonzero(sm);store=collections.defaultdict(list)
 for k in sorted(set(fold[ids])):
  te=ids[(fold[ids]==k)&cov[ids]];tr=ids[fold[ids]!=k];restricted=tr[cov[tr]]
  assert not set(H6[tr])&set(H6[te])
  sizes=list(collections.Counter(H6[restricted]).values())
  st=shared_strata(tr);y0,ys=decompose(A[common],PB[a][:,common],PB[b][:,common],st);y=ys['plant_symmetric'];X=np.column_stack((X0,y0));X1=np.column_stack((X,E0))
  design.append(dict(sample=sample,fold=k,train_H7=len(tr),train_H6=len(set(H6[tr])),restricted_H7=len(restricted),restricted_H6=len(sizes),test_H7=len(te),test_H6=len(set(H6[te])),restricted_H6_sizes=';'.join(map(str,sorted(sizes)))))
  for name,i,tt in [('full',-1,tr),('coverage_selected',-1,restricted)]+[('matched_random',i,random_training(tr,sizes)) for i in range(R)]:
   assert not set(H6[tt])&set(H6[te]);p0,e0,_=fit(tt,te,y,X,cat);p1,e1,_=fit(tt,te,y,X1,cat)
   for g in sorted(set(H6[te])):
    z=H6[te]==g;l0=float(np.mean((y[te][z]-p0[z])**2));l1=float(np.mean((y[te][z]-p1[z])**2));store[name,i].append((g,l0,l1));details.append(dict(sample=sample,training=name,draw=i,fold=k,h6=g,mse0=l0,mse1=l1,solver_diff=max(e0,e1)))
 for (name,i),ls in store.items():
  df=pd.DataFrame(ls,columns=['h6','mse0','mse1']).sort_values('h6');assert not df.h6.duplicated().any() and len(df)==len(set(H6[ids[cov[ids]]]))
  m=df[['mse0','mse1']].to_numpy();v=m.mean(0);g=100*(v[0]-v[1])/v[0];out.append(dict(sample=sample,training=name,draw=i,H6=len(df),gain=g,mse0=v[0],mse1=v[1]))
  if name!='matched_random':
   ix=np.random.default_rng(20261124).integers(0,len(m),(5000,len(m)));z=m[ix].mean(1);lo,hi=np.quantile(100*(z[:,0]-z[:,1])/z[:,0],[.025,.975]);rows.append(dict(sample=sample,training=name,H7=int((sm&cov).sum()),H6=len(df),gain=g,lo=lo,hi=hi))
res=pd.DataFrame(out)
for sample,d0 in res.groupby('sample'):
 x=d0[d0.training=='matched_random'].gain.to_numpy();s=d0[d0.training=='coverage_selected'].gain.iloc[0];f=d0[d0.training=='full'].gain.iloc[0];lo,mid,hi=np.quantile(x,[.025,.5,.975]);rows.append(dict(sample=sample,training='matched_random_distribution',gain=mid,lo=lo,hi=hi,selected_percentile=float(np.mean(x<=s)*100),full_gain=f,coverage_selected_gain=s,n_draws=len(x)))
for n,d0 in [('results',rows),('draws',out),('folds',design),('losses',details)]:pd.DataFrame(d0).to_csv(P/f'gesui_training_size_{n}.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False));print('solver max',max(x['solver_diff'] for x in details));print(pd.DataFrame(design).to_string(index=False))
