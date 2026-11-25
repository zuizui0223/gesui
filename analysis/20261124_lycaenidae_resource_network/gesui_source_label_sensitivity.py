#!/usr/bin/env python3
"""Source-label sensitivity of ORIGINAL Ferrer-Paris links, not independent catalogue.
Fixed common accepted butterfly support. HOSTS DB / BHL labels do not prove field evidence.
"""
from pathlib import Path
import sys,hashlib,json
import numpy as np,pandas as pd
P=Path(sys.argv[1])
prior=P/'gesui_strict_reference_falsifier.py'
exec(compile(prior.read_text().split('rows=[];foldrows=[];')[0],str(prior),'exec'))
# Fixed early-only geography strata, independent of graph edges.
def shared_strata(train):
 ix=a[train];hs=sorted(set(d.h6.iloc[ix]));assert len(hs)>=15
 freq=(Pcount[ix]>0).sum(0);bins=np.digitize(freq,[0,1,5,20,100,468],right=True)
 z=np.array([(Pcount[ix[d.h6.iloc[ix].to_numpy()==h]]>0).any(0) for h in hs]).T.astype(float)
 norm=np.linalg.norm(z,axis=1);lab=np.full(len(pk),-1);good=norm>0
 lab[good]=KMeans(n_clusters=6,random_state=20261114,n_init=10).fit_predict(z[good]/norm[good,None])
 return [np.flatnonzero(np.array([pf[k]+'|'+str(f)+'|'+str(g) for k,f,g in zip(pk,bins,lab)])==s) for s in sorted(set(pf[k]+'|'+str(f)+'|'+str(g) for k,f,g in zip(pk,bins,lab)))]
def exact_excess(graph,strata):
 null=np.zeros_like(graph)
 for j in strata:null[:,j]=graph[:,j].mean(1)[:,None]
 np.testing.assert_allclose(graph.sum(1),null.sum(1),atol=1e-10,rtol=0)
 return graph-null

def decompose(graph,b0,b1,strata):
 M=exact_excess(graph,strata);u=b0@M;v=b1@M;dp=PC-PA;dc=v-u
 early=np.sum(u*PA,1);late=np.sum(v*PC,1);pl=np.sum(u*dp,1);co=np.sum(dc*PA,1);inter=np.sum(dc*dp,1)
 np.testing.assert_allclose(late-early,pl+co+inter,atol=1e-12,rtol=0)
 return early,{'plant_symmetric':pl+.5*inter,'consumer_symmetric':co+.5*inter,'plant_early_consumer_fixed':pl}

gates=[];results=[];losses=[];contrasts=[]
for cutoff in ['non_wikipedia','hosts_or_bhl']:
 C=np.zeros_like(A)
 for x in edges:
  labels=set(x['source_labels'])
  keep=bool(labels-{'Wikipedia'}) if cutoff=='non_wikipedia' else bool(labels&{'HOSTS DB','BHL'})
  if keep:C[bidx[x['butterflyKey']],pidx[x['plantKey']]]=1
 # Filter only to accepted GBIF keys, with a reliable edge in both graphs.
 both=L&np.array([bmeta[k]['status']=='ACCEPTED' for k in bk])&(A.sum(1)>0)&(C.sum(1)>0)
 ao,ac=A[both],C[both];oldmapped=L&(A.sum(1)>0)
 denom=PB[:,oldmapped].sum(1);frac=np.divide(PB[:,both].sum(1),denom,out=np.zeros(len(d)),where=denom>0)
 cov=(frac[a]>=.8)&(frac[b]>=.8);overlap=int(np.sum((ao>0)&(ac>0)));union=int(np.sum((ao>0)|(ac>0)))
 for tag,sm in [('quality',np.asarray(quality)),('all',np.ones(len(q),bool))]:
  ids=np.flatnonzero(sm&cov);ng=len(set(H6[ids]));nold=int((ao.sum(0)>0).sum());nnew=int((ac.sum(0)>0).sum());gate={'source_filter':cutoff,'sample':tag,'common_butterflies':int(both.sum()),'old_plants':nold,'new_plants':nnew,'old_edges':int(ao.sum()),'new_edges':int(ac.sum()),'shared_edges':overlap,'edge_jaccard':overlap/union if union else None,'H7':len(ids),'H6':ng,'median_coverage_early':float(np.median(frac[a][sm])),'median_coverage_late':float(np.median(frac[b][sm]))}
  gate['status']='FIT' if both.sum()>=30 and min(nold,nnew)>=60 and len(ids)>=40 and ng>=20 else 'STOP_SUPPORT';gates.append(gate)
  if gate['status']!='FIT':continue
  b0,b1=PB[a][:,both],PB[b][:,both];cat=np.array([f'{c}::{m}' for c,m in zip(country,method)])
  store=collections.defaultdict(list)
  for f in sorted(set(fold[ids])):
   te=ids[fold[ids]==f];tr=ids[fold[ids]!=f];assert not set(H6[tr])&set(H6[te])
   if len(set(H6[tr]))<15:raise RuntimeError('STOP fold support')
   st=shared_strata(tr)
   for graphname,graph in [('old',ao),('filtered',ac)]:
    y0,ys=decompose(graph,b0,b1,st);x=np.column_stack((X0,y0));x1=np.column_stack((x,E0))
    for endpoint,y in ys.items():
     p0,e0,_=fit(tr,te,y,x,cat);p1,e1,coef=fit(tr,te,y,x1,cat)
     for g in sorted(set(H6[te])):
      zz=H6[te]==g;l0=float(np.mean((y[te][zz]-p0[zz])**2));l1=float(np.mean((y[te][zz]-p1[zz])**2));store[graphname,endpoint].append((g,l0,l1));losses.append(dict(source_filter=cutoff,sample=tag,catalogue=graphname,endpoint=endpoint,h6=g,mse0=l0,mse1=l1,independent_solver_diff=max(e0,e1)))
  boots={}
  for (graphname,endpoint),ls in store.items():
   dd=pd.DataFrame(ls,columns=['h6','mse0','mse1']).sort_values('h6');assert len(dd)==ng and not dd.h6.duplicated().any();m=dd[['mse0','mse1']].to_numpy();v=m[np.random.default_rng(20261124).integers(0,ng,(5000,ng))].mean(1);l0,l1=m.mean(0)
   if l0<=1e-20:results.append(dict(source_filter=cutoff,sample=tag,catalogue=graphname,endpoint=endpoint,status='STOP_CONSTANT'));continue
   bs=100*(v[:,0]-v[:,1])/v[:,0];boots[graphname,endpoint]=bs;lo,hi=np.quantile(bs,[.025,.975]);results.append(dict(source_filter=cutoff,sample=tag,catalogue=graphname,endpoint=endpoint,H7=len(ids),H6=ng,gain=100*(l0-l1)/l0,lo=lo,hi=hi,status='FIT'))
  for ep in ['plant_symmetric','consumer_symmetric','plant_early_consumer_fixed']:
   if ('old',ep) in boots and ('filtered',ep) in boots:
    diff=boots['filtered',ep]-boots['old',ep];lo,hi=np.quantile(diff,[.025,.975]);contrasts.append(dict(source_filter=cutoff,sample=tag,endpoint=ep,filtered_minus_old_lo=lo,filtered_minus_old_hi=hi))
for n,rr in [('gate',gates),('results',results),('losses',losses),('contrasts',contrasts)]:pd.DataFrame(rr).to_csv(P/f'gesui_source_label_{n}.csv',index=False)
print(pd.DataFrame(gates).to_string(index=False));print(pd.DataFrame(results).to_string(index=False))
