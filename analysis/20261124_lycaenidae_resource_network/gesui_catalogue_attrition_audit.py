#!/usr/bin/env python3
"""Posthoc attribution of old-catalogue common-support attenuation; not alternate catalogue."""
from pathlib import Path
import sys
import numpy as np,pandas as pd
P=Path(sys.argv[1]); prior=P/'gesui_source_label_sensitivity.py'
exec(compile(prior.read_text().split('gates=[];results=[];losses=[];contrasts=[]')[0],str(prior),'exec'))
C=np.zeros_like(A)
for x in edges:
 if set(x['source_labels'])-{'Wikipedia'}:C[bidx[x['butterflyKey']],pidx[x['plantKey']]]=1
old=L&(A.sum(1)>0); common=old & np.array([bmeta[k]['status']=='ACCEPTED' for k in bk]) & (C.sum(1)>0)
assert old.sum()==52 and common.sum()==40
frac=PB[:,common].sum(1)/np.maximum(PB[:,old].sum(1),1e-30)
cov=(frac[a]>=.8)&(frac[b]>=.8)
rows=[];losses=[];cat=np.array([f'{c}::{m}' for c,m in zip(country,method)])
for sample,sm in [('quality',np.asarray(quality)),('all',np.ones(len(q),bool))]:
 for sites,site_mask in [('all_eligible',sm),('common_80pct',sm&cov)]:
  ids=np.flatnonzero(site_mask); ng=len(set(H6[ids])); assert ng>=20
  store=collections.defaultdict(list)
  for k in sorted(set(fold[ids])):
   te=ids[fold[ids]==k];tr=ids[fold[ids]!=k];assert not set(H6[te])&set(H6[tr]);st=shared_strata(tr)
   for taxon,tm in [('old52',old),('common40',common)]:
    y0,ys=decompose(A[tm],PB[a][:,tm],PB[b][:,tm],st)
    X=np.column_stack((X0,y0));X1=np.column_stack((X,E0))
    for endpoint,y in ys.items():
     p0,e0,_=fit(tr,te,y,X,cat);p1,e1,coef=fit(tr,te,y,X1,cat)
     for h in sorted(set(H6[te])):
      z=H6[te]==h;l0=float(np.mean((y[te][z]-p0[z])**2));l1=float(np.mean((y[te][z]-p1[z])**2));store[taxon,endpoint].append((h,l0,l1));losses.append(dict(sample=sample,sites=sites,taxa=taxon,endpoint=endpoint,h6=h,mse0=l0,mse1=l1,independent_solver_diff=max(e0,e1)))
  for (taxon,endpoint),ls in store.items():
   dt=pd.DataFrame(ls,columns=['h6','mse0','mse1']).sort_values('h6');assert len(dt)==ng and not dt.h6.duplicated().any();m=dt[['mse0','mse1']].to_numpy();l0,l1=m.mean(0);v=m[np.random.default_rng(20261124).integers(0,ng,(5000,ng))].mean(1);bs=100*(v[:,0]-v[:,1])/v[:,0];lo,hi=np.quantile(bs,[.025,.975]);rows.append(dict(sample=sample,sites=sites,taxa=taxon,endpoint=endpoint,H7=len(ids),H6=ng,gain=100*(l0-l1)/l0,lo=lo,hi=hi))
r=pd.DataFrame(rows);r.to_csv(P/'gesui_catalogue_attrition_results.csv',index=False);pd.DataFrame(losses).to_csv(P/'gesui_catalogue_attrition_losses.csv',index=False)
print(r[r.endpoint=='plant_symmetric'].to_string(index=False));print('Independent solver maximum error',max(x['independent_solver_diff'] for x in losses))
