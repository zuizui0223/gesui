#!/usr/bin/env python3
"""Descriptive audit of fixed all-eligible training and coverage-conditioned heldout loss.
Coverage is derived from both periods, not a prospective selection covariate.
"""
from pathlib import Path
import sys
import numpy as np,pandas as pd
P=Path(sys.argv[1]);prior=P/'gesui_catalogue_attrition_audit.py'
exec(compile(prior.read_text().split('rows=[];losses=[];cat=')[0],str(prior),'exec'))
cat=np.array([f'{c}::{m}' for c,m in zip(country,method)])
loss=[];rows=[];desc=[];contrasts=[]
for sample,sm in [('quality',np.asarray(quality)),('all',np.ones(len(q),bool))]:
 ids=np.flatnonzero(sm)
 for k in sorted(set(fold[ids])):
  te=ids[fold[ids]==k];tr=ids[fold[ids]!=k];assert not set(H6[te])&set(H6[tr]);st=shared_strata(tr)
  for taxon,tm in [('old52',old),('common40',common)]:
   y0,ys=decompose(A[tm],PB[a][:,tm],PB[b][:,tm],st);X=np.column_stack((X0,y0));X1=np.column_stack((X,E0))
   y=ys['plant_symmetric'];p0,e0,_=fit(tr,te,y,X,cat);p1,e1,_=fit(tr,te,y,X1,cat)
   for j,i in enumerate(te):loss.append(dict(sample=sample,taxa=taxon,h=q.h.iloc[i],h6=H6[i],country=country[i],fold=int(k),retained=bool(cov[i]),mse0=float((y[i]-p0[j])**2),mse1=float((y[i]-p1[j])**2),independent_solver_diff=max(e0,e1)))
 for retained in [True,False]:
  z=sm&(cov==retained);i=np.flatnonzero(z)
  for c in sorted(set(country[i])):
   ii=i[country[i]==c];desc.append(dict(sample=sample,retained=retained,country=c,H7=len(ii),H6=len(set(H6[ii]))))
  desc.append(dict(sample=sample,retained=retained,country='ALL',H7=len(i),H6=len(set(H6[i])),EQR_mean=float(np.mean(E0[i])),EQR_sd=float(np.std(E0[i])),coverage_early_median=float(np.median(frac[a][i])),coverage_late_median=float(np.median(frac[b][i]))))
ll=pd.DataFrame(loss)
for (sample,taxa),t in ll.groupby(['sample','taxa']):
 assert not t.h.duplicated().any()
 for subset,tt in [('all',t),('retained',t[t.retained]),('excluded',t[~t.retained])]:
  x=tt.groupby('h6')[['mse0','mse1']].mean().to_numpy();m0,m1=x.mean(0);idx=np.random.default_rng(20261124).integers(0,len(x),(5000,len(x)));bs=x[idx].mean(1);v=100*(bs[:,0]-bs[:,1])/bs[:,0];lo,hi=np.quantile(v,[.025,.975]);rows.append(dict(sample=sample,taxa=taxa,subset=subset,H7=len(tt),H6=len(x),gain=100*(m0-m1)/m0,lo=lo,hi=hi,mse0=m0,mse1=m1))
 # Compare within shared H6 blocks only; paired reference groups, same fitted predictions.
 common_h=sorted(set(t.loc[t.retained,'h6'])&set(t.loc[~t.retained,'h6']));
 if len(common_h)>=15:
  z1=t[t.retained].groupby('h6')[['mse0','mse1']].mean().loc[common_h].to_numpy();z2=t[~t.retained].groupby('h6')[['mse0','mse1']].mean().loc[common_h].to_numpy();ix=np.random.default_rng(20261124).integers(0,len(common_h),(5000,len(common_h)));v1=z1[ix].mean(1);v2=z2[ix].mean(1);b1=100*(v1[:,0]-v1[:,1])/v1[:,0];b2=100*(v2[:,0]-v2[:,1])/v2[:,0];m1=z1.mean(0);m2=z2.mean(0);g1=100*(m1[0]-m1[1])/m1[0];g2=100*(m2[0]-m2[1])/m2[0];lo,hi=np.quantile(b2-b1,[.025,.975]);contrasts.append(dict(sample=sample,taxa=taxa,shared_H6=len(common_h),retained_gain=g1,excluded_gain=g2,difference=g2-g1,lo=lo,hi=hi,status='descriptive'))
 else:contrasts.append(dict(sample=sample,taxa=taxa,shared_H6=len(common_h),status='STOP_SUPPORT'))
for n,data in [('results',rows),('losses',loss),('coverage_description',desc),('shared_H6_contrasts',contrasts)]:pd.DataFrame(data).to_csv(P/f'gesui_catalogue_coverage_{n}.csv',index=False)
print(pd.DataFrame(rows).round(3).to_string(index=False));print(pd.DataFrame(desc).query('country=="ALL"').round(3).to_string(index=False));print(pd.DataFrame(contrasts).round(3).to_string(index=False));print('solver max',ll.independent_solver_diff.max())
