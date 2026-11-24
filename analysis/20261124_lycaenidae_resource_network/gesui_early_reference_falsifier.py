#!/usr/bin/env python3
"""Post-sync sensitivity: early-only geography null and nested genus identity.
No new field data; static host catalogue and retrospective paired GBIF record index.
"""
from pathlib import Path
import sys
P=Path(sys.argv[1] if len(sys.argv)>1 else '/mnt/data')
prior=P/'gesui_fab_resource_consumer.py'
exec(compile(prior.read_text().split('base=q[K].to_numpy(float);ee=q[E].to_numpy(float)')[0],str(prior),'exec'))
# No late source occurrences contribute to the new reference.
early=d.period.eq('2017-2019').to_numpy(); early_h6=sorted(d.loc[early,'h6'].unique())
freq=collections.Counter()
for x in raw:
 if x['group']=='Plantae' and x['period']=='2017-2019':freq.update(x['source_speciesKey_counts'].keys())
bins_early=np.digitize([freq[k] for k in pk],[0,1,5,20,100,468],right=True)
inc=np.array([(Pcount[early&(d.h6.to_numpy()==h)]>0).any(0) for h in early_h6]).T.astype(float)
norm=np.linalg.norm(inc,axis=1); labels=np.full(len(pk),-1); eligible=norm>0
labels[eligible]=KMeans(n_clusters=6,random_state=20261114,n_init=10).fit_predict(inc[eligible]/norm[eligible,None])
st=np.array([pf[k]+'|'+str(bb)+'|geo'+str(cc) for k,bb,cc in zip(pk,bins_early,labels)])
rng=np.random.default_rng(20261114); slices=[np.flatnonzero(st==s) for s in sorted(set(st))];N0=np.zeros_like(A)
for _ in range(499):
 perm=np.arange(len(pk))
 for ix in slices:perm[ix]=rng.permutation(ix)
 N0+=A[:,perm]
N0/=499
np.testing.assert_allclose(N0.sum(1),A.sum(1),atol=1e-9,rtol=0)
# Verify accepted genus keys from the original source; exact and consistent only.
from collections import defaultdict
mapping=defaultdict(set)
for x in cross:
 k=x.get('accepted_speciesKey');m=x.get('match',{})
 if k in pidx and m.get('genus') and m.get('family'):
  mapping[k].add((m['family'],m['genus']))
assert all(len(mapping[k])==1 and next(iter(mapping[k]))[0]==pf[k] for k in pk)
fam=np.array([pf[k] for k in pk]);gen=np.array(['|'.join(next(iter(mapping[k]))) for k in pk])
qout=q[['h','h6','country_holdout','method']].copy(); qout['quality']=quality
B0=q[K].to_numpy(float);ee=q[E].to_numpy(float)
# Exact symmetric decomposition, then partition plant side by nested taxa.
def split_mass(C,p0,p1,groups):
 qty=np.zeros(len(q));ident=qty.copy()
 for g in groups:
  x=p0[:,g];y=p1[:,g];c=C[:,g];m0=x.sum(1);m1=y.sum(1);t0=(x*c).sum(1);t1=(y*c).sum(1)
  u0=np.divide(t0,m0,out=np.zeros_like(t0),where=m0>1e-14);u1=np.divide(t1,m1,out=np.zeros_like(t1),where=m1>1e-14)
  u0=np.where(m0>1e-14,u0,u1);u1=np.where(m1>1e-14,u1,u0)
  qty+=(m1-m0)*.5*(u0+u1);ident+=(u1-u0)*.5*(m0+m1)
 return qty,ident
meta={'source_digests':D,'early_H7periods':int(early.sum()),'early_H6':len(early_h6),'early_only_nonzero_host_profiles':int(eligible.sum()),'early_only_strata':len(slices),'host_genus':len(set(gen)),'host_family':len(set(fam)),'paired_H7':len(q),'quality_H7':int(quality.sum()),'quality_H6':q.loc[quality,'h6'].nunique(),'variants':{},'source_note':'early-only null is retrospective/transductive on early heldout H6; no observed feeding or wastewater causality'}
country_rows=[]
def country_test(name,mask,y,X,pred):
 ids=np.where(mask & np.isfinite(y)&np.isfinite(X).all(1)&np.isfinite(pred))[0];gs=q.h6.to_numpy()[ids];ctry=q.country_holdout.to_numpy()[ids];meth=q.method.replace('','__UNKNOWN__').to_numpy()[ids];y=y[ids];X=X[ids];pred=pred[ids];out=[]
 for c in sorted(set(ctry)):
  test=ctry==c;train=(~test)&~np.isin(gs,gs[test]);
  if test.sum()<10 or len(set(gs[test]))<5 or train.sum()<30 or len(set(gs[train]))<15:continue
  assert not set(gs[test])&set(gs[train]);w=np.array([1/np.sum(gs[train]==g) for g in gs[train]]);w/=w.mean();cats=sorted(set(meth[train]));catsX=np.column_stack([meth==c0 for c0 in cats]).astype(float)
  pp=[]
  for add in (0,1):
   x=np.column_stack((X,pred)) if add else X;mu=x[train].mean(0);sd=x[train].std(0);sd[sd<1e-12]=1;z=np.column_stack(((x-mu)/sd,catsX));p=Ridge(alpha=10,fit_intercept=True).fit(z[train],y[train],sample_weight=w).predict(z[test]);pp.append(p)
  losses=np.array([[np.mean((y[test][gs[test]==g]-p[gs[test]==g])**2) for p in pp] for g in sorted(set(gs[test]))]);m0,m1=losses.mean(0);out.append(dict(name=name,country=c,n=int(test.sum()),H6=len(set(gs[test])),train_H7=int(train.sum()),train_H6=len(set(gs[train])),mse0=m0,mse1=m1,gain=100*(m0-m1)/m0))
 country_rows.extend(out)
 if out:
  m0=np.mean([v['mse0'] for v in out]);m1=np.mean([v['mse1'] for v in out]);country_rows.append(dict(name=name,country='COUNTRY_EQUAL',n=sum(v['n'] for v in out),H6=None,train_H7=None,train_H6=None,mse0=m0,mse1=m1,gain=100*(m0-m1)/m0))
for label,Nx in [('full_reference',N),('early_reference',N0)]:
 M=A[L]-Nx[L];ba=PB[a][:,L];bb=PB[b][:,L];pa=PP[a];pb=PP[b];C=.5*(ba+bb)@M
 val=lambda u,v:((u@M)*v).sum(1)
 early_ex=val(ba,pa);late_ex=val(bb,pb);e01=val(ba,pb);e10=val(bb,pa)
 plant=((pb-pa)*C).sum(1);consumer=.5*((e10-early_ex)+(late_ex-e01));delta=late_ex-early_ex
 np.testing.assert_allclose(plant+consumer,delta,atol=1e-12,rtol=0)
 if label=='full_reference':np.testing.assert_allclose(plant,q.L_plant_delta,atol=1e-12,rtol=0)
 fq,fi=split_mass(C,pa,pb,[fam==f for f in sorted(set(fam))]);gq,gi=split_mass(C,pa,pb,[gen==g for g in sorted(set(gen))]);bg=gq-fq
 np.testing.assert_allclose(fq+bg+gi,plant,atol=1e-12,rtol=0);np.testing.assert_allclose(bg+gi,fi,atol=1e-12,rtol=0)
 outcomes={'plant_delta':plant,'consumer_delta':consumer,'total_delta':delta,'family_mass':fq,'within_family_identity':fi,'genus_within_family_mass':bg,'within_genus_identity':gi}
 for key,y in outcomes.items():qout[label+'__'+key]=y
 X=np.column_stack((B0,early_ex))
 for tag,mask in [('quality',quality),('all',np.ones(len(q),bool))]:
  for key,y in outcomes.items():
   evaluate(label+'__'+tag+'__'+key,mask,y,X,ee)
  country_test(label+'__'+tag+'__plant_delta',mask,plant,X,ee)
 for key in ['plant_delta','within_genus_identity','genus_within_family_mass']:
  for pred in ['mean_original_2018_logTP','mean_original_2018_EQR']:
   evaluate(label+'__FR73__'+pred+'__'+key,frmask,outcomes[key],X,q[pred].to_numpy(),'loho')
 meta['variants'][label]={'identity_max_error':float(np.max(np.abs(fq+bg+gi-plant))),'paired_plant_sd':float(plant.std()),'paired_within_genus_sd':float(gi.std()),'paired_genus_mass_sd':float(bg.std())}
meta['fitted_models']=len(results);meta['independent_prediction_max_diff']=max(r['independent_prediction_diff'] for r in results)
pd.DataFrame(results).to_csv(P/'gesui_early_reference_results.csv',index=False)
pd.DataFrame(country_rows).to_csv(P/'gesui_early_reference_country.csv',index=False)
qout.to_csv(P/'gesui_early_reference_panel.csv',index=False)
(P/'gesui_early_reference_receipt.json').write_text(json.dumps(meta,indent=2)+'\n')
print(json.dumps(meta,indent=2));print(pd.DataFrame(results)[['name','n','H6','gain','lo','hi']].round(3).to_string(index=False));print(pd.DataFrame(country_rows).round(3).to_string(index=False))
