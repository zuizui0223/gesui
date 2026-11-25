#!/usr/bin/env python3
"""Reconstruct source-locked 205-H7 cohort, cross-country and landscape diagnostics."""
import os,sys,json,hashlib,io,zipfile,collections
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,pandas as pd
from sklearn.linear_model import Ridge
src=Path(sys.argv[1]); priorout=Path(sys.argv[2]); out=Path(sys.argv[3]);out.mkdir(exist_ok=True,parents=True)
# Replay original data sources and covariates from prior exact code; no downstream network outcomes.
prior=src/'gesui_fab_resource_consumer.py'; exec(compile(prior.read_text().split('# Recreate expected shuffled graph')[0],str(prior),'exec'))
sha=lambda b:hashlib.sha256(b).hexdigest()
ez=src/'gesui_eiv_source.zip';assert sha(ez.read_bytes())=='1d3051d4c527421bec3ec3ba917422a9162afb81cc7f38db0d8130287174a14c'
z=zipfile.ZipFile(ez);bb=z.read('EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv');assert sha(bb)=='ab45f36c0dad047442bb958bcef3135b45014bafec8a70811bcafb6109a02679'
tx=pd.read_csv(io.BytesIO(bb),dtype={'speciesKey':str});S=tx.set_index('speciesKey').EIV_NUTRI_exact.dropna().to_dict();assert len(tx)==11511 and len(S)==4011
new=pd.read_csv(priorout/'filter_clarke_edges.csv',dtype={'plantKey':str,'butterflyKey':str});assert len(new)==1974
assert hashlib.sha1(f'blob {(priorout/"filter_clarke_edges.csv").stat().st_size}\0'.encode()+(priorout/'filter_clarke_edges.csv').read_bytes()).hexdigest()=='1be721502c730e33a0249ef67b01265682550cf8'
blocked=set(pk)|set(new.plantKey);names=tx.set_index('speciesKey').gbif_canonical_name.to_dict();gen={str(names[k]).split()[0].casefold() for k in blocked if k in names and isinstance(names[k],str)};nskeys={k for k in S if str(names[k]).split()[0].casefold() not in gen};assert len(nskeys)==1876
plant=[]
for x in d.itertuples():
 c=fac[x.h,x.period,'Plantae']['source_speciesKey_counts'];v=[S[k] for k in c if k in S];u=[S[k] for k in c if k in nskeys];plant.append([np.mean(v),len(v),len(v)/len(c),np.mean(u) if u else np.nan,len(u)])
d[['plant_N','plant_score_n','plant_score_cov','nonhost_N','nonhost_n']]=np.array(plant);d['plant_log_n']=np.log1p(d.plant_score_n);d['late']=(d.period=='2020-2022').astype(float);d['EQR']=d.EQR_original_site_equal_basin_period_mean;d['fold']=d.h6.map(lambda x:int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5);d['method']=d.method.replace('','__UNKNOWN__');d['category']=d.country_holdout+'::'+d.method;d['crop_local']=[at[h]['crp_pc_sse'] for h in d.h];d['urban_local']=[at[h]['urb_pc_sse'] for h in d.h]
roster=pd.read_csv(priorout/'filter_butterfly_scores.csv',dtype={'butterflyKey':str});assert hashlib.sha1(f'blob {(priorout/"filter_butterfly_scores.csv").stat().st_size}\0'.encode()+(priorout/'filter_butterfly_scores.csv').read_bytes()).hexdigest()=='1529afe532132a88f3c5931fbb3d4c99524e6e1d'
panel=pd.read_csv(priorout/'filter_panel.csv',dtype={'h':str,'h6':str});orig=pd.read_csv(z.open('EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv'),dtype={'h':str});elig=pd.read_csv(z.open('EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv'),dtype={'h':str});assert len(orig)==451
chk=d.merge(orig[['h','period','Plant_N_incidence']],on=['h','period'],validate='one_to_one');np.testing.assert_allclose(chk.plant_N,chk.Plant_N_incidence,atol=1e-12,rtol=0)
d=d.merge(elig[['h','period','full_source_trait_eligible']],on=['h','period'],validate='one_to_one');assert (d.h.to_numpy()==prev.h.to_numpy()).all()
base=['area','riverflow','naturalout','nEEA','nPlant','nAnimal','nButterfly','plant_score_cov','plant_log_n','common_record_cov','common_log_n','late'];land=['crop','urban'];local=['crop_local','urban_local'];
results=[];losses=[];gate=[];allpanel=[]
def fit(df,y,cols,tr,te,mode):
 a=df[cols].to_numpy(float); yy=df[y].to_numpy(float); g=df.h6.to_numpy();cat=df['category' if mode=='H6' else 'method'].to_numpy();assert np.isfinite(a).all() and np.isfinite(yy).all();assert not set(g[tr])&set(g[te]);mu=a[tr].mean(0);sd=a[tr].std(0);sd[sd<1e-12]=1.;C=np.column_stack([cat==x for x in sorted(set(cat[tr]))]).astype(float);X=np.column_stack((np.ones(len(df)),(a-mu)/sd,C));cnt=collections.Counter(g[tr]);w=np.array([1/cnt[h] for h in g[tr]]);w/=w.mean();pen=np.eye(X.shape[1])*10;pen[0,0]=0.;beta=np.linalg.solve(X[tr].T@(w[:,None]*X[tr])+pen,X[tr].T@(w*yy[tr]));p=X[te]@beta;alt=Ridge(alpha=10,fit_intercept=True).fit(X[tr,1:],yy[tr],sample_weight=w).predict(X[te,1:]);err=float(np.max(np.abs(p-alt)));assert err<1e-8
 return p,err

def run(df,cutoff,mode,specs):
 df=df.reset_index(drop=True);g=df.h6.to_numpy();ids=np.arange(len(df));splits=[]
 if mode=='H6':
  for k in sorted(set(df.fold)):splits.append((str(k),ids[df.fold!=k],ids[df.fold==k]))
 else:
  counts=df.groupby('h').country.nunique();df=df[df.h.isin(counts[counts==1].index)].reset_index(drop=True);g=df.h6.to_numpy();ids=np.arange(len(df))
  for c in sorted(set(df.country)):
   te=ids[df.country==c];tr=ids[(df.country!=c)&~np.isin(g,g[te])];nh=df.iloc[te].h.nunique();ng=len(set(g[te]));th=df.iloc[tr].h.nunique();tg=len(set(g[tr]));ok=nh>=5 and ng>=3 and th>=30 and tg>=15
   gate.append(dict(cutoff=cutoff,country=c,H7=nh,H6=ng,train_H7=th,train_H6=tg,status='FIT' if ok else 'STOP_SUPPORT'))
   if ok:splits.append((c,tr,te))
  assert len(splits)>=2
 for label,y,b,extra in specs:
  rec=[];mx=0.
  for k,tr,te in splits:
   p0,e0=fit(df,y,b,tr,te,mode);p1,e1=fit(df,y,b+extra,tr,te,mode);mx=max(mx,e0,e1);yy=df[y].to_numpy(float)[te]
   for h in sorted(set(g[te])):
    z=g[te]==h;a=float(np.mean((yy[z]-p0[z])**2));b1=float(np.mean((yy[z]-p1[z])**2));rec.append((k,h,a,b1));losses.append(dict(cutoff=cutoff,mode=mode,test=label,outcome=y,fold=k,h6=h,mse0=a,mse1=b1))
  ll=pd.DataFrame(rec,columns=['fold','h6','mse0','mse1']);
  if mode=='H6':assert not ll.h6.duplicated().any();partitions=[('ALL',ll)]
  else:partitions=list(ll.groupby('fold'))
  for k,ss in partitions:
   m=ss[['mse0','mse1']].to_numpy();v=m[np.random.default_rng(20261125).integers(0,len(m),(5000,len(m)))].mean(1);a,b1=m.mean(0);q=np.quantile(100*(v[:,0]-v[:,1])/v[:,0],[.025,.975]);results.append(dict(cutoff=cutoff,mode=mode,test=label,outcome=y,fold=k,H6=ss.h6.nunique(),gain=100*(a-b1)/a,lo=q[0],hi=q[1],mse0=a,mse1=b1,max_solver_diff=mx))
  if mode=='country':
   m=ll.groupby('fold')[['mse0','mse1']].mean().to_numpy();a,b1=m.mean(0);results.append(dict(cutoff=cutoff,mode=mode,test=label,outcome=y,fold='COUNTRY_EQUAL',H6=ll.h6.nunique(),gain=100*(a-b1)/a,lo=np.nan,hi=np.nan,mse0=a,mse1=b1,max_solver_diff=mx))
for cutoff in [3,2]:
 r=roster[(roster.cutoff==cutoff)&roster.common_eligible].set_index('butterflyKey');mask=np.array([k in r.index for k in bk]);assert mask.sum()==94;x=B[:,mask];c=x>0;n=c.sum(1);df=d.copy();df['common_record_cov']=x.sum(1)/B.sum(1);df['common_n']=n;df['common_log_n']=np.log1p(n)
 for lab in ['old','clarke']:
  v=r.loc[np.array(bk)[mask],lab+'_N'].to_numpy(float);df[lab+'_HNS']=np.divide(c@v,n,out=np.full(len(df),np.nan),where=n>0)
 ok=(df.common_record_cov>=.5)&(n>0)&df.full_source_trait_eligible.astype(bool);pairs=df.loc[ok].groupby('h').period.nunique();df=df[df.h.isin(pairs[pairs==2].index)].copy();assert len(df)==414 and df.h6.nunique()==81
 ck=df.merge(panel[panel.cutoff==cutoff],on=['h','period'],suffixes=('','_saved'),validate='one_to_one');assert len(ck)==414
 for k in ['h6','country','fold']:assert (ck[k]==ck[k+'_saved']).all(),k
 for k in ['plant_N','old_HNS','clarke_HNS','common_record_cov','common_n','plant_score_cov','EQR']:np.testing.assert_allclose(ck[k],ck[k+'_saved'],atol=1e-12,rtol=0)
 valid=df[df.nonhost_n>=30].groupby('h').period.nunique();df=df[df.h.isin(valid[valid==2].index)].copy();assert len(df)==410 and df.h6.nunique()==81
 df['cutoff']=cutoff;allpanel.append(df[['cutoff','h','h6','period','country','method','fold','plant_N','nonhost_N','nonhost_n','old_HNS','clarke_HNS','W','EQR','crop','urban','crop_local','urban_local']+base])
 specs=[]
 for y in ['old_HNS','clarke_HNS']:
  specs += [('plant_after_W_land',y,base+land+['W'],['plant_N']),('nonhost_after_W_land',y,base+land+['W'],['nonhost_N']),('land_after_W',y,base+['W'],land),('W_after_land',y,base+land,['W']),('land_after_W_plant',y,base+['W','nonhost_N'],land),('W_after_land_plant',y,base+land+['nonhost_N'],['W']),('nonhost_after_W_local_land',y,base+local+['W'],['nonhost_N'])]
 for y in ['plant_N','EQR']:specs += [('land_after_W',y,base+['W'],land),('W_after_land',y,base+land,['W'])]
 for mode in ['H6','country']:run(df,cutoff,mode,specs)
for name,rr in [('results',results),('losses',losses),('country_gate',gate)]:pd.DataFrame(rr).to_csv(out/f'country_landuse_{name}.csv',index=False)
pd.concat(allpanel).to_csv(out/'country_landuse_panel.csv',index=False)
rr=pd.DataFrame(results);print(pd.DataFrame(gate).to_string(index=False));print(rr[(rr.cutoff==3)&(rr.fold.isin(['ALL','COUNTRY_EQUAL']))][['mode','test','outcome','gain','lo','hi']].round(3).to_string(index=False));print(rr[(rr.cutoff==3)&(rr['mode']=='country')&(rr.test=='nonhost_after_W_land')][['outcome','fold','H6','gain','lo','hi']].round(3).to_string(index=False));print('max independent solver difference',rr.max_solver_diff.max())
(out/'SOURCE_RECEIPT.json').write_text(json.dumps({'original_archives':D,'EIV_sha256':sha(ez.read_bytes()),'previous_207_H7_panel_scores_reproduced':True,'prior_filter_scores_git_blob':'1529afe532132a88f3c5931fbb3d4c99524e6e1d','prior_edges_git_blob':'1be721502c730e33a0249ef67b01265682550cf8','primary_205_H7':205,'models':len(rr),'independent_max_diff':float(rr.max_solver_diff.max()),'climate_available':False,'atlas_available_fields':sorted(next(iter(at.values())).keys())},indent=2))

# Fixed country-loss diagnostic: 5 countries, not confirmatory population inference.
ll=pd.DataFrame(losses); agg=[]
for (cut,test,y),gg in ll[ll['mode']=='country'].groupby(['cutoff','test','outcome']):
 m=gg.groupby('fold')[['mse0','mse1']].mean().to_numpy();v=m[np.random.default_rng(20261125).integers(0,len(m),(10000,len(m)))].mean(1);gain=lambda x:100*(x[:,0].mean()-x[:,1].mean())/x[:,0].mean();bs=100*(v[:,0]-v[:,1])/v[:,0];loo=[gain(np.delete(m,i,0)) for i in range(len(m))];agg.append(dict(cutoff=cut,test=test,outcome=y,countries=len(m),gain=gain(m),country_boot_lo=np.quantile(bs,.025),country_boot_hi=np.quantile(bs,.975),leave_one_country_min=min(loo),leave_one_country_max=max(loo),countries_improved=int((m[:,0]>m[:,1]).sum())))
pd.DataFrame(agg).to_csv(out/'country_aggregate_sensitivity.csv',index=False)
