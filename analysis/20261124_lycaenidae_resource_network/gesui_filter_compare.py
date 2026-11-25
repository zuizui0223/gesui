#!/usr/bin/env python3
"""Same-grain environmental filtering, observational predictive comparison."""
from pathlib import Path
import sys
P=Path(sys.argv[1]); prior=P/'gesui_fab_resource_consumer.py'
# Original signed source checks and original 454-row geography and counts; no network outcome tuning.
exec(compile(prior.read_text().split('# Recreate expected shuffled graph')[0],str(prior),'exec'))
from scipy.stats import spearmanr
import io
from collections import defaultdict
Ez=P/'gesui_eiv_source.zip';assert hashlib.sha256(Ez.read_bytes()).hexdigest()=='1d3051d4c527421bec3ec3ba917422a9162afb81cc7f38db0d8130287174a14c'
z=zipfile.ZipFile(Ez);raw=z.read('EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv');assert hashlib.sha256(raw).hexdigest()=='ab45f36c0dad047442bb958bcef3135b45014bafec8a70811bcafb6109a02679'
tx=pd.read_csv(io.BytesIO(raw),dtype={'speciesKey':str});assert len(tx)==11511
S=tx.set_index('speciesKey').EIV_NUTRI_exact.dropna().to_dict();assert len(S)==4011
r=json.loads((P/'filter_mapping_receipt.json').read_text());assert hashlib.sha256((P/'filter_butterfly_scores.csv').read_bytes()).hexdigest()==r['score_sha256']
roster=pd.read_csv(P/'filter_butterfly_scores.csv',dtype={'butterflyKey':str});gates=[];result=[];losses=[];contrast=[];panels=[];taxa=[]
orig=pd.read_csv(z.open('EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv'),dtype={'h':str});elig=pd.read_csv(z.open('EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv'),dtype={'h':str,'h6':str})
assert len(orig)==451 and len(elig)==454
# Rebuild all plant community scores from exact EIV species and actual original GBIF histograms.
allkeys=set();plant=[]
for x in d.itertuples():
 counts=fac[x.h,x.period,'Plantae']['source_speciesKey_counts'];allkeys.update(counts);vals=[S[k] for k in counts if k in S]
 plant.append(dict(h=x.h,period=x.period,plant_N=float(np.mean(vals)) if vals else np.nan,plant_score_n=len(vals),plant_score_cov=len(vals)/len(counts)))
assert allkeys==set(tx.speciesKey)
plant=pd.DataFrame(plant);check=plant.merge(orig,on=['h','period'],validate='one_to_one');np.testing.assert_allclose(check.plant_N,check.Plant_N_incidence,atol=1e-12,rtol=0);assert (check.plant_score_n==check.n_scored_species).all()
d=d.merge(plant,on=['h','period'],validate='one_to_one').merge(elig[['h','period','full_source_trait_eligible']],on=['h','period'],validate='one_to_one');assert (d.h.to_numpy()==prev.h.to_numpy()).all()
# Actual same physical 263 stations in 73 French H7; assign measured water variables to corresponding broad period.
member='FR73_HYBAS7_NATIVE_2018_2022_REAL_TP_I2M2_EQR_AND_GBIF_INDEX_DELTA.csv';bb=Z['FR73_original_14013312711.zip'].read(member);assert hashlib.sha256(bb).hexdigest()=='f2ff4e56f245e1aa8f84f5ac099d68ec1c7e27a786c1067439db388b186dd90b'
fr=pd.read_csv(io.BytesIO(bb),dtype={'HYBAS_ID_L07':str,'HYBAS_ID_L06':str});assert len(fr)==73 and fr.HYBAS_ID_L06.nunique()==24 and fr.original_same_physical_EEA_I2M2_and_WQ_TP_station_n.sum()==263
fr=fr.rename(columns={'HYBAS_ID_L07':'h','HYBAS_ID_L06':'frh6'});d=d.merge(fr,on='h',how='left',validate='many_to_one');m=d.frh6.notna();assert m.sum()==146 and (d.loc[m,'h6']==d.loc[m,'frh6']).all()
d['TP']=np.where(d.period=='2017-2019',d.mean_original_2018_logTP,d.mean_original_2022_logTP);d['actual_EQR']=np.where(d.period=='2017-2019',d.mean_original_2018_EQR,d.mean_original_2022_EQR);d['log_sites']=np.log1p(d.original_same_physical_EEA_I2M2_and_WQ_TP_station_n)
d['EQR']=d.EQR_original_site_equal_basin_period_mean
# Same H6 fold across all endpoints, never recompute by outcome.
d['fold']=d.h6.map(lambda x:int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5)
d['category']=d.country_holdout+'::'+d.method.replace('','__UNKNOWN__');d['late']=(d.period=='2020-2022').astype(float)
base=['area','riverflow','crop','urban','naturalout','nEEA','nPlant','nAnimal','nButterfly','plant_score_cov','plant_log_n','common_record_cov','common_log_n','late'];d['plant_log_n']=np.log1p(d.plant_score_n)
# Return paired per-H6 fixed heldout losses; independent numerical solver agreement.
store={}
def evaluate(cohort,cutoff,mode,df,y,x,basecols):
 name=f'{cutoff}__{cohort}__{mode}__{x}__{y}';g=df.h6.to_numpy();cat=df.category.to_numpy();yy=df[y].to_numpy(float);aa=df[basecols].to_numpy(float);xx=df[x].to_numpy(float);assert np.isfinite(yy).all() and np.isfinite(aa).all() and np.isfinite(xx).all(),name
 ids=np.arange(len(df)); f=df.fold.to_numpy(); pred0=np.full(len(df),np.nan);pred1=pred0.copy();coefs=[];mx=0
 for k in sorted(set(f)):
  tr=ids[f!=k];te=ids[f==k];assert not set(g[tr])&set(g[te]) and len(set(g[tr]))>=15
  cats=sorted(set(cat[tr]));C=np.column_stack([cat==c for c in cats]).astype(float);ctr=collections.Counter(g[tr]);w=np.array([1/ctr[h] for h in g[tr]]);w/=w.mean()
  for add in [False,True]:
   a1=np.column_stack((aa,xx)) if add else aa;mu=a1[tr].mean(0);sd=a1[tr].std(0);sd[sd<1e-12]=1;X=np.column_stack((np.ones(len(df)),(a1-mu)/sd,C));pen=np.eye(X.shape[1])*10;pen[0,0]=0
   beta=np.linalg.solve(X[tr].T@(w[:,None]*X[tr])+pen,X[tr].T@(w*yy[tr]));p=X[te]@beta;alt=Ridge(alpha=10,fit_intercept=True).fit(X[tr,1:],yy[tr],sample_weight=w).predict(X[te,1:]);mx=max(mx,float(np.max(np.abs(p-alt))))
   if add:pred1[te]=p;coefs.append(float(beta[aa.shape[1]+1]/sd[-1]))
   else:pred0[te]=p
 assert mx<1e-8 and np.isfinite(pred0).all() and np.isfinite(pred1).all()
 hs=sorted(set(g)); l=np.array([[np.mean((yy[g==h]-p[g==h])**2) for p in [pred0,pred1]] for h in hs]);m0,m1=l.mean(0);assert m0>1e-20
 v=l[np.random.default_rng(20261124).integers(0,len(hs),(5000,len(hs)))].mean(1);bs=100*(v[:,0]-v[:,1])/v[:,0];lo,hi=np.quantile(bs,[.025,.975]);gain=100*(m0-m1)/m0;store[name]=(hs,bs,gain)
 result.append(dict(name=name,cutoff=cutoff,cohort=cohort,mode=mode,predictor=x,outcome=y,H7=df.h.nunique(),H6=len(hs),rows=len(df),gain=gain,lo=lo,hi=hi,coef_min=min(coefs),coef_max=max(coefs),independent_solver_max_diff=mx))
 for h,(e0,e1) in zip(hs,l):losses.append(dict(name=name,h6=h,mse0=e0,mse1=e1))
 return name
for cutoff in [3,2]:
 rr=roster[(roster.cutoff==cutoff)&roster.common_eligible].set_index('butterflyKey');keys=sorted(rr.index);mask=np.array([k in keys for k in bk]);assert mask.sum()==len(keys)
 x=B[:,mask];c=x>0;frac=x.sum(1)/B.sum(1);n=c.sum(1);dd=d.copy();dd['common_record_cov']=frac;dd['common_log_n']=np.log1p(n);dd['common_n']=n
 for lab in ['old','clarke']:
  v=rr.loc[np.array(bk)[mask],lab+'_N'].to_numpy(float);dd[lab+'_HNS']=np.divide(c@v,n,out=np.full(len(dd),np.nan),where=n>0);dd[lab+'_HNS_record']=np.divide(x@v,x.sum(1),out=np.full(len(dd),np.nan),where=x.sum(1)>0)
 taxa.append(dict(cutoff=cutoff,common_species=len(keys),score_spearman=spearmanr(rr.old_N,rr.clarke_N).statistic,score_mean_abs_diff=float(np.mean(np.abs(rr.old_N-rr.clarke_N))),median_record_cov=float(np.median(frac)),n_rows_cov50=int((frac>=.5).sum()),median_species_n=float(np.median(n))))
 ok=(frac>=.5)&(n>0)&dd.full_source_trait_eligible.astype(bool).to_numpy();pairs=dd.loc[ok].groupby('h').period.nunique();valid=set(pairs[pairs==2].index);allp=dd[dd.h.isin(valid)].copy();assert len(allp)==2*len(valid)
 for cohort,df in [('EU',allp),('FR73',allp[allp.frh6.notna()].copy())]:
  gate=dict(cutoff=cutoff,cohort=cohort,common_species=len(keys),H7=df.h.nunique(),H6=df.h6.nunique(),rows=len(df),median_record_cov=float(df.common_record_cov.median()),median_plant_cov=float(df.plant_score_cov.median()),median_butterfly_n=float(df.common_n.median()))
  gate['status']='FIT' if len(keys)>=50 and gate['H7']>=40 and gate['H6']>=20 else 'STOP_SUPPORT';gates.append(gate)
  if gate['status']!='FIT':continue
  df=df.reset_index(drop=True);pp=base+(['log_sites'] if cohort=='FR73' else []);names=[]
  if cohort=='FR73':
   df['category']='FR::I2M2';base0=pp+['W'];names+=[evaluate(cohort,cutoff,'spatial_pooled',df,'actual_EQR','TP',base0),evaluate(cohort,cutoff,'spatial_pooled',df,'plant_N','TP',base0)]
   for lab in ['old','clarke']:
    names+=[evaluate(cohort,cutoff,'spatial_pooled',df,lab+'_HNS','TP',base0+['plant_N']),evaluate(cohort,cutoff,'spatial_pooled',df,lab+'_HNS','plant_N',base0+['TP']),evaluate(cohort,cutoff,'record_sensitivity',df,lab+'_HNS_record','plant_N',base0+['TP'])]
   # Fixed first differences, no pseudo replicate and no late outcome in predictors except defined differences.
   early=df[df.late==0].set_index('h').sort_index();late=df[df.late==1].set_index('h').sort_index();assert early.index.equals(late.index)
   dx=early[['h6','fold','category','area','crop','urban','riverflow','naturalout','W','log_sites']].copy().reset_index();change=['TP','actual_EQR','plant_N','old_HNS','clarke_HNS','nEEA','nPlant','nAnimal','nButterfly','plant_score_cov','plant_log_n','common_record_cov','common_log_n']
   for v in change:dx['d_'+v]=late[v].to_numpy()-early[v].to_numpy()
   ctl=['area','crop','urban','riverflow','naturalout','W','log_sites']+['d_'+v for v in change[5:]]
   names+=[evaluate(cohort,cutoff,'delta',dx,'d_actual_EQR','d_TP',ctl),evaluate(cohort,cutoff,'delta',dx,'d_plant_N','d_TP',ctl)]
   for lab in ['old','clarke']:
    names+=[evaluate(cohort,cutoff,'delta',dx,'d_'+lab+'_HNS','d_TP',ctl+['d_plant_N']),evaluate(cohort,cutoff,'delta',dx,'d_'+lab+'_HNS','d_plant_N',ctl+['d_TP'])]
  else:
   names+=[evaluate(cohort,cutoff,'spatial_pooled',df,'plant_N','W',pp),evaluate(cohort,cutoff,'spatial_pooled',df,'EQR','W',pp)]
   for lab in ['old','clarke']:
    names+=[evaluate(cohort,cutoff,'spatial_pooled',df,lab+'_HNS','plant_N',pp+['W']),evaluate(cohort,cutoff,'spatial_pooled',df,lab+'_HNS','EQR',pp+['W','plant_N']),evaluate(cohort,cutoff,'record_sensitivity',df,lab+'_HNS_record','plant_N',pp+['W'])]
  # paired catalogue contrasts on the same heldout loss groups; primary contrast no causal or equivalence interpretation.
  for oldname in names:
   if '__old_HNS' not in oldname:continue
   newname=oldname.replace('__old_HNS','__clarke_HNS');assert newname in names;ho,bo,go=store[oldname];hn,bn,gn=store[newname];assert ho==hn;lo,hi=np.quantile(bn-bo,[.025,.975]);contrast.append(dict(old=oldname,clarke=newname,difference=gn-go,lo=lo,hi=hi))
 allp['cutoff']=cutoff;panels.append(allp[['cutoff','h','h6','period','country','fold','plant_N','old_HNS','clarke_HNS','common_record_cov','common_n','plant_score_cov','EQR','TP','actual_EQR']])
for name,rows in [('gate',gates),('taxa',taxa),('results',result),('losses',losses),('contrasts',contrast)]:pd.DataFrame(rows).to_csv(P/f'filter_{name}.csv',index=False)
pd.concat(panels).to_csv(P/'filter_panel.csv',index=False)
print('GATES\n'+pd.DataFrame(gates).to_string(index=False));print('SCORE AGREEMENT\n'+pd.DataFrame(taxa).to_string(index=False));print('RESULTS\n'+pd.DataFrame(result)[['cutoff','cohort','mode','predictor','outcome','H7','H6','gain','lo','hi','coef_min','coef_max']].round(4).to_string(index=False))
