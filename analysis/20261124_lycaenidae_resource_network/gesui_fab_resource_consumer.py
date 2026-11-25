#!/usr/bin/env python3
"""Exploratory source-locked Lycaenidae catalogue host-group and temporal decomposition.
Run with original six zip archives and previous reconstructed CSV in the same folder.
GBIF records are not organism counts; catalogue exclusivity is not physiological exclusivity.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import hashlib, json, zipfile, pathlib, collections, sys
import numpy as np, pandas as pd
from scipy.stats import hypergeom
def incl(N,C,m):
 N=np.asarray(N,dtype=float);C=np.asarray(C,dtype=float)
 N,C=np.broadcast_arrays(N,C)
 assert np.all(N>=m)
 miss=np.ones_like(C)
 for i in range(m):miss*=np.maximum(0,1-C/(N-i))
 return 1-miss
from sklearn.cluster import KMeans
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from scipy.stats import t
P=pathlib.Path(sys.argv[1] if len(sys.argv)>1 else '/mnt/data')
D={
 'gesui_original_panel.zip':'b0babfa63a96f9b7eb5a83122b08b8664edab6e8c6c28d2faf4ff024ca84c39b',
 'gesui_host_source.zip':'1cf6ed4ee41282a7f1be2e5cc7e1542a186b4fef5e1a46adf5861543c0ebbd94',
 'gesui_gbif_facets.zip':'4b37c741d1bd28c178b0cef61c9457af0dd30797c173f9384ebe1d1d0036d07c',
 'gesui_atlas.zip':'99f42768ac0d215fca5662de449f97fba598bebab9b147288ef99e19db3633b7',
 'gesui_h6.zip':'1f94449184095225de2048bc916754ae2072ed17f385614ddb7c67bd5d52dbdd',
 'FR73_original_14013312711.zip':'7c3b9f1976e410526e396fabdfbd527319bdde5047ca4819ee6c29588443e5e2'}
for f,v in D.items(): assert hashlib.sha256((P/f).read_bytes()).hexdigest()==v,f
Z={f:zipfile.ZipFile(P/f) for f in D}
def j(f,m): return json.loads(Z[f].read(m))
H='gesui_host_source.zip';G='gesui_gbif_facets.zip'
raw=j(G,'GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json');assert len(raw)==2808 and all(x['source_status']=='ok' for x in raw)
fac={(x['HYBAS_ID_L07'],x['period'],x['group']):x for x in raw};assert len(fac)==2808
alias=j(H,'butterfly_usage_aliases.json');edges=j(H,'verified_candidate_host_edges.json');cross=j(H,'plant_name_crosswalk.json')
assert len(alias)==303 and len(edges)==1565
bk=sorted(x['original_speciesKey'] for x in alias);pk=sorted(set(x['plantKey'] for x in edges));assert len(pk)==990
bidx={k:i for i,k in enumerate(bk)};pidx={k:i for i,k in enumerate(pk)}
bmeta={x['original_speciesKey']:x for x in alias};pf={x['accepted_speciesKey']:x['match'].get('family','UNKNOWN') for x in cross if x['accepted_speciesKey'] is not None};pn={x['accepted_speciesKey']:x['match'].get('canonicalName',x['source_name']) for x in cross if x['accepted_speciesKey'] is not None}
A=np.zeros((303,990));
for x in edges:A[bidx[x['butterflyKey']],pidx[x['plantKey']]]=1
assert A.sum()==1565
bf=np.array([bmeta[k]['family'] for k in bk]);L=bf=='Lycaenidae';fab=np.array([pf[k]=='Fabaceae' for k in pk]);assert L.sum()==75
# Restore original geography independently and compare all governing fields to previous result.
d=pd.read_csv(Z['gesui_original_panel.zip'].open('EU234_H7_ORIGINAL_STUDY_NATIVE_BASINPERIOD_EQR_PLANT_BUTTERFLY_SOURCE_COHORT.csv'),dtype={'HYBAS_ID_L07':str},keep_default_na=False).rename(columns={'HYBAS_ID_L07':'h','PERIOD':'period','original_majority_BQE_method':'method','W_static_upstream_MODEL_flow_per_UpArea':'W','log_local_H7_area_km2':'area','log1p_original_EQR_site_n':'nEEA','log1p_GBIF_Animalia_source_records':'nAnimal','log1p_GBIF_Plantae_source_records':'nPlant','log1p_GBIF_5family_butterfly_source_records':'nButterfly','log1p_original_EEA_site_DIS_AV_CMS':'riverflow'})
mp={str(x['HYBAS_ID_L07']):str(x['official_parent_HYBAS_ID_L06']) for x in j('gesui_h6.zip','HYBAS7_232_TO_OFFICIAL_HYBAS6_SOURCE_PARENT_EXPLICIT_REPAIR_v2.json')['original_454_H7_study_basin_geographic_parent_mapping']}
at=j('gesui_atlas.zip','HYDROATLAS_V1_L07_1219_H7_BASIN_RAW_LANDUSE_ATTRIBUTES_SOURCE_ONLY.json')['original_1219_H7_covariate_by_real_ID']
d['h6']=d.h.map(mp);d['crop']=[at[h]['crp_pc_use'] for h in d.h];d['urban']=[at[h]['urb_pc_use'] for h in d.h];d['naturalout']=[np.log1p(float(at[h]['dis_m3_pyr'])) for h in d.h];d['country_holdout']=d.h.map({h:sorted(g.country)[0] for h,g in d.groupby('h')})
assert (len(d),d.h.nunique(),d.h6.nunique())==(454,232,90)
prev=pd.read_csv(P/'gesui_rebuilt_network.csv',dtype={'h':str,'h6':str},keep_default_na=False);assert len(prev)==454
prev=prev.set_index(['h','period']).loc[pd.MultiIndex.from_frame(d[['h','period']])].reset_index()
for k in ['h','h6','period','country','method','country_holdout']: assert (d[k].to_numpy()==prev[k].to_numpy()).all(),k
K=['W','area','riverflow','crop','urban','naturalout','nEEA','nAnimal','nPlant','nButterfly','plant_S50','butterfly_S20','catalogue_coverage','minimum_species_coverage']
for k in K[:10]:np.testing.assert_allclose(d[k],prev[k],atol=1e-10,rtol=0)
B=np.zeros((454,303),int);Pcount=np.zeros((454,990),int);NB=[];NP=[];PS=[];cov=[]
for i,x in enumerate(d.itertuples()):
 fs=[fac[x.h,x.period,f] for f in sorted(set(bf))];p=fac[x.h,x.period,'Plantae'];nb=0
 for f in fs:
  for k,v in f['source_speciesKey_counts'].items():B[i,bidx[k]]=v;nb+=v
 for k,v in p['source_speciesKey_counts'].items():
  if k in pidx:Pcount[i,pidx[k]]=v
 np_=sum(p['source_speciesKey_counts'].values());NB.append(nb);NP.append(np_);PS.append(incl(np_,list(p['source_speciesKey_counts'].values()),50).sum());cov.append(min(x['source_facet_coverage'] or 0 for x in fs+[p]))
PB=incl(np.array(NB)[:,None],B,20);PP=incl(np.array(NP)[:,None],Pcount,50)
known=A.sum(1)>0;d['plant_S50']=PS;d['butterfly_S20']=PB.sum(1);d['catalogue_coverage']=PB[:,known].sum(1)/PB.sum(1);d['minimum_species_coverage']=cov
for k in K[10:]:np.testing.assert_allclose(d[k],prev[k],atol=1e-8,rtol=0)
# Recreate expected shuffled graph as exact average adjacency under the original 499 permutations.
allfreq=collections.Counter();
for x in raw:
 if x['group']=='Plantae':allfreq.update(x['source_speciesKey_counts'].keys())
bins=np.digitize([allfreq[k] for k in pk],[0,1,5,20,100,468],right=True)
hs=sorted(d.h6.unique());inc=np.array([(Pcount[d.h6==h]>0).any(0) for h in hs]).T.astype(float);norm=np.linalg.norm(inc,axis=1);lab=np.full(len(pk),-1);mask=norm>0;lab[mask]=KMeans(n_clusters=6,random_state=20261114,n_init=10).fit_predict(inc[mask]/norm[mask,None]);strata=np.array([pf[k]+'|'+str(b)+'|geo'+str(c) for k,b,c in zip(pk,bins,lab)])
rng=np.random.default_rng(20261114);slices=[np.flatnonzero(strata==s) for s in sorted(set(strata))];N=np.zeros_like(A)
for _ in range(499):
 perm=np.arange(len(pk));
 for ix in slices:perm[ix]=rng.permutation(ix)
 N+=A[:,perm]
N/=499;Dmat=A-N
np.testing.assert_allclose(A.sum(1),N.sum(1),atol=1e-9)
for mask,label in [(np.ones(303,bool),'all_all0'),(L,'L_all0')]:
 o=np.sum((PB[:,mask]@A[mask])*PP,axis=1);n=np.sum((PB[:,mask]@N[mask])*PP,axis=1)
 np.testing.assert_allclose(o,prev[label+'_obs'],atol=1e-8,rtol=0);np.testing.assert_allclose(n,prev[label+'_null'],atol=1e-8,rtol=0)
print('454 network replay passed',flush=True)
# All catalogue-mapped Lycaenidae; categories fixed by entire resolved host list, not field outcomes.
fa=A[:,fab].sum(1);nf=A[:,~fab].sum(1)
GM={'fab_only':L&(fa>0)&(nf==0),'mixed':L&(fa>0)&(nf>0),'nonfab_only':L&(fa==0)&(nf>0),'unknown':L&(~known)}
assert np.sum(list(GM.values()),axis=0).sum()==75
assert np.all(np.sum(list(GM.values()),axis=0)==L)
report=[]
w=1/d.groupby('h6').h6.transform('size').to_numpy();w/=w.sum()
for name,g in GM.items():
 d[name+'_mass']=PB[:,g].sum(1)
 d[name+'_nrecord']=B[:,g].sum(1)
 for sub,pcols in [('all',np.ones(990,bool)),('fab',fab),('nonfab',~fab)]:
  o=((PB[:,g]@A[g][:,pcols])*PP[:,pcols]).sum(1);n=((PB[:,g]@N[g][:,pcols])*PP[:,pcols]).sum(1)
  d[name+'_'+sub+'_obs']=o;d[name+'_'+sub+'_null']=n;d[name+'_'+sub+'_ex']=o-n
  report.append(dict(group=name,plant_set=sub,taxa=int(g.sum()),records=int(B[:,g].sum()),record_share=float(B[:,g].sum()/B[:,L].sum()),host_links=int(A[g][:,pcols].sum()),obs=float(w@o),null=float(w@n),ex=float(w@(o-n)),excess_pct=float(100*(w@(o-n))/(w@n)) if w@n>0 else None))
for k in ['obs','null','ex']:
 for sub,ref in [('all','L_all0'),('fab','L_fab'),('nonfab','L_nonfab')]:np.testing.assert_allclose(sum(d[g+'_'+sub+'_'+k] for g in GM),prev[ref+'_'+k],atol=1e-8,rtol=0)
d['L_ex']=prev.L_all0_ex;d['Fab_ex']=prev.L_fab_ex;d['NonFab_ex']=prev.L_nonfab_ex
# Known-catalogue host plant marginals, not all legume species (only 990 mapped plants).
d['fab_plant_mass']=PP[:,fab].sum(1);d['all_host_plant_mass']=PP.sum(1);d['fab_plant_share']=d.fab_plant_mass/d.all_host_plant_mass
d['lyc_mass']=PB[:,L].sum(1);d['known_lyc_mass']=PB[:,L&known].sum(1)
for g in GM: d[g+'_share']=d[g+'_mass']/d.lyc_mass.where(d.lyc_mass>0)
print(pd.DataFrame(report).query('plant_set=="all"').to_string(index=False),flush=True)
# Same H7 pairs; calculate exact temporal decomposition of the source network excess.
a=d.index[d.period=='2017-2019'];b=d.index[d.period=='2020-2022'];aa=d.loc[a,['h']].assign(a=a).merge(d.loc[b,['h']].assign(b=b),on='h',validate='one_to_one');a=aa.a.to_numpy();b=aa.b.to_numpy();assert len(a)==222 and np.array_equal(d.h6.iloc[a].to_numpy(),d.h6.iloc[b].to_numpy())
q=d.iloc[a].copy().reset_index(drop=True);later=d.iloc[b].copy().reset_index(drop=True);q['i0']=a;q['i1']=b
quality=(q.minimum_species_coverage>=.8)&(later.minimum_species_coverage>=.8);assert (quality.sum(),q.loc[quality,'h6'].nunique())==(173,77)
# Split across consumer-group and resource-family. No full-sample reference needed: exact two-period symmetric decomposition.
decomp={}
for name,g,cols in [('L',L,np.ones(990,bool)),('Fab',L,fab),('NonFab',L,~fab),('fab_only',GM['fab_only'],np.ones(990,bool)),('mixed',GM['mixed'],np.ones(990,bool)),('nonfab_only',GM['nonfab_only'],np.ones(990,bool))]:
 m=Dmat[g][:,cols];ba=PB[a][:,g];bb=PB[b][:,g];pa=PP[a][:,cols];pb=PP[b][:,cols]
 val=lambda x,y:((x@m)*y).sum(1)
 e00,e01,e10,e11=val(ba,pa),val(ba,pb),val(bb,pa),val(bb,pb)
 dp=.5*((e01-e00)+(e11-e10));dc=.5*((e10-e00)+(e11-e01));np.testing.assert_allclose(dp+dc,e11-e00,atol=1e-12)
 for k,v in [('delta',e11-e00),('plant_delta',dp),('butterfly_delta',dc),('early',e00),('late',e11),('plant_change_only',e01),('butterfly_change_only',e10)]:q[name+'_'+k]=v
 decomp[name]={'delta_sd':float(np.std(e11-e00)),'plant_sd':float(np.std(dp)),'consumer_sd':float(np.std(dc)),'identity_error':float(np.max(np.abs(dp+dc-e11+e00)))}
# attach measured water summary (same station set) source check.
f='FR73_original_14013312711.zip';member='FR73_HYBAS7_NATIVE_2018_2022_REAL_TP_I2M2_EQR_AND_GBIF_INDEX_DELTA.csv';bb=Z[f].read(member);assert hashlib.sha256(bb).hexdigest()=='f2ff4e56f245e1aa8f84f5ac099d68ec1c7e27a786c1067439db388b186dd90b'
fr=pd.read_csv(Z[f].open(member),dtype={'HYBAS_ID_L07':str,'HYBAS_ID_L06':str});assert len(fr)==73 and fr.HYBAS_ID_L06.nunique()==24
fr=fr.rename(columns={'HYBAS_ID_L07':'h','HYBAS_ID_L06':'frh6'}).set_index('h');q=q.join(fr[['frh6','mean_original_2018_logTP','mean_original_2022_logTP','mean_original_2018_EQR','mean_original_2022_EQR','original_same_physical_EEA_I2M2_and_WQ_TP_station_n']],on='h');frmask=q.frh6.notna();assert frmask.sum()==73 and (q.loc[frmask,'h6']==q.loc[frmask,'frh6']).all()
E='EQR_original_site_equal_basin_period_mean'
# fixed prior H6 hash folds, preserve assignments when including all rows.
def hsh(s):
 h=2166136261
 for c in s:h=((h^ord(c))*16777619)&0xffffffff
 return h
hg=sorted(set(q.loc[quality,'h6']),key=lambda s:(hsh(s),s));gall=hg+sorted(set(q.h6)-set(hg),key=lambda s:(hsh(s),s));fold=q.h6.map({g:i%5 for i,g in enumerate(gall)}).to_numpy()
results=[];lossstore={}
def evaluate(name,mask,target,base,predictor,foldtype='five'):
 ids=np.where(mask & np.isfinite(target) & np.isfinite(predictor) & np.isfinite(base).all(1))[0];gs=q.h6.to_numpy()[ids];yy=target[ids];cat=(q.country_holdout+'::'+q.method.replace('','__UNKNOWN__')).to_numpy()[ids];X=base[ids];xx=predictor[ids];splits=fold[ids] if foldtype=='five' else gs
 if len(ids)<30 or len(set(gs))<15:return {'name':name,'n':len(ids),'H6':len(set(gs)),'status':'STOP_SUPPORT'}
 a0=np.full(len(ids),np.nan);a1=a0.copy();coef=[];mx=0
 for f in sorted(set(splits)):
  test=splits==f;train=~test;assert not set(gs[test])&set(gs[train]);ctr=collections.Counter(gs[train]);ww=np.array([1/ctr[g] for g in gs[train]]);ww/=ww.mean()
  cats=sorted(set(cat[train]));C=np.array([cat==c for c in cats]).T.astype(float)
  for add in [0,1]:
   x=np.column_stack((X,xx)) if add else X;mu=x[train].mean(0);sd=x[train].std(0);sd[sd<1e-12]=1;x=np.column_stack(((x-mu)/sd,C));x=np.column_stack((np.ones(len(ids)),x));pen=np.eye(x.shape[1])*10;pen[0,0]=0
   c=np.linalg.solve(x[train].T@(ww[:,None]*x[train])+pen,x[train].T@(ww*yy[train]));pred=x[test]@c
   independent=Ridge(alpha=10,fit_intercept=True).fit(x[train,1:],yy[train],sample_weight=ww).predict(x[test,1:]);mx=max(mx,float(np.max(np.abs(pred-independent))))
   if add:a1[test]=pred;coef.append(c[len(mu)]/sd[-1])
   else:a0[test]=pred
 assert mx<1e-8
 gs2=sorted(set(gs));l=np.array([[np.mean((yy[gs==g]-p[gs==g])**2) for p in [a0,a1]] for g in gs2]);e0,e1=l.mean(0);idx=np.random.default_rng(20261124).integers(0,len(l),(5000,len(l)));v=l[idx].mean(1);bs=100*(v[:,0]-v[:,1])/v[:,0];gain=100*(e0-e1)/e0;lo,hi=np.quantile(bs,[.025,.975]);out=dict(name=name,n=len(ids),H6=len(gs2),gain=float(gain),lo=float(lo),hi=float(hi),coef_min=float(min(coef)),coef_max=float(max(coef)),groups_improved=int((l[:,1]<l[:,0]).sum()),mse0=float(e0),mse1=float(e1),independent_prediction_diff=mx,foldtype=foldtype,status='FIT')
 lossstore[name]=(gs2,l,bs);results.append(out);return out
base=q[K].to_numpy(float);ee=q[E].to_numpy(float)
# fixed outcomes: host category-specific excess, resource/consumer margins, and exact delta decomposition.
allmask=np.ones(len(q),bool);nonfr=(q.country!='FR')&(later.country!='FR');frboth=(q.country=='FR')&(later.country=='FR')
for tag,mask in [('quality',quality),('all',allmask),('France',quality&frboth),('other',quality&nonfr)]:
 print('model',tag,flush=True)
 for k in ['L_ex','Fab_ex','NonFab_ex','fab_only_all_ex','mixed_all_ex','nonfab_only_all_ex','fab_plant_mass','fab_plant_share','fab_only_share','mixed_share','nonfab_only_share']:
  evaluate(tag+'__'+k,mask,later[k].to_numpy(float),np.column_stack((base,q[k].to_numpy(float))),ee)
 # Delta decomposition models share the same controls and previous excess for each network component.
 for k in (['L','Fab','NonFab','fab_only','mixed','nonfab_only'] if tag in ['quality','all'] else ['L','Fab']):
  for side in ['delta','plant_delta','butterfly_delta']:
   evaluate(tag+'__'+k+'_'+side,mask,q[k+'_'+side].to_numpy(),np.column_stack((base,q[k+'_early'].to_numpy())),ee)
# Controls on marginal plant resource and consumer composition; no later outcome-derived controls in forecast.
for name,ks in [('plant_only',['fab_plant_mass','all_host_plant_mass']),('consumer_only',['known_lyc_mass','fab_only_mass','mixed_mass','nonfab_only_mass']),('both',['fab_plant_mass','all_host_plant_mass','known_lyc_mass','fab_only_mass','mixed_mass','nonfab_only_mass'])]:
 for k in ['L_ex','Fab_ex','fab_only_all_ex','mixed_all_ex']:
  evaluate('quality_margins_'+name+'__'+k,quality,later[k].to_numpy(float),np.column_stack((base,q[k].to_numpy(),q[ks].to_numpy())),ee)
# France same physical station stricter exposure; whole 73, and 80% quality subset if gate passes.
for tag,mask in [('FR73',frmask),('FR73_quality',frmask&quality)]:
 print('model',tag,flush=True)
 for pred in ['mean_original_2018_logTP','mean_original_2018_EQR']:
  for k in ['L_ex','Fab_ex','fab_only_all_ex','mixed_all_ex','nonfab_only_all_ex','fab_plant_mass','fab_only_share']:
   evaluate(tag+'__'+pred+'__'+k,mask,later[k].to_numpy(float),np.column_stack((base,q[k].to_numpy(float))),q[pred].to_numpy(float),'loho')
# Save complete results; descriptive numbers not confirmatory p-values.
pd.DataFrame(report).to_csv(P/'gesui_fab_consumer_group_summary.csv',index=False)
pd.DataFrame(results).to_csv(P/'gesui_fab_consumer_model_results.csv',index=False)
cols=['h','h6','country','period','i0','i1']+[k for k in q if any(k.endswith('_'+s) for s in ['delta','early','late','change_only'])]
q[cols].to_csv(P/'gesui_fab_temporal_decomposition.csv',index=False)
pd.DataFrame([dict(butterfly_key=k,name=bmeta[k]['canonicalName'],status=bmeta[k]['status'],catalogue_group=next(g for g,m in GM.items() if m[i]),fab_hosts=int(fa[i]),nonfab_hosts=int(nf[i]),GBIF_records=int(B[:,i].sum())) for i,k in enumerate(bk) if L[i]]).to_csv(P/'gesui_lycaenidae_catalogue_groups.csv',index=False)
(P/'gesui_fab_consumer_receipt.json').write_text(json.dumps(dict(sources=D,primary_H7=int(quality.sum()),primary_H6=q.loc[quality,'h6'].nunique(),paired_H7=len(q),FR73=int(frmask.sum()),network_reproduction='assertion PASS for all454 totals, Lycaenidae totals and 3 group sums',delta_identity=decomp,model_count=len(results),max_independent_prediction_diff=max(x['independent_prediction_diff'] for x in results),note='source-locked exploratory posthoc; fixed-prediction H6 bootstrap not independent field replication; catalogue only, not feeding; static full-panel geographic null not prospective'),indent=2))
print(pd.DataFrame(report).query('plant_set=="all"').to_string(index=False))
r=pd.DataFrame(results);print(r[r.name.str.startswith('quality__')][['name','n','H6','gain','lo','hi','coef_min','coef_max']].round(4).to_string(index=False));print(r[r.name.str.startswith('quality_margins')][['name','gain','lo','hi']].round(3).to_string(index=False));print(r[r.name.str.startswith('FR73__')][['name','n','gain','lo','hi']].round(3).to_string(index=False))
# Paired fixed-loss bootstrap differences in relative prediction improvement.
# Posthoc diagnostics selected after seeing both model estimates; not confirmatory tests.
comparisons=[]
for tag in ['quality','all','France','other']:
 for a0,b0 in [('L_plant_delta','L_butterfly_delta'),('Fab_plant_delta','Fab_butterfly_delta'),('mixed_all_ex','fab_only_all_ex'),('mixed_all_ex','nonfab_only_all_ex')]:
  a=tag+'__'+a0;b=tag+'__'+b0
  if a not in lossstore or b not in lossstore:continue
  ga,la,ba=lossstore[a];gb,lb,bb=lossstore[b];assert ga==gb
  est=lambda l:100*(l[:,0].mean()-l[:,1].mean())/l[:,0].mean()
  diff=ba-bb;lo,hi=np.quantile(diff,[.025,.975]);lo2,hi2=np.quantile(diff,[.00625,.99375]);comparisons.append(dict(sample=tag,A=a0,B=b0,H6=len(ga),diff=est(la)-est(lb),lo=lo,hi=hi,bonferroni_4_lo=lo2,bonferroni_4_hi=hi2))
pd.DataFrame(comparisons).to_csv(P/'gesui_fab_paired_contrasts.csv',index=False)
print(pd.DataFrame(comparisons).round(4).to_string(index=False))
