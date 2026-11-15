"""Offline source-locked H7 candidate host opportunities; exploratory only."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import collections,hashlib,json,pathlib,sys,zipfile,math
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from network import inclusion,metrics,shuffle_columns
ROOT=pathlib.Path(__file__).resolve().parents[1];SRC=ROOT/'source';INPUT=pathlib.Path(os.environ.get('GESUI_INPUT',str(ROOT/'inputs')));OUT=ROOT/os.environ.get('GESUI_OUTPUT','output')
FAMS=['Lycaenidae','Nymphalidae','Pieridae','Papilionidae','Hesperiidae']

def zjson(name,member):
    with zipfile.ZipFile(INPUT/name) as z:return json.loads(z.read(member))

def load_panel():
    with zipfile.ZipFile(INPUT/'gesui_h7_454_original.zip') as z:
        d=pd.read_csv(z.open('EU234_H7_ORIGINAL_STUDY_NATIVE_BASINPERIOD_EQR_PLANT_BUTTERFLY_SOURCE_COHORT.csv'),dtype={'HYBAS_ID_L07':str},keep_default_na=False)
    ren={'HYBAS_ID_L07':'h','PERIOD':'period','original_majority_BQE_method':'method','W_static_upstream_MODEL_flow_per_UpArea':'W','log_local_H7_area_km2':'area','log1p_original_EQR_site_n':'nEEA','log1p_GBIF_Animalia_source_records':'nAnimal','log1p_GBIF_Plantae_source_records':'nPlant','log1p_GBIF_5family_butterfly_source_records':'nButterfly','log1p_original_EEA_site_DIS_AV_CMS':'riverflow'}
    d=d.rename(columns=ren)
    h6=zjson('gesui_h7_h6_parent_original.zip','HYBAS7_232_TO_OFFICIAL_HYBAS6_SOURCE_PARENT_EXPLICIT_REPAIR_v2.json')
    mp={str(x['HYBAS_ID_L07']):str(x['official_parent_HYBAS_ID_L06']) for x in h6['original_454_H7_study_basin_geographic_parent_mapping']}
    at=zjson('gesui_h7_atlas_original.zip','HYDROATLAS_V1_L07_1219_H7_BASIN_RAW_LANDUSE_ATTRIBUTES_SOURCE_ONLY.json')['original_1219_H7_covariate_by_real_ID']
    d['h6']=d.h.map(mp);d['crop']=[at[h]['crp_pc_use'] for h in d.h];d['urban']=[at[h]['urb_pc_use'] for h in d.h];d['naturalout']=[math.log1p(float(at[h]['dis_m3_pyr'])) for h in d.h]
    assert len(d)==454 and d.h.nunique()==232 and d.h6.nunique()==90
    # Check previously delivered table, but rebuild all governing fields from original archives.
    prev=pd.read_csv(INPUT/'gesui_H7_posthoc_S20_source.csv',dtype={'h':str,'h6':str},keep_default_na=False)
    for c in ['h','period','h6','country','method']:assert (d[c].to_numpy()==prev[c].to_numpy()).all(),c
    for c in ['W','crop','urban','naturalout','area','nEEA','nAnimal','nPlant','nButterfly','riverflow']:assert np.allclose(d[c],prev[c],atol=1e-10),c
    stable={h:sorted(g.country)[0] for h,g in d.groupby('h')}
    d['country_holdout']=d.h.map(stable)
    return d

def design(d,train,test,extra=(),country=True):
    num=['area','nEEA','nAnimal','nPlant','nButterfly','riverflow','crop','urban','naturalout','butterfly_S20','plant_S50','catalogue_coverage']+list(extra)
    a=d[num].to_numpy(float);mu=a[train].mean(0);sd=a[train].std(0);sd[sd<1e-12]=1
    X=[np.ones((len(d),1)),(a-mu)/sd]
    for col in ['period','method']+(['country'] if country else []):
        cats=sorted(set(d.loc[train,col]));X.append(np.array([[float(x==c) for c in cats] for x in d[col]]))
    X=np.column_stack(X)
    return X[train],X[test]

def predict(d,y,extra,mode='H6',independent=False):
    if mode=='H6':groups=np.array([int(hashlib.sha256(h.encode()).hexdigest()[:8],16)%5 for h in d.h6])
    else:
        groups=d.country_holdout.to_numpy()
    pred=np.full(len(d),np.nan);splits=[]
    for g in sorted(set(groups)):
        test=groups==g;train=~test
        if mode=='country':train &= ~d.h6.isin(set(d.loc[test,'h6'])).to_numpy()
        assert not set(d.loc[train,'h6'])&set(d.loc[test,'h6'])
        if train.sum()<30 or d.loc[train,'h6'].nunique()<5:raise RuntimeError('Insufficient training after independent geography purge')
        X,T=design(d,train,test,extra,mode=='H6')
        f=collections.Counter(d.loc[train,'h6']);w=np.array([1/f[h] for h in d.loc[train,'h6']]);w/=w.mean()
        if independent:
            from sklearn.linear_model import Ridge
            # Explicit unpenalized intercept implemented by weighted centering.
            xm=np.average(X[:,1:],axis=0,weights=w);ym=np.average(y[train],weights=w)
            reg=Ridge(alpha=10.,fit_intercept=False,solver='cholesky').fit(X[:,1:]-xm,y[train]-ym,sample_weight=w)
            pred[test]=reg.predict(T[:,1:]-xm)+ym
        else:
            pen=np.eye(X.shape[1])*10;pen[0,0]=0
            coef=np.linalg.solve(X.T@(w[:,None]*X)+pen,X.T@(w*y[train]));pred[test]=T@coef
        splits.append({'group':str(g),'train_rows':int(train.sum()),'train_H6':int(d.loc[train,'h6'].nunique()),'test_rows':int(test.sum()),'zero_shared_H6':True})
    assert np.isfinite(pred).all()
    return pred,splits

def model(d,col,log=False):
    d=d.reset_index(drop=True);y=d[col].to_numpy(float);y=np.log1p(y) if log else y
    a,s=predict(d,y,());b,_=predict(d,y,('W',));qa,_=predict(d,y,('W',),independent=True)
    err0=(y-a)**2;err1=(y-b)**2;hs=sorted(set(d.h6));e0=np.array([err0[d.h6==h].mean() for h in hs]);e1=np.array([err1[d.h6==h].mean() for h in hs])
    rng=np.random.default_rng(20261114);idx=rng.integers(0,len(hs),(3000,len(hs)))
    bs=100*(e0[idx].mean(1)-e1[idx].mean(1))/e0[idx].mean(1)
    ac,sc=predict(d,y,(),mode='country');bc,_=predict(d,y,('W',),mode='country')
    cs=[]
    cg=d.country_holdout
    for c in sorted(set(cg)):
        m=cg==c;eo=float(np.mean((y[m]-ac[m])**2));en=float(np.mean((y[m]-bc[m])**2))
        cs.append({'country':c,'n':int(m.sum()),'H6':int(d.loc[m,'h6'].nunique()),'baseline_mse':eo,'with_W_mse':en,'gain_pct':100*(eo-en)/eo})
    valid=[c for c in cs if c['n']>=10 and c['H6']>=3]
    ce0=np.mean([c['baseline_mse'] for c in valid]);ce1=np.mean([c['with_W_mse'] for c in valid])
    # Within-sample descriptive coefficient, not a causal dose response.
    X,_=design(d,np.ones(len(d),bool),np.zeros(len(d),bool),('W',),True)
    beta=np.linalg.lstsq(X,y,rcond=None)[0];res=y-X@beta;bread=np.linalg.pinv(X.T@X);meat=np.zeros_like(bread)
    for h in hs:
        z=X[d.h6==h].T@res[d.h6==h];meat+=np.outer(z,z)
    rank=np.linalg.matrix_rank(X);var=bread@meat@bread*len(hs)/(len(hs)-1)*(len(d)-1)/(len(d)-rank)
    wi=13 # intercept + 12 baseline numeric fields, W next
    se=float(np.sqrt(max(0,var[wi,wi])))
    return {'n_rows':len(d),'H7':int(d.h.nunique()),'H6':len(hs),'H6_gain_pct':float(100*(e0.mean()-e1.mean())/e0.mean()),'H6_fixed_error_bootstrap95_pct':np.quantile(bs,[.025,.975]).tolist(),'country_equal_gain_pct':float(100*(ce0-ce1)/ce0),'country_row_weighted_gain_pct':float(100*(np.mean((y-ac)**2)-np.mean((y-bc)**2))/np.mean((y-ac)**2)),'eligible_country_n':len(valid),'country_results':cs,'H6_splits':s,'country_splits':sc,'independent_sklearn_prediction_max_diff':float(np.max(np.abs(b-qa))),'W_beta_perSD_descriptive':float(beta[wi]),'W_beta_H6_cluster95_descriptive':[float(beta[wi]-1.96*se),float(beta[wi]+1.96*se)]}

def main():
    OUT.mkdir(exist_ok=True)
    d=load_panel();raw=zjson('gesui_h7_234_complete_source429_repaired.zip','GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json')
    rr={(r['HYBAS_ID_L07'],r['period'],r['group']):r for r in raw}
    taxa=json.loads((SRC/'butterfly_usage_aliases.json').read_text());edges=json.loads((SRC/'verified_candidate_host_edges.json').read_text());plants=json.loads((SRC/'plant_name_crosswalk.json').read_text())
    bkeys=sorted({t['original_speciesKey'] for t in taxa});pkeys=sorted({e['plantKey'] for e in edges});bi={k:i for i,k in enumerate(bkeys)};pi={k:i for i,k in enumerate(pkeys)}
    A=np.zeros((len(bkeys),len(pkeys)))
    taxstatus={t['original_speciesKey']:t['status'] for t in taxa}
    for e in edges:
        if os.environ.get('ACCEPTED_ONLY')=='1' and taxstatus[e['butterflyKey']]!='ACCEPTED':continue
        if os.environ.get('NO_WIKIPEDIA')=='1' and not any(x!='Wikipedia' for x in e['source_labels']):continue
        A[bi[e['butterflyKey']],pi[e['plantKey']]]=1
    B=np.zeros((len(d),len(bkeys)),int);P=np.zeros((len(d),len(pkeys)),int);NB=[];NP=[];cov=[];PS=[]
    for row,x in enumerate(d.itertuples()):
        allb={}
        for f in FAMS:
            for k,v in rr[(x.h,x.period,f)]['source_speciesKey_counts'].items():
                if k in allb:raise RuntimeError('Same speciesKey in >1 family')
                allb[k]=v;B[row,bi[k]]=v
        allp=rr[(x.h,x.period,'Plantae')]['source_speciesKey_counts']
        for k,v in allp.items():
            if k in pi:P[row,pi[k]]=v
        nb=sum(allb.values());np_=sum(allp.values());NB.append(nb);NP.append(np_)
        cov.append(min(rr[(x.h,x.period,g)]['source_facet_coverage'] or 0 for g in FAMS+['Plantae']))
        PS.append(float(inclusion(np_,list(allp.values()),50).sum()))
    NB=np.array(NB);NP=np.array(NP);assert min(NB)>=20 and min(NP)>=50
    PB=np.array([inclusion(n,b,20) for n,b in zip(NB,B)]);PP=np.array([inclusion(n,p,50) for n,p in zip(NP,P)])
    known=A.sum(1)>0;den=PB[:,known].sum(1);hc=P@A.T;Q=np.array([inclusion(n,c,50) for n,c in zip(NP,hc)])
    d['expected_links']=np.sum((PB@A)*PP,axis=1);d['host_support']=np.sum(PB[:,known]*Q[:,known],axis=1)/den
    d['catalogue_coverage']=den/PB.sum(1);d['butterfly_S20']=PB.sum(1);d['plant_S50']=PS;d['minimum_species_coverage']=cov
    d['candidate_links']=np.sum(((B>0)@A)*(P>0),axis=1);d['mean_host_degree']=PB@A.sum(1)/den
    families={t['original_speciesKey']:t['family'] for t in taxa};fa=np.array([families[k] for k in bkeys]);family_summary=[];spec=[]
    for f in FAMS:
        m=fa==f;km=m&known;fd=PB[:,km].sum(1)
        vals=np.divide(np.sum(PB[:,km]*Q[:,km],axis=1),fd,out=np.full(len(d),np.nan),where=fd>0)
        d[f+'_host_support']=vals
        family_summary.append({'family':f,'recorded_speciesKeys':int(m.sum()),'host_mapped_speciesKeys':int(km.sum()),'catalogue_species_coverage':float(km.sum()/m.sum()),'source_record_coverage':float(B[:,km].sum()/B[:,m].sum()),'mean_host_support':float(np.nanmean(vals)),'median_host_support':float(np.nanmedian(vals))})
    for i,k in enumerate(bkeys):spec.append({'butterflyKey':k,'canonicalName':next(t['canonicalName'] for t in taxa if t['original_speciesKey']==k),'family':fa[i],'catalogue_host_keys':int(A[i].sum()),'period_rows_recorded':int((B[:,i]>0).sum()),'period_rows_any_host':int(((B[:,i]>0)&(hc[:,i]>0)).sum()),'source_records':int(B[:,i].sum())})
    # Catalogue rewiring retains every butterfly degree, host family, and frequency-bin margin.
    p_fam={p['accepted_speciesKey']:p['match'].get('family','UNKNOWN') for p in plants if p['accepted_speciesKey'] is not None}
    allfreq=collections.Counter()
    for r in raw:
        if r['group']=='Plantae':allfreq.update(r['source_speciesKey_counts'].keys())
    freq=np.array([allfreq[k] for k in pkeys]);bounds=[0,1,5,20,100,468];bins=np.digitize(freq,bounds,right=True)
    strata=np.array([p_fam[k]+'|'+str(v) for k,v in zip(pkeys,bins)]);group_sizes=collections.Counter(strata)
    rng=np.random.default_rng(20261114);null=np.empty((499,len(d)));last=None;BA=PB@A
    for it in range(len(null)):
        perm=shuffle_columns(strata,rng);Ap=A[:,perm]
        assert np.array_equal(Ap.sum(1),A.sum(1))
        null[it]=np.sum(BA[:,perm]*PP,axis=1);last=perm
    mu=null.mean(0);sd=null.std(0,ddof=1);d['null_expected_links']=mu;d['null_link_excess']=d.expected_links-mu;d['null_link_z']=np.divide(d.expected_links-mu,sd,out=np.full(len(d),np.nan),where=sd>1e-10)
    h6weights=1/d.groupby('h6').h6.transform('size').to_numpy();h6weights/=h6weights.sum()
    observed=float(h6weights@d.expected_links);nullavg=null@h6weights
    nullqa={'draws':499,'movable_host_keys':sum(group_sizes[s]>1 for s in strata),'all_host_keys':len(pkeys),'observed_H6_equal_expected_links':observed,'null_H6_equal_mean':float(nullavg.mean()),'null_H6_equal95':np.quantile(nullavg,[.025,.975]).tolist(),'mc_fraction_null_at_least_observed_plus1':float((1+sum(nullavg>=observed))/500),'not_causal_or_realized_network':True}
    results={}
    for label,threshold in [('all',0),('coverage80',.8),('coverage90',.9)]:
        x=d[d.minimum_species_coverage>=threshold].copy()
        results[label]={}
        for col,log in [('host_support',False),('expected_links',True),('null_link_excess',False)]:
            results[label][col]=model(x,col,log)
    d.to_csv(OUT/'H7_candidate_host_network_panel.csv',index=False);pd.DataFrame(family_summary).to_csv(OUT/'family_source_coverage.csv',index=False);pd.DataFrame(spec).to_csv(OUT/'butterfly_catalogue_support.csv',index=False)
    np.savez_compressed(OUT/'candidate_network_matrices.npz',B=B,P=P,A=A,PB=PB,PP=PP,NB=NB,NP=NP,bkeys=np.array(bkeys),pkeys=np.array(pkeys),strata=strata,null=null)
    (OUT/'model_results.json').write_text(json.dumps(results,indent=2));(OUT/'null_audit.json').write_text(json.dumps(nullqa,indent=2))
    summary={'n_original_rows':len(d),'H7':d.h.nunique(),'H6':d.h6.nunique(),'butterfly_keys':len(bkeys),'host_keys':len(pkeys),'catalogue_edges':int(A.sum()),'edges_co_recorded_anywhere':int(np.sum(A*((B>0).T@(P>0)>0))),'total_candidate_links_across_rows':int(d.candidate_links.sum()),'median_host_support':float(d.host_support.median()),'median_catalogue_coverage':float(d.catalogue_coverage.median()),'family_summary':family_summary,'null':nullqa}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print(json.dumps({k:{c:{kk:v for kk,v in res.items() if kk in ['n_rows','H7','H6','H6_gain_pct','H6_fixed_error_bootstrap95_pct','country_equal_gain_pct','W_beta_perSD_descriptive','W_beta_H6_cluster95_descriptive','independent_sklearn_prediction_max_diff']} for c,res in rv.items()} for k,rv in results.items()},indent=2))
if __name__=='__main__':main()
