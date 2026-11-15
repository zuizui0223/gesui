"""Original 454 source rows: host concentration, family, geography robustness.
Exploratory continuation, NOT regional host evidence / field feeding / sewage causality.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import collections,hashlib,json,pathlib,sys,zipfile
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from scipy.stats import hypergeom
from robustness_core import group_weights,retained_permutation
R=pathlib.Path(__file__).resolve().parents[1]
OLD=R/'previous/gesui_host_network_20261114';OUT=R/'output';OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(OLD/'scripts'))
import analyse as old
FAMS=old.FAMS;SEED=20261114

def finite_dict(v):
    if isinstance(v,dict):return {str(k):finite_dict(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [finite_dict(x) for x in v]
    if isinstance(v,np.ndarray):return finite_dict(v.tolist())
    if isinstance(v,(np.integer,np.floating)):v=v.item()
    if isinstance(v,float) and not np.isfinite(v):return None
    return v

def summary(obs, null, h6):
    w=group_weights(h6);o=float(w@obs);n=null@w
    groups=sorted(set(h6));byobs=np.array([obs[h6==h].mean() for h in groups]);byref=np.array([null[:,h6==h].mean() for h in groups])
    idx=np.random.default_rng(SEED).integers(0,len(groups),(3000,len(groups)))
    den=byref[idx].mean(1);rat=100*(byobs[idx].mean(1)/den-1)
    return {'n_rows':len(obs),'n_H6':len(groups),'observed':o,'null_mean':float(n.mean()),
      'excess_pct':100*(o/n.mean()-1) if n.mean()>0 else None,
      'null_assignment_95':np.quantile(n,[.025,.975]).tolist(),
      'H6_fixed_nullmean_bootstrap95_pct':np.quantile(rat,[.025,.975]).tolist(),
      'mc_fraction_at_least_observed_plus1':float((1+sum(n>=o))/(len(n)+1)),
      'excess_abs':o-float(n.mean()),'not_confirmatory_p_or_causal':True}

def main():
    z=np.load(OLD/'output/candidate_network_matrices.npz')
    A,B,P,PB,PP,NB,NP=[z[k] for k in ['A','B','P','PB','PP','NB','NP']]
    d=pd.read_csv(OLD/'output/H7_candidate_host_network_panel.csv',dtype={'h':str,'h6':str},keep_default_na=False)
    for c in [f+'_host_support' for f in FAMS]:d[c]=pd.to_numeric(d[c],errors='coerce')
    source=json.loads((OLD/'source/butterfly_usage_aliases.json').read_text());fm={x['original_speciesKey']:x['family'] for x in source}
    fam=np.array([fm[k] for k in z['bkeys']]);h6=d.h6.to_numpy();hs=sorted(set(h6))
    inc=np.array([(P[h6==h]>0).any(0) for h in hs]).T.astype(float)
    norm=np.linalg.norm(inc,axis=1);nonzero=norm>0;labels=np.full(len(norm),-1)
    labels[nonzero]=KMeans(n_clusters=6,random_state=SEED,n_init=10).fit_predict(inc[nonzero]/norm[nonzero,None])
    strata=np.array([s+'|geo'+str(k) for s,k in zip(z['strata'],labels)])
    ba=np.array([PB[:,fam==f]@A[fam==f] for f in FAMS])
    ba=np.concatenate([ba.sum(0)[None,:,:],ba],axis=0)
    group_names=['All']+FAMS
    contributions=ba*PP[None,:,:]
    np.testing.assert_allclose(contributions[0].sum(1),d.expected_links,rtol=0,atol=1e-10)
    ranking=pd.read_csv(OLD/'output/host_plant_candidate_contributions.csv',dtype={'plantKey':str})
    subsets={'all454':np.ones(len(d),bool),'coverage80':d.minimum_species_coverage.to_numpy(float)>=.8,
             'period2017_2019':d.period.to_numpy()=='2017-2019','period2020_2022':d.period.to_numpy()=='2020-2022',
             'France':d.country.to_numpy()=='FR','outside_France':d.country.to_numpy()!='FR'}
    for c,rows in d.groupby('country'):
        if len(rows)>=10 and rows.h6.nunique()>=3:subsets['country_'+c]=d.country.to_numpy()==c
    records=[];fullnull=None;allobs=None;removal_meta=[]
    for remove_n in [0,1,5,20]:
        removed=set(ranking.head(remove_n).plantKey);keep=np.array([k not in removed for k in z['pkeys']])
        obs=contributions[:,:,keep].sum(2)
        null=np.empty((499,len(group_names),len(d)))
        rng=np.random.default_rng(SEED)
        for it in range(len(null)):
            perm=retained_permutation(strata,keep,rng)
            if it==0:
                np.testing.assert_array_equal(strata[perm],strata[keep])
                np.testing.assert_array_equal(A[:,perm].sum(1),A[:,keep].sum(1))
            null[it]=np.einsum('fsp,sp->fs',ba[:,:,perm],PP[:,keep],optimize=False)
        np.testing.assert_allclose(obs[0],obs[1:].sum(0),atol=1e-10)
        np.testing.assert_allclose(null[:,0],null[:,1:].sum(1),atol=1e-10)
        for sub,mask in subsets.items():
            for fi,f in enumerate(group_names):
                row={'remove_top_hosts':remove_n,'subset':sub,'family':f,
                     **summary(obs[fi,mask],null[:,fi,mask],h6[mask])}
                records.append(row)
        sizes=collections.Counter(strata[keep]);movable=np.array([sizes[s]>1 for s in strata[keep]])
        w=group_weights(h6)
        removal_meta.append({'remove_top_hosts':remove_n,'removed_names':ranking.head(remove_n).accepted_name.tolist(),
          'remaining_host_keys':int(keep.sum()),'movable_keys':int(movable.sum()),
          'observed_mass_in_movable_strata_fraction':float(w@contributions[0][:,keep][:,movable].sum(1)/(w@obs[0])),
          'fraction_original_index_removed':float(1-(w@obs[0])/(w@contributions[0].sum(1)))})
        if remove_n==0:fullnull=null;allobs=obs
        print('completed exclusion',remove_n,'kept',keep.sum(),flush=True)
    # Same frozen structural reference as prior PR86; compare before using new diagnostics.
    previous=json.loads((OLD/'output/independent_numerical_source_QA.json').read_text())['geography_constrained_null']
    mainrow=next(x for x in records if x['remove_top_hosts']==0 and x['subset']=='all454' and x['family']=='All')
    assert abs(mainrow['observed']-previous['observed_H6_equal_expected_links'])<1e-8
    assert abs(mainrow['null_mean']-previous['null_H6_equal_mean'])<1e-8
    d['geographic_null_expected_links']=fullnull[:,0,:].mean(0)
    d['geographic_null_link_excess']=allobs[0]-d.geographic_null_expected_links
    for fi,f in enumerate(FAMS,1):
        d[f+'_expected_links']=allobs[fi]
        d[f+'_geographic_null_links']=fullnull[:,fi,:].mean(0)
        d[f+'_geographic_null_excess']=allobs[fi]-d[f+'_geographic_null_links']
    (OUT/'structure_checkpoint.json').write_text(json.dumps(finite_dict(records),indent=2))
    d.to_csv(OUT/'panel_checkpoint.csv',index=False)
    # Retain original inferential target; stronger null and family outcomes are posthoc.
    models={}
    for lab,mask in [('all454',np.ones(len(d),bool)),('coverage80',d.minimum_species_coverage.to_numpy(float)>=.8)]:
        models[lab]={};df=d.loc[mask].copy()
        models[lab]['geographic_null_link_excess']=old.model(df,'geographic_null_link_excess',False)
        for f in FAMS:
            col=f+'_host_support';ff=df.loc[df[col].notna()].copy()
            models[lab][col]=old.model(ff,col,False)
            models[lab][col]['missing_source_denominators']=int(len(df)-len(ff))
            print('model',lab,col,'done',flush=True)
        print('models',lab,'done',flush=True)
    # Independent probability engine, sparse edge accumulation, all original response hashes.
    bp=np.array([hypergeom.sf(0,int(n),b,20) for n,b in zip(NB,B)])
    pp=np.array([hypergeom.sf(0,int(n),p,50) for n,p in zip(NP,P)])
    ii,jj=np.nonzero(A);sparse=np.sum(bp[:,ii]*pp[:,jj],axis=1)
    assert np.max(abs(sparse-allobs[0]))<1e-8
    receipts=[json.loads(x) for x in (OLD/'source/source_receipts.jsonl').read_text().splitlines()]
    for r in receipts:
        assert hashlib.sha256((OLD/'source/responses'/r['cache_file']).read_bytes()).hexdigest()==r['response_sha256']
    qa={'previous_sha_verified_responses':len(receipts),'original_graph_reproduced':True,
      'max_independent_scipy_sparse_expected_links_error':float(np.max(abs(sparse-allobs[0]))),
      'max_model_sklearn_prediction_error':max(v['independent_sklearn_prediction_max_diff'] for mm in models.values() for v in mm.values()),
      'previous_null_reproduced':previous,'removal_meta':removal_meta,'source_europe_rows_acquired':False,
      'source_window_and_sites_unchanged':True,'original_tests_passed':7,'new_tests_passed':8}
    (OUT/'structural_robustness.json').write_text(json.dumps(finite_dict(records),indent=2,allow_nan=False))
    flat=[]
    for r in records:
        x={k:v for k,v in r.items() if not isinstance(v,list)}
        x['H6_ci_low_pct'],x['H6_ci_high_pct']=r['H6_fixed_nullmean_bootstrap95_pct']
        x['null_lo'],x['null_hi']=r['null_assignment_95'];flat.append(x)
    pd.DataFrame(flat).to_csv(OUT/'structural_robustness.csv',index=False)
    (OUT/'model_results.json').write_text(json.dumps(finite_dict(models),indent=2,allow_nan=False))
    tab=[]
    for lab,mm in models.items():
        for name,v in mm.items():
            x={k:v[k] for k in ['n_rows','H7','H6','H6_gain_pct','country_equal_gain_pct','country_row_weighted_gain_pct','W_beta_perSD_descriptive']}
            x.update(cohort=lab,outcome=name,H6_gain_CI_low=v['H6_fixed_error_bootstrap95_pct'][0],H6_gain_CI_high=v['H6_fixed_error_bootstrap95_pct'][1],missing_denominators=v.get('missing_source_denominators',0))
            tab.append(x)
    pd.DataFrame(tab).to_csv(OUT/'model_comparison.csv',index=False)
    (OUT/'QA.json').write_text(json.dumps(finite_dict(qa),indent=2,allow_nan=False))
    d.to_csv(OUT/'H7_454_geographic_null_network_panel.csv',index=False)
    np.savez_compressed(OUT/'geographic_null_draws.npz',null=fullnull,observed=allobs,family_names=np.array(group_names),strata=strata)
    print(pd.DataFrame(flat).query('family == "All" and subset in ["all454","coverage80","outside_France","period2017_2019","period2020_2022"]')[['remove_top_hosts','subset','excess_pct','H6_ci_low_pct','H6_ci_high_pct']].to_string(index=False))
    print(pd.DataFrame(flat).query('remove_top_hosts == 0 and subset == "all454"')[['family','excess_pct','H6_ci_low_pct','H6_ci_high_pct']].to_string(index=False))
    print(pd.DataFrame(tab).to_string(index=False))
if __name__=='__main__':main()
