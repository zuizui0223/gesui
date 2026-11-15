#!/usr/bin/env python3
"""Actual reported 2018/2022 EEA upstream N/P index vs same-H7 HNS/EQR/plant.

Source reporting index only: no-match is NOT a clean control; annual mass / model
flow is NOT measured delivered nutrient concentration. Posthoc predictive audit.
"""
from __future__ import annotations
import argparse
import collections
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from source_join import load_locked_source,join_source,examine_source_eligibility,sha256,FROZEN_PANEL

PERIODS=('2017-2019','2020-2022')


def load_core(path:Path):
    spec=importlib.util.spec_from_file_location('frozen_multilayer_core',path)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def make_source_features(df:pd.DataFrame)->pd.DataFrame:
    d=df.copy()
    if not {'reported_positive','reporting_class','log_reported_N_mgL_equivalent','log_reported_P_ugL_equivalent'}.issubset(d):
        raise ValueError('Original matched source state not available')
    if d.reported_positive.isna().any():raise ValueError('Missing source reporting state')
    if not (d.reported_positive.isin([0,1])).all():raise ValueError('Reporting state not binary')
    return d


def refit_old_static_W(df,core)->dict:
    x=core.add_country_method(df)
    return {'HNS':core.cv_add_value(x,'primary_record_mean',core.COMMON_NUMERIC+core.HNS_ONLY,['W'])['H6_equal_MSE_gain_percent'],
            'EQR':core.cv_add_value(x,'EQR_original_site_equal_basin_period_mean',core.COMMON_NUMERIC,['W'])['H6_equal_MSE_gain_percent']}


def grouped_mse_comparison(e0,e1,groups,n_resamples=2000,seed=20261115):
    """H6-independent, equal-true-parent-weight heldout squared error comparison."""
    e0=np.asarray(e0,dtype=float);e1=np.asarray(e1,dtype=float);groups=np.asarray(groups,dtype=str)
    if e0.shape!=e1.shape or len(groups)!=len(e0):raise ValueError('H6 outcome-error shape drift')
    if not np.isfinite(e0).all() or not np.isfinite(e1).all() or (e0<0).any() or (e1<0).any():
        raise ValueError('Invalid true OOF squared errors')
    names=np.unique(groups)
    x=np.array([e0[groups==g].mean() for g in names]);y=np.array([e1[groups==g].mean() for g in names])
    if (x.mean()<=0):raise ValueError('Empty or zero baseline error')
    gain=100*(x.mean()-y.mean())/x.mean()
    rng=np.random.default_rng(seed)
    ix=rng.integers(len(names),size=(n_resamples,len(names)))
    ax=x[ix].mean(1);by=y[ix].mean(1)
    keep=ax>0
    boot=100*(ax[keep]-by[keep])/ax[keep]
    return {'n_H6':int(len(names)), 'n_rows':len(e0), 'baseline_H6equal_MSE':float(x.mean()),
            'added_H6equal_MSE':float(y.mean()),'H6_equal_MSE_gain_percent':float(gain),
            'H6_fixed_heldout_errors_bootstrap95':list(map(float,np.percentile(boot,[2.5,97.5]))),
            'row_weighted_gain_percent':float(100*(e0.mean()-e1.mean())/e0.mean()),
            'not_causal_effect':True}


def paired_source_change(df):
    a=df.pivot(index='h',columns='period',values='reported_positive').dropna()
    if len(a)!=222 or not set(PERIODS).issubset(a.columns):raise ValueError('Expected original 222 repeated H7')
    c=collections.Counter(f'{int(q)}_to_{int(r)}' for q,r in a[list(PERIODS)].itertuples(index=False,name=None))
    return {'n_pairs':int(len(a)),
            'report_status_transitions':{k:int(c.get(k,0)) for k in ('0_to_0','0_to_1','1_to_0','1_to_1')},
            'not_equivalent_to_WWTP_opening_or_shutdown':True}


def country_h6_disjoint_folds(df):
    countries=sorted(df.country.astype(str).unique())
    out=[]
    for c in countries:
        test=df.country.astype(str)==c
        h6_test=set(df.loc[test,'h6'].astype(str))
        train=(~test)&(~df.h6.astype(str).isin(h6_test))
        overlap=set(df.loc[train,'h6'].astype(str))&h6_test
        out.append({'target_country':c,'test_n':int(test.sum()), 'train_n':int(train.sum()),
                    'train_h6_n':int(df.loc[train,'h6'].nunique()),
                    'test_h6_n':int(df.loc[test,'h6'].nunique()),'overlap_H6':len(overlap),
                    'train_mask':train.to_numpy(dtype=bool),'test_mask':test.to_numpy(dtype=bool)})
    return out


def country_predictions(x,core,response,base,extras):
    """Strict leave-entire-country-out; any shared H6 removed from training."""
    x=core.add_country_method(x)
    y=x[response].to_numpy(dtype=float)
    h6=x.h6.astype(str).to_numpy()
    predictions=np.full(len(x),np.nan)
    info=[]
    for fold in country_h6_disjoint_folds(x):
        tr=fold['train_mask'];te=fold['test_mask']
        info.append({k:v for k,v in fold.items() if k not in ('train_mask','test_mask')})
        if tr.sum()<100 or len(set(h6[tr]))<25 or te.sum()<1:
            continue
        train=x.loc[tr];test=x.loc[te]
        X,T=core._matrix(train,test,base+extras,core.CATEGORICAL)
        freq=collections.Counter(h6[tr]);w=np.array([1.0/freq[g] for g in h6[tr]])
        w/=w.mean();pen=np.eye(X.shape[1])*10.;pen[0,0]=0
        beta=np.linalg.solve(X.T@(w[:,None]*X)+pen,X.T@(w*y[tr]))
        predictions[te]=T@beta
    return predictions,info


def country_compare(x,core,response,base,extra_a,extra_b):
    y=x[response].to_numpy(dtype=float)
    a,info=country_predictions(x,core,response,base,extra_a)
    b,_=country_predictions(x,core,response,base,extra_b)
    valid=np.isfinite(a)&np.isfinite(b)
    if valid.sum()<len(y) or not np.all(valid):
        return {'status':'unsupported_country_transfer','n_valid':int(valid.sum()), 'n_total':len(y),
                'countries':info}
    errors_a=(y-a)**2; errors_b=(y-b)**2
    summary={};supported=[]
    for f in info:
        c=f['target_country'];mask=(x.country.astype(str)==c).to_numpy()
        mse_a=errors_a[mask].mean();mse_b=errors_b[mask].mean()
        score=100*(mse_a-mse_b)/mse_a if mse_a>0 else None
        summary[c]={'test_rows':f['test_n'],'test_H6':f['test_h6_n'],'train_rows':f['train_n'],
                    'gain_percent':float(score) if score is not None else None}
        if f['test_n']>=10 and f['test_h6_n']>=5 and score is not None:supported.append(c)
    if supported:
        mse_a=np.mean([errors_a[(x.country.astype(str)==c).to_numpy()].mean() for c in supported])
        mse_b=np.mean([errors_b[(x.country.astype(str)==c).to_numpy()].mean() for c in supported])
        equal=100*(mse_a-mse_b)/mse_a
    else:equal=None
    return {'status':'complete','n_rows':len(x),'row_weighted_gain_percent':float(100*(errors_a.mean()-errors_b.mean())/errors_a.mean()),
            'supported_countries_equal_gain_percent':float(equal) if equal is not None else None,
            'supported_countries':supported,'country_details':summary, 'all_training_target_H6_disjoint':True}


def positive_spearman(x,core,seed=20261115):
    pos=x.loc[x.reported_positive==1].copy()
    groups=pos.h6.astype(str).unique()
    group_indices={g:np.flatnonzero(pos.h6.astype(str).to_numpy()==g) for g in groups}
    fields={'HNS':'primary_record_mean','EQR':'EQR_original_site_equal_basin_period_mean','PlantS50':'plant_S50'}
    sources={'N_source_equiv':'n_mg_l','P_source_equiv':'p_mg_l'}
    rng=np.random.default_rng(seed)
    res={}
    for response,outcome in fields.items():
        for exposure,col in sources.items():
            xx_all=pos[outcome].to_numpy(dtype=float);yy_all=pos[col].to_numpy(dtype=float)
            rho=float(spearmanr(xx_all,yy_all).statistic)
            boot=[]
            for i in range(1000):
                samp=rng.choice(groups,size=len(groups),replace=True)
                ix=np.concatenate([group_indices[g] for g in samp])
                xx=xx_all[ix];yy=yy_all[ix]
                if len(set(xx))<2 or len(set(yy))<2:continue
                b=float(spearmanr(xx,yy).statistic)
                if np.isfinite(b):boot.append(b)
            res[response+'_vs_'+exposure]={'n_positive':int(len(pos)),'n_H6':len(groups),
                    'spearman_r':rho,'H6_cluster_bootstrap95':list(map(float,np.percentile(boot,[2.5,97.5]))),
                    'selected_report_positive_and_not_adjusted':True}
    return res


def run(panel:Path,source_zip:Path,old_core:Path,outdir:Path):
    rawpanelsha=sha256(panel)
    if rawpanelsha!=FROZEN_PANEL:raise ValueError('Original source HNS panel SHA drift '+rawpanelsha)
    src,src_qa=load_locked_source(source_zip)
    p=pd.read_csv(panel,dtype={'h':str,'h6':str})
    merged,join_qa=join_source(p,src)
    source_qa=examine_source_eligibility(merged)
    if not source_qa['exploratory_support_gate_pass']:
        raise ValueError('Pre-source outcome population insufficient to fit reported-index prediction')
    core=load_core(old_core)
    x=core.add_country_method(make_source_features(merged))
    ret={'schema':'GESUI_20261115_REAL_REPORTED_NP_TO_HNS_EQR_PLANT_H7_SOURCE_GATED_POSTHOC_v1',
         'original_HNS_panel_SHA256':rawpanelsha,'source_reconciled_index':src_qa,
         'join_qa':join_qa,'source_qa':source_qa,'paired_report_status':paired_source_change(x),
         'original_prior_static_W_baseline_replication':refit_old_static_W(x,core),
         'N_P_source_index_not_actual_received_nutrient_dose':True,
         'source_no_matched_report_not_a_clean_control':True,
         'outcome_type':'EEA method-dependent EQR, GBIF plant rarefied S50, butterfly occurrence-weighted HNS',
         'study':'posthoc observational predictive check, no causal wastewater effect or SEM',
         'cv_group':'authentic H6 parent; original H7 periods never cross folds',
         'models':{}, 'source_positive_only_descriptive_rank':positive_spearman(x,core)}
    exposures=['reported_positive','log_reported_N_mgL_equivalent','log_reported_P_ugL_equivalent']
    outcomes={
        'HNS':('primary_record_mean',core.COMMON_NUMERIC+core.HNS_ONLY),
        'EQR':('EQR_original_site_equal_basin_period_mean',core.COMMON_NUMERIC),
        'PlantS50':('plant_S50',core.COMMON_NUMERIC),
    }
    summary=[]
    for i,(label,(target,common)) in enumerate(outcomes.items()):
        base=common+['W']
        y=x[target].to_numpy(dtype=float)
        h6=x.h6.astype(str).to_numpy()
        recipes={
            'B0_static_W_only':[],
            'B1_plus_matched_source_presence':['reported_positive'],
            'B2_plus_N_report_index':['log_reported_N_mgL_equivalent'],
            'B3_plus_P_report_index':['log_reported_P_ugL_equivalent'],
            'B4_plus_presence_and_N':['reported_positive','log_reported_N_mgL_equivalent'],
            'B5_plus_presence_and_P':['reported_positive','log_reported_P_ugL_equivalent'],
            'B6_plus_presence_and_N_and_P':exposures,
        }
        predictions={}
        for name,extra in recipes.items():
            predictions[name]=core.ridge_oof(x,target,base,core.CATEGORICAL,extra)
        errors={name:(y-pred)**2 for name,pred in predictions.items()}
        comparisons={}
        for num,(name,extra) in enumerate(recipes.items()):
            if name=='B0_static_W_only':continue
            contrast='B1_plus_matched_source_presence' if name in ('B4_plus_presence_and_N','B5_plus_presence_and_P','B6_plus_presence_and_N_and_P') else 'B0_static_W_only'
            d=grouped_mse_comparison(errors[contrast],errors[name],h6,seed=20261115+i*20+num)
            d['target']=target;d['baseline']=contrast;d['added_predictors']=extra
            comparisons[name]=d
            summary.append({'outcome':label,'model':name,'baseline':contrast,'gain_pct':d['H6_equal_MSE_gain_percent'],
                'ci_low':d['H6_fixed_heldout_errors_bootstrap95'][0],'ci_high':d['H6_fixed_heldout_errors_bootstrap95'][1],
                'n':len(x),'n_H6':int(x.h6.nunique()),'source_positive_n':source_qa['positive_rows']})
        # All adjusted coefficients are observational and estimated on the same source cohort.
        partials={}
        for name in ('B1_plus_matched_source_presence','B4_plus_presence_and_N','B5_plus_presence_and_P','B6_plus_presence_and_N_and_P'):
            partials[name]=core.cluster_paths(x,target,base,recipes[name],core.CATEGORICAL)
        ret['models'][label]={'target':target,'baseline_features':base,'baseline_MSE_H6equal':float(grouped_mse_comparison(errors['B0_static_W_only'],errors['B0_static_W_only'],h6)['baseline_H6equal_MSE']),
             'comparisons':comparisons,'observational_partials':partials}
        if label=='HNS':
            ret['models'][label]['country_transfer']={
             'report_presence_vs_staticW':country_compare(x,core,target,base,[],['reported_positive']),
             'N_gradient_after_presence':country_compare(x,core,target,base,['reported_positive'],['reported_positive','log_reported_N_mgL_equivalent'])}
    outdir.mkdir(parents=True,exist_ok=True)
    (outdir/'SOURCE_QA_AND_EXPOSURE.json').write_text(json.dumps({'source':src_qa,'join':join_qa,'eligibility':source_qa,'paired':ret['paired_report_status']},indent=2,sort_keys=True))
    # Persist source-index panel with names that explicitly distinguish reported load from delivered fraction.
    result_panel_cols=['h','h6','period','country','method','W','crop','urban','riverflow',
        'reported_positive','reporting_class','n_source','n_mg_l','p_mg_l','log_reported_N_mgL_equivalent','log_reported_P_ugL_equivalent',
        'primary_record_mean','EQR_original_site_equal_basin_period_mean','plant_S50']
    x[result_panel_cols].to_csv(outdir/'H7_454_REPORTED_SOURCE_INDEX_HNS_EQR_PLANT.csv',index=False)
    (outdir/'RESULTS_EEA_REPORTED_NP_TO_HNS_EQR_PLANT.json').write_text(json.dumps(ret,indent=2,sort_keys=True,allow_nan=False))
    pd.DataFrame(summary).to_csv(outdir/'MODEL_COMPARISONS.csv',index=False)
    print(json.dumps({'source_join':join_qa,'positive_support':source_qa,'baseline_W_replicated':ret['original_prior_static_W_baseline_replication'],
            'summary':summary,'source_positive_rank':ret['source_positive_only_descriptive_rank']},sort_keys=True,allow_nan=False),flush=True)
    return ret

if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--panel',type=Path,required=True)
    cli.add_argument('--source-zip',type=Path,required=True)
    cli.add_argument('--reference-core',type=Path,required=True)
    cli.add_argument('--out',type=Path,required=True)
    a=cli.parse_args();run(a.panel,a.source_zip,a.reference_core,a.out)
