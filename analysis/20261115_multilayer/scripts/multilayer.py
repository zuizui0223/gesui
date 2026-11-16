#!/usr/bin/env python3
"""Same-H7-period aquatic, plant, butterfly conditional co-response checks.

These observational associations are NOT causal mediation. Exposure W is
modelled basin-outlet wastewater source geography (not measured delivered dose).
GBIF butterfly entries are occurrence records, not individual abundance.
"""
from __future__ import annotations
import argparse
import collections
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

PANEL_SOURCE_SHA256 = '1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca'
PERIODS = ('2017-2019', '2020-2022')
OUTCOMES = {'EQR': 'EQR_original_site_equal_basin_period_mean',
            'plant': 'plant_S50', 'HNS': 'primary_record_mean'}
COMMON_NUMERIC = ['crop','urban','riverflow','naturalout','area','nEEA','nAnimal','nPlant','nButterfly']
HNS_ONLY = ['primary_record_coverage_raw','primary_host_completeness',
            'primary_rf_Lycaenidae','primary_rf_Pieridae',
            'primary_rf_Papilionidae','primary_rf_Hesperiidae']
CATEGORICAL = ['country_method','period']


def add_country_method(df: pd.DataFrame) -> pd.DataFrame:
    if 'country' not in df or 'method' not in df:
        raise ValueError('Original country and EQR assessment method required')
    out=df.copy()
    out['country_method']=out['country'].fillna('__UNKNOWN__').astype(str)+'::'+out['method'].fillna('__UNKNOWN__').astype(str)
    return out


def check_grain(df: pd.DataFrame) -> None:
    for c in ('h', 'h6', 'period'):
        if c not in df: raise ValueError('Missing original required grain '+c)
    if df.duplicated(['h', 'period']).any():
        raise ValueError('Repeated original H7×period key')
    if df['h6'].isna().any() or df['period'].isna().any():
        raise ValueError('Source geography/period must not be missing')
    if df.groupby('h')['h6'].nunique().max() > 1:
        raise ValueError('H7 changes true H6 hydrological parent')


def make_fold_ids(groups: pd.Series | np.ndarray) -> np.ndarray:
    return np.asarray([int(hashlib.sha256(str(k).encode()).hexdigest()[:8], 16) % 5
                       for k in groups], dtype=int)


def _matrix(train: pd.DataFrame, test: pd.DataFrame,
            numeric: list[str], cats: list[str]):
    """No test-outcome-dependent encoding: numeric scales and vocab train-only."""
    a=[np.ones(len(train))]; b=[np.ones(len(test))]
    for c in numeric:
        tr=train[c].astype(float).to_numpy(); te=test[c].astype(float).to_numpy()
        if not (np.isfinite(tr).all() and np.isfinite(te).all()):
            raise ValueError('Unobserved numeric '+c)
        mu, sd=float(tr.mean()), float(tr.std())
        if sd<1e-10:sd=1.
        a.append((tr-mu)/sd);b.append((te-mu)/sd)
    for c in cats:
        tr=train[c].fillna('__MISSING__').astype(str)
        te=test[c].fillna('__MISSING__').astype(str)
        levels=sorted(tr.unique())
        for lv in levels[1:]:
            a.append((tr==lv).to_numpy(dtype=float))
            b.append((te==lv).to_numpy(dtype=float))
    return np.column_stack(a),np.column_stack(b)


def ridge_oof(df: pd.DataFrame, response: str, base_numeric: list[str],
              categorical: list[str], extra_numeric: list[str],
              alpha: float = 10.) -> np.ndarray:
    check_grain(df)
    if set(base_numeric) & set(extra_numeric):raise ValueError('Repeated source covariate')
    h6=df['h6'].astype(str).to_numpy()
    fold=make_fold_ids(h6)
    y=df[response].astype(float).to_numpy()
    if not np.isfinite(y).all():raise ValueError('Response not observed '+response)
    out=np.full(len(df),np.nan)
    for k in range(5):
        tr=fold!=k;te=~tr
        if not np.any(te): continue
        if tr.sum()<=15:raise ValueError('Insufficient heldout source training sample')
        dt=df.iloc[np.where(tr)[0]];dv=df.iloc[np.where(te)[0]]
        X,T=_matrix(dt,dv,base_numeric+extra_numeric,categorical)
        cnt=collections.Counter(h6[tr]);w=np.array([1/cnt[g] for g in h6[tr]])
        w/=w.mean()
        gram=X.T@(w[:,None]*X)
        penalty=np.eye(X.shape[1])*alpha;penalty[0,0]=0
        coef=np.linalg.solve(gram+penalty,X.T@(w*y[tr]))
        out[te]=T@coef
    if not np.isfinite(out).all():raise ValueError('Some H6 heldouts unpredicted')
    return out


def cv_add_value(df: pd.DataFrame, response: str, base_numeric: list[str],
                 added_numeric: list[str], categorical: list[str]=CATEGORICAL,
                 seed: int=20261115) -> dict:
    y=df[response].astype(float).to_numpy(); h6=df.h6.astype(str).to_numpy()
    p0=ridge_oof(df,response,base_numeric,categorical,[])
    p1=ridge_oof(df,response,base_numeric,categorical,added_numeric)
    e0=(y-p0)**2; e1=(y-p1)**2
    groups=sorted(set(h6))
    a=np.array([e0[h6==g].mean() for g in groups]);b=np.array([e1[h6==g].mean() for g in groups])
    rng=np.random.default_rng(seed)
    ix=rng.integers(0,len(groups),size=(2000,len(groups)))
    bs=100*(a[ix].mean(1)-b[ix].mean(1))/a[ix].mean(1)
    return {'H6_equal_baseline_MSE':float(a.mean()),'H6_equal_added_MSE':float(b.mean()),
            'H6_equal_MSE_gain_percent':float(100*(a.mean()-b.mean())/a.mean()),
            'H6_heldout_fixed_errors_bootstrap95_percent':list(np.percentile(bs,[2.5,97.5])),
            'n_H6':len(groups), 'n_H7':int(df.h.nunique()),'n_rows':len(df),
            'row_weighted_MSE_gain_percent':float(100*(e0.mean()-e1.mean())/e0.mean())}


def cluster_paths(df: pd.DataFrame, response: str, base_numeric: list[str],
                  extras: list[str], categorical: list[str]=CATEGORICAL) -> dict:
    """Within-sample partial associations, H6 cluster descriptive intervals.

    We deliberately do not return any direct/indirect causal effect.
    """
    import statsmodels.api as sm
    num=base_numeric+extras
    X,_=_matrix(df,df,num,categorical)
    yy=df[response].astype(float).to_numpy()
    h6=df.h6.astype(str).to_numpy()
    freq=collections.Counter(h6);w=np.array([1/freq[g] for g in h6])
    w/=w.mean()
    model=sm.WLS(yy,X,weights=w).fit(cov_type='cluster',cov_kwds={'groups':h6}, use_t=False)
    p=X.shape[1]; rank=int(np.linalg.matrix_rank(X))
    d={}
    for key in extras:
        j=num.index(key)+1
        coef=float(model.params[j]);se=float(model.bse[j]);ci=[coef-1.96*se,coef+1.96*se]
        d[key]={'beta_per_plus_one_source_SD':coef,'H6_cluster_descriptive95':ci}
    return {'associations':d,'n_design_columns':p,'design_rank':rank,
            'observational_not_causal_mediation':True,
            'outcome_scale':'original units of '+response}


def paired_deltas(df: pd.DataFrame) -> pd.DataFrame:
    check_grain(df)
    fields={'primary_record_mean':'d_HNS','plant_S50':'d_plant_S50',
            'EQR_original_site_equal_basin_period_mean':'d_EQR',
            'nPlant':'d_nPlant','nButterfly':'d_nButterfly',
            'nAnimal':'d_nAnimal','nEEA':'d_nEEA'}
    pivots={}
    for c in list(fields)+['method','h6','country','W']:
        if c in df:pivots[c]=df.pivot(index='h',columns='period',values=c)
    x=pivots['method']
    if any(t not in x for t in PERIODS):raise ValueError('Missing frozen analysis period')
    valid=x[PERIODS[0]].notna() & x[PERIODS[1]].notna()
    out=pd.DataFrame({'h':x.index[valid]})
    out['h6']=pivots['h6'].loc[valid,PERIODS[0]].to_numpy()
    out['country']=pivots['country'].loc[valid,PERIODS[0]].to_numpy()
    out['same_method']=(x.loc[valid,PERIODS[0]].to_numpy()==x.loc[valid,PERIODS[1]].to_numpy())
    if 'W' in pivots:
        out['W_constant']=np.isclose(pivots['W'].loc[valid,PERIODS[0]].to_numpy(),
                                    pivots['W'].loc[valid,PERIODS[1]].to_numpy(),atol=1e-12)
    for source,target in fields.items():
        if source in pivots:
            out[target]=pivots[source].loc[valid,PERIODS[1]].to_numpy(dtype=float)-pivots[source].loc[valid,PERIODS[0]].to_numpy(dtype=float)
    return out.reset_index(drop=True)


def corr_bootstrap_group(df: pd.DataFrame, a: str, b: str, seed: int=20261115) -> dict:
    from scipy.stats import spearmanr
    d=df[[a,b,'h6']].replace([np.inf,-np.inf],np.nan).dropna().copy()
    groups=sorted(d.h6.astype(str).unique());group_rows={g:d[d.h6.astype(str)==g] for g in groups}
    rho=float(spearmanr(d[a],d[b]).statistic)
    rng=np.random.default_rng(seed);bs=[]
    for _ in range(1200):
        sel=rng.choice(groups,len(groups),replace=True)
        q=pd.concat([group_rows[x] for x in sel],ignore_index=True)
        rr=float(spearmanr(q[a],q[b]).statistic)
        if np.isfinite(rr):bs.append(rr)
    return {'rho':rho,'H6_bootstrap_descriptive95':list(np.percentile(bs,[2.5,97.5])),
            'n_pairs':len(d),'n_H6':len(groups),'not_temporal_wastewater_change':True}


def run(panel_path: Path, outdir: Path):
    raw=panel_path.read_bytes();sha=hashlib.sha256(raw).hexdigest()
    if sha!=PANEL_SOURCE_SHA256:raise ValueError('Frozen signed H7×period source SHA mismatch: '+sha)
    df=add_country_method(pd.read_csv(panel_path,dtype={'h':str,'h6':str}));check_grain(df)
    if not(len(df)==454 and df.h.nunique()==232 and df.h6.nunique()==90):
        raise ValueError('Original sample unit cardinality drift')
    if set(df.period)!=set(PERIODS):raise ValueError('Original time window changed')
    for col in list(OUTCOMES.values())+['W']+COMMON_NUMERIC+HNS_ONLY:
        if col not in df or not np.isfinite(df[col].to_numpy(dtype=float)).all():
            raise ValueError('Missing original source variable '+col)
    result={'schema':'GESUI_20261115_SAME_H7_MULTIRESPONSE_NONCAUSAL_V1',
      'original_panel_sha256':sha,'n_rows':len(df),'n_H7':df.h.nunique(),'n_H6':df.h6.nunique(),
      'n_countries':df.country.nunique(),
      'models_are_exploratory_posthoc':True,
      'W_is_modelled_static_basin_outlet_geography_not_measured_effluent_dose':True,
      'aquatic_EQR_is_not_insect_emergence_biomass':True,
      'HNS_is_GBIF_record_weighted_host_habitat_N_affinity_not_butterfly_count':True,
      'causal_SEM_or_mediation_fitted':False,'models':{}}
    recipes={
     'W_to_EQR':(OUTCOMES['EQR'],COMMON_NUMERIC,['W']),
     'W_to_PlantS50':(OUTCOMES['plant'],COMMON_NUMERIC,['W']),
     'W_to_HNS':(OUTCOMES['HNS'],COMMON_NUMERIC+HNS_ONLY,['W']),
     'Plant_after_W_to_HNS':(OUTCOMES['HNS'],COMMON_NUMERIC+HNS_ONLY+['W'],['plant_S50']),
     'EQR_after_W_to_HNS':(OUTCOMES['HNS'],COMMON_NUMERIC+HNS_ONLY+['W'],[OUTCOMES['EQR']]),
     'Plant_and_EQR_after_W_to_HNS':(OUTCOMES['HNS'],COMMON_NUMERIC+HNS_ONLY+['W'],['plant_S50',OUTCOMES['EQR']]),
     'W_after_Plant_and_EQR_to_HNS':(OUTCOMES['HNS'],COMMON_NUMERIC+HNS_ONLY+['plant_S50',OUTCOMES['EQR']],['W']),
     'EQRaquatic_after_W_to_PlantS50':(OUTCOMES['plant'],COMMON_NUMERIC+['W'],[OUTCOMES['EQR']]),
    }
    for name,(target,baseline,add) in recipes.items():
        print('FITTING',name,flush=True)
        result['models'][name]={
            'outcome':target,'baseline_numeric':baseline,'added_numeric':add,
            'predictive_added_value':cv_add_value(df,target,baseline,add),
            'source_observational_partial_coefficients':cluster_paths(df,target,baseline,add)}
    paired=paired_deltas(df)
    if not(len(paired)==222 and paired.h6.nunique()==89 and (paired.same_method.sum()==172)):
        raise ValueError('Repeated original H7 source or method compatibility drift')
    if not paired.W_constant.all():raise ValueError('Unexpected W change; historical source W should be static')
    result['paired_H7']={'total':222,'same_majority_method':172,'different_majority_method':50,
        'mean_delta_HNS':float(paired.d_HNS.mean()),
        'n_HNS_rise':int((paired.d_HNS>0).sum()),'n_HNS_fall':int((paired.d_HNS<0).sum()),'W_static_in_all_222':True,
        'correlations':{}}
    for label,d in [('all222',paired),('same_method172',paired[paired.same_method])]:
        for factor in ['d_plant_S50','d_EQR']:
            result['paired_H7']['correlations'][label+'_'+factor+'_vs_dHNS']=corr_bootstrap_group(d,'d_HNS',factor)
    outdir.mkdir(parents=True,exist_ok=True)
    (outdir/'H7_MULTIRESPONSE_NONCAUSAL_RESULTS.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False))
    paired.to_csv(outdir/'PAIRED_222_H7_OBSERVED_CHANGE.csv',index=False)
    rows=[]
    for k,v in result['models'].items():
        c=v['predictive_added_value'];co=v['source_observational_partial_coefficients']
        rows.append({'model':k,'outcome':v['outcome'],'added':';'.join(v['added_numeric']),
                     'n_rows':c['n_rows'],'n_H7':c['n_H7'],'n_H6':c['n_H6'],
                     'H6_MSE_gain_percent':c['H6_equal_MSE_gain_percent'],
                     'fixed_error95_low':c['H6_heldout_fixed_errors_bootstrap95_percent'][0],
                     'fixed_error95_high':c['H6_heldout_fixed_errors_bootstrap95_percent'][1],
                     'row_gain_percent':c['row_weighted_MSE_gain_percent'],
                     'design_rank':co['design_rank'],'design_cols':co['n_design_columns']})
    pd.DataFrame(rows).to_csv(outdir/'PREDICTIVE_MODEL_COMPARISON.csv',index=False)
    print(json.dumps({'input_sha256':sha,'model_summary':rows,'pairs':result['paired_H7']},sort_keys=True,allow_nan=False),flush=True)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--panel',required=True);parser.add_argument('--output',required=True)
    a=parser.parse_args();run(Path(a.panel),Path(a.output))
