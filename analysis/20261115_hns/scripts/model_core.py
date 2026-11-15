"""Fixed exploratory model; prediction gain is not biological effect direction."""
from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import collections
import numpy as np,pandas as pd
from score_core import blocked_splits,ridge_fit_predict
BASE=['area','nEEA','nAnimal','nPlant','nButterfly','riverflow','crop','urban','naturalout','butterfly_S20','plant_S50']
FAMS=sorted(['Lycaenidae','Nymphalidae','Pieridae','Papilionidae','Hesperiidae'])

def numeric_columns(variant='primary',kind='record',family_adjust=True,extras=()):
    cols=BASE+[variant+'_record_coverage_raw',variant+'_species_coverage',variant+'_host_completeness']
    if family_adjust:cols += [f'{variant}_{"rf" if kind=="record" else "if"}_{f}' for f in FAMS[:-1]]
    return cols+list(extras)

def design(d,train,cols,country=True):
    v=d[cols].to_numpy(float)
    if not np.isfinite(v).all():raise ValueError('Non-finite covariates; do not silently drop or impute')
    mu=v[train].mean(0);sd=v[train].std(0);sd=np.where(sd>1e-12,sd,1.)
    arr=[np.ones((len(d),1)),(v-mu)/sd];names=['intercept']+list(cols)
    for col in ['period','method']+(['country'] if country else []):
        cats=sorted(set(d.loc[train,col].astype(str)))
        arr.append(np.column_stack([(d[col].astype(str)==x).to_numpy(float) for x in cats]))
        names.extend([col+'='+x for x in cats])
    return np.column_stack(arr),names

def groupweights(groups):
    ct=collections.Counter(groups);w=np.array([1/ct[g] for g in groups]);return w/w.mean()

def predict(d,y,cols,mode='H6',independent=False):
    pred=np.full(len(d),np.nan);splits=[]
    for g,tr,te in blocked_splits(d.h6,d.country_holdout,mode):
        if tr.sum()<30 or d.loc[tr,'h6'].nunique()<5:
            raise RuntimeError(f'Insufficient independent training for {mode} {g}')
        X,_=design(d,tr,cols,country=mode=='H6');w=groupweights(d.loc[tr,'h6'])
        if independent:
            from sklearn.linear_model import Ridge
            xm=np.average(X[tr,1:],axis=0,weights=w);ym=np.average(y[tr],weights=w)
            fit=Ridge(alpha=10.,fit_intercept=False,solver='cholesky').fit(X[tr,1:]-xm,y[tr]-ym,sample_weight=w)
            pred[te]=fit.predict(X[te,1:]-xm)+ym
        else:pred[te]=ridge_fit_predict(X[tr],y[tr],X[te],w)
        splits.append({'group':str(g),'train_rows':int(tr.sum()),'test_rows':int(te.sum()),'train_H6':int(d.loc[tr,'h6'].nunique()),'test_H6':int(d.loc[te,'h6'].nunique()),'H6_overlap':len(set(d.loc[tr,'h6'])&set(d.loc[te,'h6']))})
    if not np.isfinite(pred).all():raise AssertionError('Incomplete predictions')
    return pred,splits

def beta_cluster(d,y,cols):
    cols=list(cols)+['W'];X,names=design(d,np.ones(len(d),bool),cols,True)
    w=groupweights(d.h6);a=X*np.sqrt(w)[:,None];b=y*np.sqrt(w)
    coef=np.linalg.lstsq(a,b,rcond=None)[0];r=y-X@coef
    bread=np.linalg.pinv(a.T@a);meat=np.zeros_like(bread)
    for g in sorted(set(d.h6)):
        sel=(d.h6==g).to_numpy();z=X[sel].T@(w[sel]*r[sel]);meat+=np.outer(z,z)
    G=d.h6.nunique();rank=np.linalg.matrix_rank(a)
    V=bread@meat@bread*(G/(G-1))*((len(d)-1)/(len(d)-rank))
    k=names.index('W');se=float(np.sqrt(max(0,V[k,k])))
    return {'W_beta_per1SD':float(coef[k]),'W_beta_cluster95':[float(coef[k]-1.96*se),float(coef[k]+1.96*se)],'W_sd':float(d.W.std(ddof=0)),'design_rank':int(rank),'n_parameters':X.shape[1],'fixed_traits_conditional_association_not_causal':True}

def summarize_errors(d,y,a,b,mode='H6',seed=20261115):
    err0=(y-a)**2;err1=(y-b)**2
    if mode=='H6':
        gs=sorted(set(d.h6));e0=np.array([err0[d.h6==g].mean() for g in gs]);e1=np.array([err1[d.h6==g].mean() for g in gs])
    else:
        gs=[g for g in sorted(set(d.country_holdout)) if sum(d.country_holdout==g)>=10 and d.loc[d.country_holdout==g,'h6'].nunique()>=3]
        e0=np.array([err0[d.country_holdout==g].mean() for g in gs]);e1=np.array([err1[d.country_holdout==g].mean() for g in gs])
    ix=np.random.default_rng(seed).integers(0,len(gs),size=(3000,len(gs)))
    boot=100*(1-e1[ix].mean(1)/e0[ix].mean(1))
    return {'gain_pct':float(100*(1-e1.mean()/e0.mean())),'fixed_error_bootstrap95':np.quantile(boot,[.025,.975]).tolist(),'baseline_MSE':float(e0.mean()),'with_W_MSE':float(e1.mean()),'row_weighted_gain_pct':float(100*(1-err1.mean()/err0.mean())),'evaluation_groups':[str(g) for g in gs]}

def model(d,ycol,variant='primary',kind='record',family_adjust=True,extras=(),save_prediction=None,country=True):
    d=d.reset_index(drop=True);y=d[ycol].to_numpy(float)
    if not np.isfinite(y).all():raise ValueError('Missing outcome; report rather than silently omit')
    cols=numeric_columns(variant,kind,family_adjust,extras)
    a,sp=predict(d,y,cols);b,_=predict(d,y,cols+['W'])
    bi,_=predict(d,y,cols+['W'],independent=True);diff=float(max(abs(b-bi)))
    if diff>1e-8:raise AssertionError('Independent Ridge fit disagrees')
    rs={'n_rows':len(d),'H7':int(d.h.nunique()),'n_H6':int(d.h6.nunique()),'country_labels':int(d.country.nunique()),'W_positive_rows':int(sum(d.W>0)),
        'baseline_numeric_fields':cols,'outcome':ycol,'variant':variant,'family_adjust':family_adjust,'kind':kind,
        'H6':summarize_errors(d,y,a,b),'beta':beta_cluster(d,y,cols),'H6_splits':sp,'independent_prediction_max_absdiff':diff}
    if country:
        ca,cs=predict(d,y,cols,'country');cb,_=predict(d,y,cols+['W'],'country')
        rs['country']=summarize_errors(d,y,ca,cb,'country');rs['country_splits']=cs
        rs['country_results']=[{'country':g,'n':int(sum(d.country_holdout==g)),'H6':int(d.loc[d.country_holdout==g,'h6'].nunique()),'baseline_MSE':float(np.mean((y[d.country_holdout==g]-ca[d.country_holdout==g])**2)),'with_W_MSE':float(np.mean((y[d.country_holdout==g]-cb[d.country_holdout==g])**2))} for g in sorted(set(d.country_holdout))]
    if save_prediction is not None:
        pp=d[['h','period','h6','country','country_holdout','W']].copy();pp['outcome']=y;pp['H6_baseline_prediction']=a;pp['H6_W_prediction']=b
        if country:pp['country_baseline_prediction']=ca;pp['country_W_prediction']=cb
        pp.to_csv(save_prediction,index=False)
    return rs
