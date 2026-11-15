"""Independent aggregate/groupby, sklearn design, WLS cluster and source verification."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
import statsmodels.api as sm
ROOT=Path(__file__).resolve().parents[1];O=ROOT/'output'

def main():
    d=pd.read_csv(O/'H7_454_HNS_panel.csv',dtype={'h':str,'h6':str},keep_default_na=False)
    scores=pd.read_csv(ROOT/'inputs/gesui_303_host_nutrient_scores.csv',dtype={'butterflyKey':str})
    long=pd.read_csv(O/'species_records_with_HNS.csv',dtype={'h':str,'h6':str,'butterflyKey':str})
    ok=(scores.taxon_status=='ACCEPTED') & ~scores.actual_used_alias_conflict.astype(bool)
    scores['verified_score']=scores.host_nutrient_mean.where(ok)
    scored=long[['h','period','butterflyKey','n_records']].merge(scores[['butterflyKey','verified_score']],on='butterflyKey',how='left',validate='many_to_one')
    scored=scored[scored.verified_score.notna()].copy();scored['weighted_score']=scored.n_records*scored.verified_score
    gg=scored.groupby(['h','period']).agg(total=('n_records','sum'),weighted=('weighted_score','sum'),incidence=('verified_score','mean'))
    gg['mean']=gg.weighted/gg.total;ix=pd.MultiIndex.from_frame(d[['h','period']]);ind=gg.loc[ix]
    ae=float(max(abs(ind['mean'].to_numpy()-d.primary_record_mean)));ai=float(max(abs(ind.incidence.to_numpy()-d.primary_incidence_mean)))
    assert ae<1e-12 and ai<1e-12
    base=json.load(open(O/'model_results.json'));diag=json.load(open(O/'diagnostic_results.json'))['models'];res={**base,**diag}
    geo=pd.read_csv(O/'source_polygon_centroids.csv',dtype={'h':str});gd=d.merge(geo,on='h',validate='many_to_one',sort=False);gd['lon2']=gd.longitude**2;gd['lat2']=gd.latitude**2;gd['lonlat']=gd.longitude*gd.latitude
    pair=pd.read_csv(O/'paired_222_H7_HNS_changes.csv',dtype={'h':str,'h6':str},keep_default_na=False)
    details=[];maxerr=0.;maxbeta=0.;maxci=0.
    for name,r in res.items():
        if 'status' in r:continue
        pred=pd.read_csv(O/(name+'_predictions.csv'),dtype={'h':str,'h6':str},keep_default_na=False)
        source=pair.copy() if name.startswith('temporal_') else gd
        if name.startswith('temporal_'):
            kind=r['kind'];initial=d[d.period=='2017-2019'].set_index('h')['primary_'+kind+'_mean'];source['initial_HNS']=source.h.map(initial)
        lookup=source.set_index(['h','period']);x=lookup.loc[pd.MultiIndex.from_frame(pred[['h','period']])].reset_index();assert len(x)==r['n_rows']
        y=x[r['outcome']].to_numpy(float);np.testing.assert_allclose(y,pred.outcome,atol=1e-12)
        cols=r['baseline_numeric_fields']
        for mode in ['H6']+(['country'] if 'country' in r else []):
            groups=np.array([int(hashlib.sha256(h.encode()).hexdigest()[:8],16)%5 for h in x.h6]) if mode=='H6' else x.country_holdout.to_numpy()
            for add in [False,True]:
                pp=np.full(len(x),np.nan)
                for g in sorted(set(groups)):
                    te=groups==g;tr=~te
                    if mode=='country':tr &= ~x.h6.isin(x.loc[te,'h6'].unique()).to_numpy()
                    assert not set(x.loc[tr,'h6'])&set(x.loc[te,'h6'])
                    numeric=cols+(['W'] if add else []);scaler=StandardScaler().fit(x.loc[tr,numeric])
                    X=[scaler.transform(x[numeric])]
                    for c in ['period','method']+(['country'] if mode=='H6' else []):
                        # Category codes with training-only vocabulary; unknown test labels all zero.
                        cats=sorted(set(x.loc[tr,c].astype(str)));categorical=pd.Categorical(x[c].astype(str),categories=cats)
                        X.append(pd.get_dummies(categorical,dtype=float).to_numpy())
                    X=np.column_stack(X);freq=x.loc[tr,'h6'].value_counts();w=x.loc[tr,'h6'].map(1/freq).to_numpy();w=w/w.mean()
                    fit=Ridge(alpha=10.,fit_intercept=True,solver='svd').fit(X[tr],y[tr],sample_weight=w);pp[te]=fit.predict(X[te])
                target=mode+('_W_prediction' if add else '_baseline_prediction');diff=float(max(abs(pp-pred[target])));maxerr=max(maxerr,diff);assert diff<1e-8,(name,mode,diff)
        # Independently construct unpenalized weighted OLS and uncorrected cluster sandwich.
        nc=cols+['W'];num=StandardScaler().fit_transform(x[nc]);X=[np.ones((len(x),1)),num]
        for c in ['period','method','country']:X.append(pd.get_dummies(x[c].astype(str),dtype=float).to_numpy())
        X=np.column_stack(X);freq=x.h6.value_counts();w=x.h6.map(1/freq).to_numpy();w=w/w.mean()
        ols=sm.WLS(y,X,weights=w).fit();cr=ols.get_robustcov_results(cov_type='cluster',groups=x.h6,use_correction=False)
        rank=np.linalg.matrix_rank(X*np.sqrt(w)[:,None]);G=x.h6.nunique();corr=G/(G-1)*(len(x)-1)/(len(x)-rank)
        j=1+nc.index('W');beta=float(ols.params[j]);se=float(np.sqrt(cr.cov_params()[j,j]*corr));ci=np.array([beta-1.96*se,beta+1.96*se])
        err=abs(beta-r['beta']['W_beta_per1SD']);ec=float(max(abs(ci-np.array(r['beta']['W_beta_cluster95']))));maxbeta=max(maxbeta,err);maxci=max(maxci,ec);assert err<1e-9 and ec<1e-8,(name,err,ec)
        assert r['n_H6']==x.h6.nunique()
        details.append({'analysis':name,'rows':len(x),'H6':G,'independent_prediction_checked':True,'WLS_cluster_checked':True})
    # Strong statement prevented: lower/upper uniform missing-score slopes are not independent tests.
    identity=d.missing_hosts_high_record_mean-d.missing_hosts_low_record_mean-8*(1-d.primary_host_completeness)
    assert max(abs(identity))<1e-12
    out={'all454_aggregate_max_error_record':ae,'all454_aggregate_max_error_incidence':ai,'independent_sklearn_design_SVD_predictions_max_error':maxerr,'independent_statsmodels_WLS_beta_max_error':maxbeta,'independent_cluster_interval_max_error':maxci,'model_specifications_checked':len(details),'uniform_missing_score_identity_max_error':float(max(abs(identity))),'no_new_independent_ecological_sampling':True,'specifications':details}
    (O/'INDEPENDENT_QA.json').write_text(json.dumps(out,indent=2,allow_nan=False));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
