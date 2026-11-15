"""Post-result sensitivity of selected Pieridae signal, not confirmation."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,sys,pathlib,hashlib
import numpy as np,pandas as pd
from robustness_core import hellinger,group_weights
R=pathlib.Path(__file__).resolve().parents[1];OLD=R/'previous/gesui_host_network_20261114';OUT=R/'output'
sys.path.insert(0,str(OLD/'scripts'))
import analyse as old
from network import inclusion

def compare(d,col,extra=(),countries=True):
    d=d.reset_index(drop=True);y=d[col].to_numpy(float)
    if not np.isfinite(y).all():raise ValueError('Missing outcome cannot be filled')
    a,split=old.predict(d,y,extra);b,_=old.predict(d,y,tuple(extra)+('W',))
    check,_=old.predict(d,y,tuple(extra)+('W',),independent=True)
    h=d.h6.to_numpy();hs=sorted(set(h));e0=(y-a)**2;e1=(y-b)**2
    g0=np.array([e0[h==k].mean() for k in hs]);g1=np.array([e1[h==k].mean() for k in hs])
    idx=np.random.default_rng(20261114).integers(0,len(hs),(3000,len(hs)))
    boost=100*(g0[idx].mean(1)-g1[idx].mean(1))/g0[idx].mean(1)
    X,_=old.design(d,np.ones(len(d),bool),np.zeros(len(d),bool),tuple(extra)+('W',),True)
    beta=np.linalg.lstsq(X,y,rcond=None)[0];res=y-X@beta;bread=np.linalg.pinv(X.T@X);meat=np.zeros_like(bread)
    for k in hs:
        v=X[h==k].T@res[h==k];meat+=np.outer(v,v)
    rank=np.linalg.matrix_rank(X);var=bread@meat@bread*len(hs)/(len(hs)-1)*(len(d)-1)/(len(d)-rank)
    wi=13+len(extra);se=np.sqrt(max(0,var[wi,wi]))
    out={'n_rows':len(d),'H7':d.h.nunique(),'H6':len(hs),'H6_gain_pct':100*(g0.mean()-g1.mean())/g0.mean(),
         'H6_ci':np.quantile(boost,[.025,.975]).tolist(),'W_beta_perSD':beta[wi],
         'W_cluster95':[beta[wi]-1.96*se,beta[wi]+1.96*se],
         'max_sklearn_diff':float(np.max(abs(b-check))),'splits':split}
    if countries:
        ac,ss=old.predict(d,y,extra,mode='country');bc,_=old.predict(d,y,tuple(extra)+('W',),mode='country')
        country=d.country_holdout.to_numpy();cs=[]
        for c in sorted(set(country)):
            m=country==c
            if m.sum()>=10 and d.loc[m,'h6'].nunique()>=3:
                cs.append({'country':c,'n':int(m.sum()),'m0':float(((y[m]-ac[m])**2).mean()),'m1':float(((y[m]-bc[m])**2).mean())})
        c0=np.array([v['m0'] for v in cs]);c1=np.array([v['m1'] for v in cs]);ix=np.random.default_rng(20261114).integers(0,len(cs),(3000,len(cs)))
        cb=100*(c0[ix].mean(1)-c1[ix].mean(1))/c0[ix].mean(1)
        out['country_equal_gain_pct']=100*(c0.mean()-c1.mean())/c0.mean()
        out['country_fixed_error_bootstrap95']=np.quantile(cb,[.025,.975]).tolist()
        out['n_supported_countries']=len(cs);out['country_results']=cs;out['country_splits']=ss
    return out

def main():
    z=np.load(OLD/'output/candidate_network_matrices.npz');A=z['A'];PB=z['PB'];B=z['B'];P=z['P'];NP=z['NP']
    d=pd.read_csv(OUT/'H7_454_geographic_null_network_panel.csv',dtype={'h':str,'h6':str},keep_default_na=False)
    taxa=json.loads((OLD/'source/butterfly_usage_aliases.json').read_text());td={t['original_speciesKey']:t for t in taxa}
    mask=np.array([td[k]['family']=='Pieridae' for k in z['bkeys']]);known=mask&(A.sum(1)>0)
    H=hellinger(B[:,mask]);cols=[]
    for j,k in enumerate(z['bkeys'][mask]):
        c='Pieridae_mix_'+k;d[c]=H[:,j];cols.append(c)
    d['log_Pieridae_records']=np.log1p(B[:,mask].sum(1))
    d['Pieridae_catalogue_coverage']=PB[:,known].sum(1)/PB[:,mask].sum(1)
    q=np.array([inclusion(n,counts,50) for n,counts in zip(NP,P@A[known].T)])
    freq=group_weights(d.h6.to_numpy())@PB[:,known];freq/=freq.sum()
    d['Pieridae_fixed_mix_host_support']=q@freq
    # This fixed-reference outcome is descriptive standardization of the whole source population.
    # No assertion of external ecological replication or estimated causal mediation.
    cases=[];subset=d.loc[d.minimum_species_coverage>=.8].copy();common=['log_Pieridae_records','Pieridae_catalogue_coverage']
    for label,extra in [('original',[]),('family_recording_coverage',common),('species_composition_adjusted',common+cols)]:
        qout=compare(subset,'Pieridae_host_support',tuple(extra),True);cases.append({'test':label,'cohort':'coverage80',**qout});print(label,round(qout['H6_gain_pct'],4),qout['H6_ci'],flush=True)
    qout=compare(subset,'Pieridae_fixed_mix_host_support',tuple(common),True);cases.append({'test':'fixed_butterfly_mix_plant_opportunity','cohort':'coverage80',**qout});print('fixed_mix',qout['H6_gain_pct'],flush=True)
    for name,ma in [('period2017_2019',subset.period=='2017-2019'),('period2020_2022',subset.period=='2020-2022'),('France',subset.country=='FR'),('outside_France',subset.country!='FR')]:
        ds=subset.loc[ma].copy();qout=compare(ds,'Pieridae_host_support',tuple(common+cols),False)
        cases.append({'test':'species_composition_adjusted','cohort':name,**qout});print(name,round(qout['H6_gain_pct'],4),qout['H6_ci'],flush=True)
    # Independent reproduction of original exploratory score.
    prior=json.loads((OUT/'model_results.json').read_text())['coverage80']['Pieridae_host_support']
    assert abs(cases[0]['H6_gain_pct']-prior['H6_gain_pct'])<1e-9
    (OUT/'pieridae_sensitivity.json').write_text(json.dumps(cases,indent=2,allow_nan=False))
    rows=[]
    for x in cases:
        rows.append({k:x[k] for k in ['test','cohort','n_rows','H7','H6','H6_gain_pct','W_beta_perSD','max_sklearn_diff']}|{
            'H6_CI_low':x['H6_ci'][0],'H6_CI_high':x['H6_ci'][1],
            'country_equal_gain_pct':x.get('country_equal_gain_pct'),
            'country_CI_low':x.get('country_fixed_error_bootstrap95',[None,None])[0],
            'country_CI_high':x.get('country_fixed_error_bootstrap95',[None,None])[1],
            'W_beta_CI_low':x['W_cluster95'][0],'W_beta_CI_high':x['W_cluster95'][1]})
    pd.DataFrame(rows).to_csv(OUT/'pieridae_sensitivity.csv',index=False)
    d.to_csv(OUT/'H7_454_Pieridae_composition_sensitivity.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__':main()
