#!/usr/bin/env python3
"""Explicit POSTHOC sensitivities AFTER first NRSA2018/19 effect signs were known.
Same original SHA-pinned sources and frozen UID; no causal effect. Program is
independently written without importing the primary model script.
"""
import csv,hashlib,json,os,pathlib
import numpy as np
import pandas as pd
p=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"nrsa1819_functional"
orig=json.loads((p/"NRSA1819_RIVER_FUNCTIONAL_PROXY_EXPLORATORY_REPORT.json").read_text())
manifest=json.loads(pathlib.Path(os.environ['NRSA_FROZEN_MANIFEST']).read_text())['primary_berw_UID_subset']
source=pd.DataFrame(manifest).set_index('uid')
A=pd.read_csv(p/'benthos.csv',dtype=str,keep_default_na=False,low_memory=False)
W=pd.read_csv(p/'water.csv',dtype=str,keep_default_na=False,low_memory=False).set_index('UID')
S=pd.read_csv(p/'sites.csv',dtype=str,keep_default_na=False,low_memory=False).set_index('UID')
a=A[(A.UID.isin(source.index))&(A.SAMPLE_TYPE=='BERW')].copy()
a['n']=pd.to_numeric(a.TOTAL300,errors='coerce').fillna(0)
a=a[a.n>0].copy()
a['is_ept']=a.ORDER.str.strip().str.lower().isin(['ephemeroptera','plecoptera','trichoptera'])
group=a.groupby('UID')
N=group.n.sum()
E=(a.n*a.is_ept).groupby(a.UID).sum()
f=pd.DataFrame({'n':N,'ept':E/N})
f=f.reindex(source.index)
W=W.reindex(f.index);S=S.reindex(f.index)
f['site']=source.site_id
f['year']=pd.to_numeric(source.year,errors='coerce')
f['date']=source.date_col
f['eco']=source.eco9
f['state']=source.state
f['huc8']=source.huc8
f['lat']=pd.to_numeric(source.lat,errors='coerce')
f['lon']=pd.to_numeric(source.lon,errors='coerce')
f['order']=pd.to_numeric(S.STRAH_ORD,errors='coerce')
f['cond']=pd.to_numeric(W.COND_RESULT,errors='coerce')
f['TN']=pd.to_numeric(W.NTL_RESULT,errors='coerce')
f['TP']=pd.to_numeric(W.PTL_RESULT,errors='coerce')/1000
f['logTN']=np.log10(f.TN)
f['logTP']=np.log10(f.TP)
f['logcond']=np.log10(f.cond.where(f.cond>0))
f['flags_clear']=(W.NTL_NARS_FLAG.str.strip()=='')&(W.PTL_NARS_FLAG.str.strip()=='')
f=f.replace([np.inf,-np.inf],np.nan).dropna(subset=['ept','logTN','logTP','lat','lon','year','order','logcond'])
def regression(d):
    if len(d)<250 or d.state.nunique()<20:return {'n':len(d),'no_fit':'support low'}
    cats=pd.get_dummies(d[['state','eco']],drop_first=True,dtype=float)
    # Avoid singular columns from one-size data by explicit independent intercept,
    # state and eco; use pseudo-inverse and HUC8 clustered sandwich.
    X=np.column_stack((np.ones(len(d)),d[['lat','lon','year','order','logcond','logTN','logTP']].to_numpy(float),cats.to_numpy(float)))
    y=d.ept.to_numpy(float)
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    u=y-X@beta
    xx=np.linalg.pinv(X.T@X)
    groups=d.huc8.astype(str).to_numpy()
    meat=np.zeros((X.shape[1],X.shape[1]))
    for z in np.unique(groups):
        idx=groups==z;v=X[idx].T@u[idx];meat+=v[:,None]@v[None,:]
    G=len(np.unique(groups));K=np.linalg.matrix_rank(X);n=len(y)
    V=xx@meat@xx * (G/(G-1) * (n-1)/(n-K))
    se=np.sqrt(np.maximum(np.diag(V),0))
    def slope(i):return {'beta':float(beta[i]),'descriptive95':[float(beta[i]-1.96*se[i]),float(beta[i]+1.96*se[i])]}
    return {'n':len(d),'physical_sites':d.site.nunique(),'states':d.state.nunique(),'eco9':d.eco.nunique(),
      'HUC8':d.huc8.nunique(),'mean_ept':float(d.ept.mean()),'TN':slope(6),'TP':slope(7)}
out={
 'nature':'posthoc_independent_manual_numpy_sensitivity_not_untested_prospective_vegetation_or_pollination',
 'source_sha256':orig['pre_source_cohort_SHAs'],
 'frozen_manifest_sha256':orig['source_freeze_manifest_sha256'],
 'scenarios':{},
}
out['scenarios']['all_eligible_original']=regression(f)
out['scenarios']['both_nutrient_nars_flags_blank']=regression(f[f.flags_clear])
out['scenarios']['benthic_subsample_at_least_250']=regression(f[f.n>=250])
out['scenarios']['flags_blank_AND_subsample_at_least_250']=regression(f[f.flags_clear&(f.n>=250)])
# Restrict to single original visit per physical site using *only* frozen pre-
# outcome metadata for deterministic selection. This is an outcome-seen sensitivity.
one=f.sort_values(['site','year','date']).drop_duplicates('site',keep='first')
out['scenarios']['first_site_visit_only']=regression(one)
# A source reporting asymmetry is not 0; composite FFG codes should not be
# treated as absence of feeding function.
ffg=a.FFG.str.strip()
out['source_qa']={'positive_benthic_individuals_with_multiple_FFG_codes':float(a.loc[ffg.str.contains(',',regex=False),'n'].sum()),
  'positive_benthic_individuals_total':float(a.n.sum()),
  'fraction_multi_FFG_codes':float(a.loc[ffg.str.contains(',',regex=False),'n'].sum()/a.n.sum()),
  'flags_clear_n':int(f.flags_clear.sum()),'n_total_eligible':len(f)}
out['independent_original_joint_slopes_match']=bool(abs(out['scenarios']['all_eligible_original']['TN']['beta']-
    orig['functional_proxy_results']['ept_fraction']['model_results']['TN+TP']['posthoc_partial_slopes']['log_tn']['linear_fraction_beta_per_10x_nutrient'])<1e-8)
if not out['independent_original_joint_slopes_match']:
    raise RuntimeError('Independent basic OLS EPT TN/P coefficient does not reproduce saved primary original')
reportpath=p/'NRSA1819_EPT_FLAGS_SUBSAMPLE_NATIVE_SITE_POSTHOC_INDEPENDENT_QA.json'
reportpath.write_text(json.dumps(out,indent=2,sort_keys=True))
print('NRSA1819_EPT_SENSITIVITY_START')
print(json.dumps(out,sort_keys=True))
print('NRSA1819_EPT_SENSITIVITY_END')
