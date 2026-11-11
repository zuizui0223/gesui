#!/usr/bin/env python3
"""Post-source-freeze exploratory river-only functional proxies against measured
EPA NRSA2018-19 TN/TP. NO WWTP attribution, emergence flux, butterflies or
causal chain. Independent HUC8-block OOS prediction & clustered regression.

Inputs: pinned source-only native UID manifest (from GH Actions run #43034107486)
and exact 3 SHA256-pinned EPA CSVs. No source-selection decisions on effect sign.
Output: per-contrast aggregated results; no raw organism/site rows are committed.
"""
import collections,csv,hashlib,json,math,os,pathlib,time,urllib.request
import numpy as np
import pandas as pd
source_dir=pathlib.Path(os.environ.get('RUNNER_TEMP','/tmp'))/'nrsa1819_functional'
source_dir.mkdir(exist_ok=True)
manifest_path=pathlib.Path(os.environ['NRSA_FROZEN_MANIFEST'])
EXPECTED_MANIFEST_SHA='fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91'
EXPECTED={
 'benthos':('https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_benthic_macroinvertebrate_count_-_data.csv','4e9b3a0798b25f67c7a96d4800b94f86544bbd623abaf119358e067b158324a7'),
 'sites':('https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_SiteInfo.csv','1b5b2119c881290191edeb47f1d80a7711f10440e9851574cfa86fd57be770e3'),
 'water':('https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_water_chemistry_chla_-_data.csv','cb9a3349ae3708834eab533e71423d74c1778bf76aba0737f0ae1bcc8f50baa1')
}
def sha(blob):return hashlib.sha256(blob).hexdigest()
if sha(manifest_path.read_bytes())!=EXPECTED_MANIFEST_SHA:
    raise RuntimeError('Pre-outcome-frozen UID site manifest SHA mismatch')
manifest=json.loads(manifest_path.read_text())
m=manifest['primary_berw_UID_subset']
if len(m)!=1221 or len(set(z['uid'] for z in m))!=1221:
    raise RuntimeError('UID pre-outcome freeze sample has changed')
uids={z['uid'] for z in m}
for name,(url,expected) in EXPECTED.items():
    path=source_dir/(name+'.csv')
    req=urllib.request.Request(url,headers={'User-Agent':'gesui-public-ecology-repro/1.0','Accept':'text/csv,*/*'})
    last=None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req,timeout=75) as r,path.open('wb') as f:
                while chunk:=r.read(1048576):f.write(chunk)
            if sha(path.read_bytes())!=expected:
                raise ValueError('Official source changed from prefreeze '+name)
            break
        except Exception as e:
            last=e;path.unlink(missing_ok=True);time.sleep(min(2**attempt,8))
    else:raise RuntimeError('Unable to reconstruct original pinned EPA source '+name+' '+str(last))
data={k:pd.read_csv(source_dir/(k+'.csv'),dtype=str,keep_default_na=False,low_memory=False)
      for k in EXPECTED}
B=data['benthos'];W=data['water'];S=data['sites']
B=B[B['UID'].isin(uids)&(B['SAMPLE_TYPE']=='BERW')].copy()
W=W[W['UID'].isin(uids)].copy()
S=S[S['UID'].isin(uids)].copy()
if len(W)!=1221 or len(S)!=1221 or len(set(W.UID))!=1221 or len(set(S.UID))!=1221:
    raise RuntimeError('Source numeric join not exactly one-to-one by pinned UID')
for x in (B,W,S):x['UID']=x['UID'].astype(str)
B['c']=pd.to_numeric(B['TOTAL300'],errors='coerce').fillna(0)
if (B.c<0).any():raise RuntimeError('Unexpected negative nonzero benthic count')
B=B[B.c>0].copy()
# One named taxon unit per uid unless full duplicate case is flagged rather than mass summed.
n_dup_taxon=int(B.duplicated(['UID','TAXA_ID']).sum())
B['order']=B.ORDER.str.strip().str.lower()
B['family']=B.FAMILY.str.strip().str.lower()
B['ffg']=B.FFG.str.strip().str.upper()
B['class']=B['CLASS'].str.strip().str.lower()
B['is_ept']=B.order.isin(['ephemeroptera','plecoptera','trichoptera'])
B['is_chironomid']=B.family.eq('chironomidae')
B['is_odonata']=B.order.eq('odonata')
B['is_aquatic_insect']=B['class'].eq('insecta')
B['is_shredder']=B.ffg.eq('SH')
B['is_gatherer']=B.ffg.eq('CG')
B['is_scraper']=B.ffg.eq('SC')
B['is_filterer']=B.ffg.eq('CF')
grp=B.groupby('UID')
frame=pd.DataFrame({'total_subsample':grp['c'].sum(),'taxon_nonzero_rows':grp.size()})
for name in ['is_ept','is_chironomid','is_odonata','is_aquatic_insect','is_shredder','is_gatherer','is_scraper','is_filterer']:
    frame[name+'_individuals']=(B['c']*B[name].astype(int)).groupby(B.UID).sum()
frame=frame.reindex(sorted(uids),fill_value=0)
if (frame.total_subsample<=0).any():
    raise RuntimeError('Frozen UID has no proper count; do not fit a negative-zero missing site')
for name in ['ept','chironomid','odonata','aquatic_insect','shredder','gatherer','scraper','filterer']:
    frame[name+'_fraction']=frame['is_'+name+'_individuals']/frame.total_subsample
W=W.set_index('UID').reindex(frame.index)
S=S.set_index('UID').reindex(frame.index)
M=pd.DataFrame(m).set_index('uid').reindex(frame.index)
if M.isna().any().any():raise RuntimeError('Manifest has missing source attributes')
for target in ('NTL_RESULT','PTL_RESULT','COND_RESULT'):
    W[target+'_num']=pd.to_numeric(W[target],errors='coerce')
if not np.isfinite(W['NTL_RESULT_num']).all() or not np.isfinite(W['PTL_RESULT_num']).all():
    raise RuntimeError('Original source-pinned nutrient values missing')
if (W['NTL_RESULT_num']<=0).any() or (W['PTL_RESULT_num']<=0).any():
    raise RuntimeError('Original source-eligible N/P must remain positive')
if set(W.NTL_UNITS) != {'mg N/L'} or set(W.PTL_UNITS) != {'ug/L'}:
    raise RuntimeError('Detected new nutrient unit or ambiguous conversion')
# TP micrograms/L -> mg P/L exactly.
frame['tn_mgL']=W['NTL_RESULT_num']
frame['tp_mgL']=W['PTL_RESULT_num']/1000
frame['log_tn']=np.log10(frame.tn_mgL)
frame['log_tp']=np.log10(frame.tp_mgL)
frame['log_cond']=np.where(W.COND_RESULT_num>0,np.log10(W.COND_RESULT_num),np.nan)
frame['state']=M.state
frame['eco9']=M.eco9
frame['year']=pd.to_numeric(M.year,errors='coerce')
frame['site_id']=M.site_id
frame['huc8']=M.huc8
frame['lat']=pd.to_numeric(M.lat,errors='coerce')
frame['lon']=pd.to_numeric(M.lon,errors='coerce')
frame['stream_order']=pd.to_numeric(S['STRAH_ORD'],errors='coerce')
# No posthoc imputation of missing covariates. Fixed complete-case subset for
# all outcomes and models including N/P and background region/landcover.
pred_cols=['log_tn','log_tp','lat','lon','year','stream_order','log_cond']
eligible=frame.dropna(subset=pred_cols).copy()
if len(eligible)<650 or eligible.state.nunique()<25 or eligible.eco9.nunique()<6:
    raise RuntimeError('Predefined numerical/river covariates collapse the source sample')
# Calculate fraction of original 300-individual subsamples in a plausible range
# to identify if duplicated taxon labels could create compositional numerator >1.
source_QA={
 'matched_pinned_source_visits':len(frame),'same_continental_US_eligible_covariate_cases':len(eligible),
 'n_PHYSICAL_site_id':eligible.site_id.nunique(),'n_HUC8':eligible.huc8.nunique(),
 'n_states':eligible.state.nunique(),'n_eco9':eligible.eco9.nunique(),
 'duplicate_UID_TAXA_ID_positive_taxon_rows':n_dup_taxon,
 'total300_per_sample_quantiles':{str(q):float(frame.total_subsample.quantile(q)) for q in [0,.05,.5,.95,1]},
 'original_nutrient_NARS_FLAG_distributions':{
     'TN':W.NTL_NARS_FLAG.value_counts(dropna=False).to_dict(),
     'TP':W.PTL_NARS_FLAG.value_counts(dropna=False).to_dict()},
 'original_nutrient_QA_FLAG_distributions':{
     'TN':W.NTL_QA_FLAG.value_counts(dropna=False).to_dict(),
     'TP':W.PTL_QA_FLAG.value_counts(dropna=False).to_dict()},
 'N_P_paired_log_pearson':float(eligible[['log_tn','log_tp']].corr().iloc[0,1]),
 'FFG_codes_positive_mass':B.groupby('ffg')['c'].sum().sort_values(ascending=False).head(14).to_dict(),
 'order_codes_positive_mass':B.groupby('order')['c'].sum().sort_values(ascending=False).head(15).to_dict(),
}
# Don't silently model with duplicate taxon units. >=0.5% duplicate rate stops.
if n_dup_taxon>len(B)*.005: raise RuntimeError('Too many duplicated UID/TAXA_ID taxon rows; no model')
# Missingness is strictly complete case, quantified.
# Continuous predictors standardized on TRAIN data in spatial CV, not on test.
ykeys=['ept_fraction','chironomid_fraction','odonata_fraction',
       'aquatic_insect_fraction','shredder_fraction','gatherer_fraction',
       'scraper_fraction','filterer_fraction']
X_cat=pd.get_dummies(eligible[['state','eco9']],drop_first=True,dtype=float)
X_cont=eligible[['lat','lon','year','stream_order','log_cond']].astype(float)
X0=pd.concat([X_cont,X_cat],axis=1)
X0.insert(0,'intercept',1.)
X0=X0.to_numpy(dtype=float)
plus={'baseline':[],'TN':['log_tn'],'TP':['log_tp'],'TN+TP':['log_tn','log_tp']}
Xdict={model:np.column_stack([X0]+[eligible[c].to_numpy(dtype=float) for c in covars])
       for model,covars in plus.items()}
def fixfold(huc):
    return int(hashlib.sha256(str(huc).encode()).hexdigest()[:8],16)%5
folds=np.array([fixfold(x) for x in eligible.huc8])
def ridge_predict(X,train,test,alpha=10):
    """Penalty on all standardized non-intercept columns, no test leakage."""
    arr=X.copy()
    mu=arr[train,1:].mean(axis=0);sd=arr[train,1:].std(axis=0)
    sd=np.where(sd<1e-12,1,sd)
    arr[:,1:]=(arr[:,1:]-mu)/sd
    z=arr[train];zt=arr[test]
    reg=np.eye(z.shape[1])*alpha;reg[0,0]=0
    return zt@np.linalg.solve(z.T@z+reg,z.T@Y[train])
def fit_ols_cluster(X,y,groups):
    # Full rank may be lower than columns (state/eco dummies). Pseudoinverse
    # and cluster sums are descriptive posthoc, not preregistered causality.
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    u=y-X@beta
    meat=np.zeros((X.shape[1],X.shape[1]))
    for g in sorted(set(groups)):
        ix=(groups==g)
        v=X[ix].T@u[ix]
        meat+=np.outer(v,v)
    bread=np.linalg.pinv(X.T@X)
    V=bread@meat@bread
    G=len(set(groups));N=len(y);K=np.linalg.matrix_rank(X)
    if G>1 and N>K:V*=G/(G-1)*(N-1)/(N-K)
    return beta,np.sqrt(np.maximum(np.diag(V),0)),N,G,K
outputs={}
blocks=eligible.huc8.astype(str).to_numpy()
for outcome in ykeys:
    Y=eligible[outcome].to_numpy(dtype=float)
    preds={m:np.empty(len(Y),dtype=float) for m in plus}
    for fold in range(5):
        test=folds==fold;train=~test
        for name,X in Xdict.items():
            preds[name][test]=ridge_predict(X,train,test)
    err={name:float(np.mean((Y-p)**2)) for name,p in preds.items()}
    res={'n':len(Y),'observed_fraction_mean':float(np.mean(Y)),
         'baseline_huc8_disjoint_5fold_MSE':err['baseline'],
         'mse':err,'model_results':{}}
    for name in ['TN','TP','TN+TP']:
        res['model_results'][name]={
         'heldout_MSE_gain_pct':100*(err['baseline']-err[name])/err['baseline'],
         'within_site_EQR_link_not_tested':True,
        }
        X=Xdict[name];beta,se,N,G,K=fit_ols_cluster(X,Y,blocks)
        slopes={}
        for j,v in enumerate(plus[name]):
            ii=X.shape[1]-len(plus[name])+j
            slopes[v]={
             'linear_fraction_beta_per_10x_nutrient':float(beta[ii]),
             'huc8_cluster_descriptive95':[float(beta[ii]-1.96*se[ii]),float(beta[ii]+1.96*se[ii])]}
        res['model_results'][name]['posthoc_partial_slopes']=slopes
    outputs[outcome]=res
report={
 'schema':'NRSA_2018_2019_AQUATIC_FUNCTIONAL_PROXY_EXPLORATORY_v1',
 'source_freeze_manifest_sha256':EXPECTED_MANIFEST_SHA,
 'pre_source_cohort_SHAs':{n:v[1] for n,v in EXPECTED.items()},
 'source_only_original_composite_gate_fail_retained':True,
 'corrected_native_UID_source_gate_pass':True,
 'not_2018_2019_measured_sewage_WWTP_mass':True,
 'not_causal':True,'not_emergence_or_pollination_service':True,
 'not_a_multitaxon_water_land_mediation_result':True,
 'dataset_qa':source_QA,
 'functional_proxy_results':outputs,
 'interpretation':[
  'All outcomes use normalized larval/benthic counted fractions of a standardized sample, not emergence flux, pollination, biodiversity maintenance or trophic transfer.',
  'Matched site/visit and both measured concentrations are EPA NRSA, not traced to municipal wastewater sources.',
  'NRSA random spatial sample has revisits: 1,221 UID visits from 1,125 distinct physical sites; HUC8-folds and river blocks limit but do not remove spatial confounding.',
  'TP reported in ug P/L; conversion to mg P/L /1000 checked, log slopes are per 10-fold exposure.',
  'Source-native UID correction and regression are posthoc exploratory; state/ecoregion, river order, conductivity, geography and year are covariates.',
  'No plants or butterflies at these exact river sites have yet been collected/verified.'
 ]}
output_path=source_dir/'NRSA1819_RIVER_FUNCTIONAL_PROXY_EXPLORATORY_REPORT.json'
output_path.write_text(json.dumps(report,indent=2,sort_keys=True))
# Keep stdout concise: input QC plus each out-of-site coefficient and predictive gain.
print('NRSA1819_FUNCTIONAL_PROXY_RESULT_START')
short={'source_manifest_sha':EXPECTED_MANIFEST_SHA,'dataset_qa':source_QA,
 'results':{k:{'mean':v['observed_fraction_mean'],
     'models':{m:{'mse_gain':v['model_results'][m]['heldout_MSE_gain_pct'],
       'slopes':v['model_results'][m]['posthoc_partial_slopes']} for m in ['TN','TP','TN+TP']}}
     for k,v in outputs.items()}}
print(json.dumps(short,sort_keys=True,default=float))
print('NRSA1819_FUNCTIONAL_PROXY_RESULT_END')
print('Original-data result JSON: '+str(output_path))
