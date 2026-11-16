"""Run frozen community-HNS analyses; all outcomes and failed sensitivities retained."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import json,numpy as np,pandas as pd
from scipy.stats import spearmanr
from model_core import model,groupweights,numeric_columns,design
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'output'

def clean(v):
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,list):return [clean(x) for x in v]
    if isinstance(v,(float,np.floating)) and not np.isfinite(v):return None
    return v

def main():
    d=pd.read_csv(OUT/'H7_454_HNS_panel.csv',dtype={'h':str,'h6':str},keep_default_na=False)
    results={};summary=[]
    specs=[]
    for kind in ['record','incidence']:
        specs.append((f'primary_{kind}',d,'primary',kind,True,True))
        specs.append((f'no_family_adjust_{kind}',d,'primary',kind,False,True))
    specs.append(('source80',d[d.minimum_species_coverage>=.8],'primary','record',True,True))
    for period in sorted(set(d.period)):specs.append(('period_'+period,d[d.period==period],'primary','record',True,True))
    specs.extend([('FR',d[d.country=='FR'],'primary','record',True,False),('nonFR',d[d.country!='FR'],'primary','record',True,True)])
    for name in ['host80','strict_names','no_wikipedia','drop_Aglais_io','missing_hosts_low','missing_hosts_high']:
        specs.append((name,d,name,'record',True,True))
    for label,x,v,k,fa,co in specs:
        outcol=f'{v}_{k}_mean'
        try:
            r=model(x,outcol,v,k,fa,save_prediction=OUT/(label+'_predictions.csv'),country=co)
        except RuntimeError as e:
            results[label]={'status':'failed_support','reason':str(e)};print(label,'FAILED SUPPORT',str(e),flush=True);continue
        results[label]=r
        ss={'analysis':label,'n':r['n_rows'],'H7':r['H7'],'H6':r['n_H6'],'H6_gain_pct':r['H6']['gain_pct'],'H6_low':r['H6']['fixed_error_bootstrap95'][0],'H6_high':r['H6']['fixed_error_bootstrap95'][1],'W_beta':r['beta']['W_beta_per1SD'],'beta_low':r['beta']['W_beta_cluster95'][0],'beta_high':r['beta']['W_beta_cluster95'][1],'country_gain':r.get('country',{}).get('gain_pct'),'country_low':r.get('country',{}).get('fixed_error_bootstrap95',[None,None])[0],'country_high':r.get('country',{}).get('fixed_error_bootstrap95',[None,None])[1]}
        summary.append(ss);print(ss,flush=True)
        (OUT/'model_results.json').write_text(json.dumps(clean(results),indent=2,allow_nan=False));pd.DataFrame(summary).to_csv(OUT/'model_summary.csv',index=False)
    # Trait-label diagnostics, source fixed, exchangeability not phylogenetically guaranteed.
    z=np.load(OUT/'trait_matrices.npz');B=z['B'];s=z['scores_primary'];known=np.isfinite(s);f=z['family'];s0=np.where(known,s,0.)
    traitdiag={};driver_tables=[]
    for kind in ['record','incidence']:
        weights=B*known if kind=='record' else (B>0)*known
        proportions=weights/weights.sum(1)[:,None]
        cols=numeric_columns('primary',kind,True);X,_=design(d,np.ones(len(d),bool),cols,True);w=groupweights(d.h6)
        W=(d.W.to_numpy()-d.W.mean())/d.W.std(ddof=0)
        cf=np.linalg.lstsq(X*np.sqrt(w)[:,None],W*np.sqrt(w),rcond=None)[0];rw=W-X@cf
        coeff=(w*rw)@proportions/np.sum(w*rw*rw)
        obs=float(coeff@s0);rng=np.random.default_rng(20261115);null=[]
        groups=[np.flatnonzero((f==fam)&known) for fam in sorted(set(f))]
        for it in range(999):
            perm=s0.copy()
            for ix in groups:perm[ix]=rng.permutation(s0[ix])
            null.append(float(coeff@perm))
        traitdiag[kind]={'observed_beta':obs,'null95':np.quantile(null,[.025,.975]).tolist(),'upper_tail_fraction_plus1':float((1+np.sum(np.array(null)>=obs))/1000),'two_sided_absolute_fraction_plus1':float((1+np.sum(abs(np.array(null))>=abs(obs)))/1000),'permutations':999,'species_permutation_within_family_not_phylogenetic_test':True}
        np.testing.assert_allclose(obs,results[f'primary_{kind}']['beta']['W_beta_per1SD'],atol=1e-9)
        means={g:np.mean(s[(f==g)&known]) for g in sorted(set(f))};centered=np.array([s[i]-means[f[i]] if known[i] else 0 for i in range(len(s))])
        dt=pd.DataFrame({'butterflyKey':z['keys'],'kind':kind,'family':f,'HNS':s,'conditional_W_record_composition_coefficient':coeff,'centered_HNS_slope_contribution':coeff*centered,'source_records':B.sum(0)})
        np.testing.assert_allclose(dt.centered_HNS_slope_contribution.sum(),obs,atol=1e-9)
        driver_tables.append(dt)
    (OUT/'trait_label_diagnostics.json').write_text(json.dumps(clean(traitdiag),indent=2,allow_nan=False));pd.concat(driver_tables).to_csv(OUT/'species_slope_decomposition.csv',index=False)
    # Descriptive distribution and paired observed change, no wastewater-change exposure.
    desc={'correlation_W_urban':float(spearmanr(d.W,d.urban).statistic),'raw_record_CWM_HNS_W_rho':float(spearmanr(d.W,d.primary_record_mean).statistic),'raw_incidence_CWM_HNS_W_rho':float(spearmanr(d.W,d.primary_incidence_mean).statistic),'H6_equal_record_CWM':float(np.average(d.primary_record_mean,weights=groupweights(d.h6))),'record_CWM_minmax':[float(d.primary_record_mean.min()),float(d.primary_record_mean.max())],'record_CWM_median':float(d.primary_record_mean.median()),'score_coverage_raw_minmax':[float(d.primary_record_coverage_raw.min()),float(d.primary_record_coverage_raw.max())],'n_raw_score_coverage80':int(sum(d.primary_record_coverage_raw>=.8))}
    (OUT/'descriptive.json').write_text(json.dumps(desc,indent=2));print('TRAIT_DIAG',json.dumps(traitdiag),flush=True);print('DESCRIPTIVE',json.dumps(desc),flush=True)
if __name__=='__main__':main()
