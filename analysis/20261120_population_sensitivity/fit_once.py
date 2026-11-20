"""Frozen final population sensitivity; paired models on exactly the same rows."""
from pathlib import Path
import importlib.util
import hashlib
import json
import subprocess
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
WORK=ROOT.parent
OUT=ROOT/'results/20261120_population_sensitivity'
FROZEN='223da1a'
POP='log1p_upstream_pop_density_2020'
source=ROOT/'analysis/20261115_eu454_plant_nutri/h6_model_if_gate.py'
spec=importlib.util.spec_from_file_location('original_h6_models',source)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)


def digest(p):
    b=p.read_bytes()
    if p.suffix in {'.py','.md'}:b=b.replace(b'\r\n',b'\n')
    return hashlib.sha256(b).hexdigest()


def main():
    specification=Path(__file__).with_name('SPECIFICATION.md')
    original=subprocess.check_output(['git','-C',str(ROOT),'show',FROZEN+':analysis/20261120_population_sensitivity/SPECIFICATION.md'])
    assert original.replace(b'\r\n',b'\n')==specification.read_bytes().replace(b'\r\n',b'\n')
    previous_path=ROOT/'results/20261120_eu451_W_confirmation/result.json'
    previous=json.loads(previous_path.read_text())
    for row in previous['inputs']:
        assert digest(WORK/row['path'])==row['sha256'],row['path']
    gate=json.loads((OUT/'source_receipt.json').read_text())
    if not gate['gate_pass']:
        result={'status':'SOURCE_STOP','source':gate,'population_model_fitted':False}
        (OUT/'result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print('SOURCE_STOP: no outcome model fitted');return
    d=pd.read_csv(WORK/'original/H7_454_HNS_panel.csv',dtype={'h':str,'h6':str})
    e=pd.read_csv(WORK/'plant/EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv',dtype={'h':str})
    e=e[e.full_source_trait_eligible.astype(str)=='True']
    d=d.merge(e[['h','period','NUTRI_coverage_all_Plantae_speciesKeys','exact_NUTRI_scored_vascular_accepted_species']],on=['h','period'],validate='one_to_one')
    s=pd.read_csv(WORK/'plant/EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv',dtype={'h':str})
    d=d.merge(s[['h','period','Plant_N_incidence','n_scored_species']],on=['h','period'],validate='one_to_one')
    assert (len(d),d.h.nunique(),d.h6.nunique())==(451,232,90)
    population=pd.read_csv(OUT/'population_by_original_H7.csv',dtype={'h':str})
    population=population[population.eligible]
    d=d.merge(population[['h',POP]],on='h',validate='many_to_one')
    assert d.h.nunique()>=186 and d.h6.nunique()>=72 and d.country_holdout.nunique()>=8
    d['log1p_eiv_scored_taxa']=np.log1p(d.exact_NUTRI_scored_vascular_accepted_species)
    d['country_method']=d.country.fillna('__UNKNOWN__').astype(str)+'::'+d.method.fillna('__UNKNOWN__').astype(str)
    controls=previous['numeric_controls']
    assert controls==[c for c in old.B_GEOG if c!='W']
    outcome=previous['target'];assert outcome=='primary_incidence_mean'
    results={}
    for name,extra in [('same_sample_original',[]),('population_adjusted',[POP])]:
        base=controls+extra
        results[name]={'association':old.source_cluster_regression(d,outcome,base+['W'],'W'),
                       'prediction':old.cv_h6(d,outcome,base,'W')}
    r=results['population_adjusted']
    positive=r['association']['H6_cluster_descriptive95'][0]>0
    predictive=r['prediction']['fixed_out_of_fold_error_H6_resampling95_pct'][0]>0
    label={(True,True):'positive_association_and_predictive_information',
           (True,False):'positive_association_only',(False,True):'predictive_information_only',
           (False,False):'positive_W_signal_not_supported'}[(positive,predictive)]
    w=old.H6_weights(d);sw=np.sqrt(w)
    X,names=old.design_full(d,controls+[POP,'W'])
    j=names.index('W');Z=np.delete(X,j,axis=1)
    def residual(v,z):return v-z@np.linalg.lstsq(z*sw[:,None],v*sw,rcond=None)[0]
    wx=residual(X[:,j],Z);y=residual(d[outcome].to_numpy(),Z)
    beta=np.sum(w*wx*y)/np.sum(w*wx*wx)
    np.testing.assert_allclose(beta,r['association']['adjusted_beta_per_1SD_focal'],atol=1e-12)
    pre,n=old.design_full(d,controls)
    Wres=residual(d.W.to_numpy(),pre);Pres=residual(d[POP].to_numpy(),pre)
    def corr(a,b):
        a=a-np.average(a,weights=w);b=b-np.average(b,weights=w)
        return float(np.sum(w*a*b)/np.sqrt(np.sum(w*a*a)*np.sum(w*b*b)))
    audit={'weighted_W_population_r':corr(d.W.to_numpy(),d[POP].to_numpy()),
           'weighted_partial_r_after_original_controls':corr(Wres,Pres),
           'fraction_W_residual_variance_removed_by_population':corr(Wres,Pres)**2,
           'full_design_condition_number':float(np.linalg.cond(X)),
           'FWL_beta_check':float(beta),'same_sample_W_SD':float(d.W.std(ddof=0))}
    if len(d)==451:
        np.testing.assert_allclose(results['same_sample_original']['association']['adjusted_beta_per_1SD_focal'],previous['association']['adjusted_beta_per_1SD_focal'],atol=1e-12)
        np.testing.assert_allclose(results['same_sample_original']['prediction']['H6_equal_added_MSE_gain_pct'],previous['prediction']['H6_equal_added_MSE_gain_pct'],atol=1e-10)
    result={'status':'FIT_COMPLETED','specification_commit':subprocess.check_output(['git','-C',str(ROOT),'rev-parse',FROZEN]).decode().strip(),
            'cohort':{'rows':len(d),'H7':int(d.h.nunique()),'H6':int(d.h6.nunique()),'country_groups':int(d.country_holdout.nunique())},
            'outcome':outcome,'original_numeric_controls':controls,'added_population_term':POP,
            'models':results,'fixed_rule_result':label,'identifiability':audit,
            'source_receipt_sha256':digest(OUT/'source_receipt.json'),
            'population_table_sha256':digest(OUT/'population_by_original_H7.csv'),
            'files':{str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),specification,Path(__file__).with_name('build_population.py'),previous_path,source]},
            'limits':['Same observed data, not independent replication.','2020 static upstream population; local human activity may remain.',
                      'Distinct population product, not fully independent census inputs.',
                      'Predictive CI resamples fixed H6 errors, not all model uncertainty.',
                      'No further ecological outcome/model search after this result.']}
    (OUT/'result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['cohort','models','fixed_rule_result','identifiability']},indent=2))


if __name__=='__main__':main()
