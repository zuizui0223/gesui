"""One frozen EU451 W increment comparison, using original model functions."""
from pathlib import Path
import importlib.util
import hashlib
import json
import subprocess

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'analysis/20261115_eu454_plant_nutri/h6_model_if_gate.py'
spec=importlib.util.spec_from_file_location('original_eu451_models',SOURCE)
old=importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
WORK=ROOT.parent
OUT=ROOT/'results/20261120_eu451_W_confirmation'
FROZEN='20b2cb6'


def digest(p):
    b=p.read_bytes()
    if p.suffix in {'.py','.md'}:b=b.replace(b'\r\n',b'\n')
    return hashlib.sha256(b).hexdigest()


def main():
    specification=Path(__file__).with_name('SPECIFICATION.md')
    committed=subprocess.check_output(['git','-C',str(ROOT),'show',FROZEN+':analysis/20261120_eu451_W_confirmation/SPECIFICATION.md'])
    assert committed.replace(b'\r\n',b'\n')==specification.read_bytes().replace(b'\r\n',b'\n')
    gate_path=WORK/'plant/EU454_FULL_EIV_SOURCE_GATE.json'
    gate=json.loads(gate_path.read_text())
    assert gate['schema']=='GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA'
    assert gate['source_facets_SHA256']==old.FACETS_SHA
    assert gate['source_EIV']['original_source_SHA256']==old.EIV_SHA
    assert gate['hard_source_support_pass']
    assert gate['original_EIV_NUTRI_scored_keys_unchanged'] and gate['original_plant_ELIGIBILITY_H7period_unchanged']
    assert gate['source_hybrid_score_assigned']==0 and gate['official_confirmed_valid_unscorable_nonbinomial_hybrids']==175
    panel=WORK/'original/H7_454_HNS_panel.csv'
    assert old.sha(panel)==old.ORIGINAL_PANEL_SHA
    d=pd.read_csv(panel,dtype={'h':str,'h6':str})
    ep=WORK/'plant/EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv'
    sp=WORK/'plant/EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv'
    e=pd.read_csv(ep,dtype={'h':str,'h6':str})
    s=pd.read_csv(sp,dtype={'h':str})
    assert len(e)==454 and not e.duplicated(['h','period']).any() and not s.duplicated(['h','period']).any()
    e=e[e.full_source_trait_eligible.astype(str)=='True']
    d=d.merge(e[['h','period','NUTRI_coverage_all_Plantae_speciesKeys','exact_NUTRI_scored_vascular_accepted_species']],
              on=['h','period'],validate='one_to_one')
    d=d.merge(s[['h','period','Plant_N_incidence','Plant_N_record','n_scored_species']],on=['h','period'],validate='one_to_one')
    assert (len(d),d.h.nunique(),d.h6.nunique())==(451,232,90)
    assert len(d)==gate['eligible_H7_period'] and d.h6.nunique()==gate['eligible_true_H6']
    assert d.country_holdout.nunique()==gate['eligible_majority_country_groups']
    assert d.n_scored_species.min()>=10 and d.NUTRI_coverage_all_Plantae_speciesKeys.min()>=.4
    assert d.groupby('h').h6.nunique().max()==1
    d['log1p_eiv_scored_taxa']=np.log1p(d.exact_NUTRI_scored_vascular_accepted_species.to_numpy(float))
    d['country_method']=d.country.fillna('__UNKNOWN__').astype(str)+'::'+d.method.fillna('__UNKNOWN__').astype(str)
    controls=[c for c in old.B_GEOG if c!='W']
    assert controls==['crop','urban','riverflow','naturalout','area','nEEA','nAnimal','nPlant','nButterfly',
                      'primary_record_coverage_raw','primary_host_completeness','NUTRI_coverage_all_Plantae_speciesKeys','log1p_eiv_scored_taxa']
    target='primary_incidence_mean'
    association=old.source_cluster_regression(d,target,controls+['W'],'W')
    prediction=old.cv_h6(d,target,controls,'W')
    a=association['H6_cluster_descriptive95'][0]>0
    b=prediction['fixed_out_of_fold_error_H6_resampling95_pct'][0]>0
    decision={ (True,True):'positive_association_and_predictive_information',
               (True,False):'positive_association_only', (False,True):'predictive_information_only',
               (False,False):'positive_W_signal_not_supported'}[(a,b)]
    # Independent weighted Frisch-Waugh-Lovell point check.
    X,labels=old.design_full(d,controls+['W'])
    w=old.H6_weights(d)
    j=labels.index('W'); Z=np.delete(X,j,axis=1)
    sw=np.sqrt(w)
    def residual(v):return v-Z@np.linalg.lstsq(Z*sw[:,None],v*sw,rcond=None)[0]
    x=residual(X[:,j]);y=residual(d[target].to_numpy())
    beta=np.sum(w*x*y)/np.sum(w*x*x)
    np.testing.assert_allclose(beta,association['adjusted_beta_per_1SD_focal'],atol=1e-12)
    files=[panel,ep,sp,gate_path,SOURCE,Path(__file__),specification]
    result={'specification_commit':subprocess.check_output(['git','-C',str(ROOT),'rev-parse',FROZEN]).decode().strip(),
            'cohort':{'H7_period':451,'H7':232,'H6':90,'country_holdout_groups':int(d.country_holdout.nunique())},
            'target':target,'numeric_controls':controls,'categorical_controls':old.CATS,
            'association':association,'prediction':prediction,'fixed_rule_result':decision,
            'checks':{'source_gate_pass':True,'fixed_spec_matches_commit':True,'focal_FWL_beta':beta,
                      'full_rank_design':True,'original_CV_checks_train_test_H6_disjoint':True},
            'inputs':[{'path':str(p.relative_to(WORK)),'sha256':digest(p),
                       'hash_mode':'LF-normalized' if p.suffix in {'.py','.md'} else 'original bytes'} for p in files],
            'limits':['Not independent replication or new-data preregistration.',
                      'Old models already adjusted for cropland and the other landscape variables.',
                      'Static W is not measured sewage dose; climate/elevation unadjusted.',
                      'Species-equal larval host affinity, not adult protein use or leaf N.',
                      'Neither failure nor success identifies a causal nutrient pathway.']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['cohort','association','prediction','fixed_rule_result','checks']},indent=2))


if __name__=='__main__':main()
