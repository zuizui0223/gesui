#!/usr/bin/env python3
"""Require same-source prediction, spatial grouping and no causal path overclaims."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import pandas as pd

PANEL_SHA='1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca'
EXPECTED={'W_to_EQR':6.724969789376491,
'W_to_PlantS50':-1.3581253244027647,
'W_to_HNS':10.014323780152122,
'Plant_after_W_to_HNS':4.306041146360169,
'EQR_after_W_to_HNS':5.078783567017654,
'Plant_and_EQR_after_W_to_HNS':7.941011050123354,
'W_after_Plant_and_EQR_to_HNS':5.696582815017332,
'EQRaquatic_after_W_to_PlantS50':0.4822559143669536}

def check(panel:Path, outdir:Path):
    if hashlib.sha256(panel.read_bytes()).hexdigest()!=PANEL_SHA:raise ValueError('Source original SHA changed')
    j=json.loads((outdir/'H7_MULTIRESPONSE_NONCAUSAL_RESULTS.json').read_text())
    assert j['causal_SEM_or_mediation_fitted'] is False
    assert (j['n_rows'],j['n_H7'],j['n_H6'])==(454,232,90)
    assert set(j['models'])==set(EXPECTED)
    for k,reference in EXPECTED.items():
        q=j['models'][k];x=q['predictive_added_value']
        assert abs(x['H6_equal_MSE_gain_percent']-reference)<1e-7,k
        assert x['n_H6']==90 and x['n_rows']==454
        v=q['source_observational_partial_coefficients']
        assert v['design_rank']==v['n_design_columns'],('rank deficient',k)
    assert (j['paired_H7']['total'],j['paired_H7']['same_majority_method'])==(222,172)
    assert j['paired_H7']['W_static_in_all_222']
    pairs=pd.read_csv(outdir/'PAIRED_222_H7_OBSERVED_CHANGE.csv')
    assert len(pairs)==222 and pairs.h.nunique()==222
    assert pairs.same_method.sum()==172 and pairs.W_constant.all()
    assert len(pd.read_csv(outdir/'PREDICTIVE_MODEL_COMPARISON.csv'))==len(EXPECTED)
    independent=json.loads((outdir/'INDEPENDENT_SKLEARN_QA.json').read_text())
    assert independent['passed'] is True and set(independent['models'])=={
      'W_to_EQR','W_to_PlantS50','W_to_HNS','Plant_and_EQR_after_W_to_HNS'}
    assert max(q['difference'] for q in independent['models'].values())<1e-7
    sensitivity=json.loads((outdir/'POSTHOC_FR_SOURCE80_SENSITIVITIES.json').read_text())
    assert set(sensitivity['results'])=={'FR_only','outside_FR','taxon_coverage80'}
    assert len(pd.read_csv(outdir/'POSTHOC_FR_SOURCE80_SENSITIVITIES.csv'))==12
    print(json.dumps({'PASS':True,'models':len(EXPECTED),'posthoc_cohort_models':12,
      'paired':len(pairs),'identical_H6_groups':90,
      'independent_regression_QA_models':len(independent['models']),
      'causal_SEM':False},sort_keys=True))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--panel',type=Path,required=True);a.add_argument('--out',type=Path,required=True)
    args=a.parse_args();check(args.panel,args.out)
