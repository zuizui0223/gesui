#!/usr/bin/env python3
"""Validate fixed exact source basins and noncausal interpretation for taxon ablation."""
import argparse,hashlib,json
from pathlib import Path
import pandas as pd

EXPECTED={0:0.10802204501635151,1:0.11149071287132045,5:0.09947635761512583,
  10:0.07117639040717524,20:0.044395412813850905,30:0.04761073240295205}

def validate(panel:Path,results:Path,independent:Path):
    assert hashlib.sha256(panel.read_bytes()).hexdigest()=='1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca'
    j=json.loads(results.read_text());audit=json.loads(independent.read_text())
    assert (j['n_H7_period'],j['n_H7'],j['n_H6'])==(454,232,90)
    assert j['all_source_covariates_and_family_share_frozen'] is True
    assert j['source_record_weighting_not_individual_abundance'] is True
    assert list(x['dropped_top_absolute_previous_contributions'] for x in j['models'])==[0,1,5,10,20,30]
    for m in j['models']:
        n=m['dropped_top_absolute_previous_contributions']
        assert (m['new_source_rows'],m['n_H7'],m['n_H6'])==(454,232,90)
        assert len(m['drop_keys'])==n and len(set(m['drop_keys']))==n
        assert m['original_source_denominator_not_zero'] is True
        assert abs(m['W_beta']-EXPECTED[n])<1e-7,(n,m['W_beta'])
        assert m['W_beta_H6_descriptive95_low']<=m['W_beta']<=m['W_beta_H6_descriptive95_high']
        assert m['W_gain_fixed_error95_low']<=m['W_added_MSE_gain_percent']<=m['W_gain_fixed_error95_high']
    assert audit['all_pass'] is True and audit['models_checked']==6 and audit['grains']==454
    print(json.dumps({'status':'PASS','source_rows':454,'n_H6':90,'ablation_models':6,'recordwise_independent_arithmetic':'PASS','causal_path_analysis':False}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--panel',type=Path,required=True);p.add_argument('--results',type=Path,required=True);p.add_argument('--independent',type=Path,required=True)
    x=p.parse_args();validate(x.panel,x.results,x.independent)
