#!/usr/bin/env python3
"""Check rerun source-locked directional results against proposed frozen numeric ledger."""
import argparse,json,math,pathlib

EXPECTED_ROWS=454
KEYS=('HNS','EQR','plant_S50')
def validate(saved:dict,rerun:dict):
    assert saved['n_original_rows']==rerun['n_original_rows']==EXPECTED_ROWS
    assert saved['original_H7']==rerun['original_H7']==232
    assert saved['original_H6']==rerun['original_H6']==90
    assert saved['source_counts']['original_upstream_positive_site_year']==rerun['source_counts']['original_upstream_positive_site_year']==132
    assert saved['source_counts']['original_downstream_positive_site_year']==rerun['source_counts']['original_downstream_positive_site_year']==249
    n=0;worst=0.
    for key in KEYS:
        a=saved['outcomes'][key]['scores'];b=rerun['outcomes'][key]['scores']
        assert set(a)==set(b)
        for name in a:
            assert math.isfinite(a[name]['gain']) and math.isfinite(b[name]['gain'])
            diff=abs(a[name]['gain']-b[name]['gain']);worst=max(diff,worst);n+=1
            if diff>1e-8:raise RuntimeError(f'Model {key} {name} drift {diff}')
    return {'pass':True,'original_H7_period':454,'true_H6':90,'n_h6_ridge_fit_gain_specifications':n,'largest_gain_difference':worst}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--saved',type=pathlib.Path,required=True);p.add_argument('--rerun',type=pathlib.Path,required=True)
    q=p.parse_args();print(json.dumps(validate(json.loads(q.saved.read_text()),json.loads(q.rerun.read_text())),sort_keys=True))
