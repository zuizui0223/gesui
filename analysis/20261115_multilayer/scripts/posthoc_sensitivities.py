#!/usr/bin/env python3
"""Post-result cohort sensitivity; no causal claims or sign-filtered estimates."""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd
from multilayer import (PANEL_SOURCE_SHA256,COMMON_NUMERIC,HNS_ONLY,OUTCOMES,
                        add_country_method,check_grain,cv_add_value,cluster_paths)

MODELS={
 'W_to_EQR':(OUTCOMES['EQR'],COMMON_NUMERIC,['W']),
 'W_to_HNS':(OUTCOMES['HNS'],COMMON_NUMERIC+HNS_ONLY,['W']),
 'W_to_PlantS50':(OUTCOMES['plant'],COMMON_NUMERIC,['W']),
 'Plant_and_EQR_after_W_to_HNS':(OUTCOMES['HNS'],COMMON_NUMERIC+HNS_ONLY+['W'],['plant_S50',OUTCOMES['EQR']])}


def main(panel:Path,out:Path):
    if hashlib.sha256(panel.read_bytes()).hexdigest()!=PANEL_SOURCE_SHA256:
        raise ValueError('Unfrozen source input')
    df=add_country_method(pd.read_csv(panel,dtype={'h':str,'h6':str}));check_grain(df)
    cohorts={'FR_only':df[df.country=='FR'].copy(),
             'outside_FR':df[df.country!='FR'].copy(),
             'taxon_coverage80':df[df.minimum_species_coverage>=.8].copy()}
    assert [(k,len(z),z.h.nunique(),z.h6.nunique()) for k,z in cohorts.items()]==[
        ('FR_only',299,152,36),('outside_FR',155,82,56),('taxon_coverage80',390,217,83)]
    output={'stage':'posthoc after main model','causal_mediation':False,
            'modelled_static_outlet_W_only':True,'results':{},'source_sha256':PANEL_SOURCE_SHA256}
    rows=[]
    for key,d in cohorts.items():
        output['results'][key]={'n':len(d),'H7':d.h.nunique(),'H6':d.h6.nunique(),'models':{}}
        for name,(y,b,a) in MODELS.items():
            cv=cv_add_value(d,y,b,a);cl=cluster_paths(d,y,b,a)
            output['results'][key]['models'][name]={'predictive':cv,'within_sample_associations':cl}
            first=cl['associations'][a[0]]
            rows.append({'cohort':key,'model':name,'n':len(d),'H7':d.h.nunique(),'H6':d.h6.nunique(),
                'H6_MSE_gain_pct':cv['H6_equal_MSE_gain_percent'],
                'H6_MSE_CI_low':cv['H6_heldout_fixed_errors_bootstrap95_percent'][0],
                'H6_MSE_CI_high':cv['H6_heldout_fixed_errors_bootstrap95_percent'][1],
                'first_added':a[0], 'first_added_beta':first['beta_per_plus_one_source_SD'],
                'beta_CI_low':first['H6_cluster_descriptive95'][0],
                'beta_CI_high':first['H6_cluster_descriptive95'][1],
                'full_rank':cl['design_rank']==cl['n_design_columns']})
    out.mkdir(parents=True,exist_ok=True)
    (out/'POSTHOC_FR_SOURCE80_SENSITIVITIES.json').write_text(json.dumps(output,indent=2,sort_keys=True))
    pd.DataFrame(rows).to_csv(out/'POSTHOC_FR_SOURCE80_SENSITIVITIES.csv',index=False)
    print(json.dumps({'n_rows':len(rows),'all_full_rank':all(q['full_rank'] for q in rows),
        'cohorts':{k:(len(v),v.h.nunique(),v.h6.nunique()) for k,v in cohorts.items()}},sort_keys=True))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--panel',required=True,type=Path);a.add_argument('--out',required=True,type=Path)
    q=a.parse_args();main(q.panel,q.out)
