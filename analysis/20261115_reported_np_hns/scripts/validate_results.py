#!/usr/bin/env python3
"""Cross-check stored original reported-index fit and scope guardrails."""
from pathlib import Path
import argparse,json
import pandas as pd

EXPECTED={'HNS':{'B1_plus_matched_source_presence':3.2761587356,
                 'B4_plus_presence_and_N':-0.3287},
          'EQR':{'B1_plus_matched_source_presence':1.9197,
                 'B4_plus_presence_and_N':-0.9740},
          'PlantS50':{'B1_plus_matched_source_presence':-0.3385}}

def verify(root:Path):
    src=json.loads((root/'SOURCE_QA_AND_EXPOSURE.json').read_text())
    full=json.loads((root/'RESULTS_EEA_REPORTED_NP_TO_HNS_EQR_PLANT.json').read_text())
    source=src['eligibility'];common=src['join']
    assert (source['rows'],source['positive_rows'],source['positive_H7'],source['positive_H6'])==(454,56,38,31)
    assert sum(source['n_positive_periods'].values())==56
    assert common['country_mode_mismatch_n']==11
    assert full['original_HNS_panel_SHA256']=='1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca'
    assert full['source_reconciled_index']['source_sha256']=='012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26'
    assert abs(full['original_prior_static_W_baseline_replication']['HNS']-10.014323780152122)<1e-8
    assert abs(full['original_prior_static_W_baseline_replication']['EQR']-6.724969789376491)<1e-8
    d=pd.read_csv(root/'H7_454_REPORTED_SOURCE_INDEX_HNS_EQR_PLANT.csv',dtype={'h':str,'h6':str})
    assert len(d)==454 and d.h.nunique()==232 and d.h6.nunique()==90
    assert d.duplicated(['h','period']).sum()==0 and d.reported_positive.sum()==56
    assert d.groupby('h').h6.nunique().max()==1
    assert set(d.loc[d.reported_positive==0,'reporting_class'])=={'no_matched_positive_report'}
    row=pd.read_csv(root/'MODEL_COMPARISONS.csv')
    assert len(row)==18 and not row.duplicated(['outcome','model']).any()
    for label,checks in EXPECTED.items():
        for model,value in checks.items():
            got=float(row[(row.outcome==label)&(row.model==model)].gain_pct.iloc[0])
            assert abs(got-value)<0.002,(label,model,got,value)
    qa=json.loads((root/'INDEPENDENT_SKLEARN_QA.json').read_text())
    assert qa['success'] and len([q for z in qa['results'].values() for q in z])==18
    assert qa['max_absolute_MSE_gain_deviation']<1e-7
    print(json.dumps({'status':'PASS','rows':454,'H7':232,'H6':90,'positive':56,
      'source_country_mode_mismatches_recorded':11,'comparisons':18,
      'max_sklearn_diff':qa['max_absolute_MSE_gain_deviation'],
      'noncausal_WWTP_source_index_only':True},sort_keys=True))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True)
    a=p.parse_args();verify(a.results)
