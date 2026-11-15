#!/usr/bin/env python3
"""Compare a fresh Actions run to signed original local research result scalars."""
from pathlib import Path
import argparse,json,math
import pandas as pd


def run(reference:Path,new:Path):
    ref=json.loads((reference/'RESULTS_EEA_REPORTED_NP_TO_HNS_EQR_PLANT.json').read_text())
    check=json.loads((new/'RESULTS_EEA_REPORTED_NP_TO_HNS_EQR_PLANT.json').read_text())
    assert check['original_HNS_panel_SHA256']==ref['original_HNS_panel_SHA256']
    assert check['source_reconciled_index']['source_sha256']==ref['source_reconciled_index']['source_sha256']
    assert check['source_qa']==ref['source_qa']
    assert check['paired_report_status']==ref['paired_report_status']
    deviations=[]
    for label in ref['models']:
        for model in ref['models'][label]['comparisons']:
            a=ref['models'][label]['comparisons'][model]['H6_equal_MSE_gain_percent']
            b=check['models'][label]['comparisons'][model]['H6_equal_MSE_gain_percent']
            deviations.append(abs(a-b))
    assert len(deviations)==18
    source_csv=pd.read_csv(reference/'MODEL_COMPARISONS.csv')
    actual_csv=pd.read_csv(new/'MODEL_COMPARISONS.csv')
    assert set(zip(source_csv.outcome,source_csv.model))==set(zip(actual_csv.outcome,actual_csv.model))
    assert max(deviations)<1e-7,max(deviations)
    independent=json.loads((new/'INDEPENDENT_SKLEARN_QA.json').read_text())
    assert independent['success'] and independent['max_absolute_MSE_gain_deviation']<1e-7
    grouped=pd.read_csv(new/'H7_454_REPORTED_SOURCE_INDEX_HNS_EQR_PLANT.csv',dtype={'h':str,'h6':str})
    assert grouped.duplicated(['h','period']).sum()==0
    print(json.dumps({'status':'PASS','comparisons_rechecked':18,
            'maximum_difference_percentage_points':max(deviations),
            'independent_sklearn_max_delta':independent['max_absolute_MSE_gain_deviation'],
            'all_454_original_rows_preserved':len(grouped)==454,'source_sha_locked':True},sort_keys=True))

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--reference',type=Path,required=True)
    ap.add_argument('--new',type=Path,required=True)
    a=ap.parse_args();run(a.reference,a.new)
