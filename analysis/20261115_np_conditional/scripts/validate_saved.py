#!/usr/bin/env python3
"""Assert GitHub rerun recreates the source-frozen conditional MSE comparisons."""
import argparse,json,math
from pathlib import Path


def validate(saved,recomputed,tolerance=1e-7):
 assert saved['schema']==recomputed['schema']=='GESUI_H7_454_DIRECTIONAL_NEGATIVE_CONTROL_CONDITIONAL_ADDED_VALUE_V1'
 for key,expected in [('original_study_rows',454),('H7',232),('H6',90),('upstream_report_positive_site_year',132)]:
  assert saved[key]==recomputed[key]==expected,(key,saved[key],recomputed[key])
 assert saved['cells']==recomputed['cells']=={'neither':372,'up_only':21,'down_only':26,'both':35}
 assert saved['country_cells']==recomputed['country_cells']
 n=0;maxgap=0
 for target in ('HNS','EQR','plant_S50'):
  for predictor in ('presence','reported_N_log','reported_P_log'):
   for edge in ('up_added_after_down','down_added_after_up','both_added_after_W','up_added_after_W','down_added_after_W'):
    a=saved['outcomes'][target][predictor][edge]
    b=recomputed['outcomes'][target][predictor][edge]
    for field in ('h6_mse_gain_pct','MSE_base','MSE_plus'):
     gap=abs(a[field]-b[field]);maxgap=max(maxgap,gap)
     if gap>tolerance:raise AssertionError((target,predictor,edge,field,gap))
    n+=1
 return {'pass':True,'original_cohort':454,'true_parent_H6':90,'nested_gain_specs_verified':n,'max_numeric_gap':maxgap,
         'no_causal_mediation_claim':True}
if __name__=='__main__':
 p=argparse.ArgumentParser()
 p.add_argument('--saved',type=Path,required=True)
 p.add_argument('--rerun',type=Path,required=True)
 a=p.parse_args()
 print(json.dumps(validate(json.loads(a.saved.read_text()),json.loads(a.rerun.read_text())),sort_keys=True))
