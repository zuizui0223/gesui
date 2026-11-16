#!/usr/bin/env python3
"""Independent algebraic audit using PR92 saved out-of-H6 model MSEs.

No shared transformation or model fit required: validates that each nested
conditional gain is computed from the corresponding original same-unit MSE.
"""
import argparse,json,math
from pathlib import Path


def audit(previous,new):
 assert new['original_study_rows']==previous['n_original_rows']==454
 assert new['H7']==previous['original_H7']==232
 assert new['H6']==previous['original_H6']==90
 assert new['cells']=={'neither':372,'up_only':21,'down_only':26,'both':35}
 n=0;error=0
 for outcome in ('HNS','EQR','plant_S50'):
  prev=previous['outcomes'][outcome]['scores']
  for typ,tag in (('presence','present'),('reported_N_log','N'),('reported_P_log','P')):
   up=prev['source_up_'+tag];down=prev['source_down_'+tag];both=prev['source_both_'+tag]
   recipes={
    'up_added_after_down':(down['MSE_added'],both['MSE_added']),
    'down_added_after_up':(up['MSE_added'],both['MSE_added']),
    'both_added_after_W':(up['MSE_base'],both['MSE_added']),
    'up_added_after_W':(up['MSE_base'],up['MSE_added']),
    'down_added_after_W':(up['MSE_base'],down['MSE_added'])}
   for contrast,(restricted,expanded) in recipes.items():
    expected=100*(restricted-expanded)/restricted
    actual=new['outcomes'][outcome][typ][contrast]['h6_mse_gain_pct']
    gap=abs(expected-actual);error=max(error,gap);n+=1
    if gap>1e-9:raise AssertionError((outcome,typ,contrast,gap,expected,actual))
 return {'pass':True,'independently_verified_nested_MSE_comparisons':n,
         'largest_percent_point_discrepancy':error,'original_identical_H6_groups':90,
         'QA_is_same_data_algebra_not_independent_observation':True}

if __name__=='__main__':
 p=argparse.ArgumentParser()
 p.add_argument('--previous',required=True,type=Path)
 p.add_argument('--current',required=True,type=Path)
 p.add_argument('--out',required=True,type=Path)
 args=p.parse_args()
 result=audit(json.loads(args.previous.read_text()),json.loads(args.current.read_text()))
 args.out.write_text(json.dumps(result,indent=2))
 print(json.dumps(result))
