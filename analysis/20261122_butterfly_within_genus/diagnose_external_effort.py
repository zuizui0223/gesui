#!/usr/bin/env python3
"""Posthoc GBIF recording-effort sensitivity, same source and strata as PR #138.

Add original Animalia-minus-butterfly family and focal-pair-excluded butterfly
recording proxies. These are *not* measured observer effort or independent
non-GBIF survey data. The observed butterfly species set and H6/period/paired
species strata are unchanged. Saved model variants are never cherry-picked.

Run:
 python analysis/20261122_butterfly_within_genus/diagnose_external_effort.py \
   --source ../gesui_source_cache --output results/20261122_butterfly_within_genus/effort

The original six --source input files used by analyse.py are required, plus
HYBAS7_1219x2_GBIF_PLANTAE_LYCAENIDAE_ANIMALIA_HYDROWASTE_SOURCE_ONLY_MASTER.csv
from original Actions #43202839926 artifact 13935433410.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logsumexp

from analyse import restore_sources, build_pair_strata, source_summary, CONTROLS, FAMS
from inference import prepare_design, report_fit, model_arrays, fit_conditional
from conditional import conditional_loglike

MASTER_NAME='HYBAS7_1219x2_GBIF_PLANTAE_LYCAENIDAE_ANIMALIA_HYDROWASTE_SOURCE_ONLY_MASTER.csv'
MASTER_SHA='3f95d4825455cfab67a41de3ec4dab8778acf97d2804f2148108b65c95225c52'
BASE=['W','crop','urban','naturalout','area','nPlant','plant_coverage','plant_scored_species']
GEO=BASE[:5]
CONTROL_VARIANTS={
 'geo5':GEO,
 'geo5_nonbf_same':GEO+['nonbf_same'],
 'plant8':BASE,
 'legacy_base9_nAnimal':BASE+['nAnimal'],
 'plant8_nonbf_average':BASE+['nonbf_average'],
 'plant8_nonbf_same':BASE+['nonbf_same'],
 'plant8_nonbf_opposite':BASE+['nonbf_opposite'],
 'plant8_butterfly_expair':BASE+['butterfly_expair'],
 'plant8_both_external_proxies':BASE+['nonbf_same','butterfly_expair'],
 'legacy_main13':CONTROLS,
}
FIXED_GATE=(353,1155,34,29,122,96)


def write_json(p,value):
 p.parent.mkdir(parents=True,exist_ok=True)
 p.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def construct_proxies(source:Path,cells:pd.DataFrame,panel:pd.DataFrame,counts:pd.DataFrame):
 master=source/MASTER_NAME
 if hashlib.sha256(master.read_bytes()).hexdigest()!=MASTER_SHA:
  raise RuntimeError('Original Animalia source master SHA256 mismatch')
 index=pd.read_csv(master,dtype={'HYBAS_ID_L07':str})
 if len(index)!=2438 or index.duplicated(['HYBAS_ID_L07','PERIOD']).any():
  raise RuntimeError('Original 1219×2 master support mismatch')
 index=index[index.HYBAS_ID_L07.isin(panel.h)]
 if len(index)!=366 or index.GBIF_INDEX_ANIMALIA_2017_2022_FIXED_PERIOD_RECORD_COUNT.isna().any():
  raise RuntimeError('Missing Animalia source for fixed 183 H7 two periods')
 # Butterfly group counts are the original raw query index totals; no taxonomy filling.
 groups={}
 for name,n in [('GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json',2808),
                 ('MERGED_ORIGINAL_30H7_TWO_PERIOD_SIX_FAMILY_360_SOURCE_FACETS.json',360)]:
  raw=json.loads((source/name).read_text())
  if len(raw)!=n:raise RuntimeError('Original facet count changed '+name)
  for r in raw:
   key=(str(r['HYBAS_ID_L07']),str(r['period']),str(r['group']))
   if key in groups or r['source_status']!='ok':raise RuntimeError('Duplicated/bad facet '+str(key))
   groups[key]=int(r['source_query_occurrence_count'])
 rows=[]
 for r in index.itertuples(index=False):
  h,t=str(r.HYBAS_ID_L07),r.PERIOD
  total=int(r.GBIF_INDEX_ANIMALIA_2017_2022_FIXED_PERIOD_RECORD_COUNT)
  butterfly=sum(groups[(h,t,f)] for f in FAMS)
  if total<butterfly:raise RuntimeError('Animalia count smaller than butterfly group source')
  other=total-butterfly
  if other<=0:raise RuntimeError('Zero non-butterfly source index count')
  rows.append({'h':h,'period':t,'nonbf_same':float(np.log1p(other)),
               'nAnimal_replayed':float(np.log1p(total)),
               'butterfly_five_count':butterfly,'animalia_count':total,
               'other_animalia_count':other,
               'butterfly_expair':float(np.log1p(butterfly))})
 q=pd.DataFrame(rows)
 q['nonbf_average']=q.groupby('h').nonbf_same.transform('mean')
 q['nonbf_opposite']=q.groupby('h').nonbf_same.transform(lambda a:a.iloc[::-1].to_numpy())
 error=float(np.max(abs(q.groupby('h').nAnimal_replayed.mean()-panel.set_index('h').nAnimal)))
 if error>1e-10:raise RuntimeError('Original Animalia two-period nAnimal is not reproduced')
 new=['nonbf_same','nonbf_average','nonbf_opposite','butterfly_expair']
 enhanced=cells.merge(q[['h','period']+new],on=['h','period'],how='left',validate='one_to_one')
 if len(enhanced)!=366 or enhanced[new].isna().any().any():raise RuntimeError('Source proxy merge drift')
 p=panel.merge(q.groupby('h')[new].mean().reset_index(),on='h',how='left',validate='one_to_one')
 if len(p)!=183 or p[new].isna().any().any():raise RuntimeError('Original H7 scale drift')
 # The selected butterfly source count excludes *both* focal species in each comparison.
 specific={(str(r.h),str(r.period),str(r.butterflyKey)):int(r.count) for r in counts.itertuples(index=False)}
 totals={(r.h,r.period):r.butterfly_five_count for r in q.itertuples(index=False)}
 info={'master_SHA256':MASTER_SHA,'master_rows_all':2438,'selected_cells':366,
       'original_animalia_replay_max_abs':error,
       'nonbutterfly_count_min':int(q.other_animalia_count.min()),
       'nonbutterfly_count_median':float(q.other_animalia_count.median()),
       'nonbutterfly_count_max':int(q.other_animalia_count.max()),
       'animal_butterfly_fraction_median':float((q.butterfly_five_count/q.animalia_count).median()),
       'nonbutterfly_zero':int((q.other_animalia_count==0).sum()),
       'nonbutterfly_negative':int((q.other_animalia_count<0).sum()),
       'observational_effort_proxy_not_observer_time':True,
       'same_GBIF_platform_for_nonbutterfly_Animalia':True}
 return enhanced,p,specific,totals,info


def create_strata(cells,counts,scores,specific,totals):
 strata=build_pair_strata(cells,counts,scores,covariates=CONTROLS+['nonbf_same','nonbf_average','nonbf_opposite','butterfly_expair'])
 support=source_summary(strata)
 check=(support['informative_strata'],support['discordant_pair_cell_rows'],
        support['informative_H6'],support['informative_genera'],
        support['unique_species_pairs'],support['species_in_pairs'])
 if check!=FIXED_GATE: raise RuntimeError('Original 353 stratum same-butterfly-genus sample changed '+str(check))
 for r in strata:
  a,b=r['species_a'],r['species_b'];ref=[]
  for h in r['h']:
   five=totals[h,r['period']];ca=specific.get((h,r['period'],a),0);cb=specific.get((h,r['period'],b),0)
   if five<ca+cb:raise RuntimeError('Focal butterfly counts exceed total family GBIF records')
   ref.append(float(np.log1p(five-ca-cb)))
  r['controls']['butterfly_expair']=ref
 return strata,support


def predict_gain(design,assignment):
 X,y,lengths,w=model_arrays(design)
 groups=design['clusters']
 rows=[]
 for fold in sorted(set(assignment.values())):
  test=np.array([assignment[h]==fold for h in groups]);train=~test
  # Every fold's geographic H6 is absent from both training fits.
  if not train.any() or not test.any():raise RuntimeError('Invalid geographic fold')
  f=fit_conditional(X[train],y[train],lengths[train],w[train])
  n=fit_conditional(X[train,:,1:],y[train],lengths[train],w[train])
  if not f['converged'] or not n['converged']:
   raise RuntimeError('H6-heldout nested conditional model fit failed')
  lf=conditional_loglike(np.asarray(f['params']),X[test],y[test],lengths[test],w[test])[0]
  ln=conditional_loglike(np.asarray(n['params']),X[test,:,1:],y[test],lengths[test],w[test])[0]
  rows.append({'fold':int(fold),'h6_n':len(set(groups[test])),
               'strata_n':int(test.sum()),'gain_of_plant_N_loglik':float(ln-lf)})
 return rows


def main(source:Path,output:Path,run_cv:bool=True):
 cells,counts,scores,panel,_=restore_sources(source)
 enhanced,p,specific,totals,source_info=construct_proxies(source,cells,panel,counts)
 strata,support=create_strata(enhanced,counts,scores,specific,totals)
 results={}
 for label,controls in CONTROL_VARIANTS.items():
  design=prepare_design(strata,enhanced,p,controls)
  result=report_fit(design)
  if result['status']!='FIT':raise RuntimeError('Model STOP '+label+': '+result['status'])
  results[label]={'beta':result['beta'],'odds_ratio':result['odds_ratio'],
    'odds_ratio95_H6_cluster_t':result['odds_ratio95_H6_cluster_t'],
    'nominal_p':result['nominal_two_sided_H6_cluster_t_p'],'H6_clusters':result['H6_clusters'],
    'condition':result['information_condition_number'],'controls':controls}
 if abs(results['legacy_main13']['beta']-0.3456186953934604)>1e-6 or abs(results['legacy_base9_nAnimal']['beta']-0.1334352523119997)>1e-6:
  raise RuntimeError('Published PR138 main/baseline does not reproduce')
 cv={}
 if run_cv:
  for key in ['plant8','plant8_nonbf_same','plant8_butterfly_expair']:
   design=prepare_design(strata,enhanced,p,CONTROL_VARIANTS[key]);h6=sorted(set(design['clusters']))
   loo=predict_gain(design,{h:i for i,h in enumerate(h6)})
   repeats=[]
   for seed in range(10):
    rng=np.random.default_rng(20261122+seed);shuffled=np.asarray(h6);rng.shuffle(shuffled)
    vv=predict_gain(design,{h:i%5 for i,h in enumerate(shuffled)})
    repeats.append(float(sum(r['gain_of_plant_N_loglik'] for r in vv)))
   cv[key]={'h6_leave_one_out_gain':sum(r['gain_of_plant_N_loglik'] for r in loo),
            'h6_leave_one_out_positive':sum(r['gain_of_plant_N_loglik']>0 for r in loo),
            'h6_leave_one_out_n':len(loo),'ten_shuffled_fivefold_gains':repeats}
 out={'status':'EXPLORATORY_POSTHOC_EFFORT_PROXY_DIAGNOSTIC','source':source_info,
      'source_pair_support':support,'original_PR138_replayed':True,'models':results,
      'H6_heldout_performance':cv,
      'caveats':['Nonbutterfly Animalia is from the same GBIF index, not independently measured observer-hours.',
          'Pair-excluded butterfly index remains ecologically/recording correlated.',
          'The response records identities of observed species, not abundance or true absences.',
          'All models are posthoc sensitivity checks on the same original cohort.']}
 output.mkdir(parents=True,exist_ok=True)
 write_json(output/'external_effort_diagnostic.json',out)
 pd.DataFrame([{'model':k,**{a:b for a,b in r.items() if a!='controls'}} for k,r in results.items()]).to_csv(output/'external_effort_models.csv',index=False)
 print(json.dumps({'source':source_info,'models':{k:r['odds_ratio'] for k,r in results.items()},
    'cv':cv},ensure_ascii=False,allow_nan=False),flush=True)

if __name__=='__main__':
 parser=argparse.ArgumentParser()
 parser.add_argument('--source',type=Path,required=True)
 parser.add_argument('--output',type=Path,required=True)
 parser.add_argument('--skip-cv',action='store_true')
 args=parser.parse_args();main(args.source,args.output,not args.skip_cv)

