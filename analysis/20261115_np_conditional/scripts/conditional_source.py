#!/usr/bin/env python3
"""Conditional upstream-after-downstream source placebo in original 454 H7 periods.

This is predictive source-report geography, not a delivered WWTP nutrient dose,
observed water chemistry, a truly sewage-free counterfactual or causal SEM.
"""
from __future__ import annotations
import argparse
import collections
import importlib.util
import os
import json
import pathlib
import sys
import numpy as np
import pandas as pd


def source_cells(up, down):
    u=np.asarray(up,dtype=int); d=np.asarray(down,dtype=int)
    if len(u)!=len(d) or not np.isin(u,[0,1]).all() or not np.isin(d,[0,1]).all():
        raise ValueError('Source-report presence must be aligned binary arrays')
    return {
        'neither':int(np.sum((u==0)&(d==0))),
        'up_only':int(np.sum((u==1)&(d==0))),
        'down_only':int(np.sum((u==0)&(d==1))),
        'both':int(np.sum((u==1)&(d==1)))
    }


def validate_equal_sample(keys,expected_n=454):
    keys=list(keys)
    if len(keys)!=expected_n or len(set(keys))!=expected_n:
        raise ValueError(f'Expected {expected_n} unique original H7×period IDs')
    return True


def paired_h6_gain(base_squared_error, added_squared_error, groups, seed=20261115, rep=2000):
    a=np.asarray(base_squared_error,dtype=float)
    b=np.asarray(added_squared_error,dtype=float)
    g=np.asarray(groups,dtype=str)
    if len(a)!=len(b) or len(a)!=len(g) or len(a)==0 or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Identical nonempty finite OOF error vectors required')
    if (a<0).any() or (b<0).any():raise ValueError('Error components must be squared errors')
    ids=sorted(set(g))
    A=np.array([a[g==s].mean() for s in ids])
    B=np.array([b[g==s].mean() for s in ids])
    if A.mean()<=0:raise ValueError('Zero denominator MSE')
    gain=100*(A.mean()-B.mean())/A.mean()
    rng=np.random.default_rng(seed)
    sample=rng.integers(0,len(ids),size=(rep,len(ids)))
    boot=100*(A[sample].mean(axis=1)-B[sample].mean(axis=1))/A[sample].mean(axis=1)
    return {'h6_mse_gain_pct':float(gain),'bootstrap95':list(map(float,np.quantile(boot,[.025,.975]))),
        'n_H6':len(ids),'MSE_base':float(A.mean()),'MSE_plus':float(B.mean())}


def load_joined(panel,upstream,placebo):
    # Reuse original source SHA-checked join, not regenerated historical GBIF data.
    src=pathlib.Path(__file__).resolve().parent
    prior=pathlib.Path(os.environ.get('GESUI_PR92_SCRIPTS', str(src.parents[1]/'gesui_downstream_negativecontrol_20261115'/'scripts')))
    sys.path.insert(0,str(prior))
    import analyze_placebo
    return analyze_placebo.join_direction(panel,upstream,placebo)


def run(panel,source,placebo,model,out,seed=20261115):
    sys.path.insert(0,str(pathlib.Path(__file__).resolve().parent))
    scripts=pathlib.Path(os.environ.get('GESUI_PR92_SCRIPTS',str(pathlib.Path(__file__).resolve().parents[2]/'gesui_downstream_negativecontrol_20261115'/'scripts')))
    if not scripts.is_dir():raise RuntimeError('Source-frozen PR92 code unavailable')
    sys.path.insert(0,str(scripts))
    from analyze_placebo import join_direction
    import analyze_placebo
    indspec=importlib.util.spec_from_file_location('pr92_sk_influence', pathlib.Path(analyze_placebo.__file__).with_name('independent_sklearn_qa.py'))
    ind=importlib.util.module_from_spec(indspec)
    indspec.loader.exec_module(ind)
    multilayer=analyze_placebo.load_core(model)
    df,_,orig=join_direction(panel,source,placebo)
    df=multilayer.add_country_method(df)
    df['country_method']=df.country.fillna('__UNKNOWN__').astype(str)+'::'+df.method.fillna('__UNKNOWN__').astype(str)
    validate_equal_sample(zip(df.h,df.period))
    assert (len(df),df.h.nunique(),df.h6.nunique())==(454,232,90)
    cells=source_cells(df.uPresence,df.dPresence)
    assert cells=={'neither':372,'up_only':21,'down_only':26,'both':35}
    reported=collections.defaultdict(lambda:collections.Counter())
    for _,r in df.iterrows():
        flag=('neither' if r.uPresence==0 and r.dPresence==0 else 'up_only' if r.uPresence==1 and r.dPresence==0 else 'down_only' if r.uPresence==0 else 'both')
        reported[str(r.country)][flag]+=1
    outcomes={
      'HNS':('primary_record_mean',multilayer.COMMON_NUMERIC+multilayer.HNS_ONLY),
      'EQR':('EQR_original_site_equal_basin_period_mean',multilayer.COMMON_NUMERIC),
      'plant_S50':('plant_S50',multilayer.COMMON_NUMERIC)}
    out=pathlib.Path(out);out.mkdir(parents=True,exist_ok=True)
    summary=[];detailed={};qa=[]
    exposure={'presence':('uPresence','dPresence'),'reported_N_log':('uN','dN'),'reported_P_log':('uP','dP')}
    groups=df.h6.astype(str).to_numpy()
    for yidx,(name,(target,cols)) in enumerate(outcomes.items()):
        y=df[target].to_numpy(float)
        base=cols+['W']
        preds={}
        for tag,vars_ in [('base',[])]+[(f'{mode}_{side}',features) for mode,(up,down) in exposure.items() for side,features in [('up',[up]),('down',[down]),('both',[up,down])]]:
            p=multilayer.ridge_oof(df,target,base,multilayer.CATEGORICAL,vars_)
            errors=(y-p)**2
            preds[tag]=errors
            # Actual second implementation from frozen PR92.
            sklearn_errors=ind.sk_oof(df,target,base+vars_)
            discrepancy=float(np.max(np.abs(errors-sklearn_errors)))
            qa.append({'outcome':name,'model':tag,'max_abs_squared_error_difference':discrepancy,
                       'independent_agreement':bool(discrepancy<1e-8)})
            if discrepancy>1e-8:raise ValueError(f'Independent numerical mismatch {name}/{tag}: {discrepancy}')
        detailed[name]={}
        for idx,(mode,(up,down)) in enumerate(exposure.items()):
            eu=preds[f'{mode}_up'];ed=preds[f'{mode}_down'];eb=preds[f'{mode}_both'];e0=preds['base']
            # These effects are truly nested: both vs down, and both vs up.
            measures={
                'up_added_after_down':paired_h6_gain(ed,eb,groups,seed+1000*yidx+10*idx),
                'down_added_after_up':paired_h6_gain(eu,eb,groups,seed+1000*yidx+10*idx+1),
                'both_added_after_W':paired_h6_gain(e0,eb,groups,seed+1000*yidx+10*idx+2),
                'up_added_after_W':paired_h6_gain(e0,eu,groups,seed+1000*yidx+10*idx+3),
                'down_added_after_W':paired_h6_gain(e0,ed,groups,seed+1000*yidx+10*idx+4)
            }
            detailed[name][mode]=measures
            for lab,res in measures.items():
                summary.append({'outcome':name,'contrast':lab,'predictor_type':mode,'n_H7_period':454,
                    'n_H7':232,'n_H6':90,'H6_MSE_gain_pct':res['h6_mse_gain_pct'],
                    'H6_bootstrap95_low':res['bootstrap95'][0], 'H6_bootstrap95_high':res['bootstrap95'][1],
                    'H6_baseline_MSE':res['MSE_base'],'H6_added_MSE':res['MSE_plus']})
    result={'schema':'GESUI_H7_454_DIRECTIONAL_NEGATIVE_CONTROL_CONDITIONAL_ADDED_VALUE_V1',
            'upstream_report_positive_site_year':orig['original_upstream_positive_site_year'],
            'original_study_rows':len(df),'H7':df.h.nunique(),'H6':df.h6.nunique(),
            'cells':cells,'country_cells':{k:dict(v) for k,v in sorted(reported.items())},
            'outcomes':detailed,'same_group_cv':True,'not_a_wastewater_causal_effect':True,
            'pre2026_water_chemistry_not_imputed':True,'after_result_posthoc':True,
            'baseline_ridge_source_frozen':True,
            'estimator':'Train-only Ridge alpha=10, five-fold heldout H6, H6-equal error, paired fixed-error cluster resampling. Nested gains relative to respective restricted models.',
            'report_missing_is_not_verified_zero_sewage':True,'biological_covariate_source_Sha_frozen':True}
    (out/'CONDITIONAL_NEGATIVE_CONTROL_RESULTS.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
    pd.DataFrame(summary).to_csv(out/'CONDITIONAL_SOURCE_H6_COMPARISONS.csv',index=False)
    (out/'INDEPENDENT_SKLEARN_QA.json').write_text(json.dumps({'n':len(qa),'pass':all(z['independent_agreement'] for z in qa),
        'max_error_diff':max(z['max_abs_squared_error_difference'] for z in qa),'models':qa},indent=2))
    pd.DataFrame([{'country':c,**dict(v)} for c,v in reported.items()]).fillna(0).to_csv(out/'SOURCE_DIRECTION_COUNTRY_COUNTS.csv',index=False)
    print(json.dumps({'original_rows':len(df),'H7':df.h.nunique(),'H6':df.h6.nunique(),
        'source_cells':cells,'independent_specs':len(qa),
        'conditional_ups_after_down':{n:{m:round(v['up_added_after_down']['h6_mse_gain_pct'],3)
                                        for m,v in o.items()} for n,o in detailed.items()}},ensure_ascii=False))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for x in ('panel','source','placebo','model','out'):p.add_argument('--'+x,type=pathlib.Path,required=True)
    a=p.parse_args();run(a.panel,a.source,a.placebo,a.model,a.out)
