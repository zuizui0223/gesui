#!/usr/bin/env python3
"""Posthoc same-source EEA annual WWTP N/P directional negative-control check.

Source-reported mass/mean river flow is not measured TN/TP, and downstream
source is a geography/reporting negative control, not randomized assignment.
"""
from __future__ import annotations
import argparse, collections, hashlib, importlib.util, json, pathlib, os, sys, zipfile
import numpy as np, pandas as pd

SRC=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(SRC))
# Original PR90 code reproduces locked H7×period source and the source EQR flow.
PR90=pathlib.Path(os.environ.get('GESUI_PR90_SCRIPTS', str(SRC.parents[1]/'20261115_reported_np_hns'/'scripts')))
if not PR90.is_dir(): PR90=pathlib.Path('/mnt/data/gesui_h7_actual_np_hns_20261115/scripts')
if PR90.is_dir():sys.path.insert(0,str(PR90))
from source_join import load_locked_source,join_source,sha256,FROZEN_PANEL

UPDOWN_SOURCE_SHA='c07815e270f7f48f6b24629490b4827e31c4448ae421afc47f9ff0bf85d65e90'
SOURCE_CSV='H7_454_ORIGINAL_UPSTREAM_DOWNSTREAM_REPORTED_NP_PLACEBO.csv'
SOURCE_JSON='SOURCE_ONLY_UPDOWN_PLACEBO_QA.json'
def load_core(path):
    spec=importlib.util.spec_from_file_location('frozen_multilayer_core_negativecontrol',path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def join_direction(panel,original_upstream_zip,updown_zip):
    if sha256(panel)!=FROZEN_PANEL:raise ValueError('HNS source panel drift')
    old,_=load_locked_source(original_upstream_zip)
    p=pd.read_csv(panel,dtype={'h':str,'h6':str})
    study,qa=join_source(p,old)
    if updown_zip.is_dir():
        data=pd.read_csv(updown_zip/SOURCE_CSV,dtype={'h':str})
        raw=json.loads((updown_zip/SOURCE_JSON).read_text())
    else:
        with zipfile.ZipFile(updown_zip) as z:
            data=pd.read_csv(z.open(SOURCE_CSV),dtype={'h':str})
            raw=json.loads(z.read(SOURCE_JSON))
    if len(data)!=454 or data.duplicated(['h','period']).any() or len(set(data.h))!=232:raise ValueError('Placebo H7 source grain changed')
    if raw['original_upstream_positive_site_year']!=132 or raw['H7_period_upstream_positive']!=56:
        raise ValueError('Placebo population wrong')
    m=study.merge(data,on=['h','period'],validate='one_to_one',indicator=True)
    if len(m)!=454 or not(m['_merge']=='both').all():raise ValueError('H7 source join not exact')
    if not np.array_equal(m.reported_positive.to_numpy(dtype=int),m.upstream_report_positive.to_numpy(dtype=int)):
        raise ValueError('Original upstream source-positive differs from current panel')
    if np.max(np.abs(m.n_mg_l-m.upstream_N_equiv_mgL))>1e-9 or np.max(np.abs(m.p_mg_l-m.upstream_P_equiv_mgL))>1e-9:
        raise ValueError('Original numerical upstream N/P index not reproduced')
    for k in ('downstream_report_positive','downstream_N_equiv_mgL','downstream_P_equiv_mgL'):
        if m[k].isna().any():raise ValueError('Missing report index '+k)
    m['dN']=np.log1p(m.downstream_N_equiv_mgL)
    m['dP']=np.log1p(1000*m.downstream_P_equiv_mgL)
    m['uN']=np.log1p(m.n_mg_l)
    m['uP']=np.log1p(1000*m.p_mg_l)
    m['uPresence']=m.reported_positive.astype(int)
    m['dPresence']=m.downstream_report_positive.astype(int)
    return m.drop(columns='_merge'),qa,raw

def group_error(e0,e1,groups,seed=20261115):
    ids=np.array(groups,dtype=str)
    g=sorted(set(ids));a=np.array([np.mean(e0[ids==z]) for z in g]);b=np.array([np.mean(e1[ids==z]) for z in g])
    if a.mean()<=0:raise ValueError('Baseline MSE missing')
    rng=np.random.default_rng(seed)
    ix=rng.integers(0,len(g),size=(2000,len(g)))
    boot=100*(a[ix].mean(axis=1)-b[ix].mean(axis=1))/a[ix].mean(axis=1)
    return {'gain':float(100*(a.mean()-b.mean())/a.mean()),'95_interval':list(map(float,np.quantile(boot,[.025,.975]))),'n_H6':len(g),'MSE_base':float(a.mean()),'MSE_added':float(b.mean())}

def test_negativecontrol(df,core,response,base,category,seed=20261115):
    y=df[response].to_numpy(float);h6=df.h6.astype(str).to_numpy()
    F={
      'static_W':[],
      'source_up_present':['uPresence'],
      'source_down_present':['dPresence'],
      'source_both_present':['uPresence','dPresence'],
      'source_up_N':['uN'],
      'source_down_N':['dN'],
      'source_both_N':['uN','dN'],
      'source_up_P':['uP'],
      'source_down_P':['dP'],
      'source_both_P':['uP','dP'],
    }
    preds={m:core.ridge_oof(df,response,base,core.CATEGORICAL,fields) for m,fields in F.items()}
    errors={m:(y-preds[m])**2 for m in F}
    scores={m:group_error(errors['static_W'],errors[m],h6,seed+i) for i,m in enumerate(F) if m!='static_W'}
    contrasts={}
    for i,(u,d,both) in enumerate((('source_up_present','source_down_present','source_both_present'),('source_up_N','source_down_N','source_both_N'),('source_up_P','source_down_P','source_both_P'))):
        first=errors[u];second=errors[d];gg=sorted(set(h6))
        ua=np.array([np.mean(first[h6==g]) for g in gg]);da=np.array([np.mean(second[h6==g]) for g in gg]);base0=np.array([np.mean(errors['static_W'][h6==g]) for g in gg]);
        rng=np.random.default_rng(seed+1000+i);ix=rng.integers(0,len(gg),size=(2000,len(gg)))
        obs=float(100*(da.mean()-ua.mean())/base0.mean())
        bs=100*(da[ix].mean(1)-ua[ix].mean(1))/base0[ix].mean(1)
        contrasts[('presence','N','P')[i]]={'upstream_minus_downstream_gain_percentage_points':obs,
          'H6_descriptive_bootstrap95':list(map(float,np.quantile(bs,[.025,.975]))),
          'both_source_gain_pct':scores[both]['gain'],
          'interpretation':'posthoc directionality falsifier; not causal proof'}
    return {'n_H6':len(set(h6)),'n_original_454':len(df),'outcome':response,'scores':scores,'directional_contrast':contrasts,'posthoc':True}

def run(panel,prior,placebo,model,out):
    df,joinqa,source=join_direction(panel,prior,placebo)
    core=load_core(model)
    d=core.add_country_method(df)
    outcomes={
      'HNS':('primary_record_mean',core.COMMON_NUMERIC+core.HNS_ONLY),
      'EQR':('EQR_original_site_equal_basin_period_mean',core.COMMON_NUMERIC),
      'plant_S50':('plant_S50',core.COMMON_NUMERIC)
    }
    result={'schema':'GESUI_EEA_2018_2022_ACTUAL_UPDOWN_PLACEBO_454_NONCAUSAL_V1',
       'source_GH_action_run':43679382449,'n_original_rows':len(d),
       'original_H7':d.h.nunique(),'original_H6':d.h6.nunique(),
       'source_counts':source,'original_source_join_qa':joinqa,'posthoc':True,
       'WWTP_location_and_report_not_received_dose':True,
       'downstream_only_valid_as_monitor_site_weak_negative_control':True,
       'full_population_not_true_untreated_controls':True,'outcomes':{}}
    summary=[]
    for i,(label,(response,fields)) in enumerate(outcomes.items()):
        fit=test_negativecontrol(d,core,response,fields+['W'],label,seed=20261115+10*i)
        result['outcomes'][label]=fit
        for name,score in fit['scores'].items():
            summary.append({'outcome':label,'model':name,'n':len(d),'H7':d.h.nunique(),'H6':d.h6.nunique(),
               'H6_gain_pct':score['gain'],'bootstrap95_low':score['95_interval'][0],
               'bootstrap95_high':score['95_interval'][1]})
    out.mkdir(parents=True,exist_ok=True)
    (out/'DIRECTIONAL_NEGATIVE_CONTROL_RESULTS.json').write_text(json.dumps(result,indent=2,sort_keys=True))
    pd.DataFrame(summary).to_csv(out/'DIRECTIONAL_NEGATIVE_CONTROL_MODEL_SUMMARY.csv',index=False)
    d[['h','h6','period','country','W','primary_record_mean','EQR_original_site_equal_basin_period_mean','plant_S50',
       'upstream_report_positive','downstream_report_positive','upstream_N_equiv_mgL','downstream_N_equiv_mgL',
       'upstream_P_equiv_mgL','downstream_P_equiv_mgL']].to_csv(out/'H7_454_DIRECTIONAL_SOURCE_AUDIT_PANEL.csv',index=False)
    print(json.dumps({'n454':len(d),'source_source_positive_up':int(d.uPresence.sum()),'source_source_positive_down':int(d.dPresence.sum()),
                      'results':{label:{dim:{k:round(v['gain'],4) for k,v in m['scores'].items()} for dim,m in [('MSE',res)]}
                        for label,res in result['outcomes'].items()}},sort_keys=True))
    return result
if __name__=='__main__':
    a=argparse.ArgumentParser()
    for z in ['panel','source','placebo','model','out']:a.add_argument('--'+z,required=True,type=pathlib.Path)
    x=a.parse_args();run(x.panel,x.source,x.placebo,x.model,x.out)
