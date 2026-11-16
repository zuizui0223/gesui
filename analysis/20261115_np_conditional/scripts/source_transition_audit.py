#!/usr/bin/env python3
"""Source-report-only original H7 two-period stability audit (no biology opened)."""
import argparse
import json
from pathlib import Path
import pandas as pd


def audit_report_transitions(df):
    required={'h','h6','period','country','upstream_report_positive','downstream_report_positive'}
    if not required.issubset(df.columns):raise ValueError('Missing original source fields')
    d=df[list(required)].copy()
    if len(d)!=454 or d.duplicated(['h','period']).any() or d.h.nunique()!=232 or d.h6.nunique()!=90:
        raise ValueError('Original study grain drift')
    if not (d.groupby('h').h6.nunique()==1).all():raise ValueError('Original hydrological parent drift')
    if set(d.period)!={'2017-2019','2020-2022'}:raise ValueError('Study chronology drift')
    for col in ('upstream_report_positive','downstream_report_positive'):
        if not d[col].isin([0,1]).all():raise ValueError('Source reporting flag must be binary')
    d['state']=d.upstream_report_positive.astype(int).astype(str)+d.downstream_report_positive.astype(int).astype(str)
    wide=d.pivot(index='h',columns='period',values='state').dropna()
    if len(wide)!=222:raise ValueError('Original repeated H7 support drift')
    a=wide['2017-2019'];b=wide['2020-2022']
    matrix={f'{x}_to_{y}':int(((a==x)&(b==y)).sum())
             for x in ('00','01','10','11') for y in ('00','01','10','11')}
    changes=sum(v for k,v in matrix.items() if k[:2]!=k[-2:])
    up_change=int(sum(x[0]!=y[0] for x,y in zip(a,b)))
    down_change=int(sum(x[1]!=y[1] for x,y in zip(a,b)))
    both_change=int(sum(x[0]!=y[0] and x[1]!=y[1] for x,y in zip(a,b)))
    return {'schema':'GESUI_20261115_UPSTREAM_DOWNSTREAM_SOURCE_REPORTING_PERIOD_TRANSITIONS_V1',
            'original_H7_period':len(d),'original_H7':d.h.nunique(),'original_H6':d.h6.nunique(),
            'n_repeated_H7':len(wide),'n_not_repeated_H7':d.h.nunique()-len(wide),
            'same_report_state_across_periods':int((a==b).sum()),
            'report_status_changed':int(changes),
            'upstream_status_changed':up_change,'downstream_status_changed':down_change,
            'both_directions_changed':both_change,'source_report_state_transition_matrix':matrix,
            'not_operation_start_or_stop':True,'not_treatment_assignment':True,
            'data_used':'only original report status, geography, country, and period; no biological outcomes'}

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--original-direction-panel',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    v=p.parse_args()
    d=pd.read_csv(v.original_direction_panel,usecols=['h','h6','period','country','upstream_report_positive','downstream_report_positive'],dtype={'h':str,'h6':str})
    data=audit_report_transitions(d)
    if (data['report_status_changed'],data['upstream_status_changed'],data['downstream_status_changed'],data['both_directions_changed'])!=(28,19,19,10):
        raise ValueError('Frozen original reporting transitions drifted')
    v.out.parent.mkdir(parents=True,exist_ok=True)
    v.out.write_text(json.dumps(data,indent=2,sort_keys=True))
    print(json.dumps({k:data[k] for k in ('n_repeated_H7','same_report_state_across_periods','report_status_changed','upstream_status_changed','downstream_status_changed','both_directions_changed')}))
