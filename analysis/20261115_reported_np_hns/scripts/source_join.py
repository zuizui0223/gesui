#!/usr/bin/env python3
"""Source-gated join of true same-H7 EEA reported WWTP N/P mass indices.

IMPORTANT: Reported N/P mass divided by modelled stream flow is a hypothetical
full-delivery concentration equivalent, NOT measured in-stream N/P or actual
butterfly-site effluent fraction. A no-match index is NOT pristine water.
"""
from __future__ import annotations
from pathlib import Path
import hashlib
import json
import zipfile
import numpy as np
import pandas as pd

FROZEN_SOURCE='012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26'
FROZEN_PANEL='1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca'
SOURCE_FILENAME='EEA2018_EEA2022_454_H7_ORIGINAL_REPORTED_NP_SITE_FLOW_EQUIVALENT_SOURCE_ONLY.json'
SOURCE_FIELDS={
    'HYBAS_ID_L07':'h','PERIOD':'period',
    'source_station_annual_upstream_NP_positive_n':'n_source',
    'source_year_reported_N_mgL_equivalent_site_equal_mean':'n_mg_l',
    'source_year_reported_P_mgL_equivalent_site_equal_mean':'p_mg_l',
    'source_country_modes':'original_EEA_source_country_mode',
    'source_year_station_DIS_AV_CMS_mean':'source_year_station_DIS_AV_CMS_mean',
    'source_year_reported_upstream_N_tpy_site_equal_mean':'n_tpy_mean',
    'source_year_reported_upstream_P_tpy_site_equal_mean':'p_tpy_mean',
}

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def load_locked_source(zipfile_path: Path):
    if zipfile_path.suffix.lower()=='.json':
        if zipfile_path.name!=SOURCE_FILENAME:raise ValueError('Unexpected source JSON filename')
        raw=zipfile_path.read_bytes()
    else:
        with zipfile.ZipFile(zipfile_path) as z:
            raw=z.read(SOURCE_FILENAME)
    digest=hashlib.sha256(raw).hexdigest()
    if digest!=FROZEN_SOURCE:raise ValueError('FROZEN source SHA mismatch')
    data=json.loads(raw)
    if not data['original_source_reconciliation_pre_ecology'] or not data['true_original_EEA_numeric_EQR_never_used_to_select']:
        raise ValueError('Original upstream source-only gate changed')
    df=pd.DataFrame(data['source_H7_period_reported_N_P_upstream_to_true_original_station'])
    if len(df)!=454 or not df.source_valid.all():raise ValueError('Source original support drift')
    df=df.rename(columns=SOURCE_FIELDS)
    df['h']=df.h.astype(str)
    if df.duplicated(['h','period']).any():raise ValueError('Duplicated original source H7×period')
    if df['h'].nunique()!=232 or set(df.period)!={'2017-2019','2020-2022'}:
        raise ValueError('Source original H7 or year support changed')
    return df,{'source_sha256':digest,'source_schema':data['schema'],'source_rows':len(df),
               'source_input_container_sha256':sha256(zipfile_path)}

def derive_exposure(df: pd.DataFrame) -> pd.DataFrame:
    res=df.copy()
    for col in ['n_source','n_mg_l','p_mg_l']:
        if not np.isfinite(res[col].to_numpy(dtype=float)).all() or (res[col].astype(float)<0).any():
            raise ValueError('Invalid source reported mass/flow value '+col)
    res['reported_positive']=(res.n_source.astype(float)>0).astype(int)
    if (res.loc[res.reported_positive==0,['n_mg_l','p_mg_l']].to_numpy(dtype=float)>0).any():
        raise ValueError('Original source report positive flag inconsistent with N/P mass')
    if (res.loc[res.reported_positive==1,['n_mg_l','p_mg_l']].sum(axis=1)<=0).any():
        raise ValueError('Original positive report but zero reported N/P mass')
    res['reporting_class']=res.reported_positive.map({0:'no_matched_positive_report',1:'matched_positive_report'})
    res['log_reported_N_mgL_equivalent']=np.log1p(res.n_mg_l.to_numpy(dtype=float))
    res['log_reported_P_ugL_equivalent']=np.log1p(1000*res.p_mg_l.to_numpy(dtype=float))
    res['reported_N_and_P_both_positive']=((res.n_mg_l>0)&(res.p_mg_l>0)).astype(int)
    return res

def join_source(panel: pd.DataFrame, source: pd.DataFrame):
    p=panel.copy()
    if not {'h','h6','period','country','source_mean_DIS_AV_CMS'}.issubset(set(p.columns)):
        raise ValueError('Missing original panel geometry/measurement fields')
    p['h']=p.h.astype(str)
    if len(p)!=454 or p.h.nunique()!=232 or p.h6.nunique()!=90 or p.duplicated(['h','period']).any():
        raise ValueError('Biological source cohort grain drift')
    cols=list(SOURCE_FIELDS.values())+['source_station_annual_upstream_NP_positive_fraction']
    source=source[cols]
    out=p.merge(source,on=['h','period'],how='left',validate='one_to_one',indicator=True)
    if len(out)!=len(p) or not(out._merge=='both').all():raise ValueError('Original H7 source join incomplete')
    if out['n_source'].isna().any():raise ValueError('Unavailable original source index')
    out=derive_exposure(out.drop(columns='_merge'))
    fd=(out.source_mean_DIS_AV_CMS-out.source_year_station_DIS_AV_CMS_mean).abs()
    if fd.max()>1e-9:raise ValueError('Source field river flow does not match locked panel')
    report={
        'original_joint_rows':len(out),'original_H7':int(out.h.nunique()),'original_H6':int(out.h6.nunique()),
        'country_mode_mismatch_n':int((out.country!=out.original_EEA_source_country_mode).sum()),
        'flow_value_max_abs_diff':float(fd.max()),
        'positive_reported_source_not_confirmed_no_source_controls':True,
        'same_field_site_EEA_EQR_but_butterfly_records_are_basinwide':True,
        'source_N_and_P_are_mass_to_modelled_flow_not_measured_water_concentration':True,
        'report_not_proof_of_actual_zero':True,
    }
    return out,report

def examine_source_eligibility(x: pd.DataFrame):
    pos=x[x.reported_positive==1]
    n_by_period={v:int(pos.period.eq(v).sum()) for v in ('2017-2019','2020-2022')}
    c=pos.groupby('country').size().sort_values(ascending=False)
    return {
        'rows':int(len(x)), 'positive_rows':int(len(pos)),
        'n_positive_periods':n_by_period,'positive_H7':int(pos.h.nunique()),
        'positive_H6':int(pos.h6.nunique()), 'positive_country_counts':c.to_dict(),
        'positive_both_NP':int(pos.reported_N_and_P_both_positive.sum()),
        'uncertain_source_no_matched_positive_report_rows':int(len(x)-len(pos)),
        'exploratory_support_gate_pass': bool(len(pos)>=30 and pos.h.nunique()>=20 and
                    min(n_by_period.values())>=10 and len(c)>=3),
        'positive_source_country_FR':int(pos.country.eq('FR').sum()),
        'N_P_pearson_source_all':float(x[['n_mg_l','p_mg_l']].corr().iloc[0,1]),
        'N_P_pearson_source_positive':float(pos[['n_mg_l','p_mg_l']].corr().iloc[0,1]),
    }
