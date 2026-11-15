#!/usr/bin/env python3
"""Posthoc cross-geography and source coverage sensitivity; original H7 support fixed.

Do not cherry-pick subcohorts as positive wastewater effect. Single source-index
reporting class is not a measured source dose and no-report isn't a control.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from analyze_reported_np import load_core,grouped_mse_comparison
from source_join import join_source,load_locked_source,sha256,FROZEN_PANEL


def run(panel:Path,source_zip:Path,ref_core:Path,outdir:Path):
    if sha256(panel)!=FROZEN_PANEL:raise ValueError('Baseline source panel SHA drift')
    source,_=load_locked_source(source_zip)
    df,_=join_source(pd.read_csv(panel,dtype={'h':str,'h6':str}),source)
    m=load_core(ref_core)
    subsets={
        'all454':np.ones(len(df),dtype=bool),
        'source_taxon_coverage_at_least80':(df.minimum_species_coverage>=.80).to_numpy(),
        'France_only':df.country.eq('FR').to_numpy(),
        'outside_France':df.country.ne('FR').to_numpy(),
        '2017_2019_only':df.period.eq('2017-2019').to_numpy(),
        '2020_2022_only':df.period.eq('2020-2022').to_numpy(),
    }
    out=[]
    for i,(key,mask) in enumerate(subsets.items()):
        sub=m.add_country_method(df.loc[mask].copy())
        n_pos=int(sub.reported_positive.sum())
        if len(sub)<100 or sub.h6.nunique()<25 or n_pos<3:
            out.append({'cohort':key,'status':'unsupported','n':len(sub),'positive':n_pos})
            continue
        target='primary_record_mean'
        base=m.COMMON_NUMERIC+m.HNS_ONLY+['W']
        H6=sub.h6.astype(str).to_numpy(); y=sub[target].to_numpy(dtype=float)
        coefs={}
        recipes={'B0':[],'B1':['reported_positive'],
            'B1_plus_N':['reported_positive','log_reported_N_mgL_equivalent'],
            'B1_plus_P':['reported_positive','log_reported_P_ugL_equivalent']}
        for name,extras in recipes.items():
            p=m.ridge_oof(sub,target,base,m.CATEGORICAL,extras)
            coefs[name]=(y-p)**2
        for j,(model,ref) in enumerate([('B1','B0'),('B1_plus_N','B1'),('B1_plus_P','B1')]):
            z=grouped_mse_comparison(coefs[ref],coefs[model],H6,seed=20261115+i*10+j)
            out.append({'cohort':key,'status':'exploratory','n':len(sub),'H7':int(sub.h.nunique()),
               'H6':int(sub.h6.nunique()),'positive':n_pos,'positive_H6':int(sub.loc[sub.reported_positive>0,'h6'].nunique()),
               'model':model,'base':ref,'gain_pct':z['H6_equal_MSE_gain_percent'],
               'fixed_errors95':z['H6_fixed_heldout_errors_bootstrap95']})
    outdir.mkdir(parents=True,exist_ok=True)
    (outdir/'POSTHOC_SOURCE_COMPLETENESS_COUNTRY_ERA.json').write_text(json.dumps(out,indent=2,sort_keys=True))
    pd.DataFrame(out).to_csv(outdir/'POSTHOC_SOURCE_COMPLETENESS_COUNTRY_ERA.csv',index=False)
    print(json.dumps(out,sort_keys=True),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--panel',type=Path,required=True)
    p.add_argument('--source-zip',type=Path,required=True)
    p.add_argument('--reference-core',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.panel,a.source_zip,a.reference_core,a.out)
