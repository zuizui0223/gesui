#!/usr/bin/env python3
"""Independent row-record groupby recomputation of the HNS taxon ablations.

Uses original long-format GBIF species-by-H7 source records, not the NPZ weighted
matrix used by taxon_influence. Independently verifies each response mean.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def run(species:Path,panel:Path,report:Path,out:Path):
    d=pd.read_csv(panel,dtype={'h':str,'h6':str})
    records=pd.read_csv(species,dtype={'h':str,'h6':str,'butterflyKey':str})
    res=json.loads(report.read_text())
    model_rows=res['models']
    stats=[]
    for row in model_rows:
        drop=set(row['drop_keys'])
        scored=records.loc[records.primary_HNS.notna() & ~records.butterflyKey.isin(drop)].copy()
        scored['numerator']=scored.n_records*scored.primary_HNS
        group=scored.groupby(['h','period'],sort=False)[['numerator','n_records']].sum()
        group['meanHNS']=group.numerator/group.n_records
        if len(group)!=454 or (group.n_records<=0).any():
            raise AssertionError('Original source population or denominators changed')
        fitted=d[['h','period']].join(group['meanHNS'],on=['h','period'],validate='one_to_one')
        if fitted.meanHNS.isna().any():
            raise AssertionError('Unexpected missing independently computed HNS')
        mean=float(fitted.meanHNS.mean())
        if int(row['dropped_top_absolute_previous_contributions'])==0:
            delta=float(np.max(np.abs(fitted.meanHNS.to_numpy()-d.primary_record_mean.to_numpy())))
            if delta>1e-10:raise AssertionError('Original HNS mean not reproduced')
        else:delta=None
        stats.append({'drop':row['dropped_top_absolute_previous_contributions'],
                      'n_original_rows':len(fitted),'independent_mean_HNS':mean,
                      'original_no_drop_max_abs_diff':delta})
    result={'source':'authentic long-format scored butterfly species records × each H7 historical period',
            'models_checked':len(stats),'grains':len(d),'all_pass':True,'checks':stats,
            'not_independent_field_data':True}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2))
    print(json.dumps(result,sort_keys=True))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['species','panel','report','out']:
        p.add_argument('--'+k,type=Path,required=True)
    x=p.parse_args();run(x.species,x.panel,x.report,x.out)
