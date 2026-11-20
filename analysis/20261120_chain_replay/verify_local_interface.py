"""Offline replay of local plant contrasts and frozen butterfly support.
No GBIF requests, new effects, sign tests or relaxed support criteria.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT/'docs/manuscript_riparian_interactions'
COMMIT = '2bedeae879a80b8af842703a29836ce099cd59f0'
PREFIX = 'results/independent_deu_tp_plantae_20261107/'

def main():
    p=argparse.ArgumentParser(); p.add_argument('--crosswalk',type=Path,required=True); a=p.parse_args()
    sources=[]
    def source(name):
        b=subprocess.check_output(['git','-C',str(ROOT),'show',COMMIT+':'+PREFIX+name])
        sources.append({'commit':COMMIT,'path':PREFIX+name,'sha256':hashlib.sha256(b).hexdigest()})
        return b
    plant=pd.read_csv(io.BytesIO(source('PLANT_site_species_records.csv')))
    tp=pd.read_csv(io.BytesIO(source('paired_monthly_TP_mg_L.csv'))).groupby('outfall_id').log_down_over_up.median()
    design=json.loads(source('FROZEN_DESIGN.json'))['original_source_plant_selections']
    cross=pd.read_csv(a.crosswalk)
    joined=plant.merge(cross[['speciesKey','EIV_NUTRI_exact']],left_on='gbif_speciesKey',right_on='speciesKey',how='left',validate='many_to_one')
    reconstructed=[]
    for site in design:
        z=joined[(joined.wwtp_id==site['wwtp_id'])&(joined.river_role==site['river_role'])]
        scored=z.EIV_NUTRI_exact.notna()
        reconstructed.append(dict(outfall=site['wwtp_id'],role=site['river_role'],plant_EIV=z.EIV_NUTRI_exact.mean(),
                                  scored_species=int(scored.sum()),score_coverage=scored.mean(),
                                  scored_records=z.loc[scored,'speciesKey_observation_count'].sum()))
    r=pd.DataFrame(reconstructed).set_index(['outfall','role']).sort_index()
    table=DOC/'evidence/OUTFALL82_SOURCE_SUPPORT_20261120.csv'
    frozen=pd.read_csv(table)
    assert len(frozen)==82 and not frozen.duplicated(['cohort','outfall','role','year']).any()
    de=frozen[frozen.cohort=='DE6'].set_index(['outfall','role']).sort_index()
    assert len(r)==len(de)==18 and r.index.equals(de.index)
    np.testing.assert_allclose(r.to_numpy(float),de[r.columns].to_numpy(float),atol=1e-12,rtol=1e-10)
    expected=pd.read_csv(DOC/'evidence/DE6_interface_contrasts.csv').set_index('outfall')
    r['supported']=(r.scored_species>=10)&(r.score_coverage>=.4)&(r.scored_records>=20)
    for k,row in expected.iterrows():
        u2,u1,d1=[r.loc[(k,t)] for t in ['u2','u1','d1']]
        np.testing.assert_allclose([np.exp(tp[k]),d1.plant_EIV-u1.plant_EIV,u1.plant_EIV-u2.plant_EIV,d1.plant_EIV-2*u1.plant_EIV+u2.plant_EIV],
                                   row[['TP_down_up_ratio','downstream_minus_upstream_plant_EIV','upstream_background_EIV','EIV_curvature']].to_numpy(float),atol=1e-12,rtol=1e-10)
        assert bool(u1.supported and d1.supported)==row.up_down_trait_supported
        assert bool(u2.supported and u1.supported and d1.supported)==row.three_sites_trait_supported
    support=(frozen.lyca_scored_records>=10)&(frozen.lyca_scored_species>=2)&(frozen.lyca_scored_records>=.7*frozen.lyca_records)
    assert (support==frozen.lyca_support).all()
    counts=[]
    for (cohort,period),z in frozen.groupby(['cohort','year']):
        q=z.groupby('outfall').agg(plant=('plant_records','min'),records=('lyca_records','min'),support=('lyca_support','all'))
        counts.append(dict(cohort=cohort,period=period,outfalls=len(q),all_sites_one_lyca=int((q.records>=1).sum()),
                           all_sites_plant20_and_lyca_support=int(((q.plant>=20)&q.support).sum())))
    eligible=expected[expected.up_down_trait_supported]
    result={'plant_sites_reconstructed':18,'outfall_contrasts_verified':6,'up_down_plant_supported':len(eligible),
            'three_site_plant_supported':int(expected.three_sites_trait_supported.sum()),
            'supported_plant_contrasts_below_zero':int((eligible.downstream_minus_upstream_plant_EIV<0).sum()),
            'frozen_butterfly_support':counts,'original_git_sources':sources,
            'crosswalk_sha256':hashlib.sha256(a.crosswalk.read_bytes()).hexdigest(),
            'frozen_support_table_sha256':hashlib.sha256(table.read_bytes()).hexdigest(),
            'limits':['Plant/TP reconstructed from fixed original source-derived tables and external trait crosswalk.',
                      'Butterfly thresholds rechecked from frozen site counts; individual GBIF records not requeried or independently recounted.',
                      'No new ecological test, causal effect or meta-analysis fitted.']}
    (DOC/'evidence/LOCAL_INTERFACE_REPLAY_20261120.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('original_git_sources','limits')}))

if __name__=='__main__':main()
