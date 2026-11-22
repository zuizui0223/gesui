#!/usr/bin/env python3
"""Posthoc source-locked diagnostic for four butterfly-derived covariates.

Requires the same six archived input files and original score table as analyse.py.
Does not refit a new ecological hypothesis or infer a causal mediator.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path

import numpy as np
import pandas as pd

from analyse import FAMS, CONTROLS, restore_sources, build_pair_strata, source_summary
from inference import prepare_design, report_fit

FOUR=('nButterfly','hns_coverage','hns_host_completeness','hns_scored_species')
BASE=[k for k in CONTROLS if k not in FOUR]
ORIGINAL_BETA=0.3456186953934604
ORIGINAL_BASE=0.1334352523119997


def fit(strata,cells,panel,controls):
    result=report_fit(prepare_design(strata,cells,panel,controls))
    if result.get('status')!='FIT':
        raise RuntimeError('Fixed diagnostic fit STOP: '+str(result.get('status')))
    return {'controls':list(controls),'beta':result['beta'],
            'odds_ratio':result['odds_ratio'],
            'odds_ratio95_H6_cluster_t':result['odds_ratio95_H6_cluster_t'],
            'nominal_p':result['nominal_two_sided_H6_cluster_t_p'],
            'H6_clusters':result['H6_clusters'],
            'condition':result['information_condition_number']}


def period_specific_scores(source:Path,cells:pd.DataFrame,score:pd.DataFrame):
    """Reconstruct four original-source period-specific indices, with no new taxa."""
    key_cov=dict(zip(score.butterflyKey.astype(str),score.trait_coverage.astype(float)))
    source_facets={}
    for fn,n in [('GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json',2808),
                 ('MERGED_ORIGINAL_30H7_TWO_PERIOD_SIX_FAMILY_360_SOURCE_FACETS.json',360)]:
        raw=json.loads((source/fn).read_text())
        if len(raw)!=n:raise RuntimeError('facet source size mismatch')
        for r in raw:
            key=(str(r['HYBAS_ID_L07']),str(r['period']),str(r['group']))
            if key in source_facets or r.get('source_status')!='ok':
                raise RuntimeError('duplicate/failed source facet')
            source_facets[key]=r
    output=[]
    for h,period in cells[['h','period']].itertuples(index=False,name=None):
        total=0;accepted=0;accepted_taxa=0;weighted_score_coverage=0.
        for fam in FAMS:
            facet=source_facets[(str(h),str(period),fam)]
            total+=int(facet['source_query_occurrence_count'])
            for key,count in facet['source_speciesKey_counts'].items():
                if int(count)>0 and str(key) in key_cov:
                    accepted+=int(count);accepted_taxa+=1
                    weighted_score_coverage+=int(count)*key_cov[str(key)]
        if not total or not accepted:raise RuntimeError('source original HNS eligible records empty')
        output.append(dict(h=str(h),period=str(period),
            nButterfly=float(np.log1p(total)),hns_coverage=float(accepted/total),
            hns_host_completeness=float(weighted_score_coverage/accepted),
            hns_scored_species=float(accepted_taxa)))
    return pd.DataFrame(output)


def main(source:Path,output:Path):
    cells,counts,scores,panel,source_receipt=restore_sources(source)
    original=build_pair_strata(cells,counts,scores)
    gate=source_summary(original)
    if not gate['support_gate_pass'] or gate['informative_strata']!=353 or gate['informative_H6']!=34:
        raise RuntimeError('original source support changed; analysis STOP')
    variants={'baseline_9':BASE,'primary_13':CONTROLS}
    variants.update({'add_'+k:BASE+[k] for k in FOUR})
    variants.update({'drop_'+k:[t for t in CONTROLS if t!=k] for k in FOUR})
    variants['external_geography_5']=['W','crop','urban','naturalout','area']
    result={k:fit(original,cells,panel,controls) for k,controls in variants.items()}
    if abs(result['primary_13']['beta']-ORIGINAL_BETA)>1e-6 or abs(result['baseline_9']['beta']-ORIGINAL_BASE)>1e-6:
        raise RuntimeError('the two published models did not reproduce (STOP)')
    period=period_specific_scores(source,cells,scores)
    averaged=period.groupby('h')[list(FOUR)].mean()
    orig=panel.set_index('h')[list(FOUR)]
    replay={k:float(np.max(np.abs(averaged[k]-orig[k]))) for k in FOUR}
    if max(replay.values())>1e-10:
        raise RuntimeError('original two-period mean control not reconstructed (STOP)')
    for variant in ['same_period','opposite_period']:
        current=period.copy()
        if variant=='opposite_period':
            swap={'2017-2019':'2020-2022','2020-2022':'2017-2019'}
            current['period']=current.period.map(swap)
        test_cells=cells.drop(columns=list(FOUR)).merge(current,on=['h','period'],validate='one_to_one')
        strata=build_pair_strata(test_cells,counts,scores)
        if source_summary(strata)!=gate:
            raise RuntimeError('comparison sample drift for '+variant)
        result[variant+'_13']=fit(strata,test_cells,panel,CONTROLS)
    # descriptive whole-H7 correlations; also H6-demeaned, never a mediation test
    informative_h7={h for r in original if abs(r['trait_difference'])>1e-12 for h in r['h']}
    sub=panel[panel.h.isin(informative_h7)].copy()
    correlations={}
    for v in FOUR:
        x=sub.plant_n_incidence-sub.groupby('h6').plant_n_incidence.transform('mean')
        y=sub[v]-sub.groupby('h6')[v].transform('mean')
        correlations[v]={'pearson_H7':float(sub[['plant_n_incidence',v]].corr().iloc[0,1]),
            'within_H6_demeaned_r':float(np.dot(x,y)/np.sqrt(np.dot(x,x)*np.dot(y,y))),
            'mean_by_plantN_quartile':sub.groupby(pd.qcut(sub.plant_n_incidence,4,labels=False))[v].mean().tolist()}
    response={'status':'EXPLORATORY_POSTHOC_CONTROL_DEPENDENCE',
        'original_source':source_receipt,'source_gate':gate,
        'original_primary_and_baseline_exactly_reproduced':True,
        'period_source_replay_max_error':replay,'models':result,'correlations':correlations,
        'caveats':['All fits use unchanged 353 comparisons and the same H6/period/paired-species identities.',
         'These model variations are posthoc diagnostics, not prespecified primary claims.',
         'The four butterfly-derived controls use the same GBIF records that define the response.',
         'Opposite-period controls are NOT independent spatial negative controls.',
         'Conditional recording OR is not true abundance, occupancy, feeding or wastewater impact.']}
    output.mkdir(parents=True,exist_ok=True)
    (output/'control_dependence.json').write_text(json.dumps(response,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    pd.DataFrame([dict(model=k,**v) for k,v in result.items()]).to_csv(output/'control_dependence_models.csv',index=False)
    print(json.dumps({'source':gate,'models':{k:dict(OR=v['odds_ratio'],beta=v['beta'],CI=v['odds_ratio95_H6_cluster_t']) for k,v in result.items()},
          'period_replay_error':replay},ensure_ascii=False),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args()
    main(args.source,args.output)

