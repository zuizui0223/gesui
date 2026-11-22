#!/usr/bin/env python3
"""New exploratory butterfly-congener recording test; immutable original sources.

The primary unit is a discordant congeneric pair comparison within H6 and
period, not an animal individual or a biological absence.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SPEC_COMMIT = '89c349c'
CONTROLS = ['W','crop','urban','naturalout','area','nPlant','nButterfly',
            'nAnimal','plant_coverage','hns_coverage','hns_host_completeness',
            'plant_scored_species','hns_scored_species']
FAMS = ['Hesperiidae','Lycaenidae','Nymphalidae','Papilionidae','Pieridae']
PERIODS = ['2017-2019','2020-2022']
FILES = {
 'GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json': '251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4',
 'MERGED_ORIGINAL_30H7_TWO_PERIOD_SIX_FAMILY_360_SOURCE_FACETS.json': 'e076f6c8e2d93282be19aea3c75ae8431302318f068b9d75e278e2a21c20aafb',
 'GENUINE_TERRESTRIAL_ONLY_183_H7_40_H6_PAIRED_H7_MEANS.csv': 'c0840b1c5c6e85342fea41d57dc40fa5ec295fabf6bcd6908a4072557db9a736',
 'GESUI_303_HOST_NUTRIENT_SCORES.csv': '00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb',
 'EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv': 'ab45f36c0dad047442bb958bcef3135b45014bafec8a70811bcafb6109a02679',
 'FIXED_183_H7_IN_40_NATIVE_H6_SOURCE_COHORT.csv': '8107761a7e943b0797ebc3e53e98baebe51ae6ebb66eece394d35ea0dbcd2f7b',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def eligible_scores(frame, coverage_min=0.0):
    d=frame.copy()
    bad=d.actual_used_alias_conflict.astype(str).str.lower().isin(['true','1','yes'])
    d=d[(d.taxon_status=='ACCEPTED') & ~bad & d.host_nutrient_mean.notna() &
        (d.trait_coverage>=coverage_min)].copy()
    d['butterflyKey']=d.butterflyKey.astype(str)
    if d.butterflyKey.duplicated().any():
        raise ValueError('duplicate butterfly species key')
    if not d.host_nutrient_mean.between(1,9).all():
        raise ValueError('HNS out of original EIV numeric range')
    d['genus']=d.accepted_name.str.split().str[0]
    return d.sort_values('butterflyKey').reset_index(drop=True)


def build_pair_strata(cells, counts, scores, covariates=CONTROLS):
    """Keep all discordant within-H6/period comparisons with both identities.

    No covariate-response slope is inspected here. A/B orientation is lexical
    speciesKey order. Source counts >1 affect neither y nor eligibility.
    """
    if cells.duplicated(['h','period']).any():
        raise ValueError('duplicate source H7-period')
    score=scores.set_index('butterflyKey')
    genera={g:sorted(d.butterflyKey.tolist()) for g,d in scores.groupby('genus')
            if len(d)>=2 and d.host_nutrient_mean.nunique()>1}
    present={}
    for (h,p), d in counts[counts['count']>0].groupby(['h','period']):
        present[(str(h),str(p))]=set(d.butterflyKey.astype(str)) & set(score.index)
    out=[]
    for (region,period), cell in cells.groupby(['h6','period'],sort=True):
        cell=cell.sort_values('h')
        region_species=set().union(*(present.get((str(h),str(period)),set()) for h in cell.h))
        for genus, members in sorted(genera.items()):
            regional=sorted(set(members)&region_species)
            for sa,sb in combinations(regional,2):
                rows=[]; ys=[]
                for row in cell.to_dict('records'):
                    seen=present.get((str(row['h']),str(period)),set())
                    a,b=sa in seen,sb in seen
                    if a!=b:
                        rows.append(row);ys.append(int(a))
                if not ys or min(ys)==max(ys):
                    continue
                out.append(dict(pair=sa+'__'+sb,species_a=sa,species_b=sb,
                    genus=genus,h6=str(region),period=str(period),
                    h=[str(r['h']) for r in rows],y=ys,
                    plant_n_incidence=[float(r['plant_n_incidence']) for r in rows],
                    controls={k:[float(r[k]) for r in rows] for k in covariates},
                    trait_difference=float(score.loc[sa,'host_nutrient_mean']-score.loc[sb,'host_nutrient_mean'])))
    pair_counts=Counter(r['pair'] for r in out)
    if out:
        norm=len(out)/len(pair_counts)
        for r in out:
            r['weight']=norm/pair_counts[r['pair']]
    return out


def restore_sources(source):
    receipts=[]
    for filename,expected in FILES.items():
        path=source/filename
        observed=sha(path)
        if observed!=expected:
            raise RuntimeError('immutable source SHA256 mismatch: '+filename)
        receipts.append(dict(file=filename,sha256=observed))
    panel=pd.read_csv(source/'GENUINE_TERRESTRIAL_ONLY_183_H7_40_H6_PAIRED_H7_MEANS.csv',dtype={'h':str,'h6':str})
    cohort=pd.read_csv(source/'FIXED_183_H7_IN_40_NATIVE_H6_SOURCE_COHORT.csv',dtype={'h':str,'h6':str})
    if len(panel)!=183 or panel.h.nunique()!=183 or panel.h6.nunique()!=40:
        raise RuntimeError('original183/40 source changed')
    if set(zip(panel.h,panel.h6))!=set(zip(cohort.h,cohort.h6)):
        raise RuntimeError('source geography and original panel disagree')
    all_scores=pd.read_csv(source/'GESUI_303_HOST_NUTRIENT_SCORES.csv',encoding='utf-8-sig',dtype={'butterflyKey':str})
    scores=eligible_scores(all_scores)
    if len(all_scores)!=303 or len(scores)!=186:
        raise RuntimeError('original303/186 butterfly score roster changed')
    plant=pd.read_csv(source/'EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv',dtype={'speciesKey':str})
    scored=plant[plant.EIV_NUTRI_exact.notna()]
    if len(plant)!=11511 or len(scored)!=4011 or plant.speciesKey.duplicated().any():
        raise RuntimeError('original11511/4011 plant crosswalk changed')
    pn=dict(zip(scored.speciesKey,scored.EIV_NUTRI_exact))
    bn=dict(zip(scores.butterflyKey,scores.host_nutrient_mean))
    bf=dict(zip(scores.butterflyKey,scores.family))
    facets={}
    for filename,n in [('GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json',2808),
                       ('MERGED_ORIGINAL_30H7_TWO_PERIOD_SIX_FAMILY_360_SOURCE_FACETS.json',360)]:
        raw=json.loads((source/filename).read_text())
        if len(raw)!=n:
            raise RuntimeError('original facet number changed')
        for r in raw:
            key=(str(r['HYBAS_ID_L07']),str(r['period']),str(r['group']))
            if key in facets or r.get('source_status')!='ok':
                raise RuntimeError('duplicated/failed source facet')
            facets[key]=r
    cells=[];counts=[]
    for row in panel.to_dict('records'):
        h,region=row['h'],row['h6']
        for period in PERIODS:
            rawplants=facets[(h,period,'Plantae')]['source_speciesKey_counts']
            pvals=[(pn[k],int(v)) for k,v in rawplants.items() if int(v)>0 and k in pn]
            nrawplant=sum(int(v)>0 for v in rawplants.values())
            if len(pvals)<10 or len(pvals)/nrawplant<.4:
                raise RuntimeError('original plant score source gate failed')
            bvals=[]; seen_keys=set()
            for family in FAMS:
                raw=facets[(h,period,family)]['source_speciesKey_counts']
                for key,count in raw.items():
                    if int(count)<=0:
                        continue
                    if key in seen_keys:
                        raise RuntimeError('butterfly counted in multiple source families')
                    seen_keys.add(key)
                    counts.append(dict(h=h,h6=region,period=period,butterflyKey=key,family=family,count=int(count)))
                    if key in bn:
                        if bf[key]!=family:
                            raise RuntimeError('butterfly family mismatch')
                        bvals.append((bn[key],int(count)))
            if not bvals:
                raise RuntimeError('empty source scored butterfly records')
            cells.append(dict(h=h,h6=region,period=period,
                **{k:float(row[k]) for k in CONTROLS},
                plant_n_incidence=float(np.mean([q[0] for q in pvals])),
                plant_n_record=float(np.average([q[0] for q in pvals],weights=[q[1] for q in pvals])),
                hns_incidence=float(np.mean([q[0] for q in bvals])),
                hns_record=float(np.average([q[0] for q in bvals],weights=[q[1] for q in bvals]))))
    cells=pd.DataFrame(cells);counts=pd.DataFrame(counts)
    if len(cells)!=366 or cells.duplicated(['h','period']).any():
        raise RuntimeError('original366 source cells not unique')
    check=cells.groupby('h')[['plant_n_incidence','plant_n_record','hns_incidence','hns_record']].mean()
    errors={k:float(np.max(np.abs(check[k]-panel.set_index('h')[k]))) for k in check.columns}
    if max(errors.values())>1e-10:
        raise RuntimeError('original community score reconstruction failed')
    for k in CONTROLS+['plant_n_incidence']:
        if not np.isfinite(cells[k]).all():
            raise RuntimeError('missing source covariate '+k)
    return cells,counts,scores,panel,dict(source_files=receipts,original_CWM_replay_max_errors=errors,
        original_H7=183,original_H6=40,original_cells=366,positive_butterfly_record_cells=len(counts),
        numeric_accepted_butterfly_species=len(scores),butterfly_genera=int(scores.genus.nunique()))


def source_summary(strata):
    variable=[r for r in strata if abs(r['trait_difference'])>1e-12]
    return dict(informative_strata=len(strata),
        discordant_pair_cell_rows=sum(len(r['h']) for r in strata),
        informative_H6=len({r['h6'] for r in variable}),
        informative_H7=len({h for r in variable for h in r['h']}),
        informative_genera=len({r['genus'] for r in variable}),
        all_pair_genera=len({r['genus'] for r in strata}),
        species_in_pairs=len({s for r in strata for s in [r['species_a'],r['species_b']]}),
        unique_species_pairs=len({r['pair'] for r in strata}),
        original_equal_HNS_strata=sum(abs(r['trait_difference'])<=1e-12 for r in strata),
        min_discordant_H7_per_stratum=min((len(r['h']) for r in strata),default=0),
        max_discordant_H7_per_stratum=max((len(r['h']) for r in strata),default=0),
        support_gate_pass=len({r['h6'] for r in variable})>=20 and len({r['genus'] for r in variable})>=10)


def run_source(source,output):
    cells,counts,scores,panel,receipt=restore_sources(source)
    strata=build_pair_strata(cells,counts,scores)
    gate=dict(specification_commit=SPEC_COMMIT,exploratory_after_known_community_results=True,
              source=receipt,primary=source_summary(strata))
    high=scores[scores.trait_coverage>=.8].copy()
    gate['coverage80']=source_summary(build_pair_strata(cells,counts,high))
    output.mkdir(parents=True,exist_ok=True)
    write_json(output/'source_gate.json',gate)
    report=[]
    for pair, ss in pd.Series(strata,dtype=object).groupby(lambda i: strata[i]['pair']):
        members=ss.tolist(); r=members[0]
        report.append(dict(pair=pair,genus=r['genus'],species_a=r['species_a'],species_b=r['species_b'],
            trait_difference=r['trait_difference'],n_strata=len(members),n_H6=len({x['h6'] for x in members}),
            n_discordant_rows=sum(len(x['h']) for x in members),total_stratum_weight=sum(x['weight'] for x in members)))
    pd.DataFrame(report).to_csv(output/'pair_support.csv',index=False)
    print('SOURCE_GATE='+json.dumps(gate,ensure_ascii=False,allow_nan=False),flush=True)
    return cells,counts,scores,panel,strata,gate


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--stage',choices=['source','all'],default='source')
    a=parser.parse_args()
    values=run_source(a.source,a.output)
    if a.stage=='all':
        from inference import run_all
        run_all(*values,a.output,CONTROLS,build_pair_strata,source_summary)


if __name__=='__main__':
    main()
