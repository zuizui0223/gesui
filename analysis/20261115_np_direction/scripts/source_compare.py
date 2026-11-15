#!/usr/bin/env python3
"""Source-only EEA original upstream vs downstream reported N/P negative control.

Not a river chemistry measurement and not treatment assignment. Missing report
is not a genuine WWTP-free field control. No ecological outcomes are loaded.
"""
from __future__ import annotations
import argparse, base64, collections, csv, hashlib, io, json, math, pathlib, urllib.request

H7_SHA='a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89'
NP_SHA='012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26'
EEA={
'EEA18': {'sha':'0f30bfbc3dac13f6a182c8a910e9899c04778322','pup':'2018_measured_upstream_P_t_per_year','nup':'2018_measured_upstream_N_t_per_year','pdown':'2018_measured_downstream_P_t_per_year_placebo','ndown':'2018_measured_downstream_N_t_per_year_placebo','expect':2634},
'EEA22': {'sha':'325f199e985d9862ffee5b0edb205c61a7be3833','pup':'upstream_2022_EEA_P_measured_t_per_yr','nup':'upstream_2022_EEA_N_measured_t_per_yr','pdown':'downstream_2022_EEA_P_t_per_yr_placebo','ndown':'downstream_2022_EEA_N_t_per_yr_placebo','expect':2475},
}
FACTOR=1e6/(365.25*86400)
def sha(b:bytes)->str:return hashlib.sha256(b).hexdigest()
def gitblob(b:bytes)->str:return hashlib.sha1(('blob '+str(len(b))+'\0').encode()+b).hexdigest()
def source_data(raw:bytes,kind:str):
    if kind=='mapping':
        if sha(raw)!=H7_SHA:raise ValueError('native EEA H7 source SHA mismatch')
        r=json.loads(raw)['candidate_source_site_period_rows']
        if len(r)!=6330:raise ValueError('native mapping rows drift')
        return r
    if kind=='reported':
        if sha(raw)!=NP_SHA:raise ValueError('original annual N/P report source SHA mismatch')
        r=json.loads(raw)['source_H7_period_reported_N_P_upstream_to_true_original_station']
        if len(r)!=454:raise ValueError('reported annual source rows drift')
        return r
    raise ValueError(kind)
def load_github_blob(commit_sha:str)->list[dict]:
    url='https://api.github.com/repos/zuizui0223/gesui/git/blobs/'+commit_sha
    request=urllib.request.Request(url,headers={'Accept':'application/vnd.github+json','User-Agent':'gesui-frozen-original-EEA-negativecontrol/1'})
    with urllib.request.urlopen(request,timeout=65) as resp: q=json.load(resp)
    b=base64.b64decode(q['content'].replace('\n',''))
    if gitblob(b)!=commit_sha:raise ValueError('GitHub source blob drift '+commit_sha)
    return list(csv.DictReader(io.StringIO(b.decode('utf-8-sig'))))
def bounded(v):
    try:
        x=float(v)
        if math.isfinite(x) and x>=0:return x
    except (ValueError,TypeError):pass
    return None
def process(cohort:list[dict], mapping:list[dict],prior:list[dict],originals:dict[str,list[dict]]):
    pairs={(str(r['HYBAS_ID_L07']),r['PERIOD']) for r in cohort}
    if len(cohort)!=454 or len(pairs)!=454 or len({x[0] for x in pairs})!=232:raise ValueError('cohort ID drift')
    look={k:{r['site_key']:r for r in arr} for k,arr in originals.items()}
    for k,arr in originals.items():
        if len(arr)!=EEA[k]['expect'] or len(look[k])!=len(arr):raise ValueError('raw EEA row drift')
    bypair=collections.defaultdict(list)
    for m in mapping:
        period=m.get('period')
        key=(str(m.get('HYBAS_ID_L07')),period)
        if m.get('panel') in EEA and key in pairs: bypair[key].append(m['native_id'])
    if sum(map(len,bypair.values()))!=3525 or set(bypair)!=pairs:raise ValueError('station-year original drift')
    prior_index={(str(r['HYBAS_ID_L07']),r['PERIOD']):r for r in prior}
    if len(prior_index)!=454 or set(prior_index)!=pairs:raise ValueError('prior source grain drift')
    rows=[]; station_n=0; up_pos_n=0;down_pos_n=0
    for h,p in sorted(pairs):
        panel='EEA18' if p=='2017-2019' else 'EEA22'
        spec=EEA[panel]
        vals=[]
        for site in bypair[(h,p)]:
            r=look[panel].get(site)
            if not r:raise ValueError('missing site from original EEA '+site)
            flow=bounded(r.get('DIS_AV_CMS'))
            inputs={k:bounded(r.get(spec[k])) for k in ('nup','pup','ndown','pdown')}
            if not flow or any(v is None for v in inputs.values()):
                continue
            vals.append({**{k:v/flow*FACTOR for k,v in inputs.items()},'site_key':site,'flow':flow})
        if not vals:raise ValueError('no valid original EEA site in H7 '+h)
        upcount=sum(x['nup']>0 or x['pup']>0 for x in vals)
        downcount=sum(x['ndown']>0 or x['pdown']>0 for x in vals)
        station_n+=len(vals);up_pos_n+=upcount;down_pos_n+=downcount
        fields={k:sum(x[k] for x in vals)/len(vals) for k in ('nup','pup','ndown','pdown')}
        item={'h':h,'period':p,'n_EEA_valid_sites':len(vals),
              'upstream_positive_sites':upcount,'downstream_positive_sites':downcount,
              'upstream_report_positive':int(upcount>0),'downstream_report_positive':int(downcount>0),
              'upstream_N_equiv_mgL':fields['nup'],'upstream_P_equiv_mgL':fields['pup'],
              'downstream_N_equiv_mgL':fields['ndown'],'downstream_P_equiv_mgL':fields['pdown'],
              'upstream_not_observed_concentration':True,'downstream_not_observed_concentration':True}
        previous=prior_index[(h,p)]
        if upcount!=int(previous['source_station_annual_upstream_NP_positive_n']):raise ValueError('upstream station reporter drift')
        for name,old in [('upstream_N_equiv_mgL','source_year_reported_N_mgL_equivalent_site_equal_mean'),('upstream_P_equiv_mgL','source_year_reported_P_mgL_equivalent_site_equal_mean')]:
            if abs(item[name]-float(previous[old]))>1e-9:raise ValueError('prior source upstream equivalence drift '+name)
        rows.append(item)
    report={'schema':'GESUI_ORIGINAL_454_H7_UPSTREAM_DOWNSTREAM_NP_PLACEBO_SOURCE_GATE_V1',
            'original_H7_periods':len(rows),'original_H7':len({x['h'] for x in rows}),
            'n_original_EEA_site_year':station_n,'original_upstream_positive_site_year':up_pos_n,
            'original_downstream_positive_site_year':down_pos_n,
            'H7_period_upstream_positive':sum(x['upstream_report_positive'] for x in rows),
            'H7_period_downstream_positive':sum(x['downstream_report_positive'] for x in rows),
            'H7_period_both_positive':sum(x['upstream_report_positive'] and x['downstream_report_positive'] for x in rows),
            'source_independent_of_biological_outcomes':True,
            'all_0_mean_no_matched_positive_report_not_untreated_control':True,
            'no_actual_receiver_water_chemistry':True,
            'no_randomized_contrast_or_SEM':True,
            'original_GitHub_EEA_source_blobs':{key:x['sha'] for key,x in EEA.items()},
            'original_source_annual_UP_GH_SHA256':NP_SHA,'original_source_native_site_H7_SHA256':H7_SHA}
    if report['original_H7_periods']!=454 or report['original_upstream_positive_site_year']!=132 or report['H7_period_upstream_positive']!=56:
        raise ValueError('frozen observed upstream source denominators changed')
    return rows,report

def main():
    p=argparse.ArgumentParser();p.add_argument('--cohort',type=pathlib.Path,required=True)
    p.add_argument('--mapping',type=pathlib.Path,required=True)
    p.add_argument('--upstream-source',type=pathlib.Path,required=True)
    p.add_argument('--out',type=pathlib.Path,required=True)
    q=p.parse_args()
    with q.cohort.open(newline='',encoding='utf-8-sig') as f:cohort=list(csv.DictReader(f))
    mapping=source_data(q.mapping.read_bytes(),'mapping');prior=source_data(q.upstream_source.read_bytes(),'reported')
    originals={k:load_github_blob(v['sha']) for k,v in EEA.items()}
    rows,report=process(cohort,mapping,prior,originals)
    q.out.mkdir(parents=True,exist_ok=True)
    with (q.out/'H7_454_ORIGINAL_UPSTREAM_DOWNSTREAM_REPORTED_NP_PLACEBO.csv').open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (q.out/'SOURCE_ONLY_UPDOWN_PLACEBO_QA.json').write_text(json.dumps(report,indent=2,sort_keys=True))
    print('SOURCE_PLACEBO_GATE',json.dumps(report,sort_keys=True),flush=True)
if __name__=='__main__':main()
