#!/usr/bin/env python3
"""Posthoc source-frozen taxon influence in butterfly HNS vs wastewater geography.

Replaces one response metric at a time on the same 454 H7 x period rows;
no species is treated as if biologically removed, and no SEM paths estimated.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd

PANEL_SHA = "1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
TRAIT_SHA = "051895bcc416ae6387ce5a16c3f0e0babf51a37fa62f5245267015fcb2fb1e03"
DECOMP_SHA = "daf4904f043eb3785142e507168f038177e799f588bc04a0b748a4ec1607b866"
GEO_SHA = "e490e110011b99830e31fff680a3482f6535d8010d989f2da4872aedd0557e24"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cwm_from_matrix(counts, scores, excluded_indices):
    """Mean HNS among currently score-known source records, excluding diagnostic taxa."""
    B = np.asarray(counts, dtype=float)
    s = np.asarray(scores, dtype=float)
    if B.ndim != 2 or s.ndim != 1 or B.shape[1] != len(s):
        raise ValueError("Original butterfly species dimensions differ")
    if not np.isfinite(B).all() or (B < 0).any():
        raise ValueError("Invalid original occurrence record count")
    mask = np.isfinite(s)
    for ix in excluded_indices:
        if ix < 0 or ix >= len(s):
            raise ValueError("Invalid excluded butterfly index")
        mask[int(ix)] = False
    denominator = B[:, mask].sum(axis=1)
    if np.any(denominator <= 0):
        raise ValueError("A basin-period has zero known HNS support after species exclusion; do not impute zero")
    return B[:, mask] @ s[mask] / denominator


def ranked_taxa(keys, contributions):
    """Order by absolute contribution in prior analysis; deterministic name tie-break."""
    if len(keys) != len(contributions):
        raise ValueError("Taxon names/contributions length mismatch")
    pairs = [(str(k), float(c)) for k,c in zip(keys, contributions)]
    if len(set(k for k,_ in pairs)) != len(pairs):
        raise ValueError("Duplicate butterfly taxon key")
    if not all(np.isfinite(c) for _,c in pairs):
        raise ValueError("Non-finite source contribution")
    return [k for k,c in sorted(pairs,key=lambda t:(-abs(t[1]),t[0]))]


def read_inputs(panel, trait, decomp, geo):
    for p, expected in [(panel,PANEL_SHA),(trait,TRAIT_SHA),(decomp,DECOMP_SHA),(geo,GEO_SHA)]:
        if sha(p) != expected:
            raise ValueError("Original source SHA drift: " + str(p))
    d=pd.read_csv(panel,dtype={'h':str,'h6':str})
    assert len(d)==454 and d.h.nunique()==232 and d.h6.nunique()==90
    assert not d.duplicated(['h','period']).any()
    assert sum(d.groupby('h').size()==2)==222
    z=np.load(trait)
    B=z['B'];s=z['scores_primary'];keys=[str(k) for k in z['keys']]
    if B.shape != (454,303) or len(keys)!=303 or len(set(keys))!=303:
        raise ValueError("Original scored species matrix drift")
    old=cwm_from_matrix(B,s,[])
    np.testing.assert_allclose(old,d.primary_record_mean.to_numpy(float),atol=1e-10,rtol=0)
    drv=pd.read_csv(decomp,dtype={'butterflyKey':str})
    dv=drv.loc[drv.kind.eq('record')].set_index('butterflyKey')
    if len(dv)!=303 or set(dv.index)!=set(keys):
        raise ValueError("Original signed species contributions do not align with source")
    prior_beta=float(dv.centered_HNS_slope_contribution.sum())
    np.testing.assert_allclose(prior_beta,0.1300624446205031,atol=1e-10)
    ranked=ranked_taxa(keys,dv.loc[keys].centered_HNS_slope_contribution.to_numpy(float))
    g=pd.read_csv(geo,dtype={'h':str})
    if len(g)!=234 or g.h.duplicated().any():
        raise ValueError('Original basin polygon centroid source invalid')
    d=d.merge(g[['h','longitude','latitude']],on='h',how='left',validate='many_to_one',sort=False)
    if d[['longitude','latitude']].isna().any().any():
        raise ValueError('Missing source original basin geometry')
    d['lon2']=d.longitude**2;d['lat2']=d.latitude**2;d['lonlat']=d.longitude*d.latitude
    if len(d)!=454:
        raise ValueError('Response join changed source population')
    return d,B,s,keys,ranked,dv,z['family']


def main(args):
    # Existing original source HNS fitting code remains separate and fixed.
    sys.path.insert(0,str(args.model_scripts))
    from model_core import model
    df,B,s,keys,ranked,driver,family=read_inputs(args.panel,args.trait,args.decomp,args.geo)
    index={k:i for i,k in enumerate(keys)}
    geo_terms=['longitude','latitude','lon2','lat2','lonlat']
    top=[k for k in ranked if np.isfinite(s[index[k]])]
    if len(top)<30:raise ValueError('Too few source scored taxa for influence sensitivity')
    output=[]
    for n in [0,1,5,10,20,30]:
        drops=top[:n]
        d=df.copy();d['taxon_drop_HNS']=cwm_from_matrix(B,s,[index[k] for k in drops])
        result=model(d,'taxon_drop_HNS',extras=geo_terms,country=False)
        b=result['beta'];cv=result['H6']
        output.append({
          'dropped_top_absolute_previous_contributions':n,
          'new_source_rows':len(d),'n_H7':int(d.h.nunique()),'n_H6':int(d.h6.nunique()),
          'original_source_denominator_not_zero':True,
          'HNS_response_min':float(d.taxon_drop_HNS.min()),
          'HNS_response_max':float(d.taxon_drop_HNS.max()),
          'W_beta':b['W_beta_per1SD'],'W_beta_H6_descriptive95_low':b['W_beta_cluster95'][0],
          'W_beta_H6_descriptive95_high':b['W_beta_cluster95'][1],
          'W_added_MSE_gain_percent':cv['gain_pct'],
          'W_gain_fixed_error95_low':cv['fixed_error_bootstrap95'][0],
          'W_gain_fixed_error95_high':cv['fixed_error_bootstrap95'][1],
          'drop_keys':drops,
        })
        print('TOP_DROP',n,'beta',round(b['W_beta_per1SD'],7),'gain',round(cv['gain_pct'],5),flush=True)
    observed=output[0]['W_beta']
    np.testing.assert_allclose(observed,0.10802204501635151,atol=1e-7)
    # The previously selected rank is posthoc; never relabel these as preregistered null experiments.
    # Each drop holds the same original 454 rows, underlying GBIF record matrix and all control covariates fixed.
    report={'schema':'GESUI_20261115_HNS_TAXON_ABLATION_POSTHOC_V1',
      'source_SHA256':{'H7_panel':PANEL_SHA,'trait_matrix':TRAIT_SHA,'original_species_contribution':DECOMP_SHA,'geography':GEO_SHA},
      'n_H7_period':454,'n_H7':232,'n_H6':90,
      'selection':'Ablate taxa ranked by |their posthoc original HNS W slope contribution|; can overstate sensitivity',
      'statistical_meaning':'Outcome-index sensitivity only; NOT ecological species removal or wastewater exposure',
      'all_source_covariates_and_family_share_frozen':True,
      'source_record_weighting_not_individual_abundance':True,
      'models':output,
      'top20_species':[{'taxonKey':k,'family':str(family[index[k]]),
         'HNS':float(s[index[k]]),
         'source_total_records':int(B[:,index[k]].sum()),
         'original_contribution':float(driver.loc[k,'centered_HNS_slope_contribution'])} for k in top[:20]]}
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/'RESULTS.json').write_text(json.dumps(report,indent=2,sort_keys=True))
    pd.DataFrame([{k:v for k,v in o.items() if k!='drop_keys'} for o in output]).to_csv(args.out/'TOP_TAXON_ABLATION.csv',index=False)
    pd.DataFrame(report['top20_species']).to_csv(args.out/'TOP20_SOURCE_TAXA.csv',index=False)
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--panel',type=Path,required=True)
    p.add_argument('--trait',type=Path,required=True)
    p.add_argument('--decomp',type=Path,required=True)
    p.add_argument('--geo',type=Path,required=True)
    p.add_argument('--model-scripts',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    main(p.parse_args())
