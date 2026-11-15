"""Fixed-trait composition metrics and geographically purged prediction primitives."""
from __future__ import annotations
import hashlib
import numpy as np


def aggregate_traits(counts: np.ndarray, scores: np.ndarray, families: np.ndarray,
                     host_coverage: np.ndarray, raw_totals: np.ndarray) -> dict:
    b=np.asarray(counts,dtype=float);s=np.asarray(scores,dtype=float)
    f=np.asarray(families,dtype=str);hc=np.asarray(host_coverage,dtype=float)
    raw=np.asarray(raw_totals,dtype=float)
    if b.ndim!=2 or b.shape[1]!=len(s) or len(f)!=len(s) or len(hc)!=len(s) or len(raw)!=len(b):
        raise ValueError('Inconsistent source matrix dimensions')
    if not np.isfinite(b).all() or (b<0).any() or (b!=np.floor(b)).any() or (raw<b.sum(1)).any():
        raise ValueError('Invalid count or raw source denominator')
    known=np.isfinite(s)
    if ((s[known]<1)|(s[known]>9)).any():raise ValueError('HNS outside 1–9')
    k=np.where(known,s,0.);weights=b*known;presence=(b>0)*known
    den=weights.sum(1);nknown=presence.sum(1);nseen=(b>0).sum(1)
    def div(a,d):return np.divide(a,d,out=np.full(np.shape(a),np.nan,dtype=float),where=d>0)
    record=div(weights@k,den);incidence=div(presence@k,nknown)
    fams=np.array(sorted(set(f)))
    shares=np.column_stack([div(weights[:,f==g].sum(1),den) for g in fams])
    ishares=np.column_stack([div(presence[:,f==g].sum(1),nknown) for g in fams])
    means=np.array([np.mean(s[(f==g)&known]) if ((f==g)&known).any() else np.nan for g in fams])
    between=shares@np.nan_to_num(means,nan=0.)
    return {'record_mean':record,'incidence_mean':incidence,
            'record_coverage_raw':div(den,raw),
            'record_coverage_identified':div(den,b.sum(1)),
            'species_coverage':div(nknown,nseen),
            'scored_records':den,'scored_species':nknown,
            'host_completeness':div(weights@np.nan_to_num(hc,nan=0.),den),
            'record_family_shares':shares,'incidence_family_shares':ishares,
            'families':fams,'family_reference_means':means,
            'between_family':between,'within_family':record-between}


def blocked_splits(h6: np.ndarray, countries: np.ndarray, mode: str='H6'):
    h=np.asarray(h6,dtype=str);c=np.asarray(countries,dtype=str)
    if len(h)!=len(c):raise ValueError('Grouping length mismatch')
    if mode=='H6':g=np.array([int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5 for x in h])
    elif mode=='country':g=c
    else:raise ValueError('Unknown validation mode')
    for value in sorted(set(g)):
        test=g==value;train=~test
        if mode=='country':train &= ~np.isin(h,h[test])
        if set(h[train])&set(h[test]):raise AssertionError('Hydrological parent leakage')
        yield value,train,test


def ridge_fit_predict(X: np.ndarray,y: np.ndarray,T: np.ndarray,weights: np.ndarray,
                      alpha: float=10.) -> np.ndarray:
    penalty=np.eye(X.shape[1])*alpha;penalty[0,0]=0.
    coef=np.linalg.solve(X.T@(weights[:,None]*X)+penalty,X.T@(weights*y))
    return T@coef
