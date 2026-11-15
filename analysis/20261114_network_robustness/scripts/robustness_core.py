"""Pure functions for posthoc source-record network robustness diagnostics."""
from __future__ import annotations
from collections import Counter
import numpy as np

def group_weights(groups: np.ndarray) -> np.ndarray:
    """Equal total weight per supplied hydrological group; no row independence."""
    groups=np.asarray(groups)
    if groups.ndim!=1 or len(groups)==0:raise ValueError('Nonempty one-dimensional groups required')
    ct=Counter(groups.tolist());w=np.array([1/ct[x] for x in groups],float)
    return w/w.sum()

def retained_permutation(strata: np.ndarray, keep: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Return original column indices; remove exclusions before shuffling."""
    strata=np.asarray(strata);keep=np.asarray(keep,dtype=bool)
    if strata.ndim!=1 or keep.shape!=strata.shape:raise ValueError('Aligned column vectors required')
    take=np.flatnonzero(keep);out=take.copy();ss=strata[take]
    for s in sorted(set(ss)):
        pos=np.flatnonzero(ss==s);out[pos]=rng.permutation(take[pos])
    return out

def edge_contributions(pb: np.ndarray, pp: np.ndarray, adjacency: np.ndarray) -> np.ndarray:
    """Expected source-record candidate links by site and host column."""
    pb=np.asarray(pb,float);pp=np.asarray(pp,float);adjacency=np.asarray(adjacency,float)
    if pb.ndim!=2 or pp.ndim!=2 or adjacency.shape!=(pb.shape[1],pp.shape[1]) or len(pb)!=len(pp):
        raise ValueError('Probability / graph axes differ')
    if not np.isin(adjacency,[0,1]).all():raise ValueError('Binary source graph required')
    return (pb@adjacency)*pp

def hellinger(counts: np.ndarray) -> np.ndarray:
    """Square-root within-family recorded composition; empty family is missing."""
    counts=np.asarray(counts,float)
    if counts.ndim!=2 or np.any(counts<0) or not np.isfinite(counts).all():
        raise ValueError('Finite nonnegative source count matrix required')
    total=counts.sum(1,keepdims=True)
    return np.sqrt(np.divide(counts,total,out=np.full_like(counts,np.nan),where=total>0))
