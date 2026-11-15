"""Exact record-rarefaction candidate-network metrics. Not observed feeding."""
import numpy as np

def inclusion(N, counts, m):
    """P(at least one of `counts` items in m of N records), stable product."""
    n=np.asarray(counts,dtype=float)
    if int(N)!=N or m<0 or int(m)!=m or N<m or np.any(n<0) or np.any(n>N) or np.any(n!=np.floor(n)):
        raise ValueError('Invalid hypergeometric record counts / depth')
    if m==0:return np.zeros_like(n)
    certain=(N-n<m)
    logp=np.zeros_like(n)
    with np.errstate(divide='ignore',invalid='ignore'):
        for j in range(m):logp+=np.log1p(-np.minimum(n,N-j)/(N-j))
    return np.where(certain,1.,-np.expm1(logp))

def metrics(b,p,A,mb=20,mp=50,NB=None,NP=None):
    b=np.asarray(b);p=np.asarray(p);A=np.asarray(A)
    if A.shape!=(len(b),len(p)) or not np.isin(A,[0,1]).all():raise ValueError('Binary adjacency shape differs')
    NB=int(b.sum()) if NB is None else int(NB);NP=int(p.sum()) if NP is None else int(NP)
    pb=inclusion(NB,b,mb);pp=inclusion(NP,p,mp)
    known=A.sum(1)>0;den=float(pb[known].sum())
    q=inclusion(NP,A@p,mp)
    local=(b>0)[:,None]&(p>0)[None,:]&(A>0)
    return {'expected_links':float(pb@A@pp),
            'host_support':float(pb[known]@q[known]/den) if den>0 else np.nan,
            'catalogue_coverage':den/float(pb.sum()) if pb.sum()>0 else np.nan,
            'candidate_links':int(local.sum()),'butterfly_S20':float(pb.sum()),'plant_S50_in_matrix':float(pp.sum()),
            'known_butterflies_present':int(((b>0)&known).sum()),
            'known_butterflies_with_local_host':int(((b>0)&known&(A@p>0)).sum()),
            'mean_host_degree':float(pb[known]@A.sum(1)[known]/den) if den>0 else np.nan}

def shuffle_columns(strata,rng):
    strata=np.asarray(strata);idx=np.arange(len(strata))
    for g in sorted(set(strata)):
        take=np.flatnonzero(strata==g);idx[take]=rng.permutation(take)
    return idx
