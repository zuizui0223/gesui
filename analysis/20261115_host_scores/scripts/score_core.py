"""Source-transparent nutrient-affinity scores. Not leaf N or wastewater response."""
from __future__ import annotations
import math, statistics
from typing import Iterable

def numeric_eiv(value):
    if value is None or str(value).strip().casefold() in {'','na','n/a','x','nan'}:
        return None
    v=float(value)
    if not math.isfinite(v): return None
    return v

def norm(name: str) -> str:
    return ' '.join(str(name).strip().split()).casefold()

def match_names(names, traits):
    hits=sorted({norm(x) for x in names if norm(x) in traits})
    if not hits: return {'status':'unmatched','names':[], 'record':None}
    if len(hits)>1: return {'status':'ambiguous','names':hits, 'record':None}
    return {'status':'exact','names':hits, 'record':traits[hits[0]]}

def summarize(host_values: Iterable[tuple[str,float|None]]) -> dict:
    """Equal weight per unique documented host concept; NA never becomes 0.
    Bounds concern only missing scores in this catalogue, NOT undiscovered hosts.
    """
    d={}
    for key, value in host_values:
        v=numeric_eiv(value)
        if v is not None and not 1<=v<=9: raise ValueError('N score outside [1,9]')
        if key in d and d[key]!=v: raise ValueError('Conflicting duplicate host concept')
        d[str(key)]=v
    n=len(d);vals=[v for v in d.values() if v is not None];k=len(vals);total=math.fsum(vals)
    mean=total/k if k else None
    return {'n_hosts':n,'n_scored':k,'sum_N':total if k else None,
            'coverage':k/n if n else None,'mean':mean,
            'median':statistics.median(vals) if k else None,
            'sd':statistics.pstdev(vals) if k else None,
            'min':min(vals) if k else None,'max':max(vals) if k else None,
            'lower_bound':(total+n-k)/n if n else None,
            'upper_bound':(total+9*(n-k))/n if n else None,
            'band':('unscored' if mean is None else 'low' if mean<=3 else 'high' if mean>=7 else 'intermediate'),
            'broad_span_low_high':bool(k and min(vals)<=3 and max(vals)>=7)}

def ambiguous_used_keys(edges):
    """Flag conflicts among actual host-join names, not unused historical aliases."""
    by_name={}
    for edge in edges:
        for name in edge['butterfly_names']:
            by_name.setdefault(norm(name),set()).add(str(edge['butterflyKey']))
    return {key for keys in by_name.values() if len(keys)>1 for key in keys}
