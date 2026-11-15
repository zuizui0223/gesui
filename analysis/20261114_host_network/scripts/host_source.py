"""Pure taxonomy gates; exact counterpart of the remote source collector."""
import re

def binomial(value):
    s=' '.join(str(value or '').strip().split())
    if not re.fullmatch(r'[A-Z][a-z-]+ [a-z][a-z-]+',s):return None
    if s.split()[1] in {'sp','spp','cf','aff','indet'}:return None
    return s.lower()

def plant_key(name,j):
    if j.get('matchType')!='EXACT' or j.get('rank')!='SPECIES' or j.get('kingdom')!='Plantae':return None
    if binomial(j.get('canonicalName'))!=binomial(name):return None
    if j.get('status') not in {'ACCEPTED','SYNONYM','HOMOTYPIC_SYNONYM','HETEROTYPIC_SYNONYM','PROPARTE_SYNONYM'}:return None
    key=j.get('acceptedUsageKey') or (j.get('usageKey') if j.get('status')=='ACCEPTED' else None)
    return str(key) if isinstance(key,int) and key>0 else None
