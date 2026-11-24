#!/usr/bin/env python3
"""Exact source-only all-family host nutrient matching, no outcomes fitted."""
from pathlib import Path
import sys
P=Path(sys.argv[1])
prior=P/'gesui_clarke_appendix_mapping.py'
exec(compile(prior.read_text().split('rows=[];mapped=[]')[0],str(prior),'exec'))
import numpy as np
# Original exact observed-flora EIV data.
zf=P/'gesui_eiv_source.zip'; assert sha(zf.read_bytes())=='1d3051d4c527421bec3ec3ba917422a9162afb81cc7f38db0d8130287174a14c'
with zipfile.ZipFile(zf) as z:
 raw=z.read('EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv');assert sha(raw)=='ab45f36c0dad047442bb958bcef3135b45014bafec8a70811bcafb6109a02679'
 tx=pd.read_csv(io.BytesIO(raw),dtype={'speciesKey':str})
assert len(tx)==11511 and tx.speciesKey.nunique()==11511 and tx.EIV_NUTRI_exact.notna().sum()==4011
meta=tx.set_index('speciesKey').to_dict('index'); namekeys=collections.defaultdict(set)
for x in tx.itertuples(index=False):
 n=norm(x.gbif_canonical_name)
 if re.fullmatch(r'[a-z][a-z-]+ [a-z][a-z-]+',n):namekeys[n].add(x.speciesKey)
def accepted(k):
 m=meta.get(k,{})
 return (m.get('GBIF_source_status')=='OK' and m.get('gbif_rank')=='SPECIES' and m.get('gbif_taxon_status')=='ACCEPTED' and m.get('gbif_phylum')=='Tracheophyta' and m.get('gbif_kingdom')=='Plantae')
oldgenfam=collections.defaultdict(set);sourcefamilies=collections.defaultdict(set)
for x in pmeta:
 m=x['match'];k=x.get('accepted_speciesKey');n=norm(m.get('canonicalName',''));fam=m.get('family','')
 if k and fam and m.get('rank')=='SPECIES' and m.get('matchType')=='EXACT' and m.get('confidence')==100 and re.fullmatch(r'[a-z][a-z-]+ [a-z][a-z-]+',n):oldgenfam[n.split()[0]].add(fam)
for x in source.itertuples(index=False):
 n=plant_binomial(x.accepted_plant)
 if n and x.plant_family:sourcefamilies[n].add(x.plant_family)
rows=[];mapped=[]
for i,x in enumerate(source.itertuples(index=False)):
 k=bmapped.get(x.accepted_butterfly);n=plant_binomial(x.accepted_plant);p=None;match=''
 if status[x.evidence]>3:reason='excluded_evidence'
 elif not x.countries.strip():reason='no_nonquestionable_country'
 elif not k:reason='unmapped_butterfly'
 elif not n or len(plant_forms[n])!=1:reason='rank_or_ambiguous_binomial'
 else:
  prev=pname.get(n,set());cand=namekeys.get(n,set())
  if len(prev)==1:
   p=next(iter(prev));match='old_verified_alias'
   if cand and cand!={p}:p=None;reason='canonical_alias_conflict'
  elif len(prev)>1:reason='ambiguous_alias'
  elif len(cand)==1:
   p=next(iter(cand));match='accepted_canonical'
   if p in pk: # old verified universe only allowed by previous alias crosswalk
    p=None;reason='unverified_old_key'
  else:reason='no_unique_key'
  if p:
   if not accepted(p):reason='nonaccepted_or_nonvascular'
   elif p in pk and pfam[p]!={x.plant_family}:reason='old_family_conflict'
   elif sourcefamilies[n]!={x.plant_family} or not x.plant_family:reason='source_family_conflict'
   elif oldgenfam.get(n.split()[0],set()) and oldgenfam[n.split()[0]]!={x.plant_family}:reason='genus_family_conflict'
   else:
    reason='mapped';mapped.append(dict(butterflyKey=k,plantKey=p,status=status[x.evidence]))
 rows.append(dict(source_row=i+2,butterflyKey=k or '',plantKey=p or '',evidence=x.evidence,match=match,disposition=reason))
new=pd.DataFrame(mapped).groupby(['butterflyKey','plantKey'],as_index=False).status.min().sort_values(['butterflyKey','plantKey']);old_df=pd.DataFrame([dict(butterflyKey=x['butterflyKey'],plantKey=x['plantKey']) for x in old]);old_df=old_df[old_df.plantKey.map(accepted)&old_df.butterflyKey.map(lambda k:bkey[k]['status']=='ACCEPTED' and bkey[k]['accepted_speciesKey']==k)].drop_duplicates()
new.to_csv(P/'filter_clarke_edges.csv',index=False);pd.DataFrame(rows).to_csv(P/'filter_mapping_audit.csv',index=False)
# Compare with archived lycaenid graph; new policy is stricter on original 990 (accepted EIV crosswalk + alias conflicts).
prev=pd.read_csv(P/'old_lycaenid_expanded.csv',dtype={'butterflyKey':str,'plantKey':str});oldly=set(zip(prev.butterflyKey,prev.plantKey));newly=set(zip(new.loc[new.butterflyKey.map(lambda k:bkey[k]['family']=='Lycaenidae'),'butterflyKey'],new.loc[new.butterflyKey.map(lambda k:bkey[k]['family']=='Lycaenidae'),'plantKey']));assert newly<=oldly,('new lycaenid keys outside archived graph',len(newly-oldly))
score=tx.set_index('speciesKey').EIV_NUTRI_exact.to_dict();out=[]
for cutoff in [3,2]:
 ns=new[new.status<=cutoff]
 for k in sorted(set(old_df.butterflyKey)|set(ns.butterflyKey)):
  rec={'cutoff':cutoff,'butterflyKey':k,'family':bkey[k]['family'],'name':bkey[k]['canonicalName']}
  for name,df in [('old',old_df),('clarke',ns)]:
   keys=sorted(set(df.loc[df.butterflyKey==k,'plantKey']));vals=np.array([score.get(p,np.nan) for p in keys]);good=np.isfinite(vals);rec.update({name+'_hosts':len(keys),name+'_scored':int(good.sum()),name+'_coverage':float(good.mean()) if len(vals) else 0.,name+'_N':float(vals[good].mean()) if good.any() else np.nan})
  rec['common_eligible']=all(rec[x+'_scored']>=2 and rec[x+'_coverage']>=.5 for x in ['old','clarke']);out.append(rec)
roster=pd.DataFrame(out);roster.to_csv(P/'filter_butterfly_scores.csv',index=False)
r=dict(status='SOURCE_MAPPING_COMPLETE',old_edges=len(old_df),clarke_edges=len(new),clarke_plants=int(new.plantKey.nunique()),clarke_butterflies=int(new.butterflyKey.nunique()),lycaenid_archived=len(oldly),lycaenid_preserved=len(newly),dispositions=pd.DataFrame(rows).disposition.value_counts().to_dict(),common_by_cutoff={str(i):int(g.common_eligible.sum()) for i,g in roster.groupby('cutoff')},eiv_sha256=sha(raw),edge_sha256=sha((P/'filter_clarke_edges.csv').read_bytes()),score_sha256=sha((P/'filter_butterfly_scores.csv').read_bytes()))
(P/'filter_mapping_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2));print(roster[roster.common_eligible].groupby(['cutoff','family']).size())
