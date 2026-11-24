#!/usr/bin/env python3
"""Source-only audit of published Appendix S1 to frozen GBIF keys. No ecological outcomes.
Exact accepted binomials and unique original verified aliases only, no fuzzy matches.
Reject infra ranks, hybrids, ambiguous source names, family mismatches, questionable-only countries.
"""
from pathlib import Path
import sys,hashlib,json,zipfile,io,re,collections
import pandas as pd
P=Path(sys.argv[1]);zpath=P/'gesui_clarke_pmc_xlsx_14916653169.zip'
sha=lambda b:hashlib.sha256(b).hexdigest()
assert sha(zpath.read_bytes())=='f50b0f9b9bd2b4d42cce5964c45b68350280868e601d0db46e1c7503176cb353'
z=zipfile.ZipFile(zpath);name='ECE3-14-e10834-s005.xlsx';b=z.read(name);assert sha(b)=='b6da07c1b62463f98761625bc2277fd1ed6e9e30ab35a833a89e0f1fc97ed043'
source=pd.read_excel(io.BytesIO(b),sheet_name='S1-data-1690376274364',dtype=str).fillna('')
assert source.shape==(5589,9) and source.accepted_butterfly.nunique()==464 and source.accepted_plant.nunique()==1972
assert set(source)=={'sp_ord','accepted_butterfly','plant_order','plant_family','accepted_plant','evidence','countries','questionable_countries','references'}
status={'Proven':1,'Larvae':2,'Eggs':3,'Claim only':4,'Breeding only':5,'Mistake':6};assert set(source.evidence)==set(status)
assert sha((P/'gesui_host_source.zip').read_bytes())=='1cf6ed4ee41282a7f1be2e5cc7e1542a186b4fef5e1a46adf5861543c0ebbd94'
h=zipfile.ZipFile(P/'gesui_host_source.zip');bmeta=json.loads(h.read('butterfly_usage_aliases.json'));pmeta=json.loads(h.read('plant_name_crosswalk.json'));old=json.loads(h.read('verified_candidate_host_edges.json'));assert len(bmeta)==303 and len(old)==1565
pk=set(x['plantKey'] for x in old);assert len(pk)==990
norm=lambda x:' '.join(str(x).casefold().split())
binomial=re.compile(r'^[A-Z][A-Za-z-]+ [a-z][a-z-]+$')
# All original GBIF names, including doubtful, enter ambiguity detection; only accepted keys may be used.
bname=collections.defaultdict(set);bcanon={};bkey={x['original_speciesKey']:x for x in bmeta}
for x in bmeta:
 k=x['original_speciesKey'];bcanon[k]=norm(x['canonicalName'])
 for n in [x['canonicalName'],*x['aliases']]:
  if re.fullmatch(r'[a-z][a-z-]+ [a-z][a-z-]+',norm(n)):bname[norm(n)].add(k)
# Clarke accepted names must not collapse two distinct source species onto one GBIF key.
bmapped={}
for s in sorted(set(source.accepted_butterfly)):
 if not binomial.fullmatch(s):continue
 cand=bname.get(norm(s),set())
 if len(cand)==1:
  k=next(iter(cand))
  if bkey[k]['status']=='ACCEPTED' and bkey[k]['accepted_speciesKey']==k:bmapped[s]=k
conflict=collections.Counter(bmapped.values());bmapped={s:k for s,k in bmapped.items() if conflict[k]==1}
# Plant source accepted names include authors. Recognize strict genus+species only; reject third epithets and rank/hybrid markers.
rank=re.compile(r'\b(?:subsp|ssp|var|subvar|forma|nothosubsp|nothovar|agg|sensu|auct|cf|aff)\b|\bs\.\s*l\.\b|\bf\.\s+[a-z][a-z-]+\b',re.I)
def plant_binomial(s):
 s=' '.join(s.split())
 if '×' in s or rank.search(s):return None
 m=re.match(r'^([A-Z][A-Za-z-]+ [a-z][a-z-]+)(?:\s+(.*))?$',s)
 if not m:return None
 name,author=m.groups()
 # A lower-case third word with no author punctuation is not a validated species name.
 if author and re.match(r'^[a-z][a-z-]+(?:\s|$)',author):return None
 return norm(name)
plant_forms=collections.defaultdict(set)
for s in sorted(set(source.accepted_plant)):
 n=plant_binomial(s)
 if n:plant_forms[n].add(' '.join(s.split()))
pname=collections.defaultdict(set);pfam=collections.defaultdict(set)
for x in pmeta:
 k=x.get('accepted_speciesKey');m=x['match']
 if k not in pk:continue
 if not x['api_success'] or m.get('rank')!='SPECIES' or m.get('matchType')!='EXACT' or m.get('confidence')!=100 or m.get('status') not in {'ACCEPTED','SYNONYM'}:continue
 for n in [x['source_name'],m.get('canonicalName','')]:
  if re.fullmatch(r'[a-z][a-z-]+ [a-z][a-z-]+',norm(n)):pname[norm(n)].add(k)
 pfam[k].add(m.get('family',''))
assert len(pfam)==990
rows=[];mapped=[]
for i,x in enumerate(source.itertuples(index=False)):
 rec={'source_excel_row':i+2,'accepted_butterfly':x.accepted_butterfly,'accepted_plant':x.accepted_plant,'plant_family':x.plant_family,'evidence':x.evidence,'status':status[x.evidence],'valid_country_present':bool(x.countries.strip()),'questionable_country_present':bool(x.questionable_countries.strip())}
 k=bmapped.get(x.accepted_butterfly);n=plant_binomial(x.accepted_plant);cand=pname.get(n,set()) if n else set();p=next(iter(cand)) if len(cand)==1 else None
 rec.update(butterflyKey=k or '',plantKey=p or '',plant_binomial=n or '',butterfly_match='canonical' if k and norm(x.accepted_butterfly)==bcanon[k] else ('verified_alias' if k else 'unmapped'))
 if status[x.evidence]>3:reason='excluded_evidence'
 elif not x.countries.strip():reason='no_nonquestionable_country'
 elif not k:reason='unmapped_or_ambiguous_butterfly'
 elif bkey[k]['family']!='Lycaenidae':reason='other_butterfly_family'
 elif not n:reason='non_species_or_ambiguous_rank'
 elif len(plant_forms[n])!=1:reason='ambiguous_source_binomial_author'
 elif not p:reason='outside_or_ambiguous_frozen_plant_universe'
 elif pfam[p]!={x.plant_family}:reason='plant_family_mismatch'
 else:
  reason='matched_reliable';mapped.append({'butterflyKey':k,'plantKey':p,'status':status[x.evidence]})
 rec['disposition']=reason;rows.append(rec)
audit=pd.DataFrame(rows);audit.to_csv(P/'gesui_clarke_appendix_mapping_audit.csv',index=False)
edge=pd.DataFrame(mapped).groupby(['butterflyKey','plantKey'],as_index=False).status.min().sort_values(['butterflyKey','plantKey']);edge.to_csv(P/'clarke_edges.csv',index=False)
receipt={'status':'PASS_SOURCE_AND_TAXONOMY','source_workflow':45386260304,'source_artifact':14916653169,'source_zip_sha256':sha(zpath.read_bytes()),'appendix_sha256':sha(b),'appendix_member':name,'source_rows':len(source),'source_butterflies':source.accepted_butterfly.nunique(),'source_plant_names':source.accepted_plant.nunique(),'evidence_counts':source.evidence.value_counts().to_dict(),'source_dispositions':audit.disposition.value_counts().to_dict(),'source_butterfly_alias_conflicts':sum(v>1 for v in conflict.values()),'common_universe':990,'mapped_rows':len(mapped),'mapped_edges':len(edge),'mapped_butterflies':edge.butterflyKey.nunique(),'mapped_plants':edge.plantKey.nunique(),'edge_csv_sha256':sha((P/'clarke_edges.csv').read_bytes()),'interpretation':'Source mapping pass only; common-support gate still required. Country field supplies at least one plausible country in Europe, not local field use. Species-level name matching does not resolve hidden taxonomic disagreements.'}
(P/'clarke_edges_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
