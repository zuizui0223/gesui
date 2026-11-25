#!/usr/bin/env python3
"""Expand only outside original 990 host keys using pinned GBIF all-flora taxonomy."""
from pathlib import Path
import sys, hashlib
P=Path(sys.argv[1]); prior=P/'gesui_clarke_appendix_mapping.py'
exec(compile(prior.read_text().split('rows=[];mapped=[]')[0],str(prior),'exec'))
# Preserve prior reliable old edges (hash-locked, verified from original source).
old_receipt=json.loads((P/'clarke_edges_receipt.json').read_text());assert sha((P/'clarke_edges.csv').read_bytes())==old_receipt['edge_csv_sha256']
base=pd.read_csv(P/'clarke_edges.csv',dtype={'butterflyKey':str,'plantKey':str});assert len(base)==256 and set(base.plantKey)<=pk
zf=P/'gesui_eiv_source.zip';assert sha(zf.read_bytes())=='1d3051d4c527421bec3ec3ba917422a9162afb81cc7f38db0d8130287174a14c'
with zipfile.ZipFile(zf) as ez:
 raw=ez.read('EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv');assert sha(raw)=='ab45f36c0dad047442bb958bcef3135b45014bafec8a70811bcafb6109a02679'
 tx=pd.read_csv(io.BytesIO(raw),dtype=str).fillna('')
assert len(tx)==11511 and tx.speciesKey.nunique()==11511
# Original 454 occurrence rows determine available keys; all 2808 facets and source bytes locked.
gz=P/'gesui_gbif_facets.zip';assert sha(gz.read_bytes())=='4b37c741d1bd28c178b0cef61c9457af0dd30797c173f9384ebe1d1d0036d07c'
pz=P/'gesui_original_panel.zip';assert sha(pz.read_bytes())=='b0babfa63a96f9b7eb5a83122b08b8664edab6e8c6c28d2faf4ff024ca84c39b'
with zipfile.ZipFile(pz) as z: cohort=pd.read_csv(z.open('EU234_H7_ORIGINAL_STUDY_NATIVE_BASINPERIOD_EQR_PLANT_BUTTERFLY_SOURCE_COHORT.csv'),dtype=str)
assert len(cohort)==454
units=set(zip(cohort.HYBAS_ID_L07,cohort.PERIOD))
with zipfile.ZipFile(gz) as z:facets=json.loads(z.read('GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json'))
counts=collections.Counter()
for x in facets:
 if x['group']=='Plantae' and (x['HYBAS_ID_L07'],x['period']) in units:counts.update(x['source_speciesKey_counts'])
assert set(counts)==set(tx.speciesKey),(len(counts),len(tx))
namekeys=collections.defaultdict(set);meta=tx.set_index('speciesKey').to_dict('index');oldgenfam=collections.defaultdict(set)
for r in tx.itertuples(index=False):
 n=norm(r.gbif_canonical_name)
 if re.fullmatch(r'[a-z][a-z-]+ [a-z][a-z-]+',n):namekeys[n].add(r.speciesKey)
for x in pmeta:
 m=x['match'];k=x.get('accepted_speciesKey');n=norm(m.get('canonicalName',''));fam=m.get('family','')
 if k and fam and m.get('rank')=='SPECIES' and m.get('matchType')=='EXACT' and m.get('confidence')==100 and re.fullmatch(r'[a-z][a-z-]+ [a-z][a-z-]+',n):oldgenfam[n.split()[0]].add(fam)
sourcefamilies=collections.defaultdict(set)
for x in source.itertuples(index=False):
 n=plant_binomial(x.accepted_plant)
 if n and x.plant_family:sourcefamilies[n].add(x.plant_family)
new=[];audit=[];newfamilies={};oldonly=set(zip(base.butterflyKey,base.plantKey))
for i,x in enumerate(source.itertuples(index=False)):
 n=plant_binomial(x.accepted_plant);bk0=bmapped.get(x.accepted_butterfly);cand=namekeys.get(n,set()) if n else set();k=next(iter(cand)) if len(cand)==1 else None
 rec=dict(source_excel_row=i+2,accepted_butterfly=x.accepted_butterfly,accepted_plant=x.accepted_plant,plant_family=x.plant_family,evidence=x.evidence,butterflyKey=bk0 or '',candidate_plantKey=k or '',candidate_binomial=n or '')
 if status[x.evidence]>3:reason='excluded_evidence'
 elif not x.countries.strip():reason='no_nonquestionable_country'
 elif not bk0:reason='unmapped_butterfly'
 elif bkey[bk0]['family']!='Lycaenidae':reason='other_butterfly_family'
 elif not n or len(plant_forms[n])!=1:reason='rank_or_ambiguous_source_binomial'
 elif not k:reason='no_unique_canonical_GBIF_key'
 elif k in pk:reason='original_990_not_changed'
 elif pname.get(n,set()) and pname[n]!={k}:reason='conflicting_previous_alias'
 elif not (meta[k]['GBIF_source_status']=='OK' and meta[k]['gbif_rank']=='SPECIES' and meta[k]['gbif_taxon_status']=='ACCEPTED' and meta[k]['gbif_phylum']=='Tracheophyta' and meta[k]['gbif_kingdom']=='Plantae'):reason='nonaccepted_or_nonvascular'
 elif counts[k]<=0:reason='no_observed_count'
 elif sourcefamilies[n]!={x.plant_family} or not x.plant_family:reason='ambiguous_family'
 elif oldgenfam.get(n.split()[0],set()) and oldgenfam[n.split()[0]]!={x.plant_family}:reason='conflicting_old_genus_family'
 else:
  reason='matched_new_outside990';new.append(dict(butterflyKey=bk0,plantKey=k,status=status[x.evidence]));newfamilies[k]=x.plant_family
 rec['disposition']=reason;audit.append(rec)
new=pd.DataFrame(new).groupby(['butterflyKey','plantKey'],as_index=False).status.min().sort_values(['butterflyKey','plantKey']);assert not set(new.plantKey)&pk
full=pd.concat([base,new],ignore_index=True).sort_values(['butterflyKey','plantKey']);assert not full.duplicated(['butterflyKey','plantKey']).any()
full.to_csv(P/'clarke_expanded_edges.csv',index=False);pd.DataFrame(audit).to_csv(P/'clarke_expanded_mapping_audit.csv',index=False);pd.DataFrame([dict(plantKey=k,canonical=meta[k]['gbif_canonical_name'],family=f,total_record_count=counts[k],family_source='Clarke; old genus cross-check when available') for k,f in sorted(newfamilies.items())]).to_csv(P/'clarke_expanded_plant_keys.csv',index=False)
r=dict(status='PASS_SOURCE_MAPPING',source_taxonomy_keys=len(tx),original_edges=len(base),original_plants=base.plantKey.nunique(),new_edges=len(new),new_plants=len(newfamilies),expanded_edges=len(full),expanded_plants=full.plantKey.nunique(),expanded_butterflies=full.butterflyKey.nunique(),source_dispositions=pd.DataFrame(audit).disposition.value_counts().to_dict(),source_sha256=sha(raw),original_edge_sha256=sha((P/'clarke_edges.csv').read_bytes()),expanded_edge_sha256=sha((P/'clarke_expanded_edges.csv').read_bytes()),plant_keys_sha256=sha((P/'clarke_expanded_plant_keys.csv').read_bytes()),notes='Original 990 preserved; newly matched exact accepted canonical keys outside 990; new family from Clarke with prior genus check if available. Not complete taxonomic reconciliation.')
(P/'clarke_expanded_mapping_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
