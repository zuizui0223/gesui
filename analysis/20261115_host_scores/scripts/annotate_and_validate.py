"""Evidence annotations and independent checks. Does not change host scores."""
from pathlib import Path
import csv,json,collections,math,statistics,hashlib
from score_core import ambiguous_used_keys
R=Path(__file__).resolve().parents[1];O=R/'output';P=R/'previous/gesui_host_network_20261114';S=R/'source'
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig')))
def write(p,rows):
 with p.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def n(v):return float(v) if v not in ('',None) else None
A=json.loads((P/'source/butterfly_usage_aliases.json').read_text());byalias=collections.defaultdict(set)
for a in A:
 for nm in a['aliases']:byalias[nm].add(a['original_speciesKey'])
ambiguous={k for keys in byalias.values() if len(keys)>1 for k in keys}
original_edges=json.loads((P/'source/verified_candidate_host_edges.json').read_text())
used_conflicts=ambiguous_used_keys(original_edges)
family_by_key={a['original_speciesKey']:a['family'] for a in A}
assert all(e['butterfly_family']==family_by_key[e['butterflyKey']] for e in original_edges)
E=[]
def add(id,species,stage,host,resource,design,region,endpoint,conclusion,url,restriction):
 E.append(dict(evidence_id=id,species_as_reported=species,life_stage=stage,host_plant=host,resource=resource,
   design=design,study_region=region,endpoint=endpoint,conclusion=conclusion,source_url=url,transfer_limitation=restriction))
nettle='https://www.sciencedirect.com/science/article/pii/S1146609X17303259'
for i,sp in enumerate(['Aglais urticae','Aglais io']):
 add('L'+str(i+1),sp,'larva','Urtica dioica','host plant N enrichment','rearing experiment','study-specific',
     'survival, development, pupal mass','Mostly beneficial under tested treatments',nettle,'Not wastewater; not all nutrients, doses or field conditions; not adult feeding')
for i,sp in enumerate(['Coenonympha pamphilus','Lycaena phlaeas','Lycaena tityrus','Pararge aegeria']):
 host='Rumex acetosella' if 'Lycaena' in sp else 'Poa pratensis'
 add('L'+str(i+3),sp,'larva',host,'host plant N fertilization','split-brood experiment','study-specific',
     'larval survival','Lower survival under tested fertilization treatments','https://pubmed.ncbi.nlm.nih.gov/30288608/',
     'Not wastewater; no universal whole-species or family response label')
add('L7','Araschnia levana','larva','Urtica dioica','N/P by temperature','rearing experiment','European laboratory',
 'growth, survival, adult heat tolerance','Endpoint- and treatment-specific responses','https://doi.org/10.1111/1365-2435.70347',
 'Not wastewater; beneficial growth does not establish field population increase')
add('A1','Pieris napi (as reported, 2004)','adult','','sodium vs nitrogenous substrates','field substrate-choice experiment','Colorado, USA',
 'feeding/visitation preference','Sodium-based substrate preferences','https://doi.org/10.1603/0046-225X-33.4.1020',
 'North American legacy taxon concept; do not automatically transfer to European Pieris napi; not sewage')
for i,sp in enumerate(['Papilio machaon','Papilio xuthus','Papilio macilentus','Papilio bianor','Papilio maackii','Papilio protenor','Papilio helenus','Papilio memnon']):
 add('A'+str(i+2),sp,'adult','','ions in puddling fluid','experimental fluid/excretion analysis','Japan',
     'Na/K exchange','Puddling physiology documented','https://doi.org/10.1371/journal.pone.0126632',
     'Species physiology evidence in Japan, not a measured sewage preference or European field frequency')
write(O/'mechanism_evidence_ledger.csv',E)
plants={p['plantKey']:p for p in read(O/'gesui_host_plant_EIV_crosswalk.csv')}
plantc={p['wcvp_accepted_id']:p for p in read(O/'chocho_host_plant_EIV_crosswalk.csv')}
links=read(O/'gesui_butterfly_host_EIV_links.csv');ll=collections.defaultdict(list)
for e in links:ll[e['butterflyKey']].append(e)
rows=read(O/'gesui_303_butterfly_host_nutrient_scores.csv')
source_alias={a['original_speciesKey']:a['aliases'] for a in A}
larv=collections.defaultdict(list)
for e in E:
 if e['life_stage']=='larva':larv[e['species_as_reported'].lower()].append(e)
adult={e['species_as_reported']:e for e in E if e['life_stage']=='adult' and e['evidence_id']!='A1'}
valid=[]
for x in rows:
 k=x['butterflyKey']; hosts=ll[k]
 vals=[float(e['N']) for e in hosts if e['N']!='']
 assert len(hosts)==int(x['n_resolved_host_keys'])
 assert len(vals)==int(x['n_scored_hosts'])
 if vals:
  assert abs(math.fsum(vals)/len(vals)-float(x['host_nutrient_mean']))<1e-12
  assert min(vals)==float(x['host_nutrient_min']) and max(vals)==float(x['host_nutrient_max'])
  valid.append(x)
 else:assert x['host_nutrient_mean']==''
 x['n_host_EIV_indifferent_x']=sum(plants[e['plantKey']]['N_original_code'].lower()=='x' for e in hosts)
 x['n_host_EIV_derived_estimates']=sum(plants[e['plantKey']]['N_estimated_from_cooccurrence']=='True' for e in hosts)
 x['shared_source_alias_review']=k in ambiguous
 x['actual_used_alias_conflict']=k in used_conflicts
 x['taxonomic_review_required']=(x['taxon_status']!='ACCEPTED' or k in used_conflicts)
 x['alias_flag_note']='shared_source_alias_review is informational; only actual used-name conflicts or non-ACCEPTED status enter review flag'
 x['source_aliases']=' | '.join(source_alias[k])
 le=[e for nm in source_alias[k] for e in larv[nm]]
 x['larval_enrichment_evidence_ids']=';'.join(sorted({e['evidence_id'] for e in le}))
 x['larval_enrichment_evidence_summary']=' | '.join(sorted({e['conclusion'] for e in le}))
 if x['accepted_name'] in adult:
  x['adult_uptake_status']='documented_other_region_Japan_no_sewage_test';x['adult_evidence_source']=adult[x['accepted_name']]['source_url']
 elif x['accepted_name']=='Pieris napi':
  x['adult_uptake_status']='legacy_North_American_concept_needs_review';x['adult_evidence_source']='https://doi.org/10.1603/0046-225X-33.4.1020'
write(O/'gesui_303_butterfly_host_nutrient_scores.csv',rows)
crows=read(O/'chocho_230_butterfly_host_nutrient_scores.csv');clinks=read(O/'chocho_butterfly_host_EIV_links.csv');cby=collections.defaultdict(list)
for e in clinks:cby[e['butterfly']].append(e)
for x in crows:
 vals=[float(e['N']) for e in cby[x['original_name']] if e['N']]
 assert len(vals)==int(x['n_scored_hosts'])
 if vals:assert abs(statistics.fmean(vals)-float(x['host_nutrient_mean']))<1e-12
 for nm in [x['original_name'],x['current_accepted_name']]:
  if nm in adult:x['adult_uptake_status']='documented_other_region_Japan_no_sewage_test'
 x['adult_evidence_source']=adult.get(x['original_name'],adult.get(x['current_accepted_name'],{})).get('source_url','')
write(O/'chocho_230_butterfly_host_nutrient_scores.csv',crows)
# Source licenses/definitions and extraction warnings travel with all tables.
sources=[
 {'id':'EIV','source':'Tichy et al. 2023, harmonized Europe EIV version 2.0','url':'https://doi.org/10.5281/zenodo.7427088','licence_or_access':'CC BY 4.0','meaning':'Habitat nutrient-affinity N; NOT foliar N. x=indifferent, NA=missing. Primary table sheet10 cached values.'},
 {'id':'FERRER','source':'Ferrer-Paris 2014 host source; original gesui verified keys','url':'https://doi.org/10.6084/m9.figshare.1168861.v1','licence_or_access':'Use original project source attribution; unmerged results','meaning':'Global catalogue, not local feeding; source overlap with HOSTS means not independent from chocho.'},
 {'id':'CHOCHO','source':'Frozen original230 WCVP host roster','url':'https://github.com/zuizui0223/chocho/actions/runs/40979986160','licence_or_access':'User-connected source; do not publish raw roster to repository','meaning':'Global host links, European indicator overlap only. Names preserve frozen taxonomy.'},
 {'id':'OLD_NUTRI','source':'gesui original 1006-key EIV crosswalk','url':'https://github.com/zuizui0223/gesui/actions/runs/42378520227','licence_or_access':'User-connected source','meaning':'624 matched numeric scores exactly reproduce Zenodo v2 values.'},
 {'id':'CNI_2017','source':'WallisDeVries and van Swaay 2017','url':'https://doi.org/10.1016/j.biocon.2016.11.029','licence_or_access':'Primary literature','meaning':'Published butterfly CNI uses species habitat nitrogen optima, NOT the host-list mean built here.'},
 {'id':'TISSUE_N','source':'Pellissier et al. 2012','url':'https://doi.org/10.1111/j.1420-9101.2012.02555.x','licence_or_access':'Primary literature','meaning':'Host leaf N / ant association is a different axis from habitat EIV N.'}
]
write(O/'sources.csv',sources)
comp=[x for x in read(O/'cross_catalogue_score_comparison.csv') if x['difference_chocho_minus_gesui']]
# Descriptive variance decomposition only, not phylogeny-corrected significance.
vv=[float(a['host_nutrient_mean']) for a in valid];grand=statistics.fmean(vv)
byfam=collections.defaultdict(list)
for a in valid:byfam[a['family']].append(float(a['host_nutrient_mean']))
ssb=sum(len(v)*(statistics.fmean(v)-grand)**2 for v in byfam.values());sst=sum((v-grand)**2 for v in vv)
quality={'all_303_rows_retained':len(rows)==303,'all_230_rows_retained':len(crows)==230,
 'independent_score_recalculation_passed':True,'documented_no_effect_groups_created':False,
 'scored_gesui_keys':len(valid),'scored_chocho_rows':sum(bool(x['host_nutrient_mean']) for x in crows),
 'gesui_any_indifferent_hosts':sum(int(x['n_host_EIV_indifferent_x'])>0 for x in rows),
 'scored_keys_shared_alias_review':sum(x['shared_source_alias_review'] for x in valid),
 'used_alias_conflicting_keys':len(used_conflicts),'all_actual_host_edges_family_consistent':True,
 'gesui_taxonomic_review_required':sum(x['taxonomic_review_required'] for x in rows),
 'gesui_higher_coverage_and_unflagged':sum(n(x['trait_coverage']) is not None and n(x['trait_coverage'])>=.8 and not x['taxonomic_review_required'] for x in rows),
 'gesui_band_counts_all_scored':dict(collections.Counter(x['nutrient_band_descriptive'] for x in valid)),
 'family_between_fraction_of_score_variance_descriptive':ssb/sst,
 'cross_catalogue_pairs_both_scored':len(comp),
 'max_abs_catalogue_score_difference':max(abs(float(x['difference_chocho_minus_gesui'])) for x in comp),
 'adult_species_other_region_evidence_gesui':sum(x['adult_uptake_status']=='documented_other_region_Japan_no_sewage_test' for x in rows),
 'adult_legacy_concept_pending_gesui':sum(x['adult_uptake_status']=='legacy_North_American_concept_needs_review' for x in rows),
 'adult_not_curated_gesui':sum(x['adult_uptake_status']=='not_yet_curated' for x in rows),
 'publisher_summary_N_count':7193,'source_file_numeric_N_count':7195,
 'count_note':'Pinned v2 file contains 7195 numeric N rows; paper summary says 7193. Source metadata discrepancy retained; not silently changed. All 624 prior numeric crosswalk overlaps agree exactly.'}
(O/'independent_validation.json').write_text(json.dumps(quality,indent=2))
print(json.dumps(quality,indent=2))
