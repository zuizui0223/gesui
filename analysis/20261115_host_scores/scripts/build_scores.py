"""Reproduce source-locked host nutrient scores; no ecological outcome models."""
from __future__ import annotations
import csv,gzip,hashlib,json,re,sys,zipfile,statistics,math,collections
from pathlib import Path
from xml.etree import ElementTree as ET
from score_core import numeric_eiv,norm,match_names,summarize
R=Path(__file__).resolve().parents[1];O=R/'output';S=R/'source'
P=R/'previous/gesui_host_network_20261114'
EIV_DOI='https://doi.org/10.5281/zenodo.7427088'
EIV_PAPER='https://doi.org/10.1111/jvs.13168'
FERRER_DOI='https://doi.org/10.6084/m9.figshare.1168861.v1'
HOST_ROSTER='https://github.com/zuizui0223/chocho/actions/runs/40979986160'
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def writecsv(p,rows):
 if not rows:return
 with p.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def js(p):return json.loads(p.read_text())
def xlsx_rows(z,name,ss):
 for row in ET.fromstring(z.read('xl/worksheets/'+name)).findall('m:sheetData/m:row',NS):
  d={}
  for cell in row:
   label=re.sub('[0-9]','',cell.get('r',''));v=cell.find('m:v',NS)
   value=v.text if v is not None else ''
   if cell.get('t')=='s' and value:value=ss[int(value)]
   elif cell.get('t')=='inlineStr':value=''.join(cell.itertext())
   d[label]=value
  yield int(row.get('r')),d

def run():
 src=S/'EIV/EIV_Europe_Tichy_v2.xlsx'
 assert hashlib.md5(src.read_bytes()).hexdigest()=='3519d41639647dde8d80042fc98b2ef5'
 with zipfile.ZipFile(src) as z:
  ss=[''.join(si.itertext()) for si in ET.fromstring(z.read('xl/sharedStrings.xml'))]
  traits={}
  for ri,d in xlsx_rows(z,'sheet10.xml',ss):
   if ri<3 or not d.get('B'):continue
   nm=norm(d['B']);assert nm not in traits
   traits[nm]={'name':d['B'],'N':numeric_eiv(d.get('G')),'N_source_code':d.get('G',''),
      'M':numeric_eiv(d.get('E')),'L':numeric_eiv(d.get('C')),'T':numeric_eiv(d.get('D')),
      'R':numeric_eiv(d.get('F')),'source_excel_row':ri}
  for ri,d in xlsx_rows(z,'sheet8.xml',ss):
   if ri<3 or norm(d.get('B','')) not in traits:continue
   t=traits[norm(d['B'])];assert t['N']==numeric_eiv(d.get('C'))
   t.update(N_source_count=numeric_eiv(d.get('G')),N_regional_min=numeric_eiv(d.get('E')),N_regional_max=numeric_eiv(d.get('F')))
 assert len(traits)==8908
 writecsv(O/'EIV_harmonized_source_extraction.csv',list(traits.values()))
 old=js(S/'old_EIV/GBIF_SPECIESKEY_EIV_EXACT_NAME_CROSSWALK.json')['species_taxonomy_and_matched_trait']
 oldmap={str(x['gbif_speciesKey']):x for x in old}
 old_cmp=[]
 for a in old:
  if a['trait_NUTRI'] is None:continue
  t=traits.get(norm(a['gbif_canonical_name']))
  if t is not None and t['N'] is not None:
   old_cmp.append({'gbif_speciesKey':a['gbif_speciesKey'],'name':a['gbif_canonical_name'],
       'old_NUTRI':a['trait_NUTRI'],'new_EIV_v2_N':t['N'],'difference':t['N']-a['trait_NUTRI']})
 writecsv(O/'old_NUTRI_vs_Zenodo_check.csv',old_cmp)
 aliases=js(P/'source/butterfly_usage_aliases.json')
 edges=js(P/'source/verified_candidate_host_edges.json')
 cross=js(P/'source/plant_name_crosswalk.json')
 byplant=collections.defaultdict(list)
 for a in cross:
  if a.get('accepted_speciesKey'):byplant[str(a['accepted_speciesKey'])].append(a)
 plant_rows=[]; plant={}
 for pk,xs in sorted(byplant.items()):
  accepted=sorted({x['match'].get('species','') for x in xs if x['match'].get('species')})
  raw=sorted({x['source_name'] for x in xs})
  primary=match_names(accepted,traits)
  used='accepted_exact'
  if primary['status']=='unmatched':
   primary=match_names(raw,traits);used='verified_source_synonym_exact'
  t=primary['record'] or {}
  family=sorted({x['match'].get('family','') for x in xs})
  p={'plantKey':pk,'accepted_name':' | '.join(accepted),'verified_source_names':' | '.join(raw),
     'family':' | '.join(family),'match_status':primary['status'],'match_basis':used,
     'EIV_name':' | '.join(primary['names']),'N':t.get('N'),'M':t.get('M'),'L':t.get('L'),
     'N_original_code':t.get('N_source_code',''),'N_source_count':t.get('N_source_count'),
     'N_estimated_from_cooccurrence':t.get('N') is not None and t.get('N_source_count') is None,
     'N_regional_min':t.get('N_regional_min'),'N_regional_max':t.get('N_regional_max'),
     'EIV_source_row':t.get('source_excel_row'),'EIV_source':EIV_DOI}
  plant[pk]=p;plant_rows.append(p)
 writecsv(O/'gesui_host_plant_EIV_crosswalk.csv',plant_rows)
 # Recover source-named unresolved hosts, without assuming they are absent.
 aliaskeys=collections.defaultdict(set)
 for a in aliases:
  for name in a['aliases']:aliaskeys[norm(name)].add(str(a['original_speciesKey']))
 cross_byname={norm(a['source_name']):a for a in cross}
 unresolved=collections.defaultdict(set)
 with (P/'source/ferrer_paris_2014_original.csv').open() as f:
  for a in csv.DictReader(f):
   n=norm(a['Lepidoptera Name']);pn=norm(a['Hostplant Name']);p=cross_byname.get(pn)
   if len(aliaskeys[n])!=1 or p is None or p.get('accepted_speciesKey'):continue
   unresolved[next(iter(aliaskeys[n]))].add(pn)
 bedges=collections.defaultdict(list)
 for e in edges:bedges[str(e['butterflyKey'])].append(e)
 score_rows=[];link_rows=[]
 for b in aliases:
  bk=str(b['original_speciesKey']);ee=bedges[bk]
  v=[(e['plantKey'],plant[e['plantKey']]['N']) for e in ee]
  v += [('unresolved:'+n,None) for n in sorted(unresolved[bk])]
  su=summarize(v)
  exact=[(e['plantKey'],plant[e['plantKey']]['N'] if plant[e['plantKey']]['match_basis']=='accepted_exact' else None) for e in ee]
  exact += [('unresolved:'+n,None) for n in sorted(unresolved[bk])]
  nonwiki=[e for e in ee if any(x!='Wikipedia' for x in e['source_labels'])]
  sw=summarize([(e['plantKey'],plant[e['plantKey']]['N']) for e in nonwiki])
  # Missingness and source status, never filter entire taxon out of exported roster.
  state='no_documented_hosts' if not v else 'no_EIV_match' if su['mean'] is None else 'higher_coverage' if su['coverage']>=.8 else 'partial_coverage'
  rec={'butterflyKey':bk,'accepted_name':b['canonicalName'],'family':b['family'],'taxon_status':b['status'],
   'n_documented_host_concepts':su['n_hosts'],'n_resolved_host_keys':len(ee),'n_unresolved_host_names':len(unresolved[bk]),
   'n_scored_hosts':su['n_scored'],'score_sum':su['sum_N'],'trait_coverage':su['coverage'],
   'host_nutrient_mean':su['mean'],'host_nutrient_median':su['median'],'host_nutrient_sd':su['sd'],
   'host_nutrient_min':su['min'],'host_nutrient_max':su['max'],
   'documented_host_missing_lower':su['lower_bound'],'documented_host_missing_upper':su['upper_bound'],
   'nutrient_band_descriptive':su['band'],'broad_low_high_host_span':su['broad_span_low_high'],
   'quality_status':state,'strict_accepted_name_only_mean':summarize(exact)['mean'],
   'no_wikipedia_only_mean':sw['mean'],'no_wikipedia_n_scored':sw['n_scored'],
   'host_moisture_mean':statistics.mean([plant[e['plantKey']]['M'] for e in ee if plant[e['plantKey']]['M'] is not None]) if any(plant[e['plantKey']]['M'] is not None for e in ee) else None,
   'host_light_mean':statistics.mean([plant[e['plantKey']]['L'] for e in ee if plant[e['plantKey']]['L'] is not None]) if any(plant[e['plantKey']]['L'] is not None for e in ee) else None,
   'leaf_N_content':'not_measured','adult_uptake_status':'not_yet_curated','adult_evidence_source':'',
   'predicted_sewage_effect':'not_assigned','host_source':FERRER_DOI,'trait_source':EIV_DOI}
  score_rows.append(rec)
  for e in ee:
   pp=plant[e['plantKey']]
   link_rows.append({'butterflyKey':bk,'butterfly':b['canonicalName'],'family':b['family'],
      'plantKey':e['plantKey'],'plant':pp['accepted_name'],'host_source_names':' | '.join(e['plant_names']),
      'N':pp['N'],'M':pp['M'],'L':pp['L'],'EIV_name':pp['EIV_name'],'match_basis':pp['match_basis'],
      'source_labels':' | '.join(e['source_labels']),'host_source':FERRER_DOI,'trait_source':EIV_DOI})
 writecsv(O/'gesui_303_butterfly_host_nutrient_scores.csv',score_rows)
 writecsv(O/'gesui_butterfly_host_EIV_links.csv',link_rows)
 # chocho independent frozen host roster (not an independent ecological replicate).
 ch=list(csv.DictReader((S/'chocho/original230_host_edges.csv').open(encoding='utf-8-sig')))
 tx=json.loads(gzip.decompress((S/'chocho_taxonomy/original230_plant_name_to_gbifkey_v1.json.gz').read_bytes()))['name_resolution']
 chocho_plants={};clinks=[];cg=collections.defaultdict(list)
 for e in ch:
  cg[e['original_name']].append(e)
  pk=e['wcvp_accepted_id']
  if pk not in chocho_plants:
   names=json.loads(e['accepted_names_json']); mt=match_names(names,traits); used='WCVP_accepted_exact'
   g=tx.get(e['query_name'],{});gk=g.get('gbif_speciesKey')
   if mt['status']=='unmatched' and gk in plant and plant[gk]['match_status']=='exact':
    q=plant[gk];mt=match_names(q['EIV_name'].split(' | '),traits);used='source_verified_GBIF_key_alias'
   t=mt['record'] or {}
   chocho_plants[pk]={'wcvp_accepted_id':pk,'accepted_name':e['accepted_name'],'family':e['family'],
    'GBIF_speciesKey':gk,'GBIF_gate':g.get('gate','not_queried'),
    'match_status':mt['status'],'match_basis':used,'EIV_name':' | '.join(mt['names']),
    'N':t.get('N'),'M':t.get('M'),'L':t.get('L'),'trait_source':EIV_DOI}
  pp=chocho_plants[pk]
  clinks.append({'butterfly':e['original_name'],'wcvp_accepted_id':pk,'plant':pp['accepted_name'],
     'N':pp['N'],'M':pp['M'],'L':pp['L'],'match_basis':pp['match_basis'],'EIV_name':pp['EIV_name'],
     'host_source':HOST_ROSTER,'trait_source':EIV_DOI})
 crows=[]
 for name,ee in sorted(cg.items()):
  su=summarize([(e['wcvp_accepted_id'],chocho_plants[e['wcvp_accepted_id']]['N']) for e in ee])
  nameset={norm(name),norm(ee[0]['current_accepted_name'])}
  bk=set().union(*(aliaskeys[n] for n in nameset))
  if len(bk)!=1:bk=set()
  fam=next((b['family'] for b in aliases if b['original_speciesKey'] in bk),'')
  crows.append({'original_name':name,'current_accepted_name':ee[0]['current_accepted_name'],
   'current_usage_status':ee[0]['current_usage_status'],'gesui_unique_speciesKey':next(iter(bk)) if bk else '',
   'family_in_gesui':fam,'n_documented_host_concepts':su['n_hosts'],'n_scored_hosts':su['n_scored'],
   'score_sum':su['sum_N'],'trait_coverage':su['coverage'],'host_nutrient_mean':su['mean'],
   'host_nutrient_median':su['median'],'host_nutrient_sd':su['sd'],'host_nutrient_min':su['min'],
   'host_nutrient_max':su['max'],'documented_host_missing_lower':su['lower_bound'],
   'documented_host_missing_upper':su['upper_bound'],'nutrient_band_descriptive':su['band'],
   'geographic_scope':'global_host_catalogue_European_EIV_overlap_only',
   'adult_uptake_status':'not_yet_curated','predicted_sewage_effect':'not_assigned',
   'host_source':HOST_ROSTER,'trait_source':EIV_DOI})
 writecsv(O/'chocho_230_butterfly_host_nutrient_scores.csv',crows)
 writecsv(O/'chocho_host_plant_EIV_crosswalk.csv',list(chocho_plants.values()))
 writecsv(O/'chocho_butterfly_host_EIV_links.csv',clinks)
 sr={a['butterflyKey']:a for a in score_rows};comp=[]
 for c in crows:
  g=sr.get(c['gesui_unique_speciesKey'])
  if g is None:continue
  comp.append({'chocho_name':c['original_name'],'gesui_name':g['accepted_name'],'speciesKey':g['butterflyKey'],
   'chocho_mean':c['host_nutrient_mean'],'gesui_mean':g['host_nutrient_mean'],
   'chocho_coverage':c['trait_coverage'],'gesui_coverage':g['trait_coverage'],
   'difference_chocho_minus_gesui':c['host_nutrient_mean']-g['host_nutrient_mean'] if c['host_nutrient_mean'] is not None and g['host_nutrient_mean'] is not None else None})
 writecsv(O/'cross_catalogue_score_comparison.csv',comp)
 famstats=[]
 for fam in sorted({a['family'] for a in score_rows}):
  allf=[a for a in score_rows if a['family']==fam];f=[a for a in allf if a['host_nutrient_mean'] is not None]
  vals=[a['host_nutrient_mean'] for a in f];hi=[a for a in f if a['trait_coverage']>=.8 and a['taxon_status']=='ACCEPTED']
  famstats.append({'family':fam,'n_roster':len(allf),'n_any_score':len(f),'n_coverage80_accepted':len(hi),
   'mean_of_species_scores':statistics.mean(vals) if vals else None,'minimum_species_score':min(vals) if vals else None,
   'maximum_species_score':max(vals) if vals else None,'n_low':sum(a['nutrient_band_descriptive']=='low' for a in f),
   'n_intermediate':sum(a['nutrient_band_descriptive']=='intermediate' for a in f),'n_high':sum(a['nutrient_band_descriptive']=='high' for a in f),
   'n_broad_host_span':sum(a['broad_low_high_host_span'] for a in f)})
 writecsv(O/'family_score_summary.csv',famstats)
 qa={'EIV_rows':len(traits),'EIV_numeric_N':sum(t['N'] is not None for t in traits.values()),
    'EIV_md5':hashlib.md5(src.read_bytes()).hexdigest(),'EIV_SHA256':sha(src),
    'original_EIV_overlap_comparisons':len(old_cmp),'original_EIV_max_abs_difference':max(abs(x['difference']) for x in old_cmp),
    'gesui_roster':len(score_rows),'gesui_with_hosts':sum(x['n_documented_host_concepts']>0 for x in score_rows),
    'gesui_with_scores':sum(x['host_nutrient_mean'] is not None for x in score_rows),
    'gesui_coverage80_accepted':sum(x['quality_status']=='higher_coverage' and x['taxon_status']=='ACCEPTED' for x in score_rows),
    'gesui_matched_host_keys':sum(x['N'] is not None for x in plant_rows),'gesui_host_keys':len(plant_rows),
    'gesui_ambiguous_plant_matches':sum(x['match_status']=='ambiguous' for x in plant_rows),
    'gesui_source_unresolved_name_links':sum(len(v) for v in unresolved.values()),
    'chocho_roster':len(crows),'chocho_with_scores':sum(x['host_nutrient_mean'] is not None for x in crows),
    'chocho_coverage80':sum(x['trait_coverage'] is not None and x['trait_coverage']>=.8 for x in crows),
    'chocho_matched_host_concepts':sum(x['N'] is not None for x in chocho_plants.values()),
    'chocho_host_concepts':len(chocho_plants),'cross_catalogue_taxa':len(comp),
    'no_outcome_models_fit':True,'all_scores_exploratory_catalogue_proxies':True,
    'source_files_SHA256':{str(p.relative_to(R)):sha(p) for p in [P/'source/verified_candidate_host_edges.json',P/'source/plant_name_crosswalk.json',P/'source/butterfly_usage_aliases.json',S/'chocho/original230_host_edges.csv',src]}}
 assert len(score_rows)==303 and len(crows)==230 and len(edges)==1565 and len(ch)==2524
 (O/'QA.json').write_text(json.dumps(qa,indent=2))
 print(json.dumps(qa,indent=2));print('FAMILIES',famstats)
 examples=['Aglais urticae','Aglais io','Inachis io','Vanessa atalanta','Araschnia levana','Pieris rapae','Pieris brassicae','Pieris napi','Anthocharis cardamines','Lycaena dispar','Lycaena phlaeas','Lycaena tityrus','Polyommatus icarus','Plebejus argus','Phengaris arion','Maniola jurtina','Pararge aegeria','Coenonympha pamphilus','Gonepteryx rhamni','Apatura iris','Apatura ilia']
 for a in score_rows:
  if a['accepted_name'] in examples:print('EXAMPLE',a['accepted_name'],a['host_nutrient_mean'],a['n_scored_hosts'],a['n_documented_host_concepts'],a['host_nutrient_min'],a['host_nutrient_max'],a['taxon_status'])
if __name__=='__main__':run()
