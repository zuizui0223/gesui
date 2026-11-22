#!/usr/bin/env python3
"""Original Springer Guariento 93-site butterfly/plant-N source eligibility.

No effect fitting. Accepts publisher XLSX as input; outputs a compact aggregate
manifest only. Species identities are *strict* accepted binomial name matches to
pre-existing Gesui HNS score table; do not silently guess synonyms.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,re,sys,zipfile
from collections import Counter,defaultdict
from pathlib import Path
from urllib.request import Request,urlopen
import xml.etree.ElementTree as ET

PUBLISHER='https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs10531-022-02498-3/MediaObjects/10531_2022_2498_MOESM3_ESM.xlsx'
ALTERNATE='https://static-content.springer.com/esm/art%3A10.1007%2Fs10531-022-02498-3/MediaObjects/10531_2022_2498_MOESM3_ESM.xlsx'
EXPECTED_XLSX_SHA256='e0774f313c5e008e7d37cef7d20d9f822f1cc02d7d0ba9f0b3301d5a9cf9202c'
EXPECTED_HNS_SHA256='00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb'
M='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
R='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'


def normal(s):
    return ' '.join(str(s).replace('\xa0',' ').split()).strip()


def colnum(ref):
    match=re.match(r'([A-Z]+)',ref or '')
    if not match:raise ValueError('BAD_CELL_REFERENCE')
    n=0
    for ch in match.group(1):n=26*n+ord(ch)-64
    return n


def xlsx_sheets(blob):
    if len(blob)>12_000_000 or blob[:2]!=b'PK':raise ValueError('NOT_VALID_BOUNDED_XLSX')
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        entries={i.filename for i in z.infolist()}
        if len(entries)>1000 or sum(i.file_size for i in z.infolist())>80_000_000:
            raise ValueError('ZIP_SAFETY_LIMIT')
        shared=[]
        if 'xl/sharedStrings.xml' in entries:
            shared=[''.join(x.text or '' for x in si.iter(M+'t'))
                    for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall(M+'si')]
        rel=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        targets={x.get('Id'):'xl/'+x.get('Target','').lstrip('/').removeprefix('xl/') for x in rel}
        workbook=ET.fromstring(z.read('xl/workbook.xml'))
        out={}
        for sh in workbook.findall('.//'+M+'sheet'):
            p=targets[sh.get(R+'id')]
            if p not in entries or not p.startswith('xl/') or '/../' in p:raise ValueError('BAD_SHEET_TARGET')
            root=ET.fromstring(z.read(p))
            matrix=[]
            for row in root.findall('.//'+M+'sheetData/'+M+'row'):
                line={}
                for c in row.findall(M+'c'):
                    v=c.find(M+'v');value='' if v is None or v.text is None else v.text
                    if c.get('t')=='s' and value:value=shared[int(value)]
                    elif c.get('t')=='inlineStr':value=''.join(t.text or '' for t in c.iter(M+'t'))
                    line[colnum(c.get('r'))]=value
                matrix.append(line)
            out[sh.get('name')]=matrix
        return out


def table(rows):
    if not rows:raise ValueError('EMPTY_SHEET')
    head={i:normal(s) for i,s in rows[0].items() if normal(s)}
    return head,[{head[i]:r.get(i,'') for i in head} for r in rows[1:]]


def fetch_xlsx():
    errors=[]
    for u in [PUBLISHER,ALTERNATE]:
        try:
            req=Request(u,headers={'User-Agent':'Mozilla/5.0 (compatible; ResearchSourceAudit/1.1)'})
            with urlopen(req,timeout=35) as r:blob=r.read(12_000_001)
            if hashlib.sha256(blob).hexdigest()!=EXPECTED_XLSX_SHA256:raise ValueError('PUBLISHER_SOURCE_SHA_CHANGED')
            return blob,u,errors
        except Exception as e:errors.append({'url':u,'error':str(e)[:180]})
    raise RuntimeError('OFFICIAL_XLSX_UNAVAILABLE: '+str(errors))


def count_source(blob, hns_data):
    d=xlsx_sheets(blob)
    if not {'Community','Variables','Traits','Species list'}.issubset(d):raise ValueError('MISSING_PUBLISHED_WORKBOOK_SHEETS')
    community_head, community=table(d['Community']); variables_head,variables=table(d['Variables'])
    hs=list(csv.DictReader(io.StringIO(hns_data.decode('utf-8-sig'))))
    if len(hs)!=303:raise ValueError('ORIGINAL_HNS_ROSTER_NOT_303')
    usable={normal(r['accepted_name']):r for r in hs if r['taxon_status']=='ACCEPTED'
        and r['actual_used_alias_conflict'].lower() not in ('1','true','yes')
        and r['host_nutrient_mean'] and math.isfinite(float(r['host_nutrient_mean']))}
    if len(usable)!=186:raise ValueError('ORIGINAL_HNS_ROSTER_NOT_186')
    if len(community)!=93 or len(variables)!=93:raise ValueError('EXPECTED_93_SITE_ROWS')
    if 'site_code' not in community_head.values() or 'site_code' not in variables_head.values():
        raise ValueError('NO_NATIVE_SITE_KEY')
    env={normal(r['site_code']):r for r in variables}
    if len(env)!=93:raise ValueError('DUPLICATE_NATIVE_VARIABLE_SITE')
    if {normal(r['site_code']) for r in community}!=set(env):raise ValueError('SITE_JOIN_MISMATCH')
    nutrient='Nutrients (Elenberg plant indicator)'
    light='Light (Elenberg plant indicator)'
    if nutrient not in variables_head.values() or light not in variables_head.values():
        raise ValueError('NO_NATIVE_NUTRI_OR_LIGHT')
    site_N={k:float(r[nutrient]) for k,r in env.items() if r.get(nutrient) not in ('','NA')}
    if not all(math.isfinite(v) for v in site_N.values()):raise ValueError('NONFINITE_NUTRI')
    all_species=[x for x in community_head.values() if x not in ('site_code','Land use type')]
    if len(all_species)!=100:raise ValueError(f'EXPECTED_100_COMMUNITY_COLUMNS_NOT_{len(all_species)}')
    if len(set(map(normal,all_species)))!=100:raise ValueError('SPECIES_NAME_NONUNIQUE_AFTER_WS_NORMALIZE')
    counts={s:defaultdict(float) for s in all_species}
    raw_total=0.;raw_noninteger=0;nonfinite=0
    for row in community:
        site=normal(row['site_code'])
        for sp in all_species:
            try:v=float(row.get(sp,'0') or '0')
            except ValueError:raise ValueError('NONNUMERIC_SPECIES_CELL')
            if v<0 or not math.isfinite(v):nonfinite+=1
            if abs(v-round(v))>1e-9:raw_noninteger+=1
            counts[sp][site]+=v
            raw_total+=v
    if nonfinite:raise ValueError('BAD_COMMUNITY_ABUNDANCE')
    variants=[]; eligible={}; matched=[]; ambiguous=[]
    for sp in all_species:
        key=normal(sp)
        if '/' in key or '\\' in key or ' cf.' in key or not re.fullmatch(r'[A-Z][a-zA-Z-]+ [a-z][a-zA-Z-]+',key):
            ambiguous.append(key);continue
        if key not in usable:continue
        trait=usable[key]
        eligible[key]={'species_key':str(trait['butterflyKey']),
            'genus':key.split()[0],'hns':float(trait['host_nutrient_mean']),
            'host_coverage':float(trait['trait_coverage']),'number_of_recorded_sites':sum(v>0 for v in counts[sp].values()),
            'total_recorded_count':sum(counts[sp].values()),'original_column':sp}
        matched.append(key)
    grouped=defaultdict(list)
    for key,r in eligible.items():grouped[r['genus']].append((key,r))
    for genus, species in grouped.items():
        candidates=[(name,r) for name,r in species if r['number_of_recorded_sites']>=2]
        if len(candidates)<2 or len(set(round(r['hns'],9) for _,r in candidates))<2:continue
        informative_sites=0; total_genus_site_counts=0
        for site in env:
            detected=sum(counts[r['original_column']][site]>0 for _,r in candidates)
            informative_sites+=detected>=2
            total_genus_site_counts+=detected>=1
        variants.append({'genus':genus,'species':len(candidates),'n_sites_any_member':total_genus_site_counts,
          'n_sites_2plus_congeners_recorded':informative_sites,
          'strict_host_coverage80_species':sum(r['host_coverage']>=.8 for _,r in candidates),
          'candidate_names':[name for name,_ in candidates]})
    variants.sort(key=lambda x:x['genus'])
    landuse=Counter(normal(r.get('Land use type','')) for r in variables)
    scored_sites=len(site_N)
    result={'status':'SOURCE_SUPPORT_AUDIT_NO_EFFECT_FIT','original_study':'Guariento et al. 2023 10.1007/s10531-022-02498-3',
      'source_sha256':hashlib.sha256(blob).hexdigest(),'HNS_sha256':hashlib.sha256(hns_data).hexdigest(),
      'sites_in_butterfly_matrix':len(community),'sites_in_variables':len(variables),
      'site_join_exact':True,'plant_ellenberg_N_nonmissing_sites':scored_sites,
      'plant_ellenberg_N_range':[min(site_N.values()),max(site_N.values())] if site_N else None,
      'nutrient_distinct_values':len(set(site_N.values())),
      'raw_butterfly_species_column_count':len(all_species),'butterfly_abundance_sum':raw_total,
      'fractional_species_count_cells':raw_noninteger,'land_use_sites':dict(landuse),
      'ambiguous_syntactic_name_columns':ambiguous,
      'accepted_HNS_numeric_exact_species':len(eligible),'eligible_species_total_records':sum(r['total_recorded_count'] for r in eligible.values()),
      'genus_HNS_variation_eligible_count':len(variants),'genus_support':variants,
      'genus_HNS_coverage80_eligible_count':sum(sum(r['host_coverage']>=.8 for _,r in grouped[v['genus']])>=2 for v in variants),
      'source_only_no_coefficient_fitted':True}
    result['gate']={'at_least_10_genera':len(variants)>=10,
      'at_least_40_sites_with_nutrient':scored_sites>=40,
      'at_least_10_distinct_N_values':len(set(site_N.values()))>=10,
      'no_source_taxa_name_imputation':True,
      'eligible_to_proceed_exploratory_model':len(variants)>=10 and scored_sites>=40 and len(set(site_N.values()))>=10}
    return result


def main():
    a=argparse.ArgumentParser();a.add_argument('--hns',type=Path,required=True)
    a.add_argument('--output',type=Path,required=True)
    x=a.parse_args();x.output.parent.mkdir(parents=True,exist_ok=True)
    try:
        blob,url,errors=fetch_xlsx()
        hns=x.hns.read_bytes()
        if hashlib.sha256(hns).hexdigest()!=EXPECTED_HNS_SHA256:raise ValueError('HNS_SOURCE_SHA_DRIFT')
        result=count_source(blob,hns);result['official_source_url']=url;result['download_errors']=errors
    except Exception as exc:
        result={'status':'SOURCE_STOP','error':str(exc)[:600],
           'source_only_no_coefficient_fitted':True}
    x.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'status':result['status'],'gate':result.get('gate'),'error':result.get('error'),
       'matched':result.get('accepted_HNS_numeric_exact_species'),
       'genera':result.get('genus_HNS_variation_eligible_count')},ensure_ascii=False))
    if result['status']=='SOURCE_STOP':sys.exit(3)

if __name__=='__main__':main()
