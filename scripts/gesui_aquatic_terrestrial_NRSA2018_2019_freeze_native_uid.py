#!/usr/bin/env python3
"""Reconcile authentic EPA NRSA 2018-19 source at native UID grain BEFORE any
nutrient–taxon effect. Execute in original source download runner after first
strict site+visit gate failure and duplicate UID diagnostics. No ecological fitting.
"""
import collections,csv,hashlib,json,math,os,pathlib
p=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_nrsa1819_source_gate"
report=json.loads((p/"NRSA_1819_SOURCE_ONLY_QA.json").read_text())
diagnosis=json.loads((p/"NRSA_1819_DUPLICATE_UID_PREOUTCOME_AUDIT.json").read_text())
EXPECTED={
 'benthos':'4e9b3a0798b25f67c7a96d4800b94f86544bbd623abaf119358e067b158324a7',
 'sites':'1b5b2119c881290191edeb47f1d80a7711f10440e9851574cfa86fd57be770e3',
 'water':'cb9a3349ae3708834eab533e71423d74c1778bf76aba0737f0ae1bcc8f50baa1'
}
for nm,sha in EXPECTED.items():
    v=(p/(nm+'.csv')).read_bytes()
    if hashlib.sha256(v).hexdigest()!=sha:raise RuntimeError("Original EPA SOURCE hash differs: "+nm)

def source_rows(nm):
    with (p/(nm+'.csv')).open(encoding="utf-8-sig",newline="") as f:
        yield from csv.DictReader(f)

def s(r,k):
    return (r.get(k) or "").strip()
def numpos(v):
    try:return math.isfinite(float(v)) and float(v)>0
    except (ValueError,TypeError):return False

site_idx=collections.defaultdict(list);water_idx=collections.defaultdict(list)
for r in source_rows('sites'):
    if s(r,'UID'):site_idx[s(r,'UID')].append(r)
for r in source_rows('water'):
    if s(r,'UID'):water_idx[s(r,'UID')].append(r)
dups={'site_uid_duplicates':sum(len(v)>1 for v in site_idx.values()),
      'water_uid_duplicates':sum(len(v)>1 for v in water_idx.values())}
# Exclude UID if either original table is duplicated; no post-outcome selection.
unique={uid for uid in site_idx.keys()&water_idx.keys()
        if len(site_idx[uid])==1 and len(water_idx[uid])==1
        and s(site_idx[uid][0],'SITE_ID')==s(water_idx[uid][0],'SITE_ID')
        and s(site_idx[uid][0],'VISIT_NO')==s(water_idx[uid][0],'VISIT_NO')}
benthic_positive_mass=collections.Counter()
benthic_ffg_mass=collections.Counter()
benthic_tax_mass=collections.Counter()
benthic_rows=collections.Counter()
benthic_composite=collections.defaultdict(set)
for r in source_rows('benthos'):
    if s(r,'SAMPLE_TYPE').upper()!='BERW':continue
    uid=s(r,'UID')
    if not uid or uid not in unique:continue
    try:v=float(s(r,'TOTAL300'))
    except ValueError:continue
    if not math.isfinite(v) or v<=0:continue
    if s(r,'SITE_ID')!=s(site_idx[uid][0],'SITE_ID'):continue
    benthic_rows[uid]+=1
    benthic_positive_mass[uid]+=v
    if s(r,'FFG'):benthic_ffg_mass[uid]+=v
    if s(r,'ORDER') or s(r,'FAMILY'):benthic_tax_mass[uid]+=v
    benthic_composite[(s(r,'SITE_ID'),s(r,'VISIT_NO'))].add(uid)

# Pre-effect availability only; fail rather than substitute nitrate for total N.
manifest=[]
chem_units=collections.Counter()
for uid in sorted(benthic_positive_mass):
    a=site_idx[uid][0];b=water_idx[uid][0]
    if not numpos(s(b,'NTL_RESULT')) or not numpos(s(b,'PTL_RESULT')):continue
    if not s(a,'PSTL_CODE') or not s(a,'AG_ECO9'):continue
    if s(a,'PSTL_CODE') != s(b,'STATE'):continue
    try:
        lat=float(s(a,'LAT_DD83'));lon=float(s(a,'LON_DD83'))
        if not math.isfinite(lat) or not math.isfinite(lon) or not 24<=lat<=50 or not -125<=lon<=-66:continue
    except ValueError:continue
    tn_unit=s(b,'NTL_UNITS');tp_unit=s(b,'PTL_UNITS')
    chem_units[(tn_unit,tp_unit)]+=1
    if tn_unit!='mg N/L' or tp_unit!='ug/L':continue
    manifest.append({
       'uid':uid,'site_id':s(a,'SITE_ID'),'site_unique_id':s(a,'UNIQUE_ID'),
       'visit_no':s(a,'VISIT_NO'),'year':s(a,'YEAR'),'date_col':s(a,'DATE_COL'),
       'state':s(a,'PSTL_CODE'),'eco9':s(a,'AG_ECO9'),
       'huc8':s(a,'HUC8'),'lat':lat,'lon':lon,
       'berw_taxon_positive_rows':benthic_rows[uid],
       # source support coverage, not response composition.
       'berw_total300_individuals':benthic_positive_mass[uid],
    })
selected={v['uid'] for v in manifest}
group=collections.defaultdict(list)
for z in manifest:group[z['site_id']].append(z)
site_date_dups=sum(len(set((r['year'],r['date_col']) for r in v))<len(v) for v in group.values())
mass=sum(benthic_positive_mass[i] for i in selected)
ffg=sum(benthic_ffg_mass[i] for i in selected)
tax=sum(benthic_tax_mass[i] for i in selected)
eco={r['eco9'] for r in manifest};states={r['state'] for r in manifest}
gate=bool(
   len(manifest)>=300 and len(states)>=25 and len(eco)>=6
   and mass>0 and ffg/mass>=.6 and tax/mass>=.8
   and not dups['site_uid_duplicates'] and not dups['water_uid_duplicates']
   and all(len(set(z['uid'] for z in values))==len(values) for values in group.values())
)
out={
 'original_source_only_composite_gate_remains':report['decision'],
 'source_fix_provenance':'After source-only audit identified distinct original UID per duplicate SITE_ID+VISIT_NO, and BEFORE checking any nutrient-versus-benthic response',
 'source_sha256':EXPECTED,
 'uid_records':{'site':len(site_idx),'water':len(water_idx),'joint_unique_concordant':len(unique),
  'berw_plus_chem_TN_TP_qualified_unique_native_UID':len(manifest),
  'distinct_physical_SITE_ID':len(group),
  'SITE_ID_with_repeated_UID':sum(len(v)>1 for v in group.values()),
  'SITE_ID_with_same_year_date_repeated_UID':site_date_dups,
  'nativeUID_duplicates':dups,
  'n_states':len(states),'n_ag_eco9':len(eco),'n_HUC8':len({z['huc8'] for z in manifest}),
  'BERW_TOTAL300_individual_count_sum':mass,
  'FFG_count_mass_coverage_fraction':ffg/mass if mass else None,
  'ORDER_FAMILY_count_mass_coverage_fraction':tax/mass if mass else None,
  'observed_N_and_P_units':[{ 'NTL':a,'PTL':b,'n_uid':n } for (a,b),n in sorted(chem_units.items())]},
 'no_ecological_effects_computed':True,
 'repaired_uid_source_gate_pass':gate,
 'source_freeze_manifest_sha256':None
}
frozen={
 'schema':'NRSA1819_EPASOURCE_UID_BERW_TNTP_POST_DUPLICATE_REPAIR_FROZEN_BEFORE_EFFECT_v1',
 'original_source_sha256':EXPECTED,
 'notes':'This is river benthic and river water chemistry; no WWTP causal assignment; no plants/butterflies or actual insect emergence included.',
 'primary_berw_UID_subset':manifest
}
fbytes=(json.dumps(frozen,sort_keys=True,ensure_ascii=False,separators=(',',':'))+'\n').encode()
fpath=p/'NRSA1819_BERW_NATIVE_UID_SITE_GEOGRAPHY_MANIFEST_SOURCE_ONLY.json'
fpath.write_bytes(fbytes)
out['source_freeze_manifest_sha256']=hashlib.sha256(fbytes).hexdigest()
rpath=p/'NRSA1819_NATIVE_UID_SOURCE_QA_RESULT.json'
rpath.write_text(json.dumps(out,sort_keys=True,indent=2))
print('NRSA1819_REPAIRED_UID_SOURCE_RESULT_START')
print(json.dumps(out,sort_keys=True))
print('NRSA1819_REPAIRED_UID_SOURCE_RESULT_END')
print('Manifest (source-only and no nutrient values): '+str(fpath))
