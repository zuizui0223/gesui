#!/usr/bin/env python3
"""SOURCE-ONLY audit of two duplicate site-visit keys in EPA NRSA 2018-19.
Must be executed AFTER the original source-gate acquisition job, same runner.
No cross-taxon effect, nutrient value, or ecological sign is calculated.
The original FAIL/STOP is retained; this emits whether a new explicitly revised
UID-source manifest could be defensibly frozen *before* numerical analysis.
"""
import collections,csv,hashlib,json,os,pathlib,math
p=pathlib.Path(os.environ.get('RUNNER_TEMP','/tmp'))/'gesui_nrsa1819_source_gate'
report=json.loads((p/'NRSA_1819_SOURCE_ONLY_QA.json').read_text())
paths={k:p/(k+'.csv') for k in ('sites','water','benthos')}
if any(not path.is_file() for path in paths.values()):raise RuntimeError('Original source data not present in same runner')
def key(r):
    site=(r.get('SITE_ID') or '').strip().upper()
    return (site,(r.get('VISIT_NO') or '1').strip().upper()) if site else None
def get(r,k):
    return (r.get(k) or '').strip()
def scan(path):
    with path.open(encoding='utf-8-sig',newline='') as fh:
        for r in csv.DictReader(fh):
            yield r
site=collections.defaultdict(list);water=collections.defaultdict(list)
for row in scan(paths['sites']):
    k=key(row)
    if k:site[k].append(row)
for row in scan(paths['water']):
    k=key(row)
    if k:water[k].append(row)
# No row ID printed, anonymous per-case diagnoses only.
dup_info={}
for nm,dct in [('sites',site),('water',water)]:
    duplicates=[v for v in dct.values() if len(v)>1]
    dup_info[nm]={
        'n_duplicate_native_visit_keys':len(duplicates),
        'extra_rows':sum(len(v)-1 for v in duplicates),
        'n_duplicates_with_only_one_native_UID':sum(len(set(get(r,'UID') for r in v))==1 for v in duplicates),
        'n_duplicates_with_all_columns_identical':sum(len({tuple(sorted(r.items())) for r in v})==1 for v in duplicates),
        'n_duplicates_with_same_geo_or_region':sum(
            len({tuple(get(r,c) for c in ('STATE','STATE_SF','AG_ECO9','LAT_DD83','LON_DD83')) for r in v})==1
            for v in duplicates),
        'uid_multisets_case_counts':sorted([len(set(get(r,'UID') for r in v)) for v in duplicates])
    }
# Cross-file diagnostic: every site+visit must refer to exactly one *same* UID
# across site and water. Repeated identical UID at the same visit is allowed
# in the derived source-only manifest ONLY if all metadata entries are bytewise
# identical; numerical nutrient QA duplication remains unresolved if rows vary.
site_to_uid={}
conflict_uid=0
for k,v in site.items():
    ss={get(r,'UID') for r in v}
    if len(ss)==1:site_to_uid[k]=next(iter(ss))
    else: conflict_uid+=1
water_to_uid={}
for k,v in water.items():
    ss={get(r,'UID') for r in v}
    if len(ss)==1:water_to_uid[k]=next(iter(ss))
    else: conflict_uid+=1
same_keys=set(site)&set(water)
cross_uid_mismatch=sum(site_to_uid.get(k)!=water_to_uid.get(k) for k in same_keys)
ambiguous_site_meta=sum(len({tuple(sorted(r.items())) for r in site[k]})>1 for k in same_keys)
ambiguous_water_records=sum(len({tuple(sorted(r.items())) for r in water[k]})>1 for k in same_keys)
# BERW actual primary site keys and total300 functional coverage can be counted
# independently of nutrient concentrations and ecological outcomes.
berw_sites=set()
ffg_count_mass=0.;total_mass=0.;class_count_mass=0.
for r in scan(paths['benthos']):
    if get(r,'SAMPLE_TYPE').upper()!='BERW':continue
    k=key(r)
    if k:berw_sites.add(k)
    try:m=float(get(r,'TOTAL300'))
    except ValueError:continue
    if m<=0 or not math.isfinite(m):continue
    total_mass+=m
    if get(r,'FFG'):ffg_count_mass+=m
    if get(r,'ORDER') or get(r,'FAMILY'):class_count_mass+=m
source_supported_site_keys=set(site)&set(water)&berw_sites
# Do NOT inspect water values; use original SOURCE-only gate's number for positive measured support.
out={
 'previous_original_gate_decision_retained':report['decision'],
 'source_sha256':{k:report['sources'][k]['sha256'] for k in ('sites','water','benthos')},
 'duplicate_audit':dup_info,
 'n_joint_site_water_berw_keys_prechemistry':len(source_supported_site_keys),
 'crossfile_same_sitevisit_UID_mismatch_count':cross_uid_mismatch,
 'site_uid_conflict_keys':conflict_uid,
 'n_site_meta_variants_for_same_native_key':ambiguous_site_meta,
 'n_water_meta_or_value_variants_for_same_native_key':ambiguous_water_records,
 'BERW_positive_TOTAL300_individual_mass':round(total_mass,1),
 'BERW_FFG_weighted_mass_coverage_fraction':round(ffg_count_mass/total_mass,6) if total_mass else None,
 'BERW_ORDER_FAMILY_weighted_mass_coverage_fraction':round(class_count_mass/total_mass,6) if total_mass else None
}
out['UID_reconciliation_feasible_source_only']=bool(
 conflict_uid==0 and cross_uid_mismatch==0 and ambiguous_site_meta==0
 and ambiguous_water_records==0
 and out['BERW_FFG_weighted_mass_coverage_fraction']>=.6
 and out['BERW_ORDER_FAMILY_weighted_mass_coverage_fraction']>=.8
)
out['status']='REVISED_SOURCE_FREEZE_POSSIBLE_ONLY_IF_FEASIBLE_TRUE' if out['UID_reconciliation_feasible_source_only'] else 'STOP_SOURCE_AMBIGUITY'
(p/'NRSA_1819_DUPLICATE_UID_PREOUTCOME_AUDIT.json').write_text(json.dumps(out,indent=2))
print('NRSA_1819_UID_SOURCE_DIAGNOSTIC_START')
print(json.dumps(out,sort_keys=True))
print('NRSA_1819_UID_SOURCE_DIAGNOSTIC_END')
