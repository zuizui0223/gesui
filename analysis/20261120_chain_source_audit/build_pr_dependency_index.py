"""Reconcile the frozen PR census with the earlier research register.

Groups are conservative review families, NOT independent studies or effect sizes.
Historical judgments are copied without recoding blank or mixed entries.
"""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / 'docs/manuscript_riparian_interactions'
GROUPS = {
    'local_europe_outfall_plants': [2,3,4,6,7,8,9,10,12],
    'eea_benthos_reported_loads': [14,16,17,19,56,57],
    'german_austrian_photo_and_river_followups': list(range(20,37)),
    'danish_standardized_butterflies': list(range(37,44)),
    'global_ports_recording': list(range(44,48)),
    'us_wetland_plants': list(range(48,56)),
    'us_river_habitat_and_emergence': [58,59,60,61],
    'europe_common_basin_panel': list(range(62,70))+[84,85,86,88,89,90,94,95,96,101,102,103,104,105,108,109,110,114,119,121,134],
    'france_station_water_benthos': [70,71,72],
    'us_nrsa_multicycle': list(range(73,84))+[97,98,99,100,133],
    'integration_archives': [87,106,112,124],
    'europe_source_direction_and_support': [91,92,93,107,111,113,115,116],
    'austria_same_station_remote_sensing': [117,118,122],
    'europe_same_station_riparian_landcover': [120,123],
    'point_soil_source_resolution': [125,126,127,128],
    'sludge_water_source_routes': [129,130,131,132],
    'english_chempop': [135],
}
# Families share these source components; this does not certify identical rows.
COMPONENTS = {
    'local_europe_outfall_plants': ['europe_wwtp', 'wise_water', 'gbif_plants'],
    'eea_benthos_reported_loads': ['europe_wwtp', 'eea_benthos', 'jrc_point_loads'],
    'german_austrian_photo_and_river_followups': ['europe_wwtp', 'wise_water', 'gbif_butterfly_photos'],
    'danish_standardized_butterflies': ['europe_wwtp', 'danish_ebms'],
    'global_ports_recording': ['global_ports', 'gbif_multitaxon'],
    'us_wetland_plants': ['epa_nwca', 'epa_potw_dmr'],
    'us_river_habitat_and_emergence': ['epa_nrsa', 'external_emergence_17streams'],
    'europe_common_basin_panel': ['europe_wwtp', 'eea_benthos', 'gbif_plants', 'gbif_butterflies', 'eiv_traits', 'host_catalogue'],
    'france_station_water_benthos': ['wise_water', 'eea_benthos', 'europe_wwtp'],
    'us_nrsa_multicycle': ['epa_nrsa', 'epa_potw_dmr'],
    'integration_archives': ['archive_only_no_new_independent_sample'],
    'europe_source_direction_and_support': ['europe_wwtp', 'wise_water', 'eea_benthos', 'gbif_plants'],
    'austria_same_station_remote_sensing': ['europe_wwtp', 'wise_water', 'eea_benthos', 'sentinel2', 'hr_vpp'],
    'europe_same_station_riparian_landcover': ['europe_wwtp', 'eea_benthos', 'copernicus_riparian_zones'],
    'point_soil_source_resolution': ['jrc_point_soil_loads', 'europe_common_basin_panel'],
    'sludge_water_source_routes': ['national_wwtp_water_sludge_records'],
    'english_chempop': ['english_chempop_water_benthos_edf'],
}

def main():
    snapshot_path = DOC / 'PR_SOURCE_SNAPSHOT_20261120.json'
    historic_path = ROOT / 'docs/PR_RESEARCH_REGISTER_20261116.csv'
    snapshot = json.loads(snapshot_path.read_text(encoding='utf-8'))['prs']
    with historic_path.open(encoding='utf-8-sig', newline='') as f:
        historical = {int(r['pr_number']): r for r in csv.DictReader(f)}
    groups = {}
    for family, numbers in GROUPS.items():
        for n in numbers:
            if n in groups:
                raise ValueError(f'Duplicate primary review assignment: {n}')
            groups[n] = family
    numbers = {p['number'] for p in snapshot}
    if numbers != set(groups):
        raise ValueError(f'Missing {numbers-set(groups)}, extra {set(groups)-numbers}')
    rows = []
    for p in sorted(snapshot, key=lambda p: p['number']):
        n = p['number']; old = historical.get(n, {})
        family = groups[n]
        rows.append({
            'pr_number': n, 'title': p['title'], 'review_family': family,
            'family_possible_source_components': '|'.join(COMPONENTS[family]),
            'historical_judgment_20261116': old.get('judgment',''),
            'historical_result_document': old.get('result_document',''),
            'judgment_status': 'historical_not_reaudited' if old else 'not_in_old_register',
            'pr_url': f'https://github.com/zuizui0223/gesui/pull/{n}',
            'frozen_body_sha256': hashlib.sha256(p.get('body','').encode()).hexdigest(),
            'independent_study_count': '',
            'effect_level_eligibility': 'requires_outcome_contrast_variance_and_sample_overlap_audit',
        })
    out = DOC / 'evidence/PR_DEPENDENCY_INDEX_20261120.csv'
    with out.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    receipt = {
        'snapshot_sha256': hashlib.sha256(snapshot_path.read_bytes()).hexdigest(),
        'historical_register_sha256': hashlib.sha256(historic_path.read_bytes()).hexdigest(),
        'snapshot_pr_count': len(snapshot), 'indexed_unique_pr_count': len(rows),
        'historical_count': len(historical),
        'new_to_historical_register': sorted(numbers-set(historical)),
        'review_families_not_independent_studies': {k:len(v) for k,v in GROUPS.items()},
        'missing_prs': [], 'new_ecological_effects': 0,
        'limits': ['Family components are possible shared sources, not verified identical rows for every PR.',
                   'Historical judgments are not new effect-level audits.',
                   'No independence, covariance, effect-size compatibility or source support is inferred from PR counts.',
                   'CELLDEX and other action-only evidence must be tracked separately.'],
    }
    (out.parent/'PR_DEPENDENCY_INDEX_RECEIPT_20261120.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'indexed_prs':len(rows),'review_families':len(GROUPS),'new_to_old_register':len(numbers-set(historical))}))

if __name__ == '__main__':
    main()
