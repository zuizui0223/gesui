#!/usr/bin/env python3
"""Independent OFFLINE arithmetic review of official CY2021 POTW DMR source gate.

Uses ONLY the two original successful FY source summary artifacts and archived
original 568 NWCA wetland/2026 NPDES source-only coordinates. Does not download
DMR raw reports nor inspect H_ALL, TN, TP at NWCA. The first main join code
is neither imported nor invoked; recomputes all site/source support by sets.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
SOURCE=D/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
ORIGINAL=D/"EPA_NWCA2021_160POTW_REAL2021_DMR_CALENDAR_FEATURE_AND_NP_SOURCE_QC.json"
SETS={
 "permit_any":"source_2021_permit_ids_with_any_DMR_record_SHA256",
 "permit_numeric":"source_2021_permit_ids_with_any_numeric_DMR_SHA256",
 "permit_exo":"source_2021_permit_ids_with_EXO_effluent_numeric_DMR_SHA256",
 "feature_exo":"source_2021_near10_registered_feature_pairs_with_EXO_effluent_numeric_DMR_SHA256",
 "feature_N":"source_2021_near10_feature_pairs_with_00600_totalN_DMR_SHA256",
 "feature_P":"source_2021_near10_feature_pairs_with_00665_totalP_DMR_SHA256",
 "permit_N":"source_2021_permit_ids_with_00600_totalN_DMR_SHA256",
 "permit_P":"source_2021_permit_ids_with_00665_totalP_DMR_SHA256",
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def audit(root,fy21,fy22):
    source=json.loads((root/SOURCE).read_text())
    result=json.loads((root/ORIGINAL).read_text())
    p=json.loads(fy21.read_text())
    q=json.loads(fy22.read_text())
    assert p["fiscal_archive"]==2021 and q["fiscal_archive"]==2022
    assert p["calendar2021_monitoring_period_subset"]==["2021-01-01","2021-09-30"]
    assert q["calendar2021_monitoring_period_subset"]==["2021-10-01","2021-12-31"]
    assert result["original_NWCA_wetland_2021_POTW_geo_source_SHA256"]==sha(root/SOURCE)
    assert result["official_EPA_ECHO_NPDES_FY2021_FY2022_DMR_2021_window_source_SHA256"]["FY2021_original_source_result_SHA256"]==sha(fy21)
    assert result["official_EPA_ECHO_NPDES_FY2021_FY2022_DMR_2021_window_source_SHA256"]["FY2022_original_source_result_SHA256"]==sha(fy22)
    assert (p["original_source_near10_geographic_wetlands"],q["original_source_near10_geographic_wetlands"])==(128,128)
    assert (p["original_source_near10_distinct_NPDES_permit_numbers"],q["original_source_near10_distinct_NPDES_permit_numbers"])==(160,160)
    originalsource=source["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]
    originalrow={x["NWCA2021_source_site_visit_SHA256"]:x for x in originalsource}
    assert len(originalrow)==568
    near={sid for sid,z in originalrow.items() if z["outfall_rows_within_10km"]>0}
    assert len(near)==128
    permits=p["original_NEAR10_NWCA_wetland_to_target_NPDES_permit_hashes"]
    features=p["original_NEAR10_NWCA_wetland_to_target_feature_pair_hashes"]
    assert permits==q["original_NEAR10_NWCA_wetland_to_target_NPDES_permit_hashes"]
    assert features==q["original_NEAR10_NWCA_wetland_to_target_feature_pair_hashes"]
    assert set(permits)==near==set(features)
    # Independent original 568 wetland SOURCE-only has exact near 10km permit
    # hash memberships and 27 state×eco original matched source opportunities.
    for sid in near:
        if sorted(permits[sid])!=originalrow[sid]["independent_10km_NPDES_permits_SHA256"]:
            raise ValueError("One near wetland switched 2026 NPDES registered source permit after 2021 DMR query")
    original_permit_hashes=set(p["original_near10_potentially_2021_sewage_permits_SHA256"])
    original_feature_hashes=set(p["original_near10_feature_pairs_SHA256"])
    assert len(original_permit_hashes)==160
    assert len(original_feature_hashes)==p["original_near10_distinct_registered_source_feature_pairs"]
    assert original_permit_hashes==set(q["original_near10_potentially_2021_sewage_permits_SHA256"])
    assert original_feature_hashes==set(q["original_near10_feature_pairs_SHA256"])
    combined={label:set(p[col])|set(q[col]) for label,col in SETS.items()}
    for label,actual in combined.items():
        expected=original_feature_hashes if label.startswith("feature") else original_permit_hashes
        if not actual.issubset(expected):raise ValueError("2021 source DMR contains nonoriginal candidate registered source "+label)
    by={
        "permit_any":"any_2021_DMR_permit",
        "permit_numeric":"any_2021_numeric_DMR_permit",
        "permit_exo":"2021_external_outfall_effluent_numeric_permit",
        "feature_exo":"2021_exact_nearby_external_outfall_effluent_feature",
        "feature_N":"2021_exact_nearby_external_outfall_total_N_feature",
        "feature_P":"2021_exact_nearby_external_outfall_total_P_feature"
    }
    independently={}
    for key,original_label in by.items():
        observed=combined[key]
        row={sid for sid in near if set((features if key.startswith("feature") else permits)[sid]) & observed}
        states={originalrow[sid]["original_state"] for sid in row}
        previous=result["original_2021_DMR_actual_numerical_reporting_availability_by_source_level"][original_label]
        if previous["source_eligible_WWTP_near10_wetland_ID_SHA256"]!=sorted(row):
            raise ValueError("Independent original near site identities differ from published CY2021 DMR "+key)
        vals=(len(row),len(observed),len(states))
        expected=(previous["original_source_wetlands_with_at_least_one_reporting_permit_or_registered_feature"],
                  previous["original_distinct_administrative_permits_or_feature_source_hashes"],
                  previous["states_with_site_report"])
        if vals!=expected:raise ValueError("Reconstructed 2021 DMR count not identical: "+key+repr((vals,expected)))
        independently[original_label]={"source_wetlands":len(row),"source_permits_or_features":len(observed),"source_states":len(states)}
    both_perm=combined["permit_N"]&combined["permit_P"]
    both_feat=combined["feature_N"]&combined["feature_P"]
    both_sites={sid for sid in near if set(features[sid]) & both_feat}
    both_states={originalrow[sid]["original_state"] for sid in both_sites}
    if (len(both_perm),len(both_feat),len(both_sites),len(both_states))!=(
        result["source_permit_DMR_near_wetland_totalN_P_00600_00665_both_count"],
        result["source_exact_feature_DMR_near_wetland_totalN_P_00600_00665_both_count"],
        result["source_wetlands_with_exact_nearby_outfall_both_TN_TP_reports_n"],
        result["source_states_with_exact_nearby_outfall_both_TN_TP_reports_n"]
    ):
        raise ValueError("Original 2021 dual TN/TP source-receiver reporting coverage not reproducible")
    gate=result["predeclared_SOURCE_reporting_gate"]
    eff=(len(combined["permit_exo"])>=30 and independently[by["permit_exo"]]["source_wetlands"]>=60
         and independently[by["permit_exo"]]["source_states"]>=10
         and independently[by["feature_exo"]]["source_wetlands"]>=40
         and len(combined["feature_exo"])>=20)
    joint=(eff and len(both_perm)>=15 and len(both_feat)>=15 and
           len(both_sites)>=30 and len(both_states)>=8)
    assert eff is True and joint is False
    assert gate["actual_2021_EFFLUENT_reporting_support_PASS"]==eff
    assert gate["actual_2021_N_and_P_measurement_support_PASS"]==joint
    assert gate["biological_plant_Shannon_WWTP_inference_permitted_without_true_RECEIVING_WATER_ROUTE"] is False
    return {
        "schema":"gesui_NWCA2021_POTW_160permit_2021_actual_EPA_DMR_FY2021FY2022_independent_SOURCE_JOIN_QA_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_source_files_SHA256":{
            "original_568_NWCA_geometric_registered_POTW_before_outcomes":sha(root/SOURCE),
            "original_FY2021_ECHO_DMR_source_receipt":sha(fy21),
            "original_FY2022_ECHO_DMR_source_receipt":sha(fy22),
            "first_2021_calendar_DMR_join_original_result":sha(root/ORIGINAL)},
        "original_568_study_sites_and_128_registered_near_source_sites_reproduced":True,
        "original_160_near_administrative_permits_and_all_regulated_features_reproduced":True,
        "original_2021_municipal_DMR_9months_plus_3months_correctly_joined":True,
        "all_six_permit_and_feature_source_site_and_US_state_counts_exactly_reproduced":independently,
        "original_both_TN_and_TP_total_numeric_00600_00665_permits":len(both_perm),
        "original_both_TN_and_TP_same_registered_feature_n":len(both_feat),
        "original_NWCA_wetlands_near_same_feature_both_nutrients":len(both_sites),
        "original_NWCA_wetland_US_states_near_same_feature_both_nutrients":len(both_states),
        "source_EFFLUENT_registration_and_historical_DMR_support_gate_PASS":eff,
        "source_TN_and_TP_both_support_gate_PASS":joint,
        "sewage_outfall_to_2021_wetland_hydrological_route_STILL_UNVERIFIED":True,
        "original_2021_NWCA_plant_diversity_and_local_water_chemistry_VALUES_NOT_READ":True,
        "no_nutrient_biodiversity_causal_WWTP_effect_regression_performed":True
    }
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--fy2021-source",type=Path,required=True)
    p.add_argument("--fy2022-source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    o=audit(root,a.fy2021_source,a.fy2022_source)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in o.items() if k in (
       "all_six_permit_and_feature_source_site_and_US_state_counts_exactly_reproduced",
       "original_both_TN_and_TP_total_numeric_00600_00665_permits",
       "original_both_TN_and_TP_same_registered_feature_n",
       "original_NWCA_wetlands_near_same_feature_both_nutrients",
       "source_EFFLUENT_registration_and_historical_DMR_support_gate_PASS",
       "source_TN_and_TP_both_support_gate_PASS")},indent=2))

if __name__=="__main__":
    main()
