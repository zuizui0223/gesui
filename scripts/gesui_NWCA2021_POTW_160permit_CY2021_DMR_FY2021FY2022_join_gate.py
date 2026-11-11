#!/usr/bin/env python3
"""SOURCE-ONLY full CALENDAR YEAR 2021 join from EPA FY2021+FY2022 DMR.

No reopening standardized NWCA vegetation H_ALL/FQAI or site water N/P
measurements. Primary output is how many original 160 strictly registered
near10 POTW permits actually report EXO outfall gross effluent (with NODI
removed) and totalN/totalP 00600/00665 during calendar2021. Permit-level
report is not exact outfall-level; both are transparently distinguished.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
SOURCE=D/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
FIELDS=(
 "source_2021_permit_ids_with_any_DMR_record_SHA256",
 "source_2021_permit_ids_with_any_numeric_DMR_SHA256",
 "source_2021_permit_ids_with_EXO_effluent_numeric_DMR_SHA256",
 "source_2021_near10_registered_feature_pairs_with_EXO_effluent_numeric_DMR_SHA256",
 "source_2021_near10_feature_pairs_with_00600_totalN_DMR_SHA256",
 "source_2021_near10_feature_pairs_with_00665_totalP_DMR_SHA256",
 "source_2021_permit_ids_with_00600_totalN_DMR_SHA256",
 "source_2021_permit_ids_with_00665_totalP_DMR_SHA256",
)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build(root,a,b):
    original=json.loads((root/SOURCE).read_text())
    first=json.loads(a.read_text());second=json.loads(b.read_text())
    for year,data in ((2021,first),(2022,second)):
        if data["schema"]!="gesui_US_NWCA2021_568source_160permit_ECHO_FY2021FY2022_actual_DMR_record_availability_SOURCE_ONLY_v01":
            raise ValueError("Original DMR source branch changed")
        if data["fiscal_archive"]!=year:
            raise ValueError("FY source not in its predeclared calendar-year slice")
        if data["original_POTW_geo_before_DMR_SHA256"]!=sha(root/SOURCE):
            raise ValueError("2021 reported original outlet source changed across first and second fiscal archive")
        if data["original_NPDES_2026_11_08_outfalls_ZIP_SHA256"]!=original["source_EPA_ECHO_NPDES_2026_live_weekly_original_ZIP_SHA256"]:
            raise ValueError("NPDES registered source list changed between DMR report archives")
        if (data["original_source_near10_geographic_wetlands"],
            data["original_source_near10_distinct_NPDES_permit_numbers"])!=(128,160):
            raise ValueError("Preexisting 128 wetlands / 160 permits source not reproduced")
        if year==2021 and data["calendar2021_monitoring_period_subset"]!=["2021-01-01","2021-09-30"]:
            raise ValueError("Original FY2021 source window incorrect")
        if year==2022 and data["calendar2021_monitoring_period_subset"]!=["2021-10-01","2021-12-31"]:
            raise ValueError("Original FY2022 source window incorrect")
    for name in ("original_NEAR10_NWCA_wetland_to_target_NPDES_permit_hashes",
                 "original_NEAR10_NWCA_wetland_to_target_feature_pair_hashes"):
        if first[name]!=second[name]:raise ValueError("Original NWCA source and NPDES feature memberships differ across 2021 quarterly DMR files")
    if first["original_near10_feature_pairs_SHA256"]!=second["original_near10_feature_pairs_SHA256"]:
        raise ValueError("Original registered feature universe changed across FY downloads")
    merged={k:set(first[k])|set(second[k]) for k in FIELDS}
    known=set(first["original_near10_potentially_2021_sewage_permits_SHA256"])
    if len(known)!=160 or known!=set(second["original_near10_potentially_2021_sewage_permits_SHA256"]):
        raise ValueError("Original registered 160 POTW sources not stable")
    allfeatures=set(first["original_near10_feature_pairs_SHA256"])
    for k,v in merged.items():
        expected=known if "permit_ids" in k else allfeatures
        if not v.issubset(expected):
            raise ValueError("A fiscal source DMR report does not match an originally registered near10 permit/feature: "+k)
    rows=original["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]
    rowmap={r["NWCA2021_source_site_visit_SHA256"]:r for r in rows}
    if len(rowmap)!=568:raise ValueError("Original 568 wetland source controls missing")
    near={r["NWCA2021_source_site_visit_SHA256"] for r in rows if r["outfall_rows_within_10km"]>0}
    if len(near)!=128:raise ValueError("Original near10 number changed")
    permit_by_site=first["original_NEAR10_NWCA_wetland_to_target_NPDES_permit_hashes"]
    feature_by_site=first["original_NEAR10_NWCA_wetland_to_target_feature_pair_hashes"]
    if set(permit_by_site)!=near or set(feature_by_site)!=near:
        raise ValueError("Source-only 128 near10 wetland features not consistently mapped")
    measures={}
    state_by={s:rowmap[s]["original_state"] for s in near}
    for name,group in (
        ("any_2021_DMR_permit","source_2021_permit_ids_with_any_DMR_record_SHA256"),
        ("any_2021_numeric_DMR_permit","source_2021_permit_ids_with_any_numeric_DMR_SHA256"),
        ("2021_external_outfall_effluent_numeric_permit","source_2021_permit_ids_with_EXO_effluent_numeric_DMR_SHA256"),
        ("2021_exact_nearby_external_outfall_effluent_feature","source_2021_near10_registered_feature_pairs_with_EXO_effluent_numeric_DMR_SHA256"),
        ("2021_exact_nearby_external_outfall_total_N_feature","source_2021_near10_feature_pairs_with_00600_totalN_DMR_SHA256"),
        ("2021_exact_nearby_external_outfall_total_P_feature","source_2021_near10_feature_pairs_with_00665_totalP_DMR_SHA256"),
    ):
        values=merged[group]
        bysite=feature_by_site if "feature" in name else permit_by_site
        eligible={sid for sid in near if set(bysite[sid])&values}
        measures[name]={
           "original_source_wetlands_with_at_least_one_reporting_permit_or_registered_feature":len(eligible),
           "original_distinct_administrative_permits_or_feature_source_hashes":len(values),
           "states_with_site_report":len({state_by[s] for s in eligible}),
           "source_eligible_WWTP_near10_wetland_ID_SHA256":sorted(eligible),
           "site_and_regional_counts_NOT_external_sewage_to_wetland_flow_paths":True,
        }
    bothpermit=merged["source_2021_permit_ids_with_00600_totalN_DMR_SHA256"]&merged["source_2021_permit_ids_with_00665_totalP_DMR_SHA256"]
    bothfeature=merged["source_2021_near10_feature_pairs_with_00600_totalN_DMR_SHA256"]&merged["source_2021_near10_feature_pairs_with_00665_totalP_DMR_SHA256"]
    withboth={sid for sid in near if set(feature_by_site[sid])&bothfeature}
    verified_EFFLUENT=(measures["2021_external_outfall_effluent_numeric_permit"]["original_distinct_administrative_permits_or_feature_source_hashes"]>=30
       and measures["2021_external_outfall_effluent_numeric_permit"]["original_source_wetlands_with_at_least_one_reporting_permit_or_registered_feature"]>=60
       and measures["2021_external_outfall_effluent_numeric_permit"]["states_with_site_report"]>=10
       and measures["2021_exact_nearby_external_outfall_effluent_feature"]["original_source_wetlands_with_at_least_one_reporting_permit_or_registered_feature"]>=40
       and measures["2021_exact_nearby_external_outfall_effluent_feature"]["original_distinct_administrative_permits_or_feature_source_hashes"]>=20)
    verified_NP=(verified_EFFLUENT and len(bothpermit)>=15 and len(bothfeature)>=15
                 and len(withboth)>=30 and len({state_by[s] for s in withboth})>=8)
    return {
       "schema":"gesui_US_NWCA2021_568wetland_160POTW_2021_calendaryear_ECHO_ACTUAL_DMR_sourceonly_gate_v01",
       "UTC":datetime.now(timezone.utc).isoformat(),
       "original_NWCA_wetland_2021_POTW_geo_source_SHA256":sha(root/SOURCE),
       "official_EPA_ECHO_NPDES_FY2021_FY2022_DMR_2021_window_source_SHA256":{
          "FY2021_original_source_result_SHA256":sha(a),
          "FY2022_original_source_result_SHA256":sha(b),
          "FY2021_download_zip_SHA256":first["original_ECHO_DMR_FY_archive_SHA256"],
          "FY2022_download_zip_SHA256":second["original_ECHO_DMR_FY_archive_SHA256"],
          "2026_11_registered_NPDES_outfall_source_SHA256":original["source_EPA_ECHO_NPDES_2026_live_weekly_original_ZIP_SHA256"]
        },
       "full_calendar_original_monitoring_period_2021":["2021-01-01","2021-12-31"],
       "FY2021_used_only_calendar_2021_Jan_Sep_and_FY2022_used_only_2021_Oct_Dec":True,
       "original_568_preplant_quality_SOURCE_n":568,
       "original_geographic_near10_potential_sites":128,
       "original_geographic_2026_sewage_permits_near10_n":160,
       "original_admin_NPDES_NAICS221320_or_SIC4952_and_permitted_feature_verified":True,
       "original_2021_DMR_actual_numerical_reporting_availability_by_source_level":measures,
       "source_permit_DMR_near_wetland_totalN_P_00600_00665_both_count":len(bothpermit),
       "source_exact_feature_DMR_near_wetland_totalN_P_00600_00665_both_count":len(bothfeature),
       "source_wetlands_with_exact_nearby_outfall_both_TN_TP_reports_n":len(withboth),
       "source_states_with_exact_nearby_outfall_both_TN_TP_reports_n":len({state_by[s] for s in withboth}),
       "predeclared_SOURCE_reporting_gate":{
          "min_original_2021_DMR_effluent_permits":30,
          "min_near10_wetland_original_sites_with_2021_any_effluent_permit_report":60,
          "min_distinct_states_with_2021_effluent_permits":10,
          "min_near10_wetlands_with_2021_exact_feature_EXO_numeric_report":40,
          "min_original_exact_feature_EXO_numeric_sources":20,
          "min_original_N_and_P_both_reporting_permits":15,
          "min_original_N_and_P_both_exact_feature_sources":15,
          "min_2021_near10_wetlands_with_both_totalN_totalP_feature_reports":30,
          "min_states_with_near10_bothNP_feature_reports":8,
          "actual_2021_EFFLUENT_reporting_support_PASS":bool(verified_EFFLUENT),
          "actual_2021_N_and_P_measurement_support_PASS":bool(verified_NP),
          "biological_plant_Shannon_WWTP_inference_permitted_without_true_RECEIVING_WATER_ROUTE":False,
       },
       "DMR_missing_or_nodi_not_called_zero_true_discharge":True,
       "actual_2021_DMR_numeric_values_not_exported_or_fitted":True,
       "DMR_2021_fiscal_available_in_2026_regulatory_system_does_not_prove_complete_historical_monitoring":True,
       "2026_registered_feature_location_NOT_actual_verified_receiving_wetland":True,
       "prior_2021_NWCA_area_weighted_positive_vs_site_equal_negative_NOT_used_to_select_sources":True
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--fy2021-source",required=True,type=Path)
    ap.add_argument("--fy2022-source",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=build(root,args.fy2021_source,args.fy2022_source)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
      "source":{k:{n:v[n] for n in ("original_source_wetlands_with_at_least_one_reporting_permit_or_registered_feature",
                 "original_distinct_administrative_permits_or_feature_source_hashes","states_with_site_report")}
                 for k,v in out["original_2021_DMR_actual_numerical_reporting_availability_by_source_level"].items()},
      "N_P_both_permit":out["source_permit_DMR_near_wetland_totalN_P_00600_00665_both_count"],
      "N_P_both_exact_feature":out["source_exact_feature_DMR_near_wetland_totalN_P_00600_00665_both_count"],
      "N_P_both_near_wetland":out["source_wetlands_with_exact_nearby_outfall_both_TN_TP_reports_n"],
      "gates":out["predeclared_SOURCE_reporting_gate"]
    },indent=2))

if __name__=="__main__":
    main()
