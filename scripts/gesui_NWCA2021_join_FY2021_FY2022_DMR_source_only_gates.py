#!/usr/bin/env python3
"""Join original FY2021/FY2022 DMR REPORTING opportunities on locked 2021 sites.

Only counts of valid historical 2021 monitoring at the SAME frozen 2026
registered municipal permit+EXTERNAL outfall feature. Date dedup across FY,
not summed effluent mass; 2021 NODI/qualifier and unit heterogeneity retained.
No 2021 wetland vegetation Shannon / TN TP outcomes are touched.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

DIR=Path("results/us_NWCA2021_plant_TNTP_20261110")
CROSS=DIR/"NWCA2021_568WETLAND_160POTW_EXACT_PERMIT_FEATURE_SOURCE_BEFORE_2021DMR.json"
GEO=DIR/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(root,inputs):
    src=json.loads((root/CROSS).read_text())
    prior=json.loads((root/GEO).read_text())
    if src["source_previous_568_site_NPDES_geo_qualifier_SHA256"]!=sha(root/GEO):
        raise ValueError("Original 2026 municipal NPDES wetland source geometry was changed")
    if len(src["registered_2026_site_to_permit_feature_pairs_geodesic_le10km_without_2021_biology"])!=392:
        raise ValueError("Original 392 wetland/source feature pairs changed")
    features={
        (x["original_public_NPDES_permit_id"],x["original_public_NPDES_feature_number"])
        for x in src["registered_2026_site_to_permit_feature_pairs_geodesic_le10km_without_2021_biology"]
    }
    if len(features)!=344:
        raise ValueError("Original administrative source features no longer 344")
    years={}
    for p in inputs:
        x=json.loads(Path(p).read_text())
        y=int(x["source_FY"])
        if y in years or y not in (2021,2022):
            raise ValueError("Expected one official US EPA DMR 2021 and one DMR 2022 source report")
        if x["original_source_568wetland_160permit_344feature_PRE_DMR_SHA256"]!=sha(root/CROSS):
            raise ValueError("Original frozen 160-source permit/feature membership mismatch")
        if x["original_total_160_registered_candidate_permits"]!=160 or x["original_total_344_candidate_permitted_features"]!=344:
            raise ValueError("Original source candidate number changed")
        years[y]={"result":x,"file_SHA256":sha(Path(p))}
    if set(years)!={2021,2022}:raise ValueError("Both original fiscal year sources mandatory for calendar 2021")
    by=defaultdict(lambda:{"reporting_fys":set(),"reports":0,"valid_numeric":0,"valid_positive":0,
                           "dates":set(),"units":Counter(),"qualifiers":Counter(),"nodis":Counter()})
    for fy,z in sorted(years.items()):
        for x in z["result"]["all_exact_original_registered_permit_FEATURE_DMR_2021_reporting_records"]:
            key=(x["original_public_NPDES_permit_id"],
                 x["original_public_NPDES_feature_number"],
                 x["DMR_2021_source_pollutant_code"])
            if key[:2] not in features:
                raise ValueError("Reconstructed 2021 discharge report for UNSELECTED source feature")
            v=by[key]
            v["reporting_fys"].add(fy)
            v["reports"]+=x["source_2021_FY"+str(fy)+"_EXO_reported_rows"]
            v["valid_numeric"]+=x["source_2021_uncensored_numeric_records"]
            v["valid_positive"]+=x["source_2021_positive_numeric_records"]
            v["dates"].update(x["source_2021_valid_numeric_monitoring_end_dates"])
            v["units"].update(x["source_2021_unit_value_type_statistic_counts"])
            v["qualifiers"].update(x["source_2021_qualifier_counts"])
            v["nodis"].update(x["source_2021_NODI_code_counts"])
    # Existence of valid recorded numeric 2021 DMR at EXACT EXO,
    # not proof effluent was received by a wetland or that mass was emitted.
    n={key[:2] for key,v in by.items() if key[2]=="00600" and v["dates"]}
    p={key[:2] for key,v in by.items() if key[2]=="00665" and v["dates"]}
    positive_n={key[:2] for key,v in by.items() if key[2]=="00600" and v["valid_positive"]>0}
    positive_p={key[:2] for key,v in by.items() if key[2]=="00665" and v["valid_positive"]>0}
    both=n&p
    site_to_feature=defaultdict(set)
    site_groups={}
    for x in src["registered_2026_site_to_permit_feature_pairs_geodesic_le10km_without_2021_biology"]:
        sid=x["source_Wetland_SITE_ID_VISIT_SHA256"]
        feature=(x["original_public_NPDES_permit_id"],x["original_public_NPDES_feature_number"])
        site_to_feature[sid].add(feature)
        site_groups[sid]={"state":x["source_original_state"],
                          "ecoregion":x["source_original_NWCA_ecoregion"]}
    # Original 568 near+far is needed to establish source-fixed 32 stratum
    # overlap, NOT by choosing a favorable plant outcome direction.
    previous=prior["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]
    stateeco_near={(x["original_state"],x["original_ecoregion"]) for x in previous if x["outfall_rows_within_10km"]>0}
    stateeco_far={(x["original_state"],x["original_ecoregion"]) for x in previous
                  if x["nearest_outfall_geodesic_km_within25km_else_null"] is None or
                     x["nearest_outfall_geodesic_km_within25km_else_null"]>=20}
    shared=stateeco_near&stateeco_far
    if len(shared)!=32:raise ValueError("Original 32 state/ecoregion overlaps changed")
    eligible_near={sid for sid in site_to_feature if (site_groups[sid]["state"],site_groups[sid]["ecoregion"]) in shared}
    if len(eligible_near)!=103:raise ValueError("Original 103 biologically comparable source candidate sites changed")
    farelig={x["NWCA2021_source_site_visit_SHA256"] for x in previous
             if ((x["original_state"],x["original_ecoregion"]) in shared and
                  (x["nearest_outfall_geodesic_km_within25km_else_null"] is None or
                   x["nearest_outfall_geodesic_km_within25km_else_null"]>=20))}
    if len(farelig)!=137:raise ValueError("Original 137 matched distance control sites changed")
    sampled_any={sid for sid,fset in site_to_feature.items() if fset&(n|p)}
    sampled_both={sid for sid,fset in site_to_feature.items() if fset&both}
    numeric_positive_both={sid for sid,fset in site_to_feature.items() if fset&(positive_n&positive_p)}
    def stats(site_ids,featuresets):
        ft={f for sid in site_ids for f in site_to_feature[sid]&featuresets}
        return {
            "n_original_NWCA_source_wetland_centroids":len(site_ids),
            "distinct_2021_measured_external_outfall_permit_IDs":len({a for a,b in ft}),
            "distinct_2021_measured_external_outfall_permit_features":len(ft),
            "n_US_states":len({site_groups[sid]["state"] for sid in site_ids}),
            "n_NWCA_ecoregions":len({site_groups[sid]["ecoregion"] for sid in site_ids}),
            "source_original_state_counts":dict(Counter(site_groups[sid]["state"] for sid in site_ids)),
            "source_original_site_visit_SHA256s":sorted(site_ids)
        }
    any_stats=stats(sampled_any,n|p)
    both_stats=stats(sampled_both,both)
    site_both_comparable=eligible_near&sampled_both
    matched_stats=stats(site_both_comparable,both)
    passed=bool(
        any_stats["n_original_NWCA_source_wetland_centroids"]>=20 and
        any_stats["distinct_2021_measured_external_outfall_permit_IDs"]>=12 and
        any_stats["n_US_states"]>=8 and
        both_stats["n_original_NWCA_source_wetland_centroids"]>=12 and
        both_stats["distinct_2021_measured_external_outfall_permit_IDs"]>=8 and
        both_stats["n_US_states"]>=5 and
        matched_stats["n_original_NWCA_source_wetland_centroids"]>=10 and
        len(farelig)>=60
    )
    result={
        "schema":"gesui_US_NWCA2021_568source_160NPDES_344FEATURE_EPA_2021_calendar_DMR_FY2021_FY2022_JOIN_SOURCE_GATE_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_568site_160permit_344feature_prior_to_DMR_SHA256":sha(root/CROSS),
        "original_568site_2026ECHO_STRICT_SEWAGE_geographic_PREOUTCOME_SHA256":sha(root/GEO),
        "source_original_US_EPA_DMR_FY_file_SHA256s":{
           str(y):{"original_fiscal_ZIP_SHA256":z["result"]["original_official_EPA_ICIS_NPDES_fiscal_DMR_zip_SHA256"],
                   "original_report_SHA256":z["file_SHA256"],
                   "original_National_FY_count_filter":z["result"]["original_fiscal_DMR_national_row_and_filter_QC"]}
           for y,z in years.items()},
        "original_source_near10_2026_registered_2021_plausible_568wetland_sites":128,
        "original_source_27state_32strata_matched_near_sites":len(eligible_near),
        "original_source_27state_32strata_matched_far_sites":len(farelig),
        "original_2021_calendar_DMR_EXO_feature_counts":{
            "any_valid_numeric_TN":len(n),
            "any_valid_numeric_TP":len(p),
            "both_TN_and_TP_same_permit_external_feature":len(both),
            "both_strict_numeric_positive_TN_and_TP":len(positive_n&positive_p),
        },
        "original_2021_DMR_sites_with_valid_TN_OR_TP_same_source_feature":any_stats,
        "original_2021_DMR_sites_with_valid_BOTH_TN_and_TP_same_source_feature":both_stats,
        "original_2021_DMR_32matched_stateecoregion_NEAR_sites_with_BOTH_TN_TP":matched_stats,
        "strict_numeric_positive_TN_and_TP_same_feature_site_n":len(numeric_positive_both),
        "precommitted_DMR_SOURCE_reporting_opportunity_gate":{
           "at_least_20_NEAR_wetlands_with_any_real_TN_or_TP_records":20,
           "at_least_12_distinct_permits_with_any_report":12,
           "at_least_8_states_with_any_report":8,
           "at_least_12_NEAR_wetlands_with_both_TN_and_TP_records":12,
           "at_least_8_distinct_permits_with_both":8,
           "at_least_5_states_with_both":5,
           "at_least_10_NEAR_sites_in_original_32_stateecoregion_matched_plots_with_both":10,
           "at_least_60_existing_far_matched_controls":60,
           "DMR_monitored_same_2021_source_feature_biological_sensitivity_allowed":passed},
        "full_source_permit_feature_parameter_monitoring_receipts":[{
            "original_public_NPDES_permit_id":key[0],
            "original_public_NPDES_feature_number":key[1],
            "original_pollutant_code":key[2],
            "original_2021_reported_event_rows_two_FY_raw_additive_only":v["reports"],
            "original_FY_rows_MAY_be_duplicated_between_archives":True,
            "original_2021_valid_numeric_monitor_period_end_dates_unique":sorted(v["dates"]),
            "original_2021_valid_numeric_rows_no_dedup":v["valid_numeric"],
            "original_2021_positive_numeric_rows_no_dedup":v["valid_positive"],
            "original_DMR_units_valueTypes_statistics":dict(v["units"]),
            "original_DMR_qualifier_codes":dict(v["qualifiers"]),
            "original_DMR_NODI_codes":dict(v["nodis"]),
            "source_fiscal_archive_years_with_reports":sorted(v["reporting_fys"])
        } for key,v in sorted(by.items())],
        "prior_2021_NWCA_plant_H_ALL_FQAI_actual_water_TNTP_measurements_NOT_READ_in_this_source":True,
        "original_calendar2021_monitoring_end_dates_not_FY_filename_surrogate":True,
        "same2021_2026registered_source_to_wetland_hydrological_connection_STILL_NOT_verified":True,
        "DMR_source_numeric_concentration_and_reported_quantity_NOT_INTERCHANGEABLE":True,
        "no_WWTP_shannon_positive_or_negative_causal_effect_estimated":True
    }
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--fy-json",type=Path,required=True,action="append")
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if len(a.fy_json)!=2:raise ValueError("Both exact FY2021/FY2022 national reporting source summaries required")
    root=Path(__file__).resolve().parents[1]
    out=run(root,a.fy_json)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "EPA_original_features_2021_NP_report":out["original_2021_calendar_DMR_EXO_feature_counts"],
        "near_any":{k:v for k,v in out["original_2021_DMR_sites_with_valid_TN_OR_TP_same_source_feature"].items()
                    if not k.endswith("SHA256s")},
        "near_both":{k:v for k,v in out["original_2021_DMR_sites_with_valid_BOTH_TN_and_TP_same_source_feature"].items()
                    if not k.endswith("SHA256s")},
        "matched_both":{k:v for k,v in out["original_2021_DMR_32matched_stateecoregion_NEAR_sites_with_BOTH_TN_TP"].items()
                    if not k.endswith("SHA256s")},
        "gate":out["precommitted_DMR_SOURCE_reporting_opportunity_gate"],
    },indent=2))

if __name__=="__main__":
    main()
