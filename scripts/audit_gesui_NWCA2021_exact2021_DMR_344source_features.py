#!/usr/bin/env python3
"""Independent source-only audit of 2021 DMR exact feature-to-wetland support.

No original NWCA vegetation or water-chemistry responses opened. Recompute
exact 344 permitted features, 392 source pairs, calendar 2021 numeric
N/P/flow by feature, candidate site support and distinct source geography.
This verifies reporting opportunity, NOT wastewater transport or a causal dose.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from collections import defaultdict
from pathlib import Path

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
S=D/"EPA_NWCA2021_160_SEWAGE_PERMITS_BEFORE_DMR2021_SOURCE.json"
F=D/"EPA_NWCA2021_128WETLAND_160POTW_EXACT_NEAR_OUTFALL_FEATURE_SOURCE.json"
J=D/"EPA_NWCA2021_ORIGINAL_160POTW_2021_DMR_TNTP_SOURCE_FEASIBILITY.json"
G=D/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
R=D/"EPA_NWCA2021_EXACT_344POTW_OUTFALL_2021_TNTP_DMR_NHDPLUS_SOURCE_GATE.json"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def verify(root):
    source=json.loads((root/S).read_text())
    f=json.loads((root/F).read_text())
    d=json.loads((root/J).read_text())
    g=json.loads((root/G).read_text())
    original=json.loads((root/R).read_text())
    if d["source_preDMR_160_permits_SHA256"]!=sha(root/S) or f["source_original160_permit_DMR_manifest_SHA256"]!=sha(root/S):
        raise ValueError("Original prior DMR source permits not frozen")
    if original["original_pre2021_DMR_exact344geo_features_SHA256"]!=sha(root/F) or original["original_FY2021andFY2022_source_DMR_join_SHA256"]!=sha(root/J):
        raise ValueError("EPA original DMR and mapped permitted features do not match original frozen source receipts")
    if f["original_128_near10_wetland_site_location_count"]!=128 or f["original_near10_exact_permit_feature_ID_count"]!=344 or f["original_site_outfall_candidate_feature_pairs"]!=392:
        raise ValueError("Original source near128 feature344 pair392 population changed")
    pairs=f["original_source_case_site_to_exact_permit_feature_sha256_and_distance"]
    if len(pairs)!=392 or len({(x["NWCA2021_site_visit_SHA256"],x["NPDES_exact_permit_feature_SHA256"]) for x in pairs})!=392:
        raise ValueError("Source features repeated or original source candidate count changed")
    for x in pairs:
        if not 0<=x["source_geodesic_outfall_to_wetland_km"]<=10:
            raise ValueError("Original outfall points not actually within 10km of public NWCA source coordinates")
    original_sites=g["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]
    near=[v for v in original_sites if v["outfall_rows_within_10km"]>0]
    far=[v for v in original_sites if v["nearest_outfall_geodesic_km_within25km_else_null"] is None or v["nearest_outfall_geodesic_km_within25km_else_null"]>=20]
    regions={v["original_state"]+"|"+v["original_ecoregion"] for v in near}&{v["original_state"]+"|"+v["original_ecoregion"] for v in far}
    originalmatched={v["NWCA2021_source_site_visit_SHA256"] for v in near if v["original_state"]+"|"+v["original_ecoregion"] in regions}
    if len(originalmatched)!=103:raise ValueError("Original state-ecoregion within-NWCA comparison no longer 103 near sites")
    features={v["NPDES_exact_permit_feature_SHA256"] for v in pairs}
    if len(features)!=344:raise ValueError("Original WPI feature source not 344")
    measurement={}
    for row in d["source_original_2021_exact_permit_feature_DMR_reporting_summary"]:
        key=row["NPDES_exact_permit_feature_SHA256"]
        if key in measurement:raise ValueError("Duplicate source permit feature summary")
        measurement[key]={
            "N":set(row["reported_original_numeric_2021_TN00600_period_end_dates"]),
            "P":set(row["reported_original_numeric_2021_TP00665_period_end_dates"]),
            "FLOW":set(row["reported_original_numeric_2021_flow50050_period_end_dates"])
        }
    selected={k:v for k,v in measurement.items() if k in features}
    props={
        "N":lambda v:bool(v["N"]),
        "P":lambda v:bool(v["P"]),
        "bothNandP":lambda v:bool(v["N"] and v["P"]),
        "eitherN_orP":lambda v:bool(v["N"] or v["P"]),
        "FLOW":lambda v:bool(v["FLOW"]),
        "NandP_same_report_period":lambda v:bool(v["N"]&v["P"])
    }
    report={}
    for name,fn in props.items():
        chosen={k for k,v in selected.items() if fn(v)}
        matching=[x for x in pairs if x["NPDES_exact_permit_feature_SHA256"] in chosen]
        month_count=set()
        for k in chosen:
            for kind in ("N","P","FLOW"):
                month_count.update(z for z in selected[k][kind] if "2021" in z)
        states={x["source_original_US_state"] for x in matching}
        report[name]={
            "matched_exact_source_features":len(chosen),
            "original_128near_wetland_site_count":len({x["NWCA2021_site_visit_SHA256"] for x in matching}),
            "original_103matched_near_wetland_site_count":len({x["NWCA2021_site_visit_SHA256"] for x in matching if x["NWCA2021_site_visit_SHA256"] in originalmatched}),
            "unique_NPDES_permits_with_genuine_feature_measurement":len({x["NPDES_registered_permit_SHA256"] for x in matching}),
            "original_states_with_feature_measurements":len(states),
            "original_2021_reported_periods":len(month_count)
        }
        old=original["exact_original_outfall_GEOMETRY_to_2021_DMR_feature_matches"][name]
        for k,oldkey in (("matched_exact_source_features","exact_nearby_source_feature_ids_with_2021_reported_numeric"),
                         ("original_128near_wetland_site_count","matched_source_NEAR_wetlands_128"),
                         ("original_103matched_near_wetland_site_count","matched_original_27state_study_NEAR_wetlands_103"),
                         ("unique_NPDES_permits_with_genuine_feature_measurement","original_2021_source_administrative_permit_IDs"),
                         ("original_states_with_feature_measurements","original_states_with_near_registered_outfall_and_2021_DMR")):
            if report[name][k]!=old[oldkey]:
                raise ValueError(f"Independent exact-feature source count mismatch {name}/{k}")
    bothNP_flow=sum(1 for k,v in selected.items() if (v["N"] or v["P"]) and v["FLOW"])
    if bothNP_flow!=original["original_reported_NP_outfall_features_with_numeric_flow_at_same_feature"]:
        raise ValueError("Independent exact nearby outfall with nutrient and flow reporting differs")
    need=original["source_precommitted_2021_nutrient_monitoring_feature_opportunity_gate"]
    observed=(report["eitherN_orP"]["original_128near_wetland_site_count"]>=need["min_near_wetlands_with_actual_2021_nutrient_DMR_at_exact_outfall_number"] and
              report["eitherN_orP"]["unique_NPDES_permits_with_genuine_feature_measurement"]>=need["min_original_independent_registered_municipal_sewage_permit_ids"] and
              report["eitherN_orP"]["original_states_with_feature_measurements"]>=need["min_original_states"] and
              bothNP_flow/max(1,report["eitherN_orP"]["matched_exact_source_features"])>=need["min_fraction_of_same_outfall_features_with_flow"])
    if observed!=need["eligible_to_proceed_to_separate_NHDPlus_downstream_flow_route"]:
        raise ValueError("Original pre-ecology macro source feasibility gate not independently reproduced")
    return {
        "schema":"gesui_NWCA2021_128registered_near_wetland_344POTW_features_CY2021_DMR_TNTP_original_independent_arithmetic_QA_v01",
        "source_file_SHA256":{"original_preDMR_160permit_source":sha(root/S),
                              "original_344_feature_geo":sha(root/F),
                              "FY2021_FY2022_nutrient_DMR":sha(root/J),
                              "original_568_source_near_far":sha(root/G),
                              "original_exact_feature_feasibility":sha(root/R)},
        "source_frozen_original_128near_344feature_392pair_reconfirmed":True,
        "source_matched_state_ecoregion_103original_near_reconfirmed":True,
        "independent_original_2021_reported_numeric_nutrient_and_flow_exact_outfall_counts":report,
        "number_same_outfall_feature_with_original_2021_NP_and_flow":bothNP_flow,
        "original_strict_feature_report_SOURCE_availability_passed":observed,
        "no_actual_2021_NWCA_Shannon_or_2021_wetland_water_chemistry_read":True,
        "no_verified_downstream_NHDPlus_hydrologic_route":True,
        "no_2021_plant_biodiversity_sewage_effect_computed":True,
    }

def main():
    a=argparse.ArgumentParser();a.add_argument("--out",type=Path,required=True);v=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    result=verify(root)
    v.out.parent.mkdir(parents=True,exist_ok=True)
    v.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps(result,indent=2))
if __name__=="__main__":main()
