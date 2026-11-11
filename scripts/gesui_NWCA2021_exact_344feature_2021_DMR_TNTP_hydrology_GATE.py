#!/usr/bin/env python3
"""Actual CY2021 EPA DMR reports at EXACT original 344 NEAR feature IDs.

This is the source-only 2021 nutrient monitoring receptor feasibility gate.
Compares frozen 2026 ECHO registered geometry SHA × actual EPA FY2021 &
FY2022 numeric final effluent reports. Still NO verified hydrologic path, dose
to NWCA wetland, botanical Shannon, site measured TN/TP, or causal model.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
FEATURE=D/"EPA_NWCA2021_128WETLAND_160POTW_EXACT_NEAR_OUTFALL_FEATURE_SOURCE.json"
FY=D/"EPA_NWCA2021_ORIGINAL_160POTW_2021_DMR_TNTP_SOURCE_FEASIBILITY.json"
GEO=D/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def quarter(date):
    s=str(date)
    m=re.search(r"^2021-(\d{2})",s)
    if m:
        month=int(m.group(1))
        return (month-1)//3+1 if 1<=month<=12 else None
    m=re.search(r"^(\d{1,2})/.*2021$",s)
    if m:
        mo=int(m.group(1))
        return (mo-1)//3+1 if 1<=mo<=12 else None
    return None

def run(root):
    f=json.loads((root/FEATURE).read_text())
    d=json.loads((root/FY).read_text())
    geo=json.loads((root/GEO).read_text())
    if f["schema"]!="gesui_NWCA2021_original128geoNear_registered_POTW2026_permitFEATURE_hash_source_only_v01":
        raise ValueError("Original exact-feature map not source frozen")
    if d["schema"]!="gesui_NWCA2021_160_originalPOTW_permit_FY2021_FY2022_calendar2021_DMR_nutrient_SOURCE_FEASIBILITY_v01":
        raise ValueError("Original EPA FY2021+FY2022 DMR source missing")
    if d["source_geo_568_NWCA_registered_outfall_SHA256"]!=sha(root/GEO):
        raise ValueError("2021 DMR permit source differs from original near/far geometry")
    if f["source_NWCA_original568_POTW_geo_result_SHA256"]!=sha(root/GEO):
        raise ValueError("Exact permitted feature coordinates not same 2026 NPDES source")
    fixed=f["original_source_case_site_to_exact_permit_feature_sha256_and_distance"]
    if f["original_128_near10_wetland_site_location_count"]!=128 or f["original_near10_exact_permit_feature_ID_count"]!=344 or len(fixed)!=392:
        raise ValueError("Original fixed 128 NWCA near sites, 344 true feature IDs, 392 site×feature pairs changed")
    ids={r["NPDES_exact_permit_feature_SHA256"] for r in fixed}
    sites={r["NWCA2021_site_visit_SHA256"] for r in fixed}
    original_sites=geo["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]
    near={r["NWCA2021_source_site_visit_SHA256"] for r in original_sites if r["outfall_rows_within_10km"]>0}
    far={r["NWCA2021_source_site_visit_SHA256"] for r in original_sites if r["nearest_outfall_geodesic_km_within25km_else_null"] is None or r["nearest_outfall_geodesic_km_within25km_else_null"]>=20}
    nearrow=[r for r in original_sites if r["NWCA2021_source_site_visit_SHA256"] in near]
    farrow=[r for r in original_sites if r["NWCA2021_source_site_visit_SHA256"] in far]
    region={r["original_state"]+"|"+r["original_ecoregion"] for r in nearrow}&{r["original_state"]+"|"+r["original_ecoregion"] for r in farrow}
    original240_near={r["NWCA2021_source_site_visit_SHA256"] for r in nearrow if r["original_state"]+"|"+r["original_ecoregion"] in region}
    if len(original240_near)!=103 or sites!=near:
        raise ValueError("Original near/far prebiology state×ecoregion study population changed")
    # DMR calendar/fiscal joined JSON currently stores per permit monthly
    # count but NOT every feature; read source original fiscal artifacts via
    # original combined if detailed source member available.
    if "source_original_2021_DMR_selected_values_no_raw_biology" in d:
        records=d["source_original_2021_DMR_selected_values_no_raw_biology"]
    elif "source_original_2021_full_source_selected_DMR_rows_no_raw_reported_values" in d:
        records=d["source_original_2021_full_source_selected_DMR_rows_no_raw_reported_values"]
    else:
        raise ValueError("Combined 160-permit FY2021+FY2022 JSON missing exact original feature-level DMR receipts; no inference from permit-level aggregates")
    feature=defaultdict(lambda:{"N":set(),"P":set(),"FLOW":set()})
    for r in records:
        if not r["strict_external_outfall_effluent_numeric_uncensored"]:continue
        key=r["source_permit_feature_ID_SHA256"]
        if key not in ids:continue
        var={"00600":"N","00665":"P","50050":"FLOW"}[r["original_parameter_code"]]
        feature[key][var].add(r["original_monitoring_period_end_date"])
    keys={
        "eitherN_orP":{k for k,v in feature.items() if v["N"] or v["P"]},
        "bothNandP":{k for k,v in feature.items() if v["N"] and v["P"]},
        "N":{k for k,v in feature.items() if v["N"]},
        "P":{k for k,v in feature.items() if v["P"]},
        "FLOW":{k for k,v in feature.items() if v["FLOW"]},
        "NandP_same_report_period":{k for k,v in feature.items() if v["N"] & v["P"]},
    }
    assignments={}
    for name,matched in keys.items():
        eligible=[r for r in fixed if r["NPDES_exact_permit_feature_SHA256"] in matched]
        dat={
           "exact_nearby_source_feature_ids_with_2021_reported_numeric":len(matched),
           "matched_source_NEAR_wetlands_128":len({r["NWCA2021_site_visit_SHA256"] for r in eligible}),
           "matched_original_27state_study_NEAR_wetlands_103":len({r["NWCA2021_site_visit_SHA256"] for r in eligible if r["NWCA2021_site_visit_SHA256"] in original240_near}),
           "original_2021_source_administrative_permit_IDs":len({r["NPDES_registered_permit_SHA256"] for r in eligible}),
           "original_states_with_near_registered_outfall_and_2021_DMR":len({r["source_original_US_state"] for r in eligible}),
           "original_distinct_calendar2021_quarters":sorted({quarter(z) for fkey in matched for ty in feature[fkey] for z in feature[fkey][ty] if quarter(z) is not None})
        }
        assignments[name]=dat
    q=assignments["eitherN_orP"]
    flowN=q["matched_source_NEAR_wetlands_128"] and sum(
       1 for r in fixed if r["NPDES_exact_permit_feature_SHA256"] in keys["eitherN_orP"] & keys["FLOW"])
    # Use distinct measured exact features, not repeated site-feature pairs.
    coflow=len(keys["eitherN_orP"]&keys["FLOW"])
    pass_measurement=(
        q["matched_source_NEAR_wetlands_128"]>=20 and
        q["original_2021_source_administrative_permit_IDs"]>=10 and
        q["original_states_with_near_registered_outfall_and_2021_DMR"]>=6 and
        len(q["original_distinct_calendar2021_quarters"])>=2 and
        coflow/max(1,len(keys["eitherN_orP"]))>=.5
    )
    return {
        "schema":"gesui_NWCA2021_2026_registered_128wetland_344outfall_FEATURE_exact_actual_CY2021_TNTP_DMR_QC_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_pre2021_DMR_exact344geo_features_SHA256":sha(root/FEATURE),
        "original_FY2021andFY2022_source_DMR_join_SHA256":sha(root/FY),
        "original_2021_POTW_source_selection_geo_SHA256":sha(root/GEO),
        "source_original_near_wetland_128_registered_outfall_feature344_pairs392":True,
        "source_reported_2021_FINAL_EFFLUENT_numeric_totalN00600_totalP00665_flow50050":True,
        "exact_original_outfall_GEOMETRY_to_2021_DMR_feature_matches":assignments,
        "original_reported_NP_outfall_features_with_numeric_flow_at_same_feature":coflow,
        "source_precommitted_2021_nutrient_monitoring_feature_opportunity_gate":{
            "min_near_wetlands_with_actual_2021_nutrient_DMR_at_exact_outfall_number":20,
            "min_original_independent_registered_municipal_sewage_permit_ids":10,
            "min_original_states":6,
            "min_original_calendar2021_reported_quarters":2,
            "min_fraction_of_same_outfall_features_with_flow":.5,
            "eligible_to_proceed_to_separate_NHDPlus_downstream_flow_route":bool(pass_measurement)
        },
        "2026_ECHO_permitted_feature_position_and_DMR2021_at_matching_feature_does_NOT_verify_receiving_wetland_flow":True,
        "NO_2021_NWCA_plant_diversity_or_wetland_tn_tp_biology_read":True,
        "no_wastewater_driven_biodiversity_coefficient_estimated":True
    }

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();root=Path(__file__).resolve().parents[1]
    r=run(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"exact_source_measured":r["exact_original_outfall_GEOMETRY_to_2021_DMR_feature_matches"],
                      "next_hydrology_gate":r["source_precommitted_2021_nutrient_monitoring_feature_opportunity_gate"]},indent=2))
if __name__=="__main__":main()
