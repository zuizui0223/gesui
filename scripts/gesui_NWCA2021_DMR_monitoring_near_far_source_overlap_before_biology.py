#!/usr/bin/env python3
"""SOURCE-ONLY state×ecoregion overlap of historical 2021 DMR monitored
municipal EXO features and NWCA2021 wetland source points.

No species/plant Shannon/FQAI nor NWCA water TN/TP observations read. Determine
whether DMR-verified (2021 CALENDAR measured TN/TP at SAME source EXO feature)
wetlands have enough strict >=20km potential comparison wetlands within the
same STATE×NWCA ECO region to support a macro-scale NEW biological sensitivity.
All 27state 240-site results remain unchanged.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

DIR=Path("results/us_NWCA2021_plant_TNTP_20261110")
DMR=DIR/"NWCA2021_160POTW_ORIGINAL_CAL2021_DMR_TN_TP_EXO_FEATURE_SOURCE_GATE.json"
GEO=DIR/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
AUD=DIR/"NWCA2021_EPA_2021DMR_160POTW_344EXO_SOURCE_OFFLINE_INDEPENDENT_QA.json"
SOURCE_LIMIT={"near_BOTH":20,"far_same_strata":40,"states":10,"state_eco_cells":10}
ANY_LIMIT={"near_N_or_P":40,"far_same_strata":60,"states":12,"state_eco_cells":15}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(root):
    d=json.loads((root/DMR).read_text())
    g=json.loads((root/GEO).read_text())
    a=json.loads((root/AUD).read_text())
    if not a["original_predeclared_biology_DMR_source_availability_gate_independently_reproduced"]:
        raise ValueError("Prior official 2021 historical DMR near-source gate not independently audited")
    if a["source_SHA256"]["DMR_2021_2022_join"]!=sha(root/DMR):
        raise ValueError("Current DMR source cohort differs from independently audited source")
    if d["original_568site_2026ECHO_STRICT_SEWAGE_geographic_PREOUTCOME_SHA256"]!=sha(root/GEO):
        raise ValueError("Source original 2026 ECHO site geometry changed")
    original=g["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]
    by={x["NWCA2021_source_site_visit_SHA256"]:x for x in original}
    if len(original)!=568 or len(by)!=568:raise ValueError("Original 568 source sites changed")
    near={sid for sid,x in by.items() if x["outfall_rows_within_10km"]>0}
    far={sid for sid,x in by.items() if x["nearest_outfall_geodesic_km_within25km_else_null"] is None
                                      or x["nearest_outfall_geodesic_km_within25km_else_null"]>=20}
    if len(near)!=128 or len(far)!=311 or near&far:
        raise ValueError("Original strict source prior 10km near / 20km far incorrect")
    ids={
      "DMR_actual2021_TN_or_TP":set(d["original_2021_DMR_sites_with_valid_TN_OR_TP_same_source_feature"]["source_original_site_visit_SHA256s"]),
      "DMR_actual2021_BOTH_TN_TP_same_EXO":set(d["original_2021_DMR_sites_with_valid_BOTH_TN_and_TP_same_source_feature"]["source_original_site_visit_SHA256s"]),
    }
    if len(ids["DMR_actual2021_TN_or_TP"])!=77 or len(ids["DMR_actual2021_BOTH_TN_TP_same_EXO"])!=26:
        raise ValueError("Original DMR monitored source near site populations changed")
    out={}
    for label,nearall in ids.items():
        if not nearall.issubset(near):
            raise ValueError("2021 exact DMR monitoring site not in original NPDES permitted proximity source")
        strata={(by[sid]["original_state"],by[sid]["original_ecoregion"]) for sid in nearall}
        far_comp={sid for sid in far if (by[sid]["original_state"],by[sid]["original_ecoregion"]) in strata}
        near_comp={sid for sid in nearall if (by[sid]["original_state"],by[sid]["original_ecoregion"]) in {
            (by[s]["original_state"],by[s]["original_ecoregion"]) for s in far_comp}}
        pairs={(by[sid]["original_state"],by[sid]["original_ecoregion"]) for sid in near_comp}
        # Drop near-only strata with no source-far comparator and count only
        # original state/eco cells with both DMR-source and fixed far.
        far_comp={sid for sid in far_comp if (by[sid]["original_state"],by[sid]["original_ecoregion"]) in pairs}
        stats=Counter()
        for sid in near_comp:stats[(by[sid]["original_state"],by[sid]["original_ecoregion"],"near")]+=1
        for sid in far_comp:stats[(by[sid]["original_state"],by[sid]["original_ecoregion"],"far")]+=1
        info={
            "original_source_candidate_near_sites":len(nearall),
            "new_2021_DMR_reported_source_near_original_NWCA_wetlands_with_same_region_far":len(near_comp),
            "original_2021_NWCA_source_far20_control_sites_same_state_ecoregion":len(far_comp),
            "states_with_NEAR_AND_FAR":len({state for state,eco in pairs}),
            "state_ecoregion_cells_with_NEAR_AND_FAR":len(pairs),
            "original_NWCA_source_near_state_counts":dict(Counter(by[sid]["original_state"] for sid in near_comp)),
            "original_NWCA_source_far_state_counts":dict(Counter(by[sid]["original_state"] for sid in far_comp)),
            "original_matched_state_ecoregion_strata":{
                state+"|"+eco:{"near":stats[(state,eco,"near")],"far":stats[(state,eco,"far")]}
                for state,eco in sorted(pairs)},
            "original_NWCA2021_NEAR_source_site_SHA256s":sorted(near_comp),
            "original_NWCA2021_FAR_source_site_SHA256s":sorted(far_comp),
        }
        if label=="DMR_actual2021_BOTH_TN_TP_same_EXO":
            limits=SOURCE_LIMIT
            pass_gate=(len(near_comp)>=limits["near_BOTH"] and len(far_comp)>=limits["far_same_strata"] and
                info["states_with_NEAR_AND_FAR"]>=limits["states"] and
                len(pairs)>=limits["state_eco_cells"])
        else:
            limits=ANY_LIMIT
            pass_gate=(len(near_comp)>=limits["near_N_or_P"] and len(far_comp)>=limits["far_same_strata"] and
                info["states_with_NEAR_AND_FAR"]>=limits["states"] and
                len(pairs)>=limits["state_eco_cells"])
        info["source_biological_coverage_thresholds_FIXED_before_new_HALL"]=limits
        info["may_descriptively_fit_2021_HALL_with_monitored_POTW_proximity_ONLY_NONCAUSAL"]=pass_gate
        out[label]=info
    return {
        "schema":"gesui_US_NWCA2021_monitored_DMR_2021_POTW_568site_geographic_sameStateEco_source_before_plant_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "source_official_2021_TNTP_DMR_GEO_SHA256":sha(root/DMR),
        "source_official_2026_POTW_568site_GEO_SHA256":sha(root/GEO),
        "source_2021_DMR_exact_feature_independent_QA_SHA256":sha(root/AUD),
        "original_source_DMR_2021_ANY_NP_and_BOTH_NP_strict_same_state_ecoregion_controls":out,
        "SOURCE_ONLY_no_2021_NWCA_plant_H_ALL_H_NAT_FQAI_or_measured_wetland_water_values_reopened":True,
        "geodesic_facility_nearness_is_not_NHDPlus_directed_WWTP_effluent_delivery":True,
        "source_monitoring_opportunity_test_not_causal_EPA_2021_outfall_loading_experiment":True
    }

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    o=run(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:{n:v[n] for n in ("original_source_candidate_near_sites",
         "new_2021_DMR_reported_source_near_original_NWCA_wetlands_with_same_region_far",
         "original_2021_NWCA_source_far20_control_sites_same_state_ecoregion",
         "states_with_NEAR_AND_FAR","state_ecoregion_cells_with_NEAR_AND_FAR",
         "may_descriptively_fit_2021_HALL_with_monitored_POTW_proximity_ONLY_NONCAUSAL")}
       for k,v in o["original_source_DMR_2021_ANY_NP_and_BOTH_NP_strict_same_state_ecoregion_controls"].items()},indent=2))

if __name__=="__main__":main()
