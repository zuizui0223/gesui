#!/usr/bin/env python3
"""Independent offline source-gate QA from original FOUR NLDI full response shards.

No calls to original indexing/routing module and no biological outcome file.
Rebuild source cohort 71 original wetlands/86 independently 2021-monitored POTW
outfall features and strict BOTH 200m direct stream snaps with distinct COMIDs.
Independently validate original 2026 ECHO and FY2021/22 2021 DMR SHA receipts.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter
from pathlib import Path
D=Path("results/us_NWCA2021_plant_TNTP_20261110")
SOURCE=D/"EPA_NWCA2021_128WETLAND_160POTW_EXACT_NEAR_OUTFALL_FEATURE_SOURCE.json"
DMR=D/"EPA_NWCA2021_ORIGINAL_160POTW_2021_DMR_TNTP_SOURCE_FEASIBILITY.json"
MAIN=D/"EPA_NWCA2021_POTW_2021_NP71WET_86OUTFALL_USGS_NLDI_STRICT200M_DIRECTION_SOURCE_GATE.json"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(root,shards):
    source=json.loads((root/SOURCE).read_text())
    dmr=json.loads((root/DMR).read_text())
    macro=json.loads((root/MAIN).read_text())
    if macro["source_original_344_EPA_registered_POTW_feature_geo_SHA256"]!=sha(root/SOURCE):
        raise ValueError("Source 344 outfall SHA mismatch")
    if macro["source_original_2021_calendaryear_NP_DMR_status_SHA256"]!=sha(root/DMR):
        raise ValueError("Source 2021 DMR actual reporting SHA mismatch")
    results={}
    original_shas={}
    original_counts=Counter()
    for p in shards:
        j=json.loads(p.read_text())
        shard=str(j["original_source_shard_id"])
        if shard in original_shas or shard not in ("0","1","2","3"):
            raise ValueError("Original USGS NLDI source shards duplicated")
        if sha(p)!=macro["source_original_USGS_NLDI_index_shards_SHA256"][shard]:
            raise ValueError("Original USGS source index shard SHA no longer matches executed source gate")
        if j["source_original_exact344_geo_SHA256"]!=sha(root/SOURCE) or j["source_original_2021_DMR160_FY2021FY2022_join_SHA256"]!=sha(root/DMR):
            raise ValueError("USGS original indexing used source-inconsistent outfalls")
        rows=j["source_original_USGS_point_source_replies"]
        if len(rows) not in (39,40):raise ValueError("Four source shards must partition 157 points")
        original_shas[shard]=sha(p)
        for r in rows:
            key=r["source_hashed_point_id"]
            if key in results:raise ValueError("Original source indexed duplicate site or outfall")
            if r["original_source_type"] not in ("NWCA_wetland","POTW_outfall"):
                raise ValueError("Unknown original source NLDI point type")
            d=r["original_distance_to_hydrolocation_m"]
            strict=(r["official_USGS_NLDI_catchment_comid"] is not None and
                    r["official_USGS_NLDI_catchment_comid"]==r["official_USGS_NLDI_hydrolocation_comid"] and
                    d is not None and 0<=d<=200.)
            sens=(r["official_USGS_NLDI_hydrolocation_comid"] is not None and
                  d is not None and 0<=d<=500.)
            if strict!=r["source_strict_200m_DIRECT_stream_candidate"] or sens!=r["source_500m_modelled_sensitivity_only"]:
                raise ValueError("Original USGS 200m/500m index classification was not reproducible")
            if (r["NLDI_position_HTTP_status"]==200 and not r["official_NLDI_position_SHA256"]) or (
                r["NLDI_hydrolocation_HTTP_status"]==200 and not r["official_NLDI_hydrolocation_SHA256"]):
                raise ValueError("Official NLDI original successful response SHA missing")
            results[key]=r
            original_counts[r["original_source_type"]]+=1
    response_QC={
        kind:{
            "NLDI_catchment_position_HTTP_status":dict(Counter(
                str(r["NLDI_position_HTTP_status"]) for r in results.values() if r["original_source_type"]==kind)),
            "NLDI_stream_hydrolocation_HTTP_status":dict(Counter(
                str(r["NLDI_hydrolocation_HTTP_status"]) for r in results.values() if r["original_source_type"]==kind)),
            "direct200m_success_n":sum(r["source_strict_200m_DIRECT_stream_candidate"]
                                      for r in results.values() if r["original_source_type"]==kind),
            "500m_modelled_position_n":sum(r["source_500m_modelled_sensitivity_only"]
                                         for r in results.values() if r["original_source_type"]==kind)
        }
        for kind in ("NWCA_wetland","POTW_outfall")
    }
    if len(results)!=157 or original_counts!=Counter({"NWCA_wetland":71,"POTW_outfall":86}):
        raise ValueError("Independent original 71/86 source cohort not reproduced")
    measured={x["NPDES_exact_permit_feature_SHA256"]
              for x in dmr["source_original_2021_exact_permit_feature_DMR_reporting_summary"]
              if x["reported_original_numeric_2021_TN00600_period_end_dates"] or
                 x["reported_original_numeric_2021_TP00665_period_end_dates"]}
    candidates=[p for p in source["original_source_case_site_to_exact_permit_feature_sha256_and_distance"]
                if p["NPDES_exact_permit_feature_SHA256"] in measured]
    if set(x["NWCA2021_site_visit_SHA256"] for x in candidates)!={k for k,v in results.items() if v["original_source_type"]=="NWCA_wetland"}:
        raise ValueError("Original measured wetland cohort changed")
    if set(x["NPDES_exact_permit_feature_SHA256"] for x in candidates)!={k for k,v in results.items() if v["original_source_type"]=="POTW_outfall"}:
        raise ValueError("Original measured POTW cohort changed")
    all200=[]
    same200=[]
    model500=[]
    for x in candidates:
        s=results[x["NPDES_exact_permit_feature_SHA256"]]
        r=results[x["NWCA2021_site_visit_SHA256"]]
        if s["source_500m_modelled_sensitivity_only"] and r["source_500m_modelled_sensitivity_only"]:
            model500.append(x)
        if not(s["source_strict_200m_DIRECT_stream_candidate"] and r["source_strict_200m_DIRECT_stream_candidate"]):continue
        if s["official_USGS_NLDI_hydrolocation_comid"]==r["official_USGS_NLDI_hydrolocation_comid"]:
            same200.append(x)
        else:
            all200.append(x)
    counts={
        "source_original_157_NLDI_indexed_points":len(results),
        "NLDI_200m_strict_wetland_source_n":sum(v["source_strict_200m_DIRECT_stream_candidate"] for v in results.values() if v["original_source_type"]=="NWCA_wetland"),
        "NLDI_200m_strict_POTW_source_n":sum(v["source_strict_200m_DIRECT_stream_candidate"] for v in results.values() if v["original_source_type"]=="POTW_outfall"),
        "source_strict200_same_comid_not_directed_pairs":len(same200),
        "source_strict200_DIFFERENT_comid_pairs":len(all200),
        "source_strict200_DIFFERENT_comid_wetlands":len({x["NWCA2021_site_visit_SHA256"] for x in all200}),
        "source_strict200_DIFFERENT_comid_POTW_features":len({x["NPDES_exact_permit_feature_SHA256"] for x in all200}),
        "source_model500m_not_directed_pairs":len(model500),
        "strict_NLDI_candidates_original_states":len({x["source_original_US_state"] for x in all200}),
    }
    old=macro["original_2021_DMR_source_2026_registered_geo_reference_counts"]
    for new_key,prev_key in (
        ("NLDI_200m_strict_wetland_source_n","USGS_200m_direct_candidate_wetlands"),
        ("NLDI_200m_strict_POTW_source_n","USGS_200m_direct_candidate_POTW"),
        ("source_strict200_same_comid_not_directed_pairs","source_site_feature_both200m_same_comid_direction_unknown_n"),
        ("source_strict200_DIFFERENT_comid_pairs","source_site_feature_both200m_DIFFERENT_comid_to_test_downstream_n"),
        ("source_strict200_DIFFERENT_comid_wetlands","source_200m_candidate_sites_different_comids_n"),
        ("source_strict200_DIFFERENT_comid_POTW_features","source_200m_candidate_outfalls_different_comids_n"),
        ("source_model500m_not_directed_pairs","source_site_feature_both500m_modelled_not_strict_n"),
        ("strict_NLDI_candidates_original_states","source_200m_candidate_states_different_comids_n")
    ):
        if counts[new_key]!=old[prev_key]:
            raise ValueError("Original source strict NLDI recorded count mismatch: "+new_key)
    directed=macro["USGS_strict200m_directed_outfall_to_wetland_potential_matches"]
    downstream=[x for x in directed if x["DM_directed_different_comid"]]
    navigated=bool(macro["USGS_downstream_DM_navigation_performed"])
    orig_gate=macro["precommitted_original_actual_DM_FLOW_SOURCE_STUDY_GATE"]
    if navigated:
        if len(downstream)!=macro["USGS_downstream_DM_50km_outcome_source"]["potential_directed_source_pairs"]:
            raise ValueError("Original strict downstream match count no longer agrees with source outcome")
        if len({x["site_sha256"] for x in downstream})!=macro["USGS_downstream_DM_50km_outcome_source"]["potential_directed_unique_NWCA_wetlands"]:
            raise ValueError("Original directed downstream wetland coverage differs")
        passed=(len({x["site_sha256"] for x in downstream})>=orig_gate["wetland_sites"] and
                len({x["feature_sha256"] for x in downstream})>=orig_gate["source_features"] and
                len({x["permit_sha256"] for x in downstream})>=orig_gate["permits"] and
                len({x["original_state"] for x in downstream})>=orig_gate["states"] and
                old["source_original_same_state_ecoregion_far_controls_n"]>=orig_gate["original_far_controls"])
    else:
        if directed or macro["USGS_downstream_DM_50km_outcome_source"]["potential_directed_source_pairs"] is not None:
            raise ValueError("NLDI navigation not performed must be represented as untested, NOT a 0 downstream result")
        passed=False
        if counts["source_strict200_DIFFERENT_comid_wetlands"]>=orig_gate["wetland_sites"] and counts["source_strict200_DIFFERENT_comid_POTW_features"]>=orig_gate["source_features"]:
            raise ValueError("NLDI navigation source prefilter was sufficient but was skipped")
    if passed!=orig_gate["sufficient_USGS_NLDI_source_supported_NWCA2021_wetlands_for_future_observational_vegetation_test"]:
        raise ValueError("Original plant-allowed NLDI source gate independently disagrees")
    return {
        "schema":"gesui_NWCA2021_71wetland_86_POTW_USGS_NLDI_NHDPlusV2_strict200m_independent_source_QA_v01",
        "source_SHA256":{"original2021_DMR_permit_and_feature":sha(root/DMR),
                         "original2026_registered_feature_geometry":sha(root/SOURCE),
                         "originalUSGS_NLDI_river_strict_source_gate":sha(root/MAIN),
                         "original4_NLDI_index_shards":original_shas},
        "independent_157_source_point_snap_and_all_2021_nutrient_DMR_eligible_pairs_QA":counts,
        "independent_USGS_NLDI_API_status_by_original_source_type":response_QC,
        "original_USGS_NLDI_navigation_performed":navigated,
        "independent_2021_measured_POTW_to_wetland_50km_directed_candidates":len(downstream) if navigated else None,
        "no_navigated_data_never_converted_to_false_ZERO_outfall_exposure":True,
        "independent_original_study_source_gate_PASS":passed,
        "NO_actual_NWCA2021_plant_biodiversity_response_ever_read":True,
        "NO_USGS_NLDI_comid_or_river_200m_means_confirmed_WWTP_effluent_delivery":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--shard",type=Path,action="append",required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if len(a.shard)!=4:raise ValueError("Exactly four original NHDPlus source shards required for independent QA")
    root=Path(__file__).resolve().parents[1]
    value=run(root,a.shard)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(value,sort_keys=True,indent=2)+"\n")
    print(json.dumps(value,indent=2))

if __name__=="__main__":main()
