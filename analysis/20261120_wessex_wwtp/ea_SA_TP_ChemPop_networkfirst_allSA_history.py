#!/usr/bin/env python3
"""All-national-SA network-first historical TP source-quality falsifier.

A prerequisite of causal wastewater source attribution, NOT an outcome fit.
Unlike 2km-geometric-preselected STWs, select 2026-registered EA native SA
treated final water-company sample points by ORIGINAL HydroRIVERS NEXT_DOWN
directed SOURCE / reverse placebo / same-reach relation ≤50km to each
original 2017/18 score-eligible ChemPop site. THEN query actual 2017/18
source 0348 total P as P mg/L using the officially validated comma-joined
pointNotation. No source-specific annual mass, no BMWP outcome, no
downstream hydraulic-delivery causality claim.

Missing small streams, STW sampling point!=official outfall, no historical
sample, and no source outflow/flow measurements are STOP/unknowns not zeroes.
Only numeric aggregate QA, hashes and site/year numbers emitted.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import sys
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis"/"20261120_wessex_wwtp"))
import ea_national_SA_registry_to_ChemPop_source_ceiling as nat
import ea_SA_TP_to_ChemPop_location_sourcegate as local
import ea_SA_TP_ChemPop_HydroRIVERS_directed_route_gate as route
import ea_WQE_national_SA_actual2017_18_TP_chemistry_bio_screen as hist

YEARS=(2017,2018)
EXPECTED_N_SA=5722
EXPECTED_BIO_YEAR={2017:389,2018:366}
EXPECTED_BIO_UNION=463
N_DISTINCT_BIO_GATE=50
N_SITEYEAR_GATE=100
N_HYDROBAS_BLOCKS_GATE=5

def hashids(ids):
    return hashlib.sha256("\n".join(sorted(map(str,ids))).encode()).hexdigest()

def relations(source_snaps,bio_snaps,year_bio,graph):
    """Report all source geometries with either-direction matched reach ≤50km.

    Current registry source selection is based ONLY on geometry, year-specific
    ecology sample PRESENCE, and river network, never BMWP or ecology values.
    """
    def paths_for(snaps):
        return {
            r: route.downstream_trace(r,graph)[0]
            for r in set(x["reach"] for x in snaps.values())
        }
    src_paths=paths_for(source_snaps)
    bio_paths=paths_for(bio_snaps)
    candidates=set()
    all_rels={}
    pre_stats={}
    for yr in YEARS:
        bio_ids=sorted(year_bio[yr])
        if len(bio_ids)!=EXPECTED_BIO_YEAR[yr]:
            raise RuntimeError("Genuine original biology annual site source count drift")
        rel_by_src={}
        missing_bio=sum(i not in bio_snaps for i in bio_ids)
        for source_id,s in source_snaps.items():
            sr=s["reach"]
            rel={"upstream":set(),"downstream_placebo":set(),"same_reach_unknown":set()}
            dsp=src_paths[sr]
            for bio_id in bio_ids:
                b=bio_snaps.get(bio_id)
                if b is None:continue
                br=b["reach"]
                if sr==br:
                    rel["same_reach_unknown"].add(bio_id)
                elif br in dsp:
                    rel["upstream"].add(bio_id)
                elif sr in bio_paths[br]:
                    rel["downstream_placebo"].add(bio_id)
            if any(rel.values()):
                rel_by_src[source_id]=rel
                candidates.add(source_id)
        pre_stats[str(yr)]={
            "original_eligible_bio_site_n":len(bio_ids),
            "bio_not_snapped_to_HydroRIVERS_500m":missing_bio,
            "national_current_SA_registered_sites_WITH_A_or_B_or_same_reach_relation":len(rel_by_src),
            "national_current_SA_with_upstream_geometry_path":
                sum(bool(v["upstream"]) for v in rel_by_src.values()),
            "national_current_SA_with_downstream_reverse_placebo_path":
                sum(bool(v["downstream_placebo"]) for v in rel_by_src.values()),
            "national_current_SA_with_same_reach_order_unknown":
                sum(bool(v["same_reach_unknown"]) for v in rel_by_src.values()),
            "bio_with_some_current_registered_STW_on_potential_upstream_reach":
                len(set().union(*(v["upstream"] for v in rel_by_src.values()))) if rel_by_src else 0,
            "bio_with_some_current_registered_STW_on_reverse_placebo_reach":
                len(set().union(*(v["downstream_placebo"] for v in rel_by_src.values()))) if rel_by_src else 0,
            "current_source_registry_not_historically_measured":True
        }
        all_rels[yr]=rel_by_src
    return candidates,all_rels,pre_stats

def validate_measured(measured,rel,graph,biosnap,year_bio):
    results={}
    union=set()
    siteyear_count=0
    upstream_stw_any=set()
    all_blocks=set()
    all_components=set()
    for yr in YEARS:
        src_codes=set(measured[yr])
        rel_year=rel[yr]
        unexpected=src_codes-set(rel_year)
        if unexpected:
            raise RuntimeError("Original WQE native source unexpectedly outside topology-selected candidate")
        groups={k:{} for k in ("upstream","downstream_placebo","same_reach_unknown")}
        n_pairs=collections.Counter()
        for s_id in sorted(src_codes):
            for kind,ids in rel_year[s_id].items():
                groups[kind][s_id]=ids
                n_pairs[kind]+=len(ids)
        sets={
            label:set().union(*v.values()) if v else set()
            for label,v in groups.items()
        }
        hit=sets["upstream"]
        union|=hit
        siteyear_count+=len(hit)
        upstream_stw_any|={s for s,v in groups["upstream"].items() if v}
        for bio_id in hit:
            r=biosnap[bio_id]["reach"]
            all_blocks.add(str(graph[r][3]))
            all_components.add(str(graph[r][2]))
        per_bio=collections.Counter()
        for stw,sites in groups["upstream"].items():
            for bio in sites:per_bio[bio]+=1
        results[str(yr)]={
            "native_WaterCompany_final_TP_samplepoint_same_year_n":len(src_codes),
            "n_routed_source_samplepoints_with_network_upstream_BIO":
                sum(bool(x) for x in groups["upstream"].values()),
            "n_distinct_source_eligible_BIO_with_upstream_SA_TP_on_network_within50km":len(hit),
            "n_distinct_source_eligible_BIO_with_reverse_downstream_placebo_STW_TP":
                len(sets["downstream_placebo"]),
            "n_distinct_source_eligible_BIO_same_reach_TP_STW_order_unknown":
                len(sets["same_reach_unknown"]),
            "n_upstream_STW_BIO_network_pairs":n_pairs["upstream"],
            "n_reverse_STW_BIO_network_pairs":n_pairs["downstream_placebo"],
            "n_same_reach_STW_BIO_pairs_order_unknown":n_pairs["same_reach_unknown"],
            "n_BIO_sites_with_two_or_more_upstream_screen_STW":sum(x>=2 for x in per_bio.values()),
            "no_actual_outfall_position_or_WWTP_P_mass_kg_year_verified":True,
            "network_candidate_2017_18_SA_codes_SHA256":hashids(src_codes),
        }
    return results,{
        "n_DISTINCT_true_BIO_site_union_network_upstream_measured_source":len(union),
        "n_SITE_YEAR_bio_network_upstream_measured_source":siteyear_count,
        "n_DISTINCT_current_SA_with_measured_TP_and_at_least_one_potential_upstream_BIO":
            len(upstream_stw_any),
        "n_distinct_HYBAS_L12_numeric_group_keys_screened":len(all_blocks),
        "n_distinct_MAIN_RIV_group_keys_screened":len(all_components),
        "source_bio_union_SHA256":hashids(union),
        "source_geographic_gate_distinct_BIO_ge_50":len(union)>=N_DISTINCT_BIO_GATE,
        "source_geographic_gate_SITE_YEAR_ge_100":siteyear_count>=N_SITEYEAR_GATE,
        "source_geographic_gate_HYBAS_L12_ge_5":
            len(all_blocks)>=N_HYDROBAS_BLOCKS_GATE,
        "original_SA_sampling_point_NOT_permitted_outfall":True,
        "river_direction_screen_not_verified_hydraulic_pollutant_delivery":True,
        "cannot_authorize_ecological_causal_regression":True
    }

def main(args):
    out={
        "schema":"GESUI_20261120_UK_WQE_ALL_SA_NETWORK_FIRST_ACTUAL2017_18_TP_BIO_SOURCE_V1",
        "posthoc_original_source_quality_only":True,
        "status":"SOURCE_STOP",
        "previous_2km_screen_falsification":"The 2km prefilter excludes 2–50km upstream sources; historical coverage must be re-evaluated on all SA network relation candidates",
        "no_BMWP_response_opened":True,
        "actual_TP_mgL_source_not_WWTP_annual_P_mass":True,
        "true_permitted_outfall_vs_WQE_sampling_point_not_verified":True,
        "river_10km2_headwater_omission_source_missingness":True
    }
    try:
        reg,regqa=nat.collect()
        if len(reg)!=EXPECTED_N_SA:
            raise RuntimeError("Original nationwide strict SA sampling point source drift")
        bio,yr,biqa=local.chempop_geo()
        if len(yr[2017]|yr[2018])!=EXPECTED_BIO_UNION:
            raise RuntimeError("Original 463 distinct sourceeligible 2017/18 BIO drift")
        with tempfile.TemporaryDirectory(prefix="gesui_uk_allSA_network_") as temp:
            rivers,graph,networkqa=route.get_rivers(Path(temp))
            bio_org={k:bio[k] for k in yr[2017]|yr[2018]}
            b_snap,bqa=route.snap(bio_org,rivers,"original ChemPop 2017/18 BIO source")
            s_snap,sqa=route.snap(reg,rivers,"all current EA WaterCompany SA source")
            candidate,rel,prestats=relations(s_snap,b_snap,yr,graph)
            if len(candidate)==0:
                raise RuntimeError("No original SA source candidates within 50km network or same reach")
            if len(candidate)>len(reg):
                raise RuntimeError("Network prefilter larger than national registered SA pool")
            # Never select by BMWP or original EDF/dose; only genuine network
            # and current official SA type + source-qualified BIO monitoring
            # presence in either biological survey year.
            actual,actualqa=hist.candidate_history(sorted(candidate))
            by_year,overall=validate_measured(actual,rel,graph,b_snap,yr)
        out.update({
            "source_registry":{k:v for k,v in regqa.items() if k not in ("source_page_SHA256",)},
            "original_ChemPop_source_frame":biqa,
            "official_HydroRIVERS_v1_source_frame":networkqa,
            "BIO_snap_QA":bqa,"all_SA_snap_QA":sqa,
            "before_2017_18_observed_TP_network_relation":prestats,
            "network_source_candidate_total_n":len(candidate),
            "network_source_candidate_code_SHA256":hashids(candidate),
            "historical_TP_source":{k:v for k,v in actualqa.items()
                                   if k not in ("source_page_SHA256",)},
            "same_year_original_TP_measured_source_vs_BIO_directed_relation":by_year,
            "network_matched_historical_source_gates":overall,
            "status":"NETWORK_FIRST_2017_18_GENUINE_TP_STW_BIO_SOURCE_SCREEN_COMPLETE_NOT_CAUSAL"
        })
    except Exception as e:
        out.update({"status":"SOURCE_OR_HYDROLOGY_GATE_STOP",
                    "error":type(e).__name__+": "+str(e)[:730]})
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_EA_ALL_SA_NETWORK_FIRST_TP_BIO_SOURCE="+json.dumps({
        "status":out["status"],"error":out.get("error"),
        "BIO_snap":out.get("BIO_snap_QA"),"all_SA_snap":out.get("all_SA_snap_QA"),
        "network_candidates":out.get("network_source_candidate_total_n"),
        "pre_measured_source":out.get("before_2017_18_observed_TP_network_relation"),
        "historical_TP_source":out.get("historical_TP_source"),
        "year_relations":out.get("same_year_original_TP_measured_source_vs_BIO_directed_relation"),
        "gates":out.get("network_matched_historical_source_gates"),
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
