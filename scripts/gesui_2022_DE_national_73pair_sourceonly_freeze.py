#!/usr/bin/env python3
"""Freeze Germany national WISE6 directed river-pair SOURCE frame before plant/butterfly data.

All inputs are original national 2022 EEA/WISE6 metadata already committed,
never using measured TP concentrations or iNaturalist butterfly photo counts.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

DIR="results/wise6_2022_national_DEAT_TP_sourceonly/"
PAIRS=DIR+"NATIONAL_2022_TP_MONITORING_PAIRS_ACTUAL_EEA_WWTP_FWINTERNAL_REACH_SOURCEONLY_GATE.json"
SITES=DIR+"NATIONAL_DE_AT_2022_WISE6_RIVER_TP_PUBLIC_STATION_SOURCE_ONLY_QC.json"
OLD=DIR+"GERMANY_2022_WISE6_TP_8SOURCE_8ZERO_SAME_RIVER_MONITOR_QUARTETS_SOURCE_ONLY_LOCK.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hav(a,b):
    rad=math.pi/180
    dlat=(a["lat"]-b["lat"])*rad
    dlon=(a["lon"]-b["lon"])*rad
    z=math.sin(dlat/2)**2 + math.cos(a["lat"]*rad)*math.cos(b["lat"]*rad)*math.sin(dlon/2)**2
    return 12742*math.asin(min(1.,math.sqrt(z)))


def build(root):
    files={name:root/path for name,path in (("national_pairs",PAIRS),("public_sites",SITES),("prior_8quartets",OLD))}
    inputs={k:json.loads(v.read_text()) for k,v in files.items()}
    source=inputs["national_pairs"]
    site=inputs["public_sites"]
    old=inputs["prior_8quartets"]
    assert source["nationwide_source_only_pair_exposure_counts"]["original_all_candidates_n"]==115
    assert source["NOT_official_river_topology_independently_verified"]
    assert site["source_only_no_iNaturalist_butterfly_or_photo_values_read"]
    registry=site["nationwide_original_source_2022_station_registry"]
    stations={}
    for s in registry:
        key=s["original_country_scheme_site_SHA256"]
        assert key not in stations
        stations[key]={
            "original_site_id_SHA256":key,"countryCode":s["countryCode"],
            "lat":float(s["original_WISE6_WGS84_lat"]),
            "lon":float(s["original_WISE6_WGS84_lon"]),
            "original_2022_distinct_TP_sample_dates_metadata":int(s["original_2022_TP_independent_sample_dates"]),
            "year2022_WISE6_TP_NUMERIC_VALUES_NOT_ACCESSED":True,
        }
    oldids={
        entry[type][side+"_station_SHA256"]
        for entry in old["selected_original_8_source_vs_clean_zero_monitor_quartets"]
        for type in ("source","registered_source_zero")
        for side in ("upstream","downstream")
    }
    assert len(oldids)==32
    final=[]
    rejected=Counter()
    for p in source["all_national_2022_same_river_source_only_pairs_with_verified_outlet_between_reaches"]:
        if p["countryCode"]!="DE":
            rejected["out_of_DE"]+=1
            continue
        x=p["real_exact_source_between"]
        if int(x["strictly_between_total_original_Fw_outfalls_n"])>0:
            typ="registered_source_interposed"
        elif int(x["unresolved_same_reach_or_first_join_into_downstream_reach_source_n"])==0:
            typ="registered_qualifying_sourcezero"
        else:
            rejected["ambiguous_reach_outlet_order"]+=1
            continue
        u=p["upstream_station_SHA256"];d=p["downstream_station_SHA256"]
        assert u!=d and u in stations and d in stations
        assert stations[u]["countryCode"]==stations[d]["countryCode"]=="DE"
        assert stations[u]["original_2022_distinct_TP_sample_dates_metadata"]>=4
        assert stations[d]["original_2022_distinct_TP_sample_dates_metadata"]>=4
        if typ=="registered_source_interposed":
            assert int(x["strictly_between_total_original_Fw_outfalls_n"])>=1
            assert float(x["strictly_between_original_P_reported_t_year"])>0
        else:
            assert x["strictly_between_total_original_Fw_outfalls_n"]==0
            assert x["unresolved_same_reach_or_first_join_into_downstream_reach_source_n"]==0
        geo=hav(stations[u],stations[d])
        rec={
            "pair_id_SHA256":hashlib.sha256(("DE|"+u+"|"+d).encode()).hexdigest(),
            "countryCode":"DE",
            "exposure_class":typ,
            "upstream_original_site_id_SHA256":u,"downstream_original_site_id_SHA256":d,
            "network_full_reach_path_km":float(p["full_reach_path_to_downstream_reach_start_km"]),
            "original_geodesic_station_pair_km":geo,
            "same_2deg_upstream_region":p["upstream_country2degree_block"],
            "upstream_3deg_region":str(math.floor(stations[u]["lat"]/3))+"_"+str(math.floor(stations[u]["lon"]/3)),
            "2022_registered_interposed_outfall_count":int(x["strictly_between_total_original_Fw_outfalls_n"]),
            "2022_registered_interposed_reported_P_t_year":float(x["strictly_between_original_P_reported_t_year"]),
            "reuses_any_prior_32WISE6_quartet_station":bool(u in oldids or d in oldids),
            "all_registration_plus_directed_hydrology_labels_fixed_before_butterfly":True,
            "original_converted_HydroRIVERS_NOT_OFFICIAL_TOPOLOGY_VALIDATED":True
        }
        assert rec["network_full_reach_path_km"]>0
        final.append(rec)
    assert len(final)==73 and len({p["pair_id_SHA256"] for p in final})==73
    assert Counter(x["exposure_class"] for x in final)=={
        "registered_qualifying_sourcezero":40,"registered_source_interposed":33
    }
    subset={k:stations[k] for k in sorted(set(v for p in final for v in (p["upstream_original_site_id_SHA256"],p["downstream_original_site_id_SHA256"])))}
    assert len(subset)==100
    assert len([k for k in subset if k in oldids])==32
    assert all(v["original_2022_distinct_TP_sample_dates_metadata"]>=4 for v in subset.values())
    frame={}
    for deg in ("2deg","3deg"):
        groups=defaultdict(lambda:Counter())
        for p in final:
            region=p["same_2deg_upstream_region"] if deg=="2deg" else p["upstream_3deg_region"]
            groups[region][p["exposure_class"]]+=1
        typelist=sorted(groups)
        matched=[min(groups[g]["registered_source_interposed"],groups[g]["registered_qualifying_sourcezero"]) for g in groups]
        frame[deg]={
            "all_geographic_groups":len(groups),
            "groups_both_source_and_qualifying_registeredzero":sum(v>0 for v in matched),
            "mathematical_sourcezero_matched_quartet_upper_bound_before_site_disjointness_and_plant_effort":sum(matched),
            "counts_per_group":{g:dict(groups[g]) for g in typelist},
        }
    assert frame["2deg"]["mathematical_sourcezero_matched_quartet_upper_bound_before_site_disjointness_and_plant_effort"]==7
    return {
        "schema":"gesui_DE_national_2022_directedWISE6_73_sourcezero_pair_pre_plantae_pre_butterfly_frame_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_frozen_source_hashes_SHA256":{k:digest(v) for k,v in files.items()},
        "new_public_plant_photo_API_requests_performed":False,
        "new_butterfly_observation_API_requests_performed":False,
        "new_numerical_WISE6_TP_analyzed":False,
        "2022_EEA_year_is_not_known_butterfly_treatment_year":True,
        "zero_class_means_no_QUALIFYING_registered_interposed_WWTP_not_sewage_free":True,
        "original_geographic_candidates_total":115,
        "excluded_candidate_source_counts":dict(rejected),
        "Germany_eligible_2022_sourceonly_directed_candidate_pairs":len(final),
        "country_by_class_counts":dict(Counter(p["exposure_class"] for p in final)),
        "original_unique_monitoring_station_locations":len(subset),
        "original_prior_32_WISE6_sites_reused_in_geographic_frame":32,
        "original_additional_distinct_monitoring_sites_unobserved_by_prior_quartets":68,
        "spatial_positive_control_overlap_diagnostic":frame,
        "original_unique_public_WISE6_monitor_site_metadata":list(subset.values()),
        "complete_original_fixed_candidate_pairs":sorted(final,key=lambda r:(r["exposure_class"],r["same_2deg_upstream_region"],r["network_full_reach_path_km"],r["pair_id_SHA256"]))
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=build(root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in out.items() if k in (
        "Germany_eligible_2022_sourceonly_directed_candidate_pairs",
        "country_by_class_counts","original_unique_monitoring_station_locations",
        "original_additional_distinct_monitoring_sites_unobserved_by_prior_quartets",
        "spatial_positive_control_overlap_diagnostic")},indent=2))
