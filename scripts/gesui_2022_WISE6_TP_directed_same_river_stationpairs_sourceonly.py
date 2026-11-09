#!/usr/bin/env python3
"""Source-only 2022 WISE6 TP monitoring-site PAIR feasibility on SAME directed river.

Uses original 187-source-station graph assignments, pinned third-party EU
HydroRIVERS reaches. Does NOT load WISE6 measured TP concentrations, site
censor/median values, iNaturalist butterfly photo outcomes, or decide effects.
Original ≤50km graph path is reach-length, NOT exact between-station distance.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import geopandas as gpd
import pandas as pd

RIVER_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
EXPECTED_STATIONS=187
EXPECTED_SNAPPED=174
EXPECTED_RIVER_REACHES=199545
MAX_PATH_KM=50.0
MIN_N_DISJOINT=10
MIN_EACH_COUNTRY=3
MIN_TWO_DEGREE_BLOCKS=4
BBOX=(-12,44,26,64)
COUNTRIES=("DE","AT")

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()

def main():
    q=argparse.ArgumentParser()
    q.add_argument("--original-directed-sourceonly",required=True,type=Path)
    q.add_argument("--original-hydrorivers-parquet",required=True,type=Path)
    q.add_argument("--out",required=True,type=Path)
    a=q.parse_args()
    prior=json.loads(a.original_directed_sourceonly.read_text())
    if prior.get("schema")!="gesui_2022_EEA_WISE6_original_TP_187_stations_HydroRIVERS_directed_WWTP2022_SOURCE_ONLY_v01":
        raise ValueError("Original 187 WISE stations source-only contract mismatch")
    if prior.get("original_2022_WISE6_numeric_TP_concentration_not_loaded_or_queried") is not True:
        raise ValueError("Original site cohort was not source-only")
    if prior.get("hydrorivers_thirdparty_EU_v1_GeoParquet_SHA256")!=RIVER_SHA:
        raise ValueError("Original 2022 directed source network does not match pinned graph")
    if sha(a.original_hydrorivers_parquet)!=RIVER_SHA:
        raise ValueError("HydroRIVERS original third-party mirror SHA drift")
    original=prior["original_187_WISE6_TP_stations_with_public_source_reach_exposures"]
    if len(original)!=EXPECTED_STATIONS:
        raise ValueError("Unexpected original WISE6 stations n")
    snapped=[x for x in original if isinstance(x.get("HYRIV_ID"),int) and
             x.get("original_receptor_snap_m") is not None and
             x["original_receptor_snap_m"]<=500]
    if len(snapped)!=EXPECTED_SNAPPED:
        raise ValueError("Original 2022 source-only reach-snap 174 station count drift")
    # Mirror limitations and same reach topology as original actual WWTP source.
    g=gpd.read_parquet(a.original_hydrorivers_parquet).to_crs("EPSG:4326")
    g=g.rename(columns={c:c.upper() for c in g.columns if c!=g.geometry.name})
    g=g.cx[BBOX[0]:BBOX[2],BBOX[1]:BBOX[3]].copy()
    needed=("HYRIV_ID","NEXT_DOWN","LENGTH_KM")
    for c in needed:g[c]=pd.to_numeric(g[c],errors="coerce")
    g=g[g.HYRIV_ID.notna()].copy()
    if g.HYRIV_ID.duplicated().any() or len(g)!=EXPECTED_RIVER_REACHES:
        raise ValueError("Original HydroRIVERS bounded study reaches changed")
    river={}
    for rec in g[["HYRIV_ID","NEXT_DOWN","LENGTH_KM"]].itertuples(index=False):
        rid=int(rec.HYRIV_ID)
        down=int(rec.NEXT_DOWN) if pd.notna(rec.NEXT_DOWN) else 0
        length=float(rec.LENGTH_KM)
        if not math.isfinite(length) or length<0:continue
        river[rid]=(down,length)
    byreach=defaultdict(list)
    for x in snapped:
        if x["HYRIV_ID"] not in river:
            raise ValueError("Pre-frozen station reach not present in original pinned graph")
        byreach[int(x["HYRIV_ID"])].append(x)
    all_pairs=[];diagnostics=Counter()
    # Each upstream source station explores only one original flow directed
    # graph path. Source WWTP original station exposure remains DIAGNOSTIC,
    # no TP magnitude or outcomes are accessed for station pair selection.
    for upstream in snapped:
        origin=int(upstream["HYRIV_ID"])
        cursor=river[origin][0]
        d_km=float(river[origin][1])
        traversed=1
        seen={origin}
        while cursor in river and cursor>0 and cursor not in seen and d_km<=MAX_PATH_KM:
            seen.add(cursor)
            for downstream in byreach.get(cursor,[]):
                diagnostics["raw_candidate_same_directed_river_pairs_total_including_other_countries"]+=1
                if downstream["countryCode"]!=upstream["countryCode"]:
                    diagnostics["excluded_country_crossing_receptor_pairs"]+=1
                    continue
                pair={
                  "countryCode":upstream["countryCode"],
                  "upstream_station_SHA256":upstream["original_WISE6_country_scheme_site_SHA256"],
                  "downstream_station_SHA256":downstream["original_WISE6_country_scheme_site_SHA256"],
                  "upstream_original_river_reach_id":int(upstream["HYRIV_ID"]),
                  "downstream_original_river_reach_id":int(downstream["HYRIV_ID"]),
                  "upstream_country2degree_block":upstream["original_country2degree_block"],
                  "downstream_country2degree_block":downstream["original_country2degree_block"],
                  "full_reach_path_to_downstream_reach_start_km":round(d_km,5),
                  "original_directed_reach_steps":traversed,
                  "original_upstream_TP_2022_sample_date_n":upstream["source_WISE6_TP_2022_distinct_sample_days"],
                  "original_downstream_TP_2022_sample_date_n":downstream["source_WISE6_TP_2022_distinct_sample_days"],
                  "registered_upstream_WWTP_count_sourcehint_upstation":
                       upstream["strict_upstream_qualified_WWTP_count_50km"],
                  "registered_upstream_WWTP_count_sourcehint_downstation":
                       downstream["strict_upstream_qualified_WWTP_count_50km"],
                  "SOURCE_HINT_IS_NOT_VERIFIED_INTERVENING_WWTP":True,
                  "no_actual_TP_numeric_concentration_opened":True,
                }
                all_pairs.append(pair)
                diagnostics["raw_candidate_same_country_pair"]+=1
            down,length=river[cursor]
            d_km+=length
            cursor=down
            traversed+=1
        if cursor in seen:diagnostics["river_directed_graph_cycle_encountered"]+=1
    # Deterministic distance-first, hashed site-ID tie-breaker, independent
    # of TP concentration AND independent of 2022 original EEA source P mass.
    all_pairs.sort(key=lambda x:(x["full_reach_path_to_downstream_reach_start_km"],
           hashlib.sha256((x["upstream_station_SHA256"]+"|"+
                 x["downstream_station_SHA256"]).encode()).hexdigest()))
    used=set();disjoint=[]
    for p in all_pairs:
        if p["upstream_station_SHA256"] in used or p["downstream_station_SHA256"] in used:continue
        used.add(p["upstream_station_SHA256"]);used.add(p["downstream_station_SHA256"])
        disjoint.append(p)
    by=Counter(p["countryCode"] for p in disjoint)
    blocs=set(p["countryCode"]+"|"+p["upstream_country2degree_block"] for p in disjoint)
    # Stage B can discover WWTPs on exact interior reaches ONLY if raw panel
    # supports a new physical same-path study, not the previous 171-site CV.
    eligible=(len(disjoint)>=MIN_N_DISJOINT and
              by.get("DE",0)>=MIN_EACH_COUNTRY and
              by.get("AT",0)>=MIN_EACH_COUNTRY and
              len(blocs)>=MIN_TWO_DEGREE_BLOCKS)
    hint=Counter()
    for p in disjoint:
        delta=p["registered_upstream_WWTP_count_sourcehint_downstation"]-p["registered_upstream_WWTP_count_sourcehint_upstation"]
        hint["upstream_source_count_increases_at_downstation"]+=int(delta>0)
        hint["upstream_source_count_equal"]+=int(delta==0)
        hint["upstream_source_count_decreases_due_to_50km_window"]+=int(delta<0)
    result={
      "schema":"gesui_2022_WISE6_original_TP_same_river_strict_directed_stationpair_SOURCE_ONLY_STAGEA_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "original_2022_187_TP_WWTP_directed_sourceonly_SHA256":sha(a.original_directed_sourceonly),
      "hydrorivers_Europe_v1_THIRDPARTY_parquet_SHA256":RIVER_SHA,
      "hydrorivers_thirdparty_not_verified_against_official":True,
      "original_stations_2022_with_sourceonly_snap":len(snapped),
      "original_directed_study_river_reaches":len(river),
      "max_full_reach_path_km":MAX_PATH_KM,
      "no_same_reach_station_pair":True,
      "both_monitor_sites_same_country_required":True,
      "pair_independence_no_repeat_monitor_station":True,
      "pair_sorting_excludes_all_TP_measured_values_and_original_source_P_mass":True,
      "path_distance_is_sum_of_whole_reach_lengths_not_exact_station_to_station_distance":True,
      "original_all_same_river_source_only_pairs":all_pairs,
      "original_source_only_disjoint_shortest_first_pairs":disjoint,
      "original_same_river_source_only_pair_QC":dict(diagnostics),
      "disjoint_same_river_pairs_n":len(disjoint),
      "disjoint_by_country":dict(by),
      "disjoint_original_2degree_upstream_blocks":len(blocs),
      "predeclared_stageA_same_river_disjoint_pair_gate":eligible,
      "source_only_disjoint_cumulative_upstream_count_DIAGNOSTIC_NOT_between_reach":dict(hint),
      "NEVER_use_50km_upstream_count_difference_as_proof_WWTP_intervenes":True,
      "no_TP_concentrations_or_butterfly_photo_data_read":True,
      "next_if_stageA_passes":"Retrieve original WWTP origin reach IDs and classify exact strict BETWEEN-source reach with original NEXT_DOWN topology and freeze pair-exposure sample before any TP numerical difference",
      "next_if_stageA_fails":"STOP, report that original 2022 WISE6 spatial monitoring design does not support 50km within-river before/after WWTP contrasts",
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_WISE6_2022_SAME_DIRECTED_RIVER_50KM_TP_MONITORING_SOURCEPAIR_GATE",
          json.dumps({"source_n":len(snapped),"raw_same_country_pairs":len(all_pairs),
            "disjoint_pairs":len(disjoint),"by_country":dict(by),
            "2degree_blocks":len(blocs),"source_count_hints_NOT_exposure":dict(hint),
            "next_real_intervening_WWTP_source_stage_supported":eligible},sort_keys=True),flush=True)

if __name__=="__main__":main()
