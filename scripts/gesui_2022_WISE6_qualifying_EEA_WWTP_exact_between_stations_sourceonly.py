#!/usr/bin/env python3
"""Original 2022 EEA freshwater WWTP outlet SOURCE REACH strictly between two
W I S E 6 river TP monitoring source sites, fixed 50km HydroRIVERS route.

Only source metadata; NO measured TP or butterfly photo outcomes opened.
Stage-A same-river candidate pairs were independently frozen and SHA-locked.
"""
from __future__ import annotations
import argparse,hashlib,json,math,sys
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import geopandas as gpd,pandas as pd,numpy as np

RIVER_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
OLD_BLOB="9c07ce9de7f5104b47052c7598a943df405a4b83"
BBOX=(-12,44,26,64)
MAX_KM=50.
ORIG_WWTP=5492;ORIG_DCP=24416;ORIG_FW=4057
MIN_POS=8;MIN_ZERO=8;MIN_BLOCKS=4

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def gitblob(p):
    raw=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def topo(g):
    t={}
    for row in g[["HYRIV_ID","NEXT_DOWN","LENGTH_KM"]].itertuples(index=False):
        a=int(row.HYRIV_ID)
        b=int(row.NEXT_DOWN) if pd.notna(row.NEXT_DOWN) else 0
        length=float(row.LENGTH_KM)
        if length<0 or not math.isfinite(length):raise ValueError("Original network invalid segment length")
        if a in t:raise ValueError("Network duplicate HYRIV_ID")
        t[a]=(b,length)
    return t

def path_from(origin,graph):
    visited={origin};routes=[origin];length=0.
    cursor=origin
    while cursor in graph:
        nxt,segmentlen=graph[cursor]
        length+=segmentlen
        if length>MAX_KM or nxt not in graph or nxt<=0 or nxt in visited:break
        visited.add(nxt);routes.append(nxt);cursor=nxt
    return routes

def station_path(origin,dest,graph):
    if origin==dest:return None
    route=[origin];seen={origin};total=0.
    curr=origin
    while curr in graph:
        nxt,length=graph[curr]
        total+=length
        if nxt<=0 or nxt in seen or nxt not in graph or total>MAX_KM:return None
        route.append(nxt);seen.add(nxt)
        if nxt==dest:return route
        curr=nxt
    return None

def source_appears_between(record,src_routes,source_index,graph):
    origin=int(record["upstream_original_river_reach_id"])
    dest=int(record["downstream_original_river_reach_id"])
    path=station_path(origin,dest,graph)
    if path is None or len(path)<2:
        raise ValueError("Stage A candidate is not a genuine original directed path")
    internal=set(path[1:-1])
    strictly_intervening_trunk=[]
    strictly_intervening_joining_tributary=[]
    potential_uncertain=[]
    already_upstream_before_first_station=0
    # Each outlet path contains original outlet reach plus directed reach
    # successors up to 50km. Qualifying source must meet INTERNAL reach
    # and be able to reach downstream station reach within the 50km horizon,
    # and must NEVER reach the upstream station reach.
    for idx,(origin_reach,massP,massN,codehash) in enumerate(source_index):
        route=src_routes[idx]
        if dest not in route:continue
        if origin in route:
            already_upstream_before_first_station+=1
            # Original outlet reach could equal upstream monitoring reach,
            # where within-reach source-vs-station order is UNKNOWN.
            if origin_reach==origin:
                potential_uncertain.append(codehash)
            continue
        if origin_reach==dest:
            potential_uncertain.append(codehash)
            continue
        if not internal.intersection(route):
            # Source could join the downstream monitoring reach itself
            # from a tributary. Within this reach, ordering is unknown.
            potential_uncertain.append(codehash)
            continue
        if origin_reach in internal:
            strictly_intervening_trunk.append((codehash,massP,massN))
        else:
            strictly_intervening_joining_tributary.append((codehash,massP,massN))
    included=strictly_intervening_trunk+strictly_intervening_joining_tributary
    counts={
      "strictly_between_trunk_original_Fw_outfall_n":len(strictly_intervening_trunk),
      "strictly_between_intersecting_tributary_original_Fw_outfall_n":len(strictly_intervening_joining_tributary),
      "strictly_between_total_original_Fw_outfalls_n":len(included),
      "strictly_between_original_P_reported_t_year":sum(x[1] for x in included),
      "strictly_between_original_N_reported_t_year":sum(x[2] for x in included),
      "unresolved_same_reach_or_first_join_into_downstream_reach_source_n":len(potential_uncertain),
      "source_already_strictly_upstream_of_first_monitor_n":already_upstream_before_first_station,
      "source_route_match_uses_original_50km_horizon_not_true_WISE2022_flow":True,
      "source_interposition_confirmed_only_in_strict_internal_reaches":True,
      "NO_TP_WATER_CONCENTRATION_OPENED":True,
    }
    return counts

def pick_strata(records):
    # Select source-free controls first (scarcer), then strictly verified
    # source-between pairs. Fixed shortest-distance then station-ID hash,
    # all selection strictly exposure-only, no TP outcome magnitudes used.
    sorted_records=sorted(records,key=lambda p:(
         p["full_reach_path_to_downstream_reach_start_km"],
         hashlib.sha256((p["upstream_station_SHA256"]+"|"+
                         p["downstream_station_SHA256"]).encode()).hexdigest()))
    controls=[x for x in sorted_records if x["real_exact_source_between"]["strictly_between_total_original_Fw_outfalls_n"]==0 and
              x["real_exact_source_between"]["unresolved_same_reach_or_first_join_into_downstream_reach_source_n"]==0]
    exposed=[x for x in sorted_records if x["real_exact_source_between"]["strictly_between_total_original_Fw_outfalls_n"]>0]
    selected={"control":[], "intervening_source": []}
    used=set()
    for key,options in (("control",controls),("intervening_source",exposed)):
        for r in options:
            a=r["upstream_station_SHA256"];b=r["downstream_station_SHA256"]
            if a in used or b in used:continue
            if len(selected[key])>=MIN_ZERO if key=="control" else len(selected[key])>=MIN_POS:continue
            used.update([a,b]);selected[key].append(r)
    return selected,{"original_all_candidates_n":len(records),
                     "source_zero_without_reach_order_ambiguity_candidate_pairs":len(controls),
                     "strict_source_between_candidate_pairs":len(exposed),
                     "unresolved_source_between_or_downstream_reach_n":
                        sum(x["real_exact_source_between"]["strictly_between_total_original_Fw_outfalls_n"]==0 and
                            x["real_exact_source_between"]["unresolved_same_reach_or_first_join_into_downstream_reach_source_n"]>0 for x in records)}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--stageA-source-only-pairs",required=True,type=Path)
    p.add_argument("--original-187-directed-source-only",required=True,type=Path)
    p.add_argument("--pinned-HydroRIVERS-eu-parquet",required=True,type=Path)
    p.add_argument("--pinned-original-UWWTD-source-script",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    first=json.loads(a.stageA_source_only_pairs.read_text())
    if not first.get("predeclared_stageA_same_river_disjoint_pair_gate"):
        raise ValueError("Stage A before/after site-pair source-only feasibility failed")
    if first["hydrorivers_Europe_v1_THIRDPARTY_parquet_SHA256"]!=RIVER_SHA or sha(a.pinned_HydroRIVERS_eu_parquet)!=RIVER_SHA:
        raise ValueError("Not same original HydroRIVERS EU v1 converted mirror")
    old=json.loads(a.original_187_directed_source_only.read_text())
    if first["original_2022_187_TP_WWTP_directed_sourceonly_SHA256"]!=sha(a.original_187_directed_source_only):
        raise ValueError("Original existing 187 TP WISE source-only paper source changed")
    if gitblob(a.pinned_original_UWWTD_source_script)!=OLD_BLOB:
        raise ValueError("Original 2022 EEA REST source script Git blob drift")
    sys.path.insert(0,str(a.pinned_original_UWWTD_source_script.parent))
    import gesui_eea_nutrient_load_river_crosswalk as prior
    import gesui_2022_WISE6_TP_original_station_vs_directed_WWTP_sourceonly as linked
    WW,prec=prior.download_layer(prior.PLANT,prior.BOTH_WHERE,prior.PLANT_FIELDS,ORIG_WWTP)
    DP,drec=prior.download_layer(prior.OUTFALL,"1=1",prior.OUT_FIELDS,ORIG_DCP)
    oldrec=old["original_source_API_page_response_SHA_receipts"]
    if prec+drec!=oldrec:
        raise ValueError("Original EEA 2022 measured N/P source REST page digests differ from already frozen source-only original")
    WW=prior.normalize_key(WW);DP=prior.normalize_key(DP)
    WW=WW[WW.compound_join_key.notna()].copy()
    DP=DP[DP.compound_join_key.notna()].copy()
    n=WW.compound_join_key.value_counts();d=DP.compound_join_key.value_counts()
    WW=WW[WW.compound_join_key.map(n).eq(1)].copy()
    WW["original_one_dcp"]=WW.compound_join_key.map(d).fillna(0).astype(int)
    unique=WW[WW.original_one_dcp.eq(1)].merge(
      DP.drop_duplicates("compound_join_key"),on="compound_join_key",
      how="inner",suffixes=("_WWTP","_DCP"),validate="one_to_one")
    for col in ("dcpLatitude","dcpLongitude","uwwPDischargeMeasured","uwwNDischargeMeasured"):
        unique[col]=pd.to_numeric(unique[col],errors="coerce")
    unique=unique[unique.dcpLatitude.between(44,64)&unique.dcpLongitude.between(-12,26)]
    unique=unique[unique.dcpWaterBodyType.astype("string").str.strip().str.upper().eq("FW")].copy()
    if len(unique)!=ORIG_FW:raise ValueError("Original reported measured 2022 one-FW source original 4057 outfalls mismatch")
    river=gpd.read_parquet(a.pinned_HydroRIVERS_eu_parquet).to_crs("EPSG:4326")
    river=river.rename(columns={x:x.upper() for x in river.columns if x!=river.geometry.name})
    river=river.cx[-12:26,44:64].copy()
    for field in ("HYRIV_ID","NEXT_DOWN","LENGTH_KM"):
        river[field]=pd.to_numeric(river[field],errors="coerce")
    river=river[river.HYRIV_ID.notna()].copy()
    sourceframe=unique.copy()
    sourceframe["source_row_id"]=np.arange(len(unique))
    sourceframe["lat"]=sourceframe.dcpLatitude
    sourceframe["lon"]=sourceframe.dcpLongitude
    snapped,sourceq=linked.snap_unambiguous(sourceframe,river,"original_public_WWTP")
    if len(snapped):
        sourceq["original_source_eligible_FW_outfalls"] = ORIG_FW
    graph=topo(river)
    if len(graph)!=199545:raise ValueError("Original HydroRIVERS 199545 study river reach count changed")
    source_index=[];route_cache=[]
    for r in snapped.itertuples(index=False):
        orig=int(r.HYRIV_ID)
        if orig not in graph:continue
        P=float(r.uwwPDischargeMeasured);N=float(r.uwwNDischargeMeasured)
        if not (math.isfinite(P) and math.isfinite(N) and P>0 and N>0):continue
        srcid=hashlib.sha256(str(r.compound_join_key).encode()).hexdigest()
        source_index.append((orig,P,N,srcid))
        route_cache.append(path_from(orig,graph))
    entries=[]
    for row in first["original_all_same_river_source_only_pairs"]:
        metrics=source_appears_between(row,route_cache,source_index,graph)
        entries.append({**row,"real_exact_source_between":metrics})
    selected,stats=pick_strata(entries)
    chosen=selected["control"]+selected["intervening_source"]
    country=Counter(x["countryCode"] for x in chosen)
    blocks=set(x["countryCode"]+"|"+x["upstream_country2degree_block"] for x in chosen)
    gate=(len(selected["control"])>=MIN_ZERO and len(selected["intervening_source"])>=MIN_POS and
          country.get("DE",0)>0 and country.get("AT",0)>0 and len(blocks)>=MIN_BLOCKS)
    result={
      "schema":"gesui_2022_WISE6_same_river_real_EEA2022_FW_intervening_outlet_SOURCE_ONLY_STAGEB_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "original_stageA_original_50km_riverpair_source_SHA256":sha(a.stageA_source_only_pairs),
      "original_frozen_187_source_SHA256":sha(a.original_187_directed_source_only),
      "original_EEA_2022_REST_31pages_MATCH_EARLIER_SHA_receipts":True,
      "original_EEA_2022_source_total_bothpositive_facilities":ORIG_WWTP,
      "original_EEA_2022_unique_public_FW_outfalls":ORIG_FW,
      "original_EEA2022_FW_outfalls_unambiguous_river_reach_snaps":len(source_index),
      "original_EEA2022_FW_outfalls_graph_snap_quality":sourceq,
      "original_EEA2022_HydroRIVERS_mirror_SHA256":RIVER_SHA,
      "source_topology_NOT_verified_against_official_HydroRIVERS_original":True,
      "strict_trunk_or_tributary_WWTP_origin_joins_internal_monitor_reaches":True,
      "same_upstream_station_or_downstream_station_reach_outfalls_ambiguous_and_excluded":True,
      "source_annual_P_tons_not_water_concentration":True,
      "original_all_sourceonly_stageA_candidate_pairs_with_exact_outlet_interposition":entries,
      "exact_original_WWTP_intervening_source_candidate_counts":stats,
      "sourceonly_independent_selected_pairs":selected,
      "sourceonly_disjoint_by_country":dict(country),
      "sourceonly_disjoint_2deg_blocks":len(blocks),
      "predeclared_stageB_real_WWTP_intervening_vs_sourcefree_pair_gate_pass":gate,
      "NO_WISE6_RIVER_TP_NUMERICAL_VALUES_ACCESSED":True,
      "NO_LYCAENIDAE_BUTTERFLY_PHOTO_OR_ABUNDANCE_QUERIED":True,
      "no_phosphorus_biological_or_pollution_causation_claim":True
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_WISE6_2022_STRICT_INTERVENING_EEA_WWTP_MONITORING_PAIR_SOURCE_GATE",
       json.dumps({"candidate":len(entries),"source_strata":stats,
          "selected_approved_control":len(selected["control"]),
          "selected_actual_source_intervening":len(selected["intervening_source"]),
          "countries":dict(country),"blocks":len(blocks),"gate":gate,
          "no_TP_concentration_values":True},sort_keys=True),flush=True)

if __name__=="__main__":main()
