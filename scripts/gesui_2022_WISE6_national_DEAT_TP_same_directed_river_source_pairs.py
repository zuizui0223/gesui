#!/usr/bin/env python3
"""National source-only 2022 DE/AT EEA WISE6 TP monitoring pairs on SAME HydroRIVERS
directed original river paths ≤50km, independently of butterfly photo regions.

NO river TP concentration, NO water quality outcome, NO butterfly/photo.
Original national 342 TP site sampling dates/source geometry were fixed
from full 102M original source rows first. Exact old 187 TP site graph
assignment must reproduce prior source-only river reach snaps.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import geopandas as gpd,pandas as pd,numpy as np
import gesui_2022_WISE6_TP_original_station_vs_directed_WWTP_sourceonly as routing

GRAPH_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
NATIONAL_ORIGINAL_N=342
MAX_KM=50.0
EXPECTED_OLD_187=187
MIN_DISJOINT=10
MIN_DE=3
MIN_AT=3
MIN_BLOCKS=4

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--national-site-source",required=True,type=Path)
    p.add_argument("--original-187-source-snap",required=True,type=Path)
    p.add_argument("--pinned-HydroRIVERS-eu",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    national=json.loads(a.national_site_source.read_text())
    prior=json.loads(a.original_187_source_snap.read_text())
    if not national["source_only_no_river_TP_numeric_values_parsed_or_selected"]:
        raise ValueError("National EEA source was not original TP-response-blind")
    if national["source_original_national_date_eligible_2022_TP_site_total"]!=NATIONAL_ORIGINAL_N:
        raise ValueError("Original national 2022 TP source 342 counts shifted")
    records=national["nationwide_original_source_2022_station_registry"]
    # Original 342 sites pass the date/month gate, but only 301 have
    # unambiguous published WISE coordinates. The original 41 without usable
    # coordinates MUST NOT be forced onto a river or called source-zero.
    if len(records)!=national["source_original_country_geocoded_valid_stations"] or len(records)<300:
        raise ValueError("Original WISE6 geocoded national site population drift")
    if sha(a.pinned_HydroRIVERS_eu)!=GRAPH_SHA:
        raise ValueError("Pinned HydroRIVERS mirror SHA drift")
    river=gpd.read_parquet(a.pinned_HydroRIVERS_eu).to_crs("EPSG:4326")
    river=river.rename(columns={x:x.upper() for x in river.columns if x!=river.geometry.name})
    river=river.cx[-12:26,44:64].copy()
    for field in ("HYRIV_ID","NEXT_DOWN","LENGTH_KM","DIS_AV_CMS","UPLAND_SKM"):
        river[field]=pd.to_numeric(river[field],errors="coerce")
    river=river[river.HYRIV_ID.notna()].copy()
    if len(river)!=199545:raise ValueError("Original pinned directed river reach count changed")
    topology=routing.river_topology(river)
    original_frame=pd.DataFrame(records)
    original_frame["source_row_id"]=np.arange(len(original_frame))
    original_frame["lat"]=original_frame.original_WISE6_WGS84_lat
    original_frame["lon"]=original_frame.original_WISE6_WGS84_lon
    snapped,QC=routing.snap_unambiguous(original_frame,river,"original_WISE6_national_TP")
    lookup={x.original_country_scheme_site_SHA256:x for x in snapped.itertuples(index=False)}
    original_old=prior["original_187_WISE6_TP_stations_with_public_source_reach_exposures"]
    if len(original_old)!=EXPECTED_OLD_187:raise ValueError("Original photo-geometry 187 site count shifted")
    old_same=0;old_unknown=0
    for old in original_old:
        r=lookup.get(old["original_WISE6_country_scheme_site_SHA256"])
        if old["HYRIV_ID"] is None:
            if r is not None:raise ValueError("Previously not snapped WISE6 site unexpectedly snapped in national scope")
            old_unknown+=1
        else:
            if r is None or int(r.HYRIV_ID)!=int(old["HYRIV_ID"]):
                raise ValueError("Original independent 187 source monitoring reach assignment failed direct replay")
            old_same+=1
    if old_same!=174 or old_unknown!=13:raise ValueError("Original independent reach mapping old 174/13 changed")
    byreach=defaultdict(list)
    for r in snapped.itertuples(index=False):
        if int(r.HYRIV_ID) not in topology:raise ValueError("Snapped WISE6 station reach missing")
        byreach[int(r.HYRIV_ID)].append(r)
    options=[];issues=Counter()
    for up in snapped.itertuples(index=False):
        origin=int(up.HYRIV_ID)
        seen={origin};nxt=topology[origin][0]
        dist=float(topology[origin][1]);steps=1
        while nxt in topology and nxt>0 and nxt not in seen and dist<=MAX_KM:
            seen.add(nxt)
            for down in byreach.get(nxt,[]):
                issues["original_same_river_possible_pairs_total_with_country_crossing"]+=1
                if down.countryCode!=up.countryCode:
                    issues["excluded_monitor_pair_crosses_national_country"]+=1
                    continue
                options.append({
                    "countryCode":up.countryCode,
                    "upstream_station_SHA256":up.original_country_scheme_site_SHA256,
                    "downstream_station_SHA256":down.original_country_scheme_site_SHA256,
                    "upstream_original_reach_HYRIV_ID":origin,
                    "downstream_original_reach_HYRIV_ID":int(down.HYRIV_ID),
                    "upstream_monitor_2022_TP_source_sample_dates":int(up.original_2022_TP_independent_sample_dates),
                    "downstream_monitor_2022_TP_source_sample_dates":int(down.original_2022_TP_independent_sample_dates),
                    "upstream_2degree_region":up.original_2degree_block,
                    "downstream_2degree_region":down.original_2degree_block,
                    "full_reach_path_approx_km":round(dist,5),
                    "directed_hyriv_segments":steps,
                    "NO_numeric_river_TP_or_original_WWTP_source_outcome_used":True,
                })
                issues["source_only_same_country_directed_pair_options"]+=1
            downid,length,_,_=topology[nxt]
            dist+=length
            nxt=downid;steps+=1
        if nxt in seen:issues["original_river_topology_cycle"]+=1
    options.sort(key=lambda x:(x["full_reach_path_approx_km"],
          hashlib.sha256((x["upstream_station_SHA256"]+"|"+
                           x["downstream_station_SHA256"]).encode()).hexdigest()))
    used=set();chosen=[]
    for x in options:
        if x["upstream_station_SHA256"] in used or x["downstream_station_SHA256"] in used:continue
        used.update((x["upstream_station_SHA256"],x["downstream_station_SHA256"]))
        chosen.append(x)
    by=Counter(x["countryCode"] for x in chosen)
    blocks=set(x["countryCode"]+"|"+x["upstream_2degree_region"] for x in chosen)
    gate=(len(chosen)>=MIN_DISJOINT and by["DE"]>=MIN_DE and by["AT"]>=MIN_AT and len(blocks)>=MIN_BLOCKS)
    out={
      "schema":"gesui_2022_original_WISE6_NATIONAL_DE_AT_RIVER_TP_SAME_REACHPATH_50km_stationpair_SOURCE_ONLY_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "independent_new_all_national_source_metadata_sha256":sha(a.national_site_source),
      "original_prior_187_source_directed_WWTP_noTP_QC_SHA256":sha(a.original_187_source_snap),
      "original_thirdparty_HydroRIVERS_EU_SHA256":GRAPH_SHA,
      "original_thirdparty_HydroRIVERS_NOT_official_topology_verified":True,
      "original_national_date_eligible_sites":len(records),
      "national_stations_uniquely_snapped_500m":len(snapped),
      "original_187_prior_station_map_exact_reproduced":{"old_snapped":old_same,"old_unsnapped":old_unknown},
      "original_2022_national_WISE6_river_snap_source_QC":QC,
      "original_same_country_same_river_50km_pair_QC":dict(issues),
      "all_original_national_source_only_directed_same_river_station_pair_candidates":options,
      "source_only_disjoint_shortest_first_independent_monitor_pairs":chosen,
      "source_only_disjoint_pairs_n":len(chosen),
      "source_only_disjoint_pair_country_counts":dict(by),
      "source_only_disjoint_country2degree_blocks":len(blocks),
      "source_only_predeclared_10_pairs_DE3_AT3_fourblocks_gate":gate,
      "source_only_NO_riverTP_concentration_magnitudes_or_photographs_READ":True,
      "next_if_supported":"Locate exact original direct FW WWTP reported source outlet reach origins BETWEEN pairs, minimum 8 disjoint source-positive and 8 clean source-zero controls; no TP numerical outcome until source-only gate",
      "do_not_automatically_rebrand_national_TP_cohort_as_Lycaenidae_photo_effect":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_2022_NATIONAL_DE_AT_RIVER_TP_SAME_DIRECTED_RIVER_SOURCEONLY_PAIRS",
      json.dumps({"original_national":len(records),"snapped":len(snapped),
        "raw_directed_pair_source_count":len(options),"disjoint":len(chosen),
        "country":dict(by),"2degree_blocks":len(blocks),"stageA_gate":gate,
        "old_187_reachIDs_reproduced":old_same},sort_keys=True),flush=True)

if __name__=="__main__":main()
