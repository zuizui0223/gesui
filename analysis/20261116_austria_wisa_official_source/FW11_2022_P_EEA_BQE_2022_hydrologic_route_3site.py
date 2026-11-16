#!/usr/bin/env python3
"""Posthoc *source-only* river-flow topology check for the actual BMLUK
2022 annual-TP × GZÜV original-225 target PKINT×EEA2022 BQE overlap.

Three genuine BMLUK FW1.1 2022 TP site IDs joined to original 225 GZÜV
station IDs and 117 EEA annual measured-NP 2022 outlets (3 outfalls).
Exactly two of three have original 2022 EEA benthos site IDs, one with a
955m government-vs-EEA coordinate conflict. Determine whether the nominal
benthic reach is truly DOWNSTREAM of the original outfall in authentic
HydroRIVERS NEXT_DOWN, not just within 2km in straight-line distance.

No actual numeric 2022 TP / benthic EQR / NDVI / butterflies used in the
route test, and BMLUK 2022 FW1 means are NOT evidence of wastewater-derived
nutrient delivery to any particular riparian plant.
"""
from __future__ import annotations

import argparse,collections,csv,hashlib,json,math
from pathlib import Path
import geopandas as gpd
import numpy as np

SCHEMA="GESUI_20261116_BMLUK_FW11_2022_TP_225_GZUEV_TARGET_2022_BENTHIC_TRUE_RIVER_FLOW_V1"
RIVER_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
OUTFALL_GEO_SHA="225_FIXED_117_AT_TRUE_NP_YEAR_SOURCE_USING_SIGNED_PR116"
BQE_BLOB="325f199e985d9862ffee5b0edb205c61a7be3833"
AT_SITE_SOURCE={
    "ATFW54110117":{"eea2022_benthic":True, "gov_EEA_benthic_station_coordinate_discrepancy_m":954.614},
    "ATFW73190907":{"eea2022_benthic":False,"gov_EEA_benthic_station_coordinate_discrepancy_m":None},
    "ATFW73290907":{"eea2022_benthic":True, "gov_EEA_benthic_station_coordinate_discrepancy_m":0.5}
}
ORIG_STATION_TO_OUTFALL={
 "ATFW54110117":"AT_UWWTD_2022|ATTP_5-A1699622R154",
 "ATFW73190907":"AT_UWWTD_2022|ATTP_7-7218603",
 "ATFW73290907":"AT_UWWTD_2022|ATTP_7-7930111"
}
MAX_HOPS=500

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def gitblob(path):
    b=path.read_bytes()
    return hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()

def source_outfalls(folder):
    gate=json.loads((folder/"EEA2022_NP_SOURCEFIRST_LOCAL_DIRECTED_RIVER_GEOMETRY_GATE.json").read_text())
    file=folder/"EEA2022_ALL_SOURCE_QUALIFIED_UWWTD_NP_LOCAL_RIVER_SITES.csv"
    if (gate.get("original_all_true_receiving_outfall_geo_rows")!=282 or
        gate.get("original_all_source_valid_geography_CSV_SHA256")!=sha(file)):
        raise RuntimeError("Original EEA UWWTD 282 genuine directed river outfall source SHA/cardinality drift")
    with file.open(newline="",encoding="utf-8-sig") as f:
        o={q["original_EEA2022_plant_id"]:q for q in csv.DictReader(f)}
    if len(o)!=282:
        raise RuntimeError("Original actual reported NP outfall source ID count changed")
    return o,gate

def benthic_sites(original):
    if gitblob(original)!=BQE_BLOB:
        raise RuntimeError("Original EEA2022 actual benthic source site file has changed")
    with original.open(encoding="utf-8-sig",newline="") as f:
        r=csv.DictReader(f)
        benthic={q["site_key"].split("|")[1]:q for q in r
            if q["countryCode"]=="AT" and q["site_key"].split("|")[2]=="euMonitoringSiteCode"}
    if len(benthic)!=68:
        raise RuntimeError("Original authentic 68 EEA2022 Austrian benthic site IDs changed")
    return benthic

def directed_graph(path):
    if sha(path)!=RIVER_SHA:
        raise RuntimeError("Original third-party HydroRIVERS Europe v1 GeoParquet public SHA changed")
    river=gpd.read_parquet(path)
    river=river.rename(columns={k:k.upper() for k in river.columns if k!=river.geometry.name})
    if "HYRIV_ID" not in river or "NEXT_DOWN" not in river:
        raise RuntimeError("Original HydroRIVERS native ID / flow direction schema changed")
    ids=river["HYRIV_ID"].astype(np.int64).to_numpy()
    if len(set(ids))!=len(ids):
        raise RuntimeError("Original river hydrology source has duplicated native reach IDs")
    n=river["NEXT_DOWN"].to_numpy()
    graph={}
    for k,v in zip(ids,n):
        try:child=int(v)
        except (ValueError,TypeError):child=-1
        graph[int(k)]=child
    return graph,{"original_HydroRIVERS_native_reach_ID_count":len(graph),
       "original_river_SHA256":RIVER_SHA}

def trace(start,target,network):
    seen=set()
    path=[]
    current=int(start)
    for hop in range(MAX_HOPS+1):
        if current==target:return {"within_reaches":True,
                "hops":hop,"path_head":path[:12],
                "source_river_line_PATH_tail":path[-8:]}
        if current<=0 or current not in network:return {"within_reaches":False,
                "outcome":"downstream_next_reach_unavailable","hops_to_exit":hop}
        if current in seen:return {"within_reaches":False,
                "outcome":"invalid_source_river_graph_loop","hops_to_exit":hop}
        seen.add(current)
        path.append(current)
        current=network[current]
    return {"within_reaches":False,"outcome":"source_max500_hops_no_directed_connection"}

def hav(lat,lon,lat2,lon2):
    a,b,c,d=map(math.radians,(float(lat),float(lon),float(lat2),float(lon2)))
    p=math.sin((a-c)/2)**2+math.cos(a)*math.cos(c)*math.sin((b-d)/2)**2
    return 12742017.6*math.asin(math.sqrt(min(1,max(0,p))))

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    source,receipt=source_outfalls(a.outfalls)
    benthic=benthic_sites(Path("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv"))
    directed,riverqa=directed_graph(a.rivers)
    records=[]
    for site,meta in AT_SITE_SOURCE.items():
        out=source[ORIG_STATION_TO_OUTFALL[site]]
        outlet_reach=int(out["source_GIS_original_river_snapped_reach"])
        record={
            "original_national_station_PKINT":site,
            "2022_total_unfiltered_phosphorus_P_mean_in_official_FW11":True,
            "official_BMLUK_2020_2022_FW11_source_sha256":"17de22544e9d043f83a928154847fb3fcd938894ce7007cfd06c000444b6eeb2",
            "original_EEA2022_reported_NP_outfall":ORIG_STATION_TO_OUTFALL[site],
            "original_EEA2022_outfall_river_native_HYRIV_ID":outlet_reach,
            "outfall_original_immediate_downstream_D1":int(out["source_true_graph_downstream_reach_id"]),
            "original_true_2022_benthic_EEA_same_PKINT_ID":meta["eea2022_benthic"],
            "original_2022_EEA_BQE_GZUEV_coordinate_discrepancy_m":meta["gov_EEA_benthic_station_coordinate_discrepancy_m"],
            "original_EEA_2022_source_annual_NP_report_not_river_concentration":True,
            "same_river_reach_network_path_NOT_inferrable_from_euclidean_distance":True}
        if site in benthic:
            actual=benthic[site]
            bio_reach=int(actual["HYRIV_ID"])
            record.update({
                "EEA2022_benthic_original_reach_native_HYRIV_ID":bio_reach,
                "original_WWTP_to_EEA2022_benthic_coordinate_distance_m":hav(
                    out["original_source_outfall_lat"],out["original_source_outfall_lon"],
                    actual["lat"],actual["lon"]),
                "true_benthic_2022_observation_coordinate_consistent_with_gov_geodata":
                   meta["gov_EEA_benthic_station_coordinate_discrepancy_m"]<=5})
            b_down=trace(outlet_reach,bio_reach,directed)
            b_up=trace(bio_reach,outlet_reach,directed)
            if outlet_reach==bio_reach:relationship="benthic_and_outfall_in_same_original_river_reach"
            elif b_down["within_reaches"]:relationship="benthic_SITE_is_hydrologically_DOWNSTREAM_of_outfall_source_reach"
            elif b_up["within_reaches"]:relationship="benthic_SITE_is_hydrologically_UPSTREAM_of_outfall_source_reach"
            else:relationship="different_source_tributary_or_river_no_direct_reach_in_500_hops"
            record["original_reach_level_hydrologic_relation"]=relationship
            record["source_outfall_to_benthic_reach_trace"]=b_down
            record["source_benthic_to_outfall_reach_trace"]=b_up
            # This is a reach-to-reach test, not within-reach along-line
            # movement or direct measured P at the outfall itself.
            record["qualifies_true_2022_receiver_reach_downstream_BQE"]=(
                relationship=="benthic_SITE_is_hydrologically_DOWNSTREAM_of_outfall_source_reach"
                and meta["gov_EEA_benthic_station_coordinate_discrepancy_m"]<=5)
        else:
            record["original_reach_level_hydrologic_relation"]="NO_EEA_2022_BENTHIC_SITE_ID_IN_ORIGINAL_FIELD_SOURCE"
            record["qualifies_true_2022_receiver_reach_downstream_BQE"]=False
        records.append(record)
    final={
      "schema":SCHEMA,
      "source_original_national_2022_annual_P_225_GZUEV_matched_sites":len(records),
      "source_original_annual_national_TN_not_present_in_FW11":True,
      "source_original_same_site_2022_benthic_ID_matched":sum(q["original_true_2022_benthic_EEA_same_PKINT_ID"] for q in records),
      "source_original_BQE_at_verified_downstream_river_reach":sum(q["qualifies_true_2022_receiver_reach_downstream_BQE"] for q in records),
      "source_2022_TP_from_official_national_yearbook_not_local_WWTP_P_apportionment":True,
      "original_benthic_EQR_values_NOT_READ":True,
      "new_local_Sentinel_NDVI_or_butterfly_values_NOT_READ":True,
      "same_reach_station_position_and_reference_unpolluted_reach_NOT_verified":True,
      "reaches_are_based_on_source_SHA_pinned_third_party_HydroRIVERS_v1":riverqa,
      "source_results":records,
      "not_eco_treatment_or_cross_boundary_mediation":True}
    (a.output/"BMLUK_FW11_2022_P_EEA_BENTHIC_2022_NP_OUTFALL_3SITE_DIRECTED_RIVER_QA.json").write_text(
       json.dumps(final,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_BMLUK_2022_P_BENTHIC_REAL_NETWORK_SOURCE="+json.dumps({
       "TP_2022_original_site":len(records),
       "source_benthic_same2022_PKINT":final["source_original_same_site_2022_benthic_ID_matched"],
       "source_downstream_benthic_reaches":final["source_original_BQE_at_verified_downstream_river_reach"],
       "classification":[{"site":q["original_national_station_PKINT"],
           "relation":q["original_reach_level_hydrologic_relation"],
           "coordinate_consistent":q.get("true_benthic_2022_observation_coordinate_consistent_with_gov_geodata"),
           "benthic_outfall_straight_m":q.get("original_WWTP_to_EEA2022_benthic_coordinate_distance_m")}
           for q in records],
       "no_causal_sewage_terrestrial_effect":True},sort_keys=True,ensure_ascii=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--outfalls",required=True,type=Path)
    p.add_argument("--rivers",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
