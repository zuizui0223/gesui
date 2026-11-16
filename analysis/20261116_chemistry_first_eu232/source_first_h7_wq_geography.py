#!/usr/bin/env python3
"""SOURCE-ONLY new water-chemistry-first geography, no Ecology outcomes.

Map ALL genuine official 2026 EEA Water Quality stations with valid WGS84
coordinates onto the ORIGINAL genuine 232 HydroBASINS Level7 polygons
whose GBIF plant/butterfly original monitoring opportunities are already
frozen. Unlike historical biology-station-aligned analyses, select WQ
stations first: NO EEA benthic EQR site match, proximity radius, positive
WWTP report, NUTRI/HNS score, or observed TN/TP value selects them.

Target years 2018/2022 are original GBIF 2017-19/2020-22 cross-sectional
windows, NOT true year-matched observed plant leaves or butterfly feeding.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,sys,collections
from pathlib import Path
from shapely.geometry import Point
from shapely.strtree import STRtree

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_reported_np_source"))
from audit_biology_wq_spatial_source_v1 import wq_sites,WQ_ZIP_SHA
sys.path.insert(0,str(ROOT/"analysis/20261115_h6_source_expansion"))
from map_official_eu501_to_h6 import fetch,geom_records,L7_ZIP_SHA,PANEL_SHA

PROTOCOL="GESUI_EEA2026_CHEMISTRY_FIRST_ORIGINAL_EU232_H7_SOURCE_ONLY_V1"
YEARS={2018:"2017-2019",2022:"2020-2022"}
MIN_SUPPORT={"H7period":30,"distinct_H7":20,"true_H6":20,"country_groups":3}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def panel_geography():
    p=ROOT/"results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"
    if sha(p)!=PANEL_SHA:
        raise RuntimeError("Original 454 H7 study source SHA mismatch")
    with p.open(newline="",encoding="utf-8-sig") as f:
        records=list(csv.DictReader(f))
    if len(records)!=454 or len({r["h"] for r in records})!=232 or len({r["h6"] for r in records})!=90:
        raise RuntimeError("Original 454 source panel geography cardinality changed")
    idx={}
    for r in records:
        k=(r["h"],r["period"])
        if k in idx:raise RuntimeError("Original H7 period duplicated")
        idx[k]={"h6":r["h6"],"country":r["country"]}
    all_parents=collections.defaultdict(set)
    for r in records:
        all_parents[r["h"]].add(r["h6"])
    if any(len(v)!=1 for v in all_parents.values()):
        raise RuntimeError("Genuine original H6 parent changed across original periods")
    return idx,set(all_parents)

def main(args):
    args.output.mkdir(parents=True,exist_ok=True)
    index,original_h7=panel_geography()
    water,original_wq=wq_sites(args.output)
    raw,geometry_sha=fetch(7)
    if geometry_sha["sha256"]!=L7_ZIP_SHA:
        raise RuntimeError("Original official HydroBASINS source ZIP changed")
    h7,polys,geoqa=geom_records(raw,7,original_h7)
    if set(h7)!=original_h7 or len(h7)!=232 or geoqa["n_unusable"]:
        raise RuntimeError("Missing source-native H7 geometry")
    tree=STRtree(polys)
    per_station=[]
    invalid_border=0
    outside=0
    for wqkey,coord in sorted(water.items()):
        if len(wqkey)!=3 or not all(wqkey):
            raise RuntimeError("Water Quality source ID missing country/site/scheme")
        # STRtree predicate operates as Point.covered_by(polygon), includes boundaries.
        p=Point(coord[1],coord[0])
        candidates=[int(i) for i in tree.query(p,predicate="covered_by")]
        if len(candidates)>1:
            invalid_border+=1
            continue
        if len(candidates)==0:
            outside+=1
            continue
        h=h7[candidates[0]]
        if not polys[candidates[0]].covers(p):
            raise RuntimeError("Native polygon gate self-disagreement")
        per_station.append((wqkey,h))
    selections=[]
    for (country,station,scheme),h in per_station:
        for year,period in YEARS.items():
            g=index.get((h,period))
            if g is None:continue
            selections.append({
                "source_WQ_country":country,
                "source_WQ_monitoringSiteIdentifier":station,
                "source_WQ_monitoringSiteIdentifierScheme":scheme,
                "source_WQ_calendar_year":year,
                "original_GBIF_source_H7_period":period,
                "official_original_HYBAS_ID_L07":h,
                "original_true_HYBAS_ID_L06":g["h6"],
                "source_original_454_majority_country":g["country"],
                "WQ_station_country_equals_original_GBIF_basin_majority_country":country==g["country"],
                "same_official_H7_polygon":True,
                "same_EEA_Biology_EQR_station_NOT_REQUIRED":True,
                "original_2026_Part2_TNTP_values_NOT_YET_OPENED":True,
                "no_actual_terrestrial_N_uptake_or_sewage_attribution":True
            })
    ids=[(q["source_WQ_country"],q["source_WQ_monitoringSiteIdentifier"],
          q["source_WQ_monitoringSiteIdentifierScheme"],q["source_WQ_calendar_year"]) for q in selections]
    if len(ids)!=len(set(ids)):raise RuntimeError("WQ station-year assigned twice to original H7")
    selected=sorted(selections,key=lambda r:(r["source_WQ_country"],
        r["source_WQ_monitoringSiteIdentifier"],r["source_WQ_monitoringSiteIdentifierScheme"],
        r["source_WQ_calendar_year"]))
    primary=[r for r in selected if r["WQ_station_country_equals_original_GBIF_basin_majority_country"]]
    def counts(rows):
        return {"original_WQ_station_year_candidates":len(rows),
                "distinct_original_WQ_stations":len({
                  (r["source_WQ_country"],r["source_WQ_monitoringSiteIdentifier"],
                   r["source_WQ_monitoringSiteIdentifierScheme"]) for r in rows}),
                "source_H7period_opportunities":len({(r["official_original_HYBAS_ID_L07"],
                                                     r["original_GBIF_source_H7_period"]) for r in rows}),
                "source_H7_distinct":len({r["official_original_HYBAS_ID_L07"] for r in rows}),
                "source_true_H6_distinct":len({r["original_true_HYBAS_ID_L06"] for r in rows}),
                "source_WQ_country_groups":len({r["source_WQ_country"] for r in rows})}
    path=args.output/"ORIGINAL_2026_EEA_WQ_STATIONS_IN_ORIGINAL_232_H7_SOURCE_ONLY.csv"
    if not selected:raise RuntimeError("Official water-quality source has no overlap with original H7 geography")
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(selected[0]));w.writeheader();w.writerows(selected)
    result={
        "schema":PROTOCOL,
        "source_original_H7panel_SHA256":PANEL_SHA,
        "source_original_WQ_spatial_ZIP_SHA256":WQ_ZIP_SHA,
        "source_original_HydroBASINS_L7_ZIP_SHA256":L7_ZIP_SHA,
        "official_spatial_H7_QA":geoqa,
        "source_original_2026_water_sites_QA":{
            "original_WQ_spatial_site_keys_with_unambiguous_coords":len(water),
            "original_WQ_geographic_source_rows":original_wq["EEA_WQ_spatial_raw_rows"]},
        "station_point_outside_original_selected_H7":outside,
        "station_on_ambiguous_polygon_border_excluded":invalid_border,
        "station_inside_one_original_H7_polygon":len(per_station),
        "full_source_any_country":counts(selected),
        "PRIMARY_WQ_country_matches_original_GBIF_basin_country":counts(primary),
        "primary_country_match_only":True,
        "source_candidate_CSV_sha256":sha(path),
        "source_quality_threshold_for_measured_TN_and_TP_both":MIN_SUPPORT,
        "source_stage": "GEOMETRY_ONLY_OUTCOMES_NOT_OPENED",
        "not_same_biological_EQR_station":True,
        "not_sewage_exposure_or_river_reach_verified":True}
    (args.output/"CHEMISTRY_FIRST_232H7_WQ_GEOGRAPHY_GATE.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_CHEMISTRY_FIRST_H7_GEO_GATE="+json.dumps(result,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
