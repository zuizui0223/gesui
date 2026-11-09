#!/usr/bin/env python3
"""2022 observed TP MONITORING SITE SUPPORT within original DE/AT Lycaenidae 42 cells.

Peters 2022 actual measured TP station list was independently published in
gesui source-aligned river-chemistry PR9, fixed Git blob. The input CSV does
contain numeric TP measurements, but **this script never accesses, parses or
compares the numeric TP concentration/log-TP fields**. Original sample
days/months & station coordinates only. Country membership is *exact public
Natural Earth50m polygon*, not an inaccurate bounding box.

No ecological effect, adult butterfly puddling/abundance, or WWTP-caused TP
is fitted; declare support before reading numeric TP or doing bridge models.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,urllib.request
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from shapely.geometry import Point,shape

ORIGINAL_PETERS_COMMIT="6f037f65e9b948d6f8e324ad4eec9d8c03424b35"
ORIGINAL_PETERS_BLOB_SHA1="a75a4cf140de64fd6038079a3fac6756a461f3f8"
PETERS_TP="https://raw.githubusercontent.com/zuizui0223/gesui/"+ORIGINAL_PETERS_COMMIT+"/results/eea2022_measured_nutrient_load_river_20261107/TP_2022_reported_wwtp_load_and_measured_river_station.csv"
NATURAL_EARTH_URL="https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_0_countries.geojson"
NATURAL_EARTH_SHA256="3e458fc036ad0a66411f2c1e6cac49c5d7bfb81cb1123bc513b22511a2b7fdeb"
SOURCE_COUNTRY=("DE","AT")
MIN_STATIONS_PER_TILE=2
MIN_TILES_TOTAL=20
MIN_DE_TILES=10
MIN_AT_TILES=6
MIN_BLOCKS=10
MIN_TILES_WITH_REGISTERED_UPSTREAM_SOURCE=10
TP_YEARMONTH_MIN=3
TP_SAMPLES_MIN=4
OBS_EXPECTED=1291

def get(url):
    req=urllib.request.Request(url,headers={"User-Agent":"gesui-measured-river-TP-2022-source-only-geocoverage-QC/1.0"})
    with urllib.request.urlopen(req,timeout=95) as r:
        if r.status!=200:raise ValueError("Original public source download HTTP "+str(r.status))
        raw=r.read(4500000)
    return raw

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(1048576),b""):h.update(chunk)
    return h.hexdigest()

def gitblob_sha(b):
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def geom(raw):
    parsed=json.loads(raw.decode("utf8"))
    countries={}
    for c,a3 in (("DE","DEU"),("AT","AUT")):
        eligible=[f for f in parsed["features"] if
            (f.get("properties",{}).get("ISO_A2")==c or f.get("properties",{}).get("ADM0_A3")==a3)
             and f.get("geometry")]
        if len(eligible)!=1:raise ValueError("Original Natural Earth country mask ambiguous: "+c)
        countries[c]=shape(eligible[0]["geometry"])
    return countries

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--2022-frozen-worldcover-43cell",type=Path,required=True,dest="world")
    p.add_argument("--2023-same43-photo-source",type=Path,required=True,dest="previous_year")
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    world=json.loads(a.world.read_text())
    previous=json.loads(a.previous_year.read_text())
    if world["status"]!="POSTHOC_LANDCOVER_OUTOFREGION_PREDICTION_AND_PLACEBO_COMPLETED":
        raise ValueError("Original 2022 WorldCover source sample not frozen")
    if not previous.get("2023_outcome_year_source_gate_pass"):
        raise ValueError("Original same 2021–2023 country×1degree photographic cell source not fixed")
    if previous["original_2022_WorldCover_land_photo_source_sha256"]!=sha(a.world):
        raise ValueError("2021–2023 comparison doesn't match original 2022 land source")
    orig=world["all_original_same_sample_cells_2022_photo_counts_and_WorldCover_categorical_2021_covariates"]
    photo2023=previous["source_2023_all_original43_butterfly_photo_count_and_receipts"]
    prev={(r["countryCode"],int(float(r["lat_floor"])),int(float(r["lon_floor"]))):r for r in photo2023}
    strict=[r for r in orig if prev[(r["countryCode"],int(float(r["lat_floor"])),int(float(r["lon_floor"])))]["2023_photo_support_ge20"]]
    if len(strict)!=42 or Counter(r["countryCode"] for r in strict)!=Counter({"DE":27,"AT":15}):
        raise ValueError("Posthoc strict three-year source geography changed")
    keys={(r["countryCode"],int(float(r["lat_floor"])),int(float(r["lon_floor"]))):r for r in strict}
    if len(keys)!=42:raise ValueError("Duplicate fixed country-1degree photo source cells")

    origbytes=get(PETERS_TP)
    if gitblob_sha(origbytes)!=ORIGINAL_PETERS_BLOB_SHA1:
        raise ValueError("Peters TP 2022 original measured-station manifest Git SHA1 not verified")
    natural=get(NATURAL_EARTH_URL)
    if hashlib.sha256(natural).hexdigest()!=NATURAL_EARTH_SHA256:
        raise ValueError("Natural Earth50m country polygon original prior WorldCover source SHA changed")
    boundaries=geom(natural)
    lines=origbytes.decode("utf8-sig").splitlines()
    headers=lines[0].split(",")
    expected=["nutrient","WQ.ID","WQ.lat","WQ.lon",
        "HYRIV_ID","spatial_3deg_block",
        "actual_2022_observation_count","distinct_2022_months","distinct_2022_dates",
        "eea2022_upstream_wwtp_both_positive_count"]
    if any(x not in headers for x in expected):raise ValueError("Original TP 2022 metadata schema drift")
    # STRICTLY DO NOT PARSE OR ACCESS:
    # real_2022_station_median_mgL, log_real_2022_median_mgL, and P source mass fields.
    stationrows=[];seen=set();status=Counter()
    reader=csv.DictReader(lines)
    for s in reader:
        status["source_rows"]+=1
        if s["nutrient"]!="TP":raise ValueError("Original 2022 stations include non-TP nutrient")
        sid=s["WQ.ID"]
        if sid in seen:raise ValueError("Duplicate original TP station ID")
        seen.add(sid)
        try:
            lat=float(s["WQ.lat"]);lon=float(s["WQ.lon"])
            months=int(s["distinct_2022_months"])
            dates=int(s["distinct_2022_dates"])
            count=int(s["actual_2022_observation_count"])
            source_count=int(s["eea2022_upstream_wwtp_both_positive_count"])
        except (TypeError,ValueError):raise ValueError("Original station-only TP metadata NA")
        if count<TP_SAMPLES_MIN or months<TP_YEARMONTH_MIN or dates<TP_SAMPLES_MIN:
            raise ValueError("Original previously-qualified 2022 station no longer meets measurement SUPPORT")
        matched=[c for c,g in boundaries.items() if g.covers(Point(lon,lat))]
        if len(matched)>1:status["station_on_ambiguous_country_boundary"]+=1;continue
        if not matched:status["station_outside_DE_AT_exact_border"]+=1;continue
        c=matched[0];status["exact_country_"+c+"_all_2022_TP_stations"]+=1
        tile=(c,math.floor(lat),math.floor(lon))
        if tile not in keys:
            status["station_outside_original_42_photo_source_cells"]+=1
            continue
        status["station_in_original_42_country_masked_photo_cells"]+=1
        stationrows.append({
          "WQ_station_ID":sid,"countryCode":c,"lat_WGS84":lat,"lon_WGS84":lon,
          "original_country_1deg_lat_floor":tile[1],
          "original_country_1deg_lon_floor":tile[2],
          "source_2022_original_TP_sample_dates_n":dates,
          "source_2022_original_TP_sample_months_n":months,
          "source_2022_original_TP_valid_sample_count":count,
          "original_hydro_river_reach_ID":s["HYRIV_ID"],
          "registered_source_upstream_bothNP_2022_count_original_50km":source_count,
          "numeric_total_P_measurement_not_read":True,
        })
    if status["source_rows"]!=OBS_EXPECTED:
        raise ValueError("Original Peters real 2022 TP monitored site source count unexpectedly changed")
    siteby=defaultdict(list)
    for q in stationrows:
        siteby[(q["countryCode"],q["original_country_1deg_lat_floor"],q["original_country_1deg_lon_floor"])].append(q)
    tilelist=[]
    for k,r in sorted(keys.items()):
        entries=siteby[k]
        tilelist.append({
            "countryCode":k[0],"lat_floor":k[1],"lon_floor":k[2],
            "original_photo_years":["2021","2022","2023"],
            "original_country_2degree_block":r["country_2deg_block"],
            "country_masked_distinct_actual_2022_TP_station_n":len(entries),
            "country_masked_TP_stations_with_reported_upstream_EEA_WWTP_source_n":sum(e["registered_source_upstream_bothNP_2022_count_original_50km"]>0 for e in entries),
            "station_metadata_support_ge2":len(entries)>=MIN_STATIONS_PER_TILE,
            "NO_2022_TP_NUMERICAL_VALUES_READ":True,
        })
    supported=[t for t in tilelist if t["station_metadata_support_ge2"]]
    byc=Counter(t["countryCode"] for t in supported)
    blocks=set(t["countryCode"]+"|"+str(t["original_country_2degree_block"]) for t in supported)
    exposed_tiles=sum(t["country_masked_TP_stations_with_reported_upstream_EEA_WWTP_source_n"]>0 for t in supported)
    gate=(len(supported)>=MIN_TILES_TOTAL and byc.get("DE",0)>=MIN_DE_TILES and
          byc.get("AT",0)>=MIN_AT_TILES and len(blocks)>=MIN_BLOCKS and
          exposed_tiles>=MIN_TILES_WITH_REGISTERED_UPSTREAM_SOURCE)
    output={
      "schema":"gesui_original_Peters2022_TP_metadata_site_support_for_2021_2023_Lycaenidae_DE_AT_country_cells_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "original_2022_Peters_TP_qualified_station_file":PETERS_TP,
      "original_Peters_TP_2022_Gitblob_SHA1":ORIGINAL_PETERS_BLOB_SHA1,
      "original_NaturalEarth50m_country_polygon_URL":NATURAL_EARTH_URL,
      "original_NaturalEarth50m_country_polygon_SHA256":NATURAL_EARTH_SHA256,
      "original_2022_WorldCover_full_43cell_source_SHA256":sha(a.world),
      "original_2023_photo_sample_full_source_SHA256":sha(a.previous_year),
      "original_2021_2022_2023_common_country1degree_photo_source_cells":len(keys),
      "original_TP_sites_country_all_2022_prequalified":dict(status),
      "country_masked_TP_stations_inside_frozen_photo_cells_n":len(stationrows),
      "exact_country_masked_original_2022_TP_station_site_metadata":stationrows,
      "all_original42_country1degree_cell_TP_monitoring_site_support":tilelist,
      "min_original_distinct_TP_monitoring_stations_per_cell":MIN_STATIONS_PER_TILE,
      "min_supported_cells_total":MIN_TILES_TOTAL,
      "min_supported_DE_cells":MIN_DE_TILES,
      "min_supported_AT_cells":MIN_AT_TILES,
      "min_country2degree_blocks":MIN_BLOCKS,
      "min_supported_tiles_with_original_upstream_registered_WWTP_2022_source":MIN_TILES_WITH_REGISTERED_UPSTREAM_SOURCE,
      "supported_cells_ge2_original_2022_TP_stations_n":len(supported),
      "supported_cells_by_country":dict(byc),
      "supported_distinct_country2degree_blocks":len(blocks),
      "supported_tiles_with_observed_TP_and_reported_upstream_WWTP_source":exposed_tiles,
      "three_layer_river_TP_chemistry_plus_Lycaenidae_photo_source_gate":gate,
      "strictly_no_numeric_TP_or_other_chemical_outcomes_accessed_in_source_selection":True,
      "original_2022_TP_station_file_does_contain_numeric_outcomes_but_they_are_not_used":True,
      "no_new_photo_data_or_observer_personal_data_read":True,
      "if_failed_do_not_fit_photorecord_vs_2022_TP_model":True,
      "no_wastewater_ecological_effect_or_pollution_mechanism_inferred":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_PETERS2022_ORIGINAL_TP_STATION_EXACT_COUNTRY_3LAYER_SUPPORT_GATE",
      json.dumps({"2022_TP_station_country_source_counts":dict(status),
        "source_only_country1degree_matched_stations":len(stationrows),
        "source_phototiles":len(keys),"tiles_with_2_stations":len(supported),
        "by_country":dict(byc),"2deg_blocks":len(blocks),
        "river_site_exposed_tiles":exposed_tiles,"gate_pass":gate,
        "no_2022_TP_values_used":True},sort_keys=True),flush=True)

if __name__=="__main__":main()
