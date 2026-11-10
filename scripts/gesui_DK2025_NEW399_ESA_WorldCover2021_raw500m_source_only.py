#!/usr/bin/env python3
"""SOURCE-ONLY ESA WorldCover2021 genuine 10m pixel habitat around fixed Denmark eBMS NEW399.

No butterfly individualCount/taxonomy/outcome file read. Source cohort and site
ID membership immutable from archived 2025 pre-butterfly, pre-2025-new399 source gate.
"""
from __future__ import annotations

import argparse,hashlib,importlib.util,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import rasterio

BASE="results/ebms_denmark_2023_2024_WWTP_source_only/"
SOURCE=BASE+"DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json"
EVENT=BASE+"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
LANDLIB="archive/chocho_pr122/scripts/exploratory_port_worldcover_urban.py"
RADIUS=500
KNOWN={10:"tree",20:"shrub",30:"grass",40:"crop",50:"built",60:"bare",70:"snow_ice",
       80:"permanent_water",90:"wetland",95:"mangrove",100:"moss_lichen"}
NONWATER=set(KNOWN)-{80}

def sha(b):return hashlib.sha256(b).hexdigest()

def load_module(path):
    spec=importlib.util.spec_from_file_location("pinned_worldcover_lib",path)
    if spec is None or spec.loader is None:raise ValueError("Missing WorldCover v200 source helper")
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    if m.CLASSES!=KNOWN or m.RADIUS_M!=RADIUS:
        raise ValueError("WorldCover raw categorical codebook or source version differs")
    return m

def run(root):
    source_file=root/SOURCE
    event_file=root/EVENT
    libfile=root/LANDLIB
    src=json.loads(source_file.read_text(encoding="utf-8"))
    evt=json.loads(event_file.read_text(encoding="utf-8"))
    if src["schema"]!="gesui_DK_eBMS_2025_NEW399_original_event_and_any_occurrence_link_source_only_v01":
        raise ValueError("Expected unchanged 2025 source-first survey frame")
    if not src["frozen_before_2025_NEW399_butterfly_counts_gate"]["2025_NEW399_butterfly_family_count_outcome_permitted"]:
        raise ValueError("Original source-qualified 2025 eBMS frame does not permit analysis")
    if not src["NEW2025_Lycaenidae_or_any_butterfly_family_individualCount_NOT_READ"]:
        raise ValueError("Earlier source-eligibility cohort contaminated by outcomes")
    if src["original_2023_24_Event_source_only_735site_SHA256"]!=sha(event_file.read_bytes()):
        raise ValueError("Original public Event-only source has changed")
    if len(evt["original_2023_2024_repeated_event_site_source_metadata"])!=735:
        raise ValueError("Expected 735 original no-species Event-only source site locations")
    all_sites={s["site_id_SHA256"]:s for s in evt["original_2023_2024_repeated_event_site_source_metadata"]}
    records=src["all_original_399_source_only_new_site_2025_documented_visit_metadata"]
    if len(records)!=399 or len({s["site_id_SHA256"] for s in records})!=399:
        raise ValueError("Original source cohort not 399 distinct public sites")
    subset=sorted((r for r in records if r["pre_2025_family_any_occurrence_linked_site_source_gate_pass"]),key=lambda r:r["site_id_SHA256"])
    if len(subset)!=286 or Counter(r["source_band"] for r in subset)!=Counter({"near5":70,"remote10":92,"intermediate":124}):
        raise ValueError("Source-only documented site cohort changed")
    world=load_module(libfile)
    rows=[]
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                      CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                      GDAL_HTTP_MAX_RETRY="3",GDAL_HTTP_RETRY_DELAY="2",
                      GDAL_HTTP_CONNECTTIMEOUT="30",GDAL_HTTP_TIMEOUT="90",
                      GDAL_CACHEMAX=256):
        with world.SourceTiles() as tilelib:
            for i,s in enumerate(subset,1):
                base=all_sites.get(s["site_id_SHA256"])
                if base is None:raise ValueError("Original public source point absent")
                lat,lon=float(base["lat"]),float(base["lon"])
                if sha(f"{lat:.5f}|{lon:.5f}".encode())!=s["site_id_SHA256"]:
                    raise ValueError("Original public coordinate SHA mismatch")
                if abs(float(base["nearest_qualified_EEA2022_DK_FW_outfall_km"])-float(s["source_distance_km"]))>1e-7:
                    raise ValueError("Original source location/proximity changed")
                expected=world.window_tiles(lon,lat,RADIUS)
                counts=Counter()
                missing=[]
                total=0
                for tile in expected:
                    raster=tilelib.get(tile)
                    if raster is None:
                        missing.append({"tile":tile,"error":tilelib.errors.get(tile,"unavailable")[:250]})
                        continue
                    pixels=world.disk_stat(raster,lon,lat,RADIUS)
                    if pixels["tiles_overlap"]:
                        counts.update(pixels["counts"])
                        total+=pixels["circle_pixels"]
                nland=sum(counts[v] for v in NONWATER)
                unknown=sum(n for c,n in counts.items() if c not in KNOWN)
                water=counts[80]
                complete=not missing and total>=500 and nland>=200 and unknown==0
                frac={name:(counts[code]/nland if complete else None) for code,name in KNOWN.items() if code!=80}
                frac["water_fraction_of_all_circle_pixels"]=water/total if complete and total else None
                frac["grass_plus_tree_fraction_of_land"]=((counts[30]+counts[10])/nland if complete else None)
                result={
                    "source_original_public_site_SHA256":s["site_id_SHA256"],
                    "source_band":s["source_band"],
                    "original_1degree_region":s["original_1degree_region"],
                    "source_near_WWTP_id_SHA256":s["nearest_original_2022_EEA_outfall_id_SHA256"] if s["source_band"]=="near5" else None,
                    "source_near_WWTP_5km_component_SHA256":s.get("5km_original_WWTP_spatial_component_SHA256"),
                    "source_2022_qualifying_EEA_outfall_distance_km":s["source_distance_km"],
                    "original_2025_any_published_Occurrence_linked_dates":s["original_linked_occurrence_dates"],
                    "original_2025_all_Event_survey_dates":s["original_event_date_count"],
                    "WorldCover2021_raw_categorical_10m_radius_m":RADIUS,
                    "WorldCover2021_original_tiles":expected,
                    "original_COG_missing_source_tiles":missing,
                    "original_circle_pixel_centers_total":total,
                    "original_nonwater_known_terrestrial_pixel_centers":nland,
                    "original_water_class_80_pixel_centers":water,
                    "original_unknown_or_nodata_pixels":unknown,
                    "categorical_pixel_count":{str(c):int(counts[c]) for c in sorted(KNOWN)},
                    "fractions":frac,
                    "source_QC":"COMPLETE" if complete else "INCOMPLETE",
                    "no_butterfly_taxon_or_counts_accessed_for_landcover":True,
                    "year2021_landcover_not_2025_larval_host_abundance":True,
                }
                rows.append(result)
                if i%25==0 or i==len(subset):
                    print("WORLD_COVER_SOURCE_SITE",i,"/",len(subset),
                          "nCompleted",sum(x["source_QC"]=="COMPLETE" for x in rows),
                          "lastTile",",".join(expected),
                          "status",result["source_QC"],flush=True)
    if len(rows)!=286:raise ValueError("Unexpected WorldCover sample")
    complete=[x for x in rows if x["source_QC"]=="COMPLETE"]
    bands=Counter(x["source_band"] for x in complete)
    geoclass=defaultdict(set)
    for r in complete:geoclass[r["original_1degree_region"]].add(r["source_band"])
    shared=sum("near5" in v and "remote10" in v for v in geoclass.values())
    good=bool(bands["near5"]>=40 and bands["remote10"]>=60 and shared>=4)
    out={
        "schema":"gesui_DK2025_new399_source286_worldcover2021_raw500m_categorical_QC_BEFORE_habitat_family_adjustment_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "all_286_original_source_selected_prior_to_this_external_habitat_measurement":True,
        "original_source_only_new399_sha256":sha(source_file.read_bytes()),
        "original_Event_only_public_735site_sha256":sha(event_file.read_bytes()),
        "original_ESA_WorldCover2021_classification_helper_sha256":sha(libfile.read_bytes()),
        "source_original_ESA_WorldCover_raster_provider":world.SOURCE,
        "ESA_WorldCover2021_version":world.VERSION,
        "ESA_WorldCover2021_original_pixel_class_codes":KNOWN,
        "source_radius_m":RADIUS,
        "original_286_source_eligible_monitoring_site_count":len(rows),
        "2021_worldcover_valid_source_sites":len(complete),
        "2021_worldcover_valid_source_sites_by_original_source_band":dict(bands),
        "one_degree_regions_with_both_near_remote_and_valid_2021_worldcover":shared,
        "fixed_before_habitat_adjustment_quality_gate":{
            "minimum_near_source_complete_raster_sites":40,
            "minimum_remote_complete_raster_sites":60,
            "minimum_1degree_regions_with_both_complete":4,
            "habitat_adjusted_2025_butterfly_analysis_may_proceed":good},
        "original_site_geo_worldcover_source_counts":rows,
        "raster_tiles_NOT_redistributed":True,
        "no_2025_butterfly_families_or_numbers_in_this_source":True,
        "2021_static_landcover_NOT_2025_host_plants_or_actual_WWTP_nutrient_receptor":True,
    }
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    v=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    x=run(root)
    v.out.parent.mkdir(parents=True,exist_ok=True)
    v.out.write_text(json.dumps(x,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"source_sites":x["original_286_source_eligible_monitoring_site_count"],
                      "landcover_complete":x["2021_worldcover_valid_source_sites"],
                      "complete_by_band":x["2021_worldcover_valid_source_sites_by_original_source_band"],
                      "shared_1degree_regions":x["one_degree_regions_with_both_near_remote_and_valid_2021_worldcover"],
                      "habitat_family_analysis_eligible":x["fixed_before_habitat_adjustment_quality_gate"]},indent=2))

if __name__=="__main__":
    main()
