#!/usr/bin/env python3
"""Freeze ALL original 1055 WPI 42country WWp ports, original 247 pilot subset.

Uses no new 808 GBIF plant outcomes. Original old candidate-host source model
must exactly reproduce old beta and country-robust CI. Global ALL-Plantae
species inventories at 808 source-frozen ports are collected only later.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from datetime import datetime,timezone
from pathlib import Path

from gesui_global_WWp2015_original_1055port_country_transfer import build_original_sample

ROOT_SOURCE=Path("results/global_all_Plantae_WWp2015_20261110")
PILOT=ROOT_SOURCE/"ORIGINAL_42COUNTRY_WWp_STRATIFIED_PORTS_PRE_GBIF_ALLFLORA_SOURCE.json"
PILOT_RESULT=ROOT_SOURCE/"ORIGINAL_247PORT_42COUNTRY_GBIF_ALL_PLANTAE_2011_24_RAREFACTION_SOURCE.json"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(root,jones,world):
    data,old_beta=build_original_sample(root,jones,world)
    pilot=json.loads((root/PILOT).read_text())
    plant=json.loads((root/PILOT_RESULT).read_text())
    if plant["source_only_study_manifest_sha256"]!=sha(root/PILOT):
        raise ValueError("Old 247 pilot all-Plantae outcome not on same source freeze")
    selected={x["original_port_id"]:x for x in pilot["all_pilot_port_georeferenced_original_sites"]}
    previous={x["original_port_id"]:x for x in plant["original_site_GBIF_facet_summary"]}
    if len(selected)!=247 or set(selected)!=set(previous):
        raise ValueError("Original pilot 247 sites changed before all1055 expansion")
    rows=[]
    for r in data.sort_values(["country","port_id"]).itertuples():
        row={
            "original_port_id":int(r.port_id),
            "original_source_country":str(r.country),
            "original_port_lat":float(r.latitude),
            "original_port_lon":float(r.longitude),
            "source_WWp2015_million_m3_per_year_within25km":float(r.produced_25km_millionm3yr),
            "source_WWp_log1p":float(r.ln_produced),
            "source_WWTP_count_25km":float(r.WWTP_count_25km),
            "source_ESA2021_built_fraction":float(r.built_fraction_2021),
            "source_ESA2021_crop_fraction":float(r.crop_fraction_2021),
            "previously_extracted_247_AllPlantae_pilot_source_port":int(r.port_id) in selected,
            "all_808_additional_Plantae_2011_24_GBIF_queries_NOT_PERFORMED_in_stage0":True
        }
        if row["original_port_id"] in selected:
            s=selected[row["original_port_id"]]
            fields=("original_source_country","original_port_lat","original_port_lon",
                    "source_WWp_log1p","source_WWp2015_million_m3_per_year_within25km",
                    "source_WWTP_count_25km","source_ESA2021_built_fraction","source_ESA2021_crop_fraction")
            for k in fields:
                if row[k]!=s[k]:
                    raise ValueError("Original pilot source WPI2019 port/2015 WWp/2021 WorldCover changed "+k)
        rows.append(row)
    if len(rows)!=1055 or len({x["original_port_id"] for x in rows})!=1055 or len({x["original_source_country"] for x in rows})!=42:
        raise ValueError("Source census not 1055 sites in 42 original countries")
    new=[x["original_port_id"] for x in rows if not x["previously_extracted_247_AllPlantae_pilot_source_port"]]
    if len(new)!=808:
        raise ValueError("Expected exactly 808 previously unqueried original country port IDs")
    out={
        "schema":"gesui_original_1055_WPI_ports42country_808NEW_2011_24_ALLPlantae_PREOUTCOME_frame_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_1596_WPI_source_jones2015_csv_SHA256":sha(jones),
        "original_1596_WPI_source_ESA2021_worldcover_csv_SHA256":sha(world),
        "original_previous_247_pilot_source_SHA256":sha(root/PILOT),
        "original_previous_247_pilot_GBIF_all_Plantae_facet_summary_SHA256":sha(root/PILOT_RESULT),
        "old_positive_candidate_host_countryFE_cohort_and_coefficient_exact_replay":old_beta,
        "all_original_country_WWp_source_sites_n":len(rows),
        "original_country_count":42,
        "previous_247_fixed_pilot_site_n":247,
        "new_allPlantae_2011_24_GBIF_port_IDs_original_not_yet_queried_n":808,
        "source_all_1055_public_country_port_site_metadata_no_new_GBIF_response":rows,
        "source_GBIF_query_all_Plantae":{"taxonKey":6,"year":"2011,2024",
               "geoDistance_radius_km":5,"hasCoordinate":"true",
               "hasGeospatialIssue":"false","occurrenceStatus":"PRESENT",
               "facet":"speciesKey","facetLimit":1000,"limit":0},
        "selection_determined_by_original_1055_port_candidate_host_model_source_not_any_new808_all_Plantae_results":True,
        "this_scale_up_is_POSTHOC_after_positive_raw_all_Plantae_247pilot_and_null_rarefied50":True,
        "GBIF_new808_species_records_NOT_OBSERVED_in_stage0":True,
        "2015_modelled_WWp_generated_NOT_effluent_discharge_or_NP_nutrient_flux":True,
    }
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--jones",type=Path,required=True)
    p.add_argument("--worldcover",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    o=run(root,a.jones,a.worldcover)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(o,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"original_old_positive":o["old_positive_candidate_host_countryFE_cohort_and_coefficient_exact_replay"],
                      "original_country_ports":o["all_original_country_WWp_source_sites_n"],
                      "new_allflora_port_queries_allowed":o["new_allPlantae_2011_24_GBIF_port_IDs_original_not_yet_queried_n"]},indent=2))
if __name__=="__main__":
    main()
