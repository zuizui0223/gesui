#!/usr/bin/env python3
"""Source-ONLY freeze: fixed all-Plantae GBIF pilot across original 42 WWp countries.

Original port table is NOT a new all-plant outcome. Choose sites only from
original 2015 PANGAEA WW-production distribution within each original
country, never from GBIF plant counts, candidate-host outcome or coefficient.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from gesui_global_WWp2015_original_1055port_country_transfer import build_original_sample

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def freeze(root,jones,world):
    data,previous=build_original_sample(root,jones,world)
    rows=[]
    summary={}
    for country,group in data.groupby("country",sort=True):
        group=group.sort_values(["ln_produced","port_id"]).reset_index(drop=True)
        n=len(group)
        if n<5:raise ValueError("Original 42country eligible sample violated pre-existing n>=5")
        take=min(6,n)
        indices=np.linspace(0,n-1,take).round().astype(int).tolist()
        if len(indices)!=len(set(indices)):
            raise ValueError("Original WWp quantile indices not unique for frozen site choice")
        chosen=group.iloc[indices]
        summary[str(country)]={
            "original_matched_source_port_n":n,
            "fixed_chosen_ports_n":take,
            "source_WWp_min":float(group.produced_25km_millionm3yr.min()),
            "source_WWp_max":float(group.produced_25km_millionm3yr.max()),
            "original_within_country_exposure_variation_n":int(group.ln_produced.nunique()),
        }
        for original_rank,i in enumerate(indices):
            site=group.iloc[i]
            lat=float(site.latitude);lon=float(site.longitude)
            if not(-90<=lat<=90 and -180<=lon<=180 and math.isfinite(lat) and math.isfinite(lon)):
                raise ValueError("WPI source coordinate missing or invalid")
            rows.append({
                "original_port_id":int(site.port_id),
                "original_source_country":str(country),
                "original_port_lat":lat,
                "original_port_lon":lon,
                "source_WWp2015_million_m3_per_year_within25km":float(site.produced_25km_millionm3yr),
                "source_WWp_log1p":float(site.ln_produced),
                "source_WWTP_count_25km":float(site.WWTP_count_25km),
                "source_ESA2021_built_fraction":float(site.built_fraction_2021),
                "source_ESA2021_crop_fraction":float(site.crop_fraction_2021),
                "source_country_WWp_order_rank":original_rank,
                "source_original_country_port_count":n,
                "GBIF_all_Plantae_OUTCOME_NOT_FETCHED":True
            })
    if len(rows)!=len({x["original_port_id"] for x in rows}):
        raise ValueError("Duplicate original public port IDs after source-only sampling")
    if len(summary)!=42 or not(220<=len(rows)<=252):
        raise ValueError("Source country-supported full 42country pilot not fixed")
    sites=sorted(rows,key=lambda x:(x["original_source_country"],x["source_country_WWp_order_rank"],x["original_port_id"]))
    return {
        "schema":"gesui_WPI2019_42country_WWp2015_selected6percountry_AllPlantae_2011_24_PREOUTCOME_FRAME_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_source_2015_PANGAEA_WWp_1596port_jonesCSV_SHA256":sha(jones),
        "original_source_HydroWASTE_WPI1596_WorldCover2021_portCSV_SHA256":sha(world),
        "original_1055port_WWp_countryFE_positive_model_verified_before_selecting_new_allflora":previous,
        "source_selection_order":"IN_EACH_COUNTRY sort original Jones2015 log1p WWproduction25km then original port_ID; retain six equally spaced indices including both original ends, or all five for countries with five eligible ports",
        "source_selection_is_independent_of_all_new_GBIF_Plantae_records":True,
        "original_WPI2019_full_site_n":1596,
        "original_WWp2015_city_countryFE_site_n":1055,
        "original_WWp2015_countryFE_eligible_country_n":42,
        "new_all_Plantae_GBIF_pilot_site_n":len(sites),
        "countries_in_source_pilot":len(summary),
        "original_countries_and_WWp_range":summary,
        "all_pilot_port_georeferenced_original_sites":sites,
        "GBIF_new_all_Plantae_requests_in_this_stage":0,
        "precommitted_new_GBIF_taxonKingdomKey":6,
        "precommitted_GBIF_window":"2011,2024",
        "precommitted_original_port_search_radius_km":5,
        "precommitted_GBIF_filters":{"hasCoordinate":"true","hasGeospatialIssue":"false","occurrenceStatus":"PRESENT"},
        "global_source_wastewater_is_2015_MODELED_production_not_measured_effluent":True,
        "original_plant_response_is_2001_24_candidate_butterfly_HOSTS_not_all_plants":True,
        "no_new_allplant_positive_coefficient_claimed":True
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--jones",type=Path,required=True)
    a.add_argument("--worldcover",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    x=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    result=freeze(root,x.jones,x.worldcover)
    x.out.parent.mkdir(parents=True,exist_ok=True)
    x.out.write_text(json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"cities":result["new_all_Plantae_GBIF_pilot_site_n"],"countries":result["countries_in_source_pilot"],"original_positive":result["original_1055port_WWp_countryFE_positive_model_verified_before_selecting_new_allflora"]},indent=2))
if __name__=="__main__":
    main()
