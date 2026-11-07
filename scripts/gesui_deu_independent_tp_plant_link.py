#!/usr/bin/env python3
"""New 18-site Plantae GBIF panel and TP–plant observational join.

The six outfalls and 12 TP station IDs were fixed from ORIGINAL hydro/monitoring
metadata before any new plant GBIF response or measured TP concentration.
No source reranking or outcome-optimized thresholds are permitted.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import csv
import time
import math
import numpy as np
from scipy.stats import spearmanr
import gesui_directed_river_plantae as source_plants

SOURCE_RUN=42371319376
SOURCE_FILE="STRICT_STATION_INDEPENDENCE_COHORTS.json"
YEAR="2011,2022"
RADIUS=1

def source_locked_sites(path):
    src=json.loads(path.read_text())
    assert src["selection_made_without_plant_responses_or_nutrient_concentrations"] is True
    a=src["independent_source_only_cohorts"]["TP"]
    assert len(a)==6
    assert len({e["wwtp_id"] for e in a})==6
    assert len({e["hydrobasin_l12"] for e in a})==6
    assert sum(e["independent_source_fraction_stratum"]=="high" for e in a)==3
    ids=[v["station_id"] for e in a for v in e["measured_source_stations"]["TP"].values()]
    assert len(ids)==len(set(ids))==12
    modified=[{**e,"exposure_stratum":e["independent_source_fraction_stratum"]} for e in a]
    source_plants.SOURCE_RUN_ID=SOURCE_RUN
    source_plants.YEAR=YEAR
    source_plants.RADIUS_KM=RADIUS
    sites=source_plants.site_registry(modified)
    assert len(sites)==18 and len({s["site_id"] for s in sites})==18
    return modified,sites

def joined_stats(plants,tp):
    a={int(x["wwtp_id"]):x for x in plants["wwtp_outcomes"]}
    b={int(x["outfall_id"]):x for x in tp["all_results"]}
    assert set(a)==set(b) and len(a)==6
    records=[]
    for wid in sorted(a):
        x=a[wid]
        y=b[wid]
        primary=bool(x["strict_primary_eligible"] and y["primary_source_months_gate_met"])
        sens=bool(x["depth10_sensitivity"]["eligible_by_records_and_species"]
                  and y["primary_source_months_gate_met"])
        records.append({
            "wwtp_id":wid,"source_country":x["country"],"source_fraction":float(x["source_fraction"]),
            "source_fraction_stratum":x["source_fraction_stratum"],
            "TP_months":y["calendar_months_with_measured_up_down"],
            "source_TP_primary_gate":y["primary_source_months_gate_met"],
            "observed_median_log_TP_down_up":y["median_monthly_log_down_up"],
            "observed_median_TP_down_up_ratio":y["median_monthly_down_up_ratio"],
            "plant_richness_strict_record_support":x["strict_primary_eligible"],
            "plant_richness_depth10_support":x["depth10_sensitivity"]["eligible_by_records_and_species"],
            "primary_joint_support":primary,
            "depth10_joint_support":sens,
            "plant_Y_log_rarefied30":x["depth30"]["extra_change_downstream_minus_prior_upstream"]
                if primary else None,
            "plant_Y_log_rarefied10":x["depth10_sensitivity"]["extra_change_downstream_minus_prior_upstream"]
                if sens else None,
            "TP_gradient_is_sewage_P_load":False,
        })
    def check(xs,field):
        use=[v for v in xs if v[field] is not None and v["observed_median_log_TP_down_up"] is not None]
        z={"eligible_outfalls":len(use),"high":sum(v["source_fraction_stratum"]=="high" for v in use),
           "low":sum(v["source_fraction_stratum"]=="low" for v in use),
           "country_count":len({v["source_country"] for v in use})}
        if len(use)>=4 and len({v[field] for v in use})>=3:
            rho,p=spearmanr([v[field] for v in use],
                            [v["observed_median_log_TP_down_up"] for v in use])
            z.update({"rho_observed_TP_gradient_vs_GBIF_richness_curvature":float(rho),
                      "two_sided_p_descriptive_only":float(p),
                      "no_causal_or_population_significance":True})
        else:z["status"]="INSUFFICIENT_JOINED_RESPONSE_SUPPORT"
        return z
    def TP_vs_flow(xs):
        arr=[v for v in xs if v["observed_median_log_TP_down_up"] is not None]
        if len(arr)<4:return {"n":len(arr),"status":"INSUFFICIENT"}
        rho,p=spearmanr([v["source_fraction"] for v in arr],
                        [v["observed_median_log_TP_down_up"] for v in arr])
        return {"n":len(arr),"rho":float(rho),"p_descriptive":float(p)}
    return {
        "status":"NEW_SOURCE_BLIND_DEU_TP_PLUS_ALL_PLANTAE_SPATIAL_PILOT",
        "source_selected_outfalls":6,
        "plant_sites":18,
        "new_outcomes_selected_based_on_data":False,
        "all_six_country":"DEU",
        "joined":records,
        "strict_richness30_TP_gradient":check(records,"plant_Y_log_rarefied30"),
        "depth10_sensitivity_TP_gradient":check(records,"plant_Y_log_rarefied10"),
        "observed_TP_gradient_vs_model_wwflow":TP_vs_flow(records),
        "interpretation_limit":[
           "Six independent source pairs from one country, unequal station-to-outfall and site segment lengths",
           "Measured TP at river monitoring stations not WWTP effluent concentration or P loads",
           "Plantae GBIF observation index cannot measure community abundance or colonization",
           "No urban agricultural/catchment landcover or other wastewater source controls yet",
           "Even a positive rho is descriptive and not evidence of P causing plant change"
        ],
    }

def run(args):
    selected,sites=source_locked_sites(args.source)
    tp=json.loads(args.tp_results.read_text())
    assert tp["source_selected_outfall_count"]==6 and tp["source_doi"]=="10.25400/lincolnuninz.30546371.v5"
    args.out.mkdir(parents=True,exist_ok=True)
    frozen={
        "original_source_run":SOURCE_RUN,
        "source_cohort_sha256":hashlib.sha256(args.source.read_bytes()).hexdigest(),
        "TP_results_sha256":hashlib.sha256(args.tp_results.read_bytes()).hexdigest(),
        "source_selected_without_biology_or_TP_values":True,
        "original_source_plant_selections":sites,
        "GBIF_kingdom":6,
        "GBIF_year":YEAR,
        "GBIF_radius_km":RADIUS,
        "primary_rarefaction_records":30,
        "primary_min_species":6,
        "primary_record_max_min_ratio":10,
        "sensitivity_rarefaction_records":10,
        "sensitivity_min_species":4,
    }
    (args.out/"FROZEN_DESIGN.json").write_text(json.dumps(frozen,indent=2,sort_keys=True)+"\n")
    print("GESUI_DEU_TP_NEW_PLANT_SOURCE_LOCK",
        json.dumps({"wwtp":len(selected),"sites":len(sites),"period":YEAR},sort_keys=True),flush=True)
    observed=[]
    for i,site in enumerate(sites,1):
        item=source_plants.query_all_plantae(site)
        observed.append(item)
        print("GESUI_DEU_TP_PLANTAE_SITE",json.dumps({
            "index":i,"outfall_id":item["wwtp_id"],"role":item["river_role"],
            "indexed_species":item["distinct_species_keys"],
            "records_specieskeyed":item["indexed_species_record_count"]
        },sort_keys=True),flush=True)
        time.sleep(.4)
    assert len(observed)==18
    raw={
        "schema":"gesui_new_independent_deu_tponly_gbif_plantae_facets_v01",
        "source_outfall_selection_run":SOURCE_RUN,
        "retrieved_utc":datetime.now(timezone.utc).isoformat(),
        "gbif_window":YEAR,
        "sites":observed,
    }
    (args.out/"PLANT_original_GBIF_species_facets.json").write_text(
        json.dumps(raw,ensure_ascii=False,separators=(",",":"))+"\n")
    with (args.out/"PLANT_site_species_records.csv").open("w",newline="") as f:
        fieldnames=["site_id","wwtp_id","river_role","river_reach_id",
            "gbif_speciesKey","speciesKey_observation_count","gbif_window"]
        writer=csv.DictWriter(f,fieldnames=fieldnames)
        writer.writeheader()
        for site in observed:
            for key,value in sorted(site["species_key_record_counts"].items(),key=lambda x:int(x[0])):
                writer.writerow({"site_id":site["site_id"],"wwtp_id":site["wwtp_id"],
                    "river_role":site["river_role"],"river_reach_id":site["river_reach_id"],
                    "gbif_speciesKey":key,"speciesKey_observation_count":value,"gbif_window":YEAR})
    plant=source_plants.summarise(selected,observed)
    plant.update({
        "status":"FIRST_ALL_PLANTAE_OBSERVATIONS_AT_SIX_NEW_SOURCE_BLIND_TP_MONITORED_OUTFALLS",
        "source_selected_outfalls":len(selected),"biological_sites":len(observed),
        "GBIF_year":YEAR,"source_cohort":SOURCE_RUN,
    })
    (args.out/"PLANT_RICHNESS_AND_COMPOSITION.json").write_text(
        json.dumps(plant,indent=2,sort_keys=True)+"\n")
    joint=joined_stats(plant,tp)
    (args.out/"JOINT_TP_PLANT_DESCRIPTIVE.json").write_text(
        json.dumps(joint,indent=2,sort_keys=True)+"\n")
    print("GESUI_DEU_TP_PLANT_OBSERVATIONAL_RESULT",
          json.dumps({k:v for k,v in plant.items() if k not in ("wwtp_outcomes","limits")},
                     sort_keys=True),flush=True)
    print("GESUI_DEU_TP_PLANTS_SOURCE_LINKAGE",
          json.dumps({k:v for k,v in joint.items() if k not in ("joined","interpretation_limit")},
                     sort_keys=True),flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--tp-results",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    run(ap.parse_args())
if __name__=="__main__":
    main()
