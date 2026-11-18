#!/usr/bin/env python3
"""Fixed source-check: eight WWTPs, one river. Compare observed normalized
up/down spider response against aquatic insect taxa on same facility.
This is a ratio of source relative responses, NOT a true trophic transfer
factor, original ingestion path, causal drug dose, or independent watershed.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,statistics
from pathlib import Path
SCHEMA="GESUI_20261118_GRAND8_PUBLISHED_MP_SPIDER_AQUATIC_RELATIVE_SOURCE_RESPONSE_V1"
PIN_SOURCE_SHA="1b69eb88fe38f601b11ff466d6a9b3177546a9b489a74c61eaf2a64d770dec66"
PIN_PARENT_SCHEMA="GESUI_20261117_GRAND8_GENUINE_WWTP_AQUATIC_ADULT_RIPARIAN_SPIDER_MP_WITHIN_FACILITY_V1"
TAXA=("Hydropsychidae_larvae","Hydropsychidae_adults","Heptageniidae_larvae","Chironomidae_larvae")
SITES=("Ayr","Caledonia","Guelph","Hespeler","Kitchener","Paris","Preston","Waterloo")
def run(a):
    raw=a.original.read_bytes()
    src=json.loads(raw)
    if src["schema"]!=PIN_PARENT_SCHEMA or src["source_original_SHA256"]!=PIN_SOURCE_SHA or src["original_source_animal_specimens_n"]!=906:
        raise RuntimeError("Original verified source tissue inventory changed")
    if src["n_independent_original_hydrologic_basins"]!=1:
        raise RuntimeError("Number of independent river basins is not one")
    m=src["numerical_original_per_material_stage_contrasts"]["all_original_taxonomic_stages"]
    sp=m["Tetragnathidae_spider"]["original_8WWTP_facility_site_summaries"]
    if set(sp)!=set(SITES):raise RuntimeError("Original eight source WWTP IDs drift")
    rows=[]
    for site in SITES:
        ref=sp[site]
        if ref["n_source_matched_visits"]<2:raise RuntimeError("Original spider source gate failed")
        a_sp=ref["downstream_to_upstream_ratio_if_both_positive"]
        if a_sp is None or a_sp<=0:raise RuntimeError("Original spider source ratio not identifiable")
        for tax in TAXA:
            orig=m[tax]["original_8WWTP_facility_site_summaries"][site]
            if orig["n_source_matched_visits"]<2:raise RuntimeError("Original insect stage source gate failed")
            a_in=orig["downstream_to_upstream_ratio_if_both_positive"]
            if a_in is None or a_in<=0:raise RuntimeError("Upstream zero or invalid concentration: fail closed")
            rows.append({"WWTP":site,"aquatic_taxon_stage":tax,
              "n_original_insect_matched_visits":orig["n_source_matched_visits"],
              "n_original_spider_matched_visits":ref["n_source_matched_visits"],
              "spider_concentration_down_up_ratio":a_sp,
              "insect_concentration_down_up_ratio":a_in,
              "ratio_of_down_up_response_ratios_spider_over_insect":a_sp/a_in,
              "spider_relative_down_up_response_exceeds_insect":a_sp>a_in})
    summary=[]
    for tax in TAXA:
        r=[x for x in rows if x["aquatic_taxon_stage"]==tax]
        if len(r)!=8:raise RuntimeError("Original 8 source stages not fully available")
        vals=[x["ratio_of_down_up_response_ratios_spider_over_insect"] for x in r]
        summary.append({"original_taxon_stage":tax,
          "n_wwtp_source_facilities_with_both_stage_ratios":8,
          "n_spider_response_greater_than_insect_response":sum(v>1 for v in vals),
          "median_relative_spider_over_insect_down_up_ratio_of_ratios":statistics.median(vals),
          "mean_log_relative_spider_over_insect_down_up_ratio_of_ratios":statistics.mean(math.log(v) for v in vals),
          "min_relative_ratio":min(vals),"max_relative_ratio":max(vals),
          "not_a_trophic_transfer_factor":True})
    result={"schema":SCHEMA,"status":"SOURCE_VERIFIED_EIGHT_FACILITY_RELATIVE_MP_TAXON_RESPONSE_DESCRIPTIVE_ONLY",
       "source_published_original_906_assay_SHA256":PIN_SOURCE_SHA,
       "source_parent_actions_artifact_byte_SHA256":hashlib.sha256(raw).hexdigest(),
       "original_independent_wwtp_facilities":8,
       "original_independent_hydrologic_basins":1,
       "animal_unit_all_concentrations_are_particle_per_mg_wet":True,
       "ratio_of_ratios_not_trophic_transfer_factor_or_energy_flux":True,
       "no_plant_nitrogen_phosphorus_nutrient_flux_computed":True,
       "not_independent_publication_experiment":True,
       "exploratory_not_multiple_testing_corrected_or_confirmatory":True,
       "fixed_taxon_stage_summary":summary,
       "original_facility_x_stage_source_values":rows}
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"GRAND8_SPIDER_VS_AQUATIC_TAXA_RELATIVE_DOWNUP_STAGE_RESPONSE.json").write_text(
       json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    with (a.out/"GRAND8_SPIDER_VS_AQUATIC_TAXA_RELATIVE_RESPONSE_PER_WWTP.csv").open(
        "w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print("GRAND8_RELATIVE_STAGE_RESPONSE="+json.dumps(
       {"status":result["status"],"taxa":summary},ensure_ascii=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--original",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    run(p.parse_args())
