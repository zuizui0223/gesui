#!/usr/bin/env python3
"""SOURCE-ONLY ceiling audit for true European H7-period plant trait analysis.

Read genuine 234 H7 x 2 period x 6 group GBIF speciesKey source,
restrict to the *unaltered* original 454 H7-period panel, and report plant
record/identified species opportunity. Does NOT compute EIV, plant nitrogen
index, environmental regression, wastewater dose, or biological mechanisms.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

EXPECTED_RAW_JSON_SHA256 = "251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
RAW_JSON_NAME = "GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json"

def sha256(file: Path) -> str:
    return hashlib.sha256(file.read_bytes()).hexdigest()

def source_gate(panel_path: Path, archived_dir: Path, out_dir: Path):
    matches=list(archived_dir.rglob(RAW_JSON_NAME))
    if len(matches)!=1:
        raise RuntimeError("Expected exactly one original publisher-derived H7 source JSON: "+repr([str(p) for p in matches]))
    raw_path=matches[0]
    if sha256(raw_path)!=EXPECTED_RAW_JSON_SHA256:
        raise RuntimeError("Original H7 source JSON SHA256 mismatch")
    original=json.loads(raw_path.read_text())
    if len(original)!=2808 or any(x.get("source_status")!="ok" for x in original):
        raise RuntimeError("Genuine 2808 facet source/status condition failed")
    plants=[x for x in original if x.get("group")=="Plantae"]
    if len(plants)!=468:
        raise RuntimeError("Expected 234 x 2 original complete Plantae source facets")
    src={(str(x["HYBAS_ID_L07"]),str(x["period"])):x for x in plants}
    if len(src)!=468:raise RuntimeError("Duplicate plant H7 source keys")
    with panel_path.open(newline="",encoding="utf-8-sig") as f:
        panel=list(csv.DictReader(f))
    if len(panel)!=454:raise RuntimeError("Original 454 H7-period cohort was changed")
    for name in ("h","period","h6","country_holdout"):
        if name not in panel[0]:raise RuntimeError("Missing original panel field "+name)
    if len({(str(x["h"]),str(x["period"])) for x in panel})!=454:
        raise RuntimeError("Duplicate original original H7-period source keys")
    rows=[]
    for x in panel:
        key=(str(x["h"]),str(x["period"]))
        if key not in src:raise RuntimeError("Plant source missing for H7-period "+repr(key))
        e=src[key];cts=e.get("source_speciesKey_counts")
        if not isinstance(cts,dict):raise RuntimeError("Original key counts unavailable")
        if any(int(v)<0 or str(v)!=str(int(v)) for v in cts.values()):
            raise RuntimeError("Non-integer or negative source species counts")
        raw_n=int(e["source_query_occurrence_count"])
        counted=sum(int(v) for v in cts.values())
        if raw_n<counted:
            raise RuntimeError("Original identified records exceed original query records")
        facetcov=e.get("source_facet_coverage")
        if facetcov is not None and not (0<=float(facetcov)<=1):
            raise RuntimeError("Facet coverage outside [0,1]")
        rows.append({
            "h":key[0],"period":key[1],"h6":str(x["h6"]),
            "country_holdout":str(x["country_holdout"]),
            "raw_Plantae_query_records":raw_n,
            "indexed_speciesKey_records":counted,
            "indexed_unique_plant_speciesKeys":sum(int(v)>0 for v in cts.values()),
            "speciesKey_record_coverage":counted/raw_n if raw_n else None,
            "original_source_facet_coverage":facetcov,
            "pretrait_only_record50_species10":raw_n>=50 and len(cts)>=10,
            "EIV_NUTRI_score_join": "NOT_PERFORMED",
        })
    pre=[x for x in rows if x["pretrait_only_record50_species10"]]
    summary={
        "schema":"GESUI_EU454_PLANT_SOURCE_ONLY_OPPORTUNITY_V1",
        "analysis_date":"2026-11-15",
        "original_H7_source_SHA256":EXPECTED_RAW_JSON_SHA256,
        "source_json_n_facets":len(original),
        "all_original_plantae_facets":len(plants),
        "matched_original_H7_periods":len(rows),
        "original_H7":len({x["h"] for x in rows}),
        "original_H6":len({x["h6"] for x in rows}),
        "original_country_labels":len({x["country_holdout"] for x in rows}),
        "pretrait_only_record50_species10_H7_periods":len(pre),
        "pretrait_only_record50_species10_H6":len({x["h6"] for x in pre}),
        "pretrait_only_record50_species10_countries":len({x["country_holdout"] for x in pre}),
        "upper_bound_passes_before_EIV_link":len(pre)>=100 and len({x["h6"] for x in pre})>=40 and len({x["country_holdout"] for x in pre})>=5,
        "min_scored_EIV_taxa_requirement":10,
        "min_scored_EIV_species_coverage_requirement":0.40,
        "EIV_NUTRI_taxon_join_and_actual_trait_eligibility":"NOT_PERFORMED",
        "site_vegetation_census":"NOT_AVAILABLE",
        "soil_or_river_actual_N_concentrations":"NOT_AVAILABLE",
        "wastewater_causal_mediation":"NOT_ESTIMATED",
        "no_posthoc_cohort_replacement":True,
        "original_butterfly_HNS_used_to_choose_plants":False,
        "summary_raw_records":dict(collections.Counter("50_or_more" if x["raw_Plantae_query_records"]>=50 else "below50" for x in rows)),
        "source_traits_need_new_taxonomic_join":True,
        "interpretation":"Plant faceted GBIF opportunity and taxon-ID completeness only; cannot assert NUTRI trait support, selection, wastewater influence, temporal mediation or ecosystem services.",
    }
    out_dir.mkdir(parents=True,exist_ok=True)
    with (out_dir/"EU454_PLANT_SOURCE_ONLY_OPPORTUNITY.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (out_dir/"EU454_PLANT_SOURCE_ONLY_GATE.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_EU454_PLANT_SOURCE_GATE="+json.dumps(summary,ensure_ascii=False,sort_keys=True),flush=True)
    return summary

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--panel",type=Path,default=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"))
    ap.add_argument("--original",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    source_gate(args.panel,args.original,args.output)

