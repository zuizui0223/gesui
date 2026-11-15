#!/usr/bin/env python3
"""SOURCE ONLY: original GBIF speciesKey facets for all 30 new H7 that pass
the independent prior 1219-index count source gate in 10 distinct H6.

Uses EXACT original 234 H7 source facet task logic; original family GBIF
backbone and original 1219 native WKT, plus frozen 40-candidate source
ranking. No EIV NUTRI, HNS response or wastewater models.
"""
from __future__ import annotations
import argparse
import concurrent.futures
import hashlib
import json
import csv
import sys
from pathlib import Path
from collections import Counter

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from gesui_H7_original_fivefamily_plants_speciesKey_retention import task,FAMILIES,PERIODS

WKT_SHA="47f597e0a798cd6155ee9fb69ce8916ed5629951d046d114feb9f667963310e4"
TAXON_SHA="afc524ebfa353e4d1c6ea51d208ba529a85fedb57c4d0f28d83545010bb76189"
SOURCE_SCHEMA="GESUI_ORIGINAL_1219_GBIF_COUNT_SOURCE_OPPORTUNITY_FOR_12_NEW_SINGLETON_H6_V1"
GROUPS=("Plantae",)+FAMILIES
SHARDS=3
EXPECTED_H7=30
EXPECTED_NEW_H6=10

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def one(folder,name):
    f=list(folder.rglob(name))
    if len(f)!=1:raise RuntimeError("Exactly one source file required: "+name+":"+repr([str(x) for x in f]))
    return f[0]
def source_inputs(args):
    scope=json.loads(one(args.prior,"ORIGINAL1219_GBIF_NEW_H7_RECORD_OPPORTUNITY_GATE.json").read_text())
    if (scope.get("schema")!=SOURCE_SCHEMA
        or not scope.get("original_fixed_40_H6_target_source_index_count_gate_pass")
        or scope.get("source_target_12_H6_with_any_new_H7_both_era_index_tier_pass")!=EXPECTED_NEW_H6
        or scope.get("source_original_40_new_H7_both_era_tier_pass")!=EXPECTED_H7
        or scope.get("original_1219_WKT_source_sha")!=WKT_SHA):
        raise RuntimeError("Original frozen 30-H7/10-H6 source-only count opportunity gate invalid")
    p=one(args.prior,"ORIGINAL1219_EEA_NEW_H7_PLANT_LYCA_ANIMALIA_RECORD_OPPORTUNITIES.csv")
    with p.open(newline="",encoding="utf-8-sig") as handle:
        candidates=list(csv.DictReader(handle))
    qualified=[q for q in candidates if q["both_periods_original_three_taxa_index_tiers_pass"]=="True"]
    if len(candidates)!=40 or len(qualified)!=EXPECTED_H7 or len({q["native_H7"] for q in qualified})!=EXPECTED_H7:
        raise RuntimeError("Original candidates no longer 40 source H7, 30 count-eligible")
    if len({q["official_native_H6"] for q in qualified})!=EXPECTED_NEW_H6:
        raise RuntimeError("Unexpected original 10 independent singleton-H6 geographic targets")
    path=one(args.wkt,"HYBAS7_1219_ORIGINAL_NATIVE_POLYGONS_WKT_BEFORE_NEW_GBIF_COUNTS.json")
    if digest(path)!=WKT_SHA:raise RuntimeError("Previously frozen native 1219 WKT manifest SHA mismatch")
    j=json.loads(path.read_text())
    if len(j["basins"])!=1219 or not j.get("geography_source_gate_pass"):
        raise RuntimeError("Genuine 1219 native WKT source dataset changed")
    by={str(x["HYBAS_ID_L07"]):x for x in j["basins"]}
    if len(by)!=1219:raise RuntimeError("Duplicate original 1219 Level7 source IDs")
    taxa_path=one(args.taxonomy,"HYBAS7_234_REPEATED_EU_EQR_ORIGINAL_WKT_TAXON_FAMILY_PRE_NEW_OUTCOME.json")
    if digest(taxa_path)!=TAXON_SHA:raise RuntimeError("Original EEA five butterfly family taxon manifest SHA mismatch")
    fixed=json.loads(taxa_path.read_text())
    families=fixed.get("verified_butterfly_family_taxa",{})
    if (set(families)!=set(FAMILIES) or families["Lycaenidae"].get("taxonKey")!=5473
        or not all(families[x].get("rank")=="FAMILY" and families[x].get("scientificName")==x for x in FAMILIES)):
        raise RuntimeError("GBIF original family keys taxonomy drift")
    selected={q["native_H7"] for q in qualified}
    if not selected.issubset(by):
        raise RuntimeError("Original source 1219 WKT does not include selected 30 H7")
    if any(by[k]["status"]!="all_original_polygons_WKT_queryable" for k in selected):
        raise RuntimeError("Some geographic candidate WKT cannot be queried exactly")
    return by,families,selected,qualified

def main(args):
    by,families,selected,qualified=source_inputs(args)
    if args.mode=="pilot":
        # Source-only, deterministic two IDs before querying any new speciesKey.
        ids=sorted(selected,key=int)[:2]
    else:
        ids=[h for h in selected if int(hashlib.sha256(h.encode()).hexdigest()[:8],16)%SHARDS==args.shard]
        ids.sort(key=int)
    if not ids:raise RuntimeError("No valid source H7 in this deterministic shard")
    tasks=[(by[h],period,group,families) for h in ids
           for period in PERIODS for group in GROUPS]
    result=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futs=[pool.submit(task,*item) for item in tasks]
        for f in concurrent.futures.as_completed(futs):
            result.append(f.result())
    result.sort(key=lambda r:(str(r["HYBAS_ID_L07"]),r["period"],r["group"]))
    if len(result)!=len(tasks):raise RuntimeError("Original 30-H7 source task input/return count changed")
    keyset={(str(r["HYBAS_ID_L07"]),r["period"],r["group"]) for r in result}
    if len(keyset)!=len(tasks):raise RuntimeError("Duplicated original source GBIF H7-period-facet result")
    for r in result:
        if r["source_status"]=="ok":
            if not isinstance(r["source_query_occurrence_count"],int) or not isinstance(r["source_speciesKey_counts"],dict):
                raise RuntimeError("Original indexed speciesKey source success malformed")
        elif r["source_status"]=="missing":
            if r.get("source_query_occurrence_count") is not None or r.get("source_speciesKey_counts") is not None:
                raise RuntimeError("Source retrieval failure was silently zero-imputed")
        else:raise RuntimeError("Unrecognized genuine GBIF source query status")
    out=args.output
    out.mkdir(parents=True,exist_ok=True)
    stem=f"EXTRA_EEA_H7_REAL_SPECIESKEY_{args.mode}_SHARD_{args.shard:02d}"
    filepath=out/(stem+".json")
    filepath.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,separators=(",",":"))+"\n")
    missing=[r for r in result if r["source_status"]!="ok"]
    report={
      "schema":"GESUI_SOURCE_ONLY_30_NEW_H7_2ERA_SIX_ORIGINAL_FACET_ACQUISITION_V1",
      "source_only_mode":args.mode,"source_shard":args.shard,
      "original_source_1219WKT_SHA256":WKT_SHA,
      "original_234_fixed_family_taxa_source_SHA256":TAXON_SHA,
      "original_count_supported_30H7_10H6":True,
      "source_H7_in_shard":len(ids),
      "source_H7_native_ids":ids,
      "source_6_group_2period_facets_requested":len(tasks),
      "source_facet_ok":len(tasks)-len(missing),
      "source_facet_missing":len(missing),
      "source_missing_which":[[str(r["HYBAS_ID_L07"]),r["period"],r["group"]] for r in missing],
      "original_request_method_reused_from_original_234_H7_source":True,
      "original_output_source_JSON_SHA256":digest(filepath),
      "site_EIV_NUTRI_and_butterfly_HNS_quality_not_yet_calculated":True,
      "NUTRI_HNS_EEA_W_numeric_or_causal_effect_not_computed":True}
    (out/(stem+"_QA.json")).write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_NEW_H7_SOURCE_SPECIESKEY_FACETS="+json.dumps(report,ensure_ascii=False,sort_keys=True),flush=True)
    if missing:
        print("SOURCE_ONLY_MISSING_OR_API_RATE_LIMIT_NOT_ECOLOGICAL_ZERO",flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--prior",type=Path,required=True)
    p.add_argument("--wkt",type=Path,required=True)
    p.add_argument("--taxonomy",type=Path,required=True)
    p.add_argument("--mode",choices=["pilot","full"],required=True)
    p.add_argument("--shard",type=int,default=0)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    if a.shard<0 or a.shard>=SHARDS:raise RuntimeError("Source shard number must be 0,1,2")
    main(a)
