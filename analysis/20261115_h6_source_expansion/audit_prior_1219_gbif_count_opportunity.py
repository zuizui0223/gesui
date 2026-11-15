#!/usr/bin/env python3
"""Reuse PREVIOUS genuine 1219 H7×two era×3 taxon GBIF record source to audit
the 40 geography-only new H7 candidates in 12 source-original singleton H6.

No new GBIF calls, no species trait scores, no EQR, no W dose or biological
coefficients. Results are bounded RECORDING OPPORTUNITY, never true absence.
"""
from __future__ import annotations
import argparse
import csv
import json
import hashlib
from collections import Counter,defaultdict
from pathlib import Path

GEO_SCHEMA="GESUI_EEA_SINGLEERA_H7_OFFICIAL_H6_GEO_SOURCE_OPPORTUNITY_V1"
SHARD_SCHEMA_PREFIX="HYBAS7_1219_GBIF_PLANTAE_LYCAENIDAE_ANIMALIA_TWO_PERIOD_INDEX_SOURCE_"
WKT_SHA="47f597e0a798cd6155ee9fb69ce8916ed5629951d046d114feb9f667963310e4"
SOURCE_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
H7_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
ALL_TAXA={"Plantae","Lycaenidae","Animalia"}
ALL_PERIODS={"2017-2019","2020-2022"}
SOURCE_TIERS={"Plantae":50,"Lycaenidae":10,"Animalia":50}
ORIGINAL_SOURCE_H7=1219
REQUIRED_NEW_H6=10

def unique(folder,name):
    paths=list(folder.rglob(name))
    if len(paths)!=1:raise RuntimeError("Exactly one original artifact file expected: "+name+" "+repr([str(x) for x in paths]))
    return paths[0]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def jsonfile(path):return json.loads(path.read_text())
def main(args):
    geompath=unique(args.geometry,"EEA_SINGLEERA_H7_OFFICIAL_H6_SOURCE_GEO_GATE.json")
    geom=jsonfile(geompath)
    if geom.get("schema")!=GEO_SCHEMA or not geom.get("source_target_geo_gate_min_10_distinct_singleton_H6"):
        raise RuntimeError("Official source-only 12-H6 geometry opportunity gate invalid")
    if geom.get("original_6330_source_SHA256")!=SOURCE_SHA or geom.get("original_454_panel_SHA256")!=H7_SHA:
        raise RuntimeError("Original source station/H7 panel SHA changed")
    if (geom.get("target_existing_singleton_H6_with_geographic_additional_EEA_H7")!=12 or
        geom.get("unique_extra_EEA_H7_candidates_in_target_singleton_H6")!=40):
        raise RuntimeError("Previously source-frozen 12 H6/40 H7 geography opportunity drift")
    registry=unique(args.geometry,"EEA_SINGLEERA_NEW_H7_IN_ORIGINAL_SINGLETON_H6_SOURCE_REGISTRY.csv")
    with registry.open(encoding="utf-8-sig",newline="") as handle:
        candidates=list(csv.DictReader(handle))
    if len(candidates)!=40 or len({z["native_H7"] for z in candidates})!=40:
        raise RuntimeError("Original 40 new candidate H7 source registry missing or duplicate")
    if len({z["official_native_H6"] for z in candidates})!=12:
        raise RuntimeError("Original 12 singleton H6 source registry missing")
    if sum(z["source_selected_1_per_parent"]=="True" for z in candidates)!=12:
        raise RuntimeError("Original deterministic 1-per-H6 ranked geo source selection changed")
    for c in candidates:
        if c["not_new_EIV_or_HNS_eligible"]!="True" or c["EEA_years"] not in ("both","2018_only","2022_only"):
            raise RuntimeError("Source geo record was misidentified as ecological trait-eligible")
    originals={}
    subsummary=[]
    for i in range(4):
        d=args.shards/str(i)
        path=unique(d,"HYBAS7_GBIF_ORIGINAL_SOURCE_COUNTS_URL_SHA.json")
        qa=jsonfile(unique(d,"HYBAS7_GBIF_ORIGINAL_SOURCE_SUPPORT_QA.json"))
        if (qa.get("source_wkt_manifest_SHA256")!=WKT_SHA or
            qa.get("fixed_original_H7_eligible_basins")!=ORIGINAL_SOURCE_H7 or
            qa.get("mode")!="full" or int(qa.get("shard",-1))!=i or
            qa.get("fixed_source_requests")!=qa.get("fixed_source_shard_basin_count")*2*3):
            raise RuntimeError("Prior archived 1219 original GBIF 3-taxon shard QA invalid")
        blob=jsonfile(path)
        if len(blob)!=qa["fixed_source_requests"]:
            raise RuntimeError("Original 1219 species group count shard incomplete")
        subsummary.append({"shard":i,"source_original_rows":len(blob),
            "original_audit_SHA256":sha(d/"HYBAS7_GBIF_ORIGINAL_SOURCE_SUPPORT_QA.json"),
            "original_numeric_source_SHA256":sha(path),
            "successful_source_requests":qa.get("GBIF_source_successful_requests"),
            "missing_source_requests":qa.get("GBIF_source_failed_requests")})
        for x in blob:
            k=(str(x["HYBAS_ID_L07"]),str(x["period"]),str(x["taxon"]))
            if k in originals or k[1] not in ALL_PERIODS or k[2] not in ALL_TAXA:
                raise RuntimeError("Duplicate/invalid original true H7-period-taxon source identity")
            count=x.get("GBIF_index_record_count")
            if x.get("status")=="ok":
                if not isinstance(count,int) or count<0:
                    raise RuntimeError("Source-OK numerical GBIF index not a nonnegative integer")
            else:
                if count is not None:raise RuntimeError("Failed source was filled with a biological count")
            originals[k]=x
    expected=ORIGINAL_SOURCE_H7*len(ALL_PERIODS)*len(ALL_TAXA)
    if len(originals)!=expected:
        raise RuntimeError(f"Original 1219 * 2 * 3 GBIF source table mismatch: {len(originals)} vs {expected}")
    unique_source_h7={k[0] for k in originals}
    if len(unique_source_h7)!=ORIGINAL_SOURCE_H7:
        raise RuntimeError("Original source frame lost unique native H7 basins")
    missing_h7=sorted(set(r["native_H7"] for r in candidates)-unique_source_h7,key=int)
    if missing_h7:
        raise RuntimeError("Original preexisting 1219 H7 source lacks frozen EEA singleton source candidates "+str(missing_h7))
    records=[]
    for c in candidates:
        hid=c["native_H7"];d={**c}
        full_ok=True; count_ok=True; missing=[]
        for period in sorted(ALL_PERIODS):
            for taxon in sorted(ALL_TAXA):
                raw=originals[(hid,period,taxon)]
                status=str(raw.get("status"))
                n=raw.get("GBIF_index_record_count")
                d[f"{period}_{taxon}_original_source_status"]=status
                d[f"{period}_{taxon}_original_GBIF_index_records"]=n
                if status!="ok":
                    full_ok=False;count_ok=False
                    missing.append(period+":"+taxon+":SOURCE_MISSING")
                elif n<SOURCE_TIERS[taxon]:
                    count_ok=False
                    missing.append(period+":"+taxon+":INDEX_BELOW_"+str(SOURCE_TIERS[taxon]))
        d["both_periods_all_original_three_taxa_source_retrieved"]=full_ok
        d["both_periods_original_three_taxa_index_tiers_pass"]=count_ok
        d["opportunity_source_fail_or_low_types"]=";".join(missing)
        d["not_original_454_study_EEA_two_period_member"]=True
        d["original_speciesKey_trait_matrix_NOT_available_in_this_count_table"]=True
        records.append(d)
    by_h6=defaultdict(list)
    for d in records:by_h6[d["official_native_H6"]].append(d)
    selected=[]
    byH6=[]
    for parent in sorted(by_h6,key=int):
        arr=sorted(by_h6[parent],key=lambda a:int(a["source_rank_within_original_H6"]))
        # Fall through original *pre-frozen* geography rank ONLY upon source
        # retrieval/record-effort failure; no NUTRI/HNS score or coefficient read.
        valid=[r for r in arr if r["both_periods_original_three_taxa_index_tiers_pass"]]
        if valid:
            chosen=valid[0]
            selected.append(chosen)
        byH6.append({"original_singleton_H6":parent,
             "n_original_geo_candidates":len(arr),
             "n_source_retrieved_both_periods":sum(r["both_periods_all_original_three_taxa_source_retrieved"] for r in arr),
             "n_count_tier_supported_both_periods":len(valid),
             "top_source_rank_candidate_H7":arr[0]["native_H7"],
             "first_original_geographic_rank_with_GBIF_source_tiers":chosen["native_H7"] if valid else "",
             "eligible_GBIF_counts_only_NOT_traits":bool(valid)})
    count_supported_h6=len(selected)
    result={
      "schema":"GESUI_ORIGINAL_1219_GBIF_COUNT_SOURCE_OPPORTUNITY_FOR_12_NEW_SINGLETON_H6_V1",
      "original_1219_WKT_source_sha":WKT_SHA,
      "original_previous_source_registry_H7":len(candidates),
      "original_genuine_1219_H7_source_total_ids":len(unique_source_h7),
      "source_4_archived_shard_provenance":subsummary,
      "original_H7period_x_3_taxa_source_rows":len(originals),
      "original_geography_target_singleton_H6":len(by_h6),
      "frozen_2period_3taxa_count_thresholds":SOURCE_TIERS,
      "original_two_fixed_periods":sorted(ALL_PERIODS),
      "source_original_40_new_H7_all_six_queries_succeed":sum(r["both_periods_all_original_three_taxa_source_retrieved"] for r in records),
      "source_original_40_new_H7_both_era_tier_pass":sum(r["both_periods_original_three_taxa_index_tiers_pass"] for r in records),
      "source_target_12_H6_with_any_new_H7_both_era_index_tier_pass":count_supported_h6,
      "original_prior_30_multiH7_H6":30,
      "count_only_possible_new_geographic_H6_total":30+count_supported_h6,
      "original_fixed_40_H6_target_source_index_count_gate_pass":count_supported_h6>=REQUIRED_NEW_H6,
      "source_ranked_count_eligible_selected_H7_ids":[r["native_H7"] for r in selected],
      "original_source_same_H6_repetition_not_yet_trait_or_covariate_suitable":True,
      "original_prior_study_454_plant_speciesKey_facet_not_automatically_extended":True,
      "no_new_species_traits_NUTRI_HNS_EEA_EQR_W_or_causal_effect_fitted":True,
      "original_GBIF_index_count_not_true_species_abundance_or_absence":True
    }
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"ORIGINAL1219_GBIF_NEW_H7_RECORD_OPPORTUNITY_GATE.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    def dump(name,rows):
        if not rows:return
        with (args.output/name).open("w",newline="",encoding="utf-8") as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]))
            writer.writeheader();writer.writerows(rows)
    dump("ORIGINAL1219_EEA_NEW_H7_PLANT_LYCA_ANIMALIA_RECORD_OPPORTUNITIES.csv",records)
    dump("ORIGINAL1219_SINGLETON_H6_INDEX_OPPORTUNITY_BY_PARENT.csv",byH6)
    print("GESUI_1219_GBIF_NEW_H7_H6_SOURCE="+json.dumps(result,sort_keys=True,ensure_ascii=False),flush=True)
    print("SOURCE_NEXT_STAGE_GATE="+("PASS_COUNTS_ONLY_REQUIRE_NEW_SPECIES_FACET" if result["original_fixed_40_H6_target_source_index_count_gate_pass"]
        else "STOP_STRICT_40_H6_GBIF_OPPORTUNITY"),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--geometry",type=Path,required=True)
    p.add_argument("--shards",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
