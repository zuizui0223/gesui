#!/usr/bin/env python3
"""Merge original1055 frozen GBIF Animalia+Fungi occurrence COUNT comparators.

Source-only strict pre-outcome gate: source 1055/42, original all-Plantae
record-support 752/38 frozen, all 4 source shards and no post-hoc radius or
biological species selection; do not interpret missing as genuine 0.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

D=Path("results/global_all_Plantae_WWp2015_20261110")
S=D/"ORIGINAL_ALL1055_42COUNTRY_WWp2015_BEFORE_NEW808_GBIF_SOURCE.json"
P=D/"ORIGINAL_ALL1055_42COUNTRY_GBIF_PLANTAE_2011_24_50RECORD_RAREFACTION_SOURCE.json"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def join(root,paths):
    source=json.loads((root/S).read_text())
    plant=json.loads((root/P).read_text())
    if source["schema"]!="gesui_original_1055_WPI_ports42country_808NEW_2011_24_ALLPlantae_PREOUTCOME_frame_v01":
        raise ValueError("Original source cohort changed")
    if plant["original_all1055_WWp_WPI_country_source_frame_SHA256"]!=sha(root/S):
        raise ValueError("Prior full all-Plantae source failed recorded/source SHA")
    frozen={x["original_port_id"]:x for x in source["source_all_1055_public_country_port_site_metadata_no_new_GBIF_response"]}
    plantrows={x["original_port_id"]:x for x in plant["original_1055_port_all_Plantae_source_GBIF_facet_summaries"]}
    assert len(frozen)==len(plantrows)==1055 and set(frozen)==set(plantrows)
    original_plantae_support={
        pid for pid,z in plantrows.items()
        if z["source_quality_min50_faceted_records_min10_species_coverage70pct"]
    }
    if len(original_plantae_support)!=752:raise ValueError("Original plant quality subpopulation changed")
    batches={}
    grouped={}
    for path in paths:
        p=Path(path)
        row=json.loads(p.read_text())
        shard=row["shard_id"]
        if shard in batches or shard not in range(4):raise ValueError("Original 4-shard order duplicated")
        if row["original_WPI2019_1055_country2015_WWp_source_SHA256"]!=sha(root/S):
            raise ValueError("Animalia/Fungi source WPI port SHA differs")
        if row["original_2011_24_allPlantae_1055_response_source_SHA256"]!=sha(root/P):
            raise ValueError("Animalia/Fungi source differs from prior Flora file SHA")
        if row["original_query_taxonomy_ids"]!={"Animalia":1,"Fungi":5}:
            raise ValueError("Precommitted two kingdom taxonKeys changed")
        expected=(264 if shard<3 else 263)
        records=row["per_original_source_port_taxon_indexed_record_counts"]
        if len(records)!=expected:raise ValueError("Unexpected number source ports in cross kingdom shard")
        for x in records:
            pid=x["original_port_id"]
            if pid in grouped or pid not in frozen:
                raise ValueError("Duplicated/non-source recording control WPI port")
            if x["original_source_country"]!=frozen[pid]["original_source_country"]:
                raise ValueError("Source country changed in taxon recording comparator")
            if x["source_original_all_Plantae_quality_752_flag"]!=(pid in original_plantae_support):
                raise ValueError("Original plant quality was reselected after Animalia/Fungi observations")
            for taxon,key in (("Animalia",1),("Fungi",5)):
                z=x[taxon]
                if z["GBIF_status"]=="COMPLETE":
                    n=z["original_indexed_record_count"]
                    if type(n) is not int or n<0 or not z["response_sha256"] or len(z["response_sha256"])!=64:
                        raise ValueError("Original successful GBIF count is invalid")
                elif z["GBIF_status"]=="API_ERROR_NOT_ZERO":
                    if z["original_indexed_record_count"] is not None or z["response_sha256"] is not None:
                        raise ValueError("API failed count was falsely set to biological zero")
                else:raise ValueError("Unknown source API status")
            grouped[pid]=x
        batches[shard]={
           "original_shard_ID":shard,"source_receipt_SHA256":sha(p),
           "raw_source_count_cache_SHA256":row["original_GBIF_crosskingdom_cache_SHA256"],
           "source_site_rows":len(records)
        }
    if set(batches)!=set(range(4)) or len(grouped)!=1055 or set(grouped)!=set(frozen):
        raise ValueError("Full original1055 globally locked WPI sites not independently sampled")
    allgood=[
        x for x in grouped.values()
        if all(x[t]["GBIF_status"]=="COMPLETE" for t in ("Animalia","Fungi"))
    ]
    bothids={x["original_port_id"] for x in allgood}
    plantgood=original_plantae_support & bothids
    countries=defaultdict(list)
    for pid in plantgood:
        countries[frozen[pid]["original_source_country"]].append(pid)
    n_source_varied=sum(
        len(ids)>=2 and len({frozen[p]["source_WWp_log1p"] for p in ids})>=2
        for country,ids in countries.items()
    )
    admissible=bool(
        len(allgood)>=1002 and len(plantgood)>=700 and
        len(countries)>=30 and n_source_varied>=25
    )
    per_status={taxon:dict(Counter(x[taxon]["GBIF_status"] for x in grouped.values()))
                for taxon in ("Animalia","Fungi")}
    return {
       "schema":"gesui_global1055_Animalia_Fungi_2011_24_recording_count_SOURCE_BEFORE_fitting_v01",
       "UTC":datetime.now(timezone.utc).isoformat(),
       "original_source1055_WPI2019_WWp2015_manifest_SHA256":sha(root/S),
       "previously_frozen_original_1055_allPlantae_source_SHA256":sha(root/P),
       "original_previously_known_Plant_sample752_not_selected_with_animal_fungi":True,
       "original_port_source_N":1055,"original_source_countries_N":42,
       "original_2011_24_5km_GBIF_taxonKey_controls":{"Animalia":1,"Fungi":5},
       "four_original_1055_port_id_shards_source_sha256":{str(k):v for k,v in sorted(batches.items())},
       "official_GBIF_indexed_total_record_counts_only_NOT_taxon_richness":True,
       "per_taxon_GBIF_count_API_status":per_status,
       "source_port_count_with_both_kingdom_queries_successful":len(allgood),
       "previous_plant_quality752_ports_with_both_kingdom_queries_successful":len(plantgood),
       "original_source_countries_with_plant_quality_and_both_record_count":len(countries),
       "original_countries_with_2plus_source_Varying_WWp_and_two_record_controls":n_source_varied,
       "precommitted_SOURCE_QC_before_effect":{
         "min_both_kingdom_success_in_1055":1002,
         "min_existing_plant_quality_and_two_kingdom_count":700,
         "min_original_countries_with_plant_and_both":30,
         "min_countries_with_2plus_source_WWp_and_both":25,
         "cross_kingdom_recording_control_statistical_models_PERMITTED":admissible
       },
       "all_original_source_port_animal_and_fungi_indexed_counts":[grouped[i] for i in sorted(grouped)],
       "API_missing_not_biological_zero":True,
       "Animalia_and_Fungi_NOT_proven_WWp_unaffected_control_taxa":True,
       "original_GBIF_live_index_year_2026_not_synced_exactly_to_older_Plantae":True,
       "no_floristic_species_diversity_or_sewage_causal_effect_tested_here":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--shard",action="append",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if len(a.shard)!=4:raise ValueError("Four exact original 1055-port source shards required")
    root=Path(__file__).resolve().parents[1]
    out=join(root,a.shard)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in out.items() if k in (
        "source_port_count_with_both_kingdom_queries_successful",
        "previous_plant_quality752_ports_with_both_kingdom_queries_successful",
        "original_source_countries_with_plant_quality_and_both_record_count",
        "original_countries_with_2plus_source_Varying_WWp_and_two_record_controls",
        "per_taxon_GBIF_count_API_status","precommitted_SOURCE_QC_before_effect")},indent=2))

if __name__=="__main__":
    main()
