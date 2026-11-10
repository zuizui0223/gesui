#!/usr/bin/env python3
"""Four fixed shards: GBIF all Plantae 2011-24 species facets at 808 NEW WWp ports.

Source original 1055 ports is SHA-frozen in Git before these requests. Reuse
previous 247 responses (not queried here). Four shards have independent
caches and strictly deterministic port IDs. GBIF failure is missing not zero.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
from gesui_global_42country_2011_24_allPlantae_GBIF_facet_census import one_port

D=Path("results/global_all_Plantae_WWp2015_20261110")
FROZEN=D/"ORIGINAL_ALL1055_42COUNTRY_WWp2015_BEFORE_NEW808_GBIF_SOURCE.json"
OLD=D/"ORIGINAL_247PORT_42COUNTRY_GBIF_ALL_PLANTAE_2011_24_RAREFACTION_SOURCE.json"
SHARDS=4

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(root,shard,cache):
    source=json.loads((root/FROZEN).read_text())
    old=json.loads((root/OLD).read_text())
    if source["schema"]!="gesui_original_1055_WPI_ports42country_808NEW_2011_24_ALLPlantae_PREOUTCOME_frame_v01":
        raise ValueError("Full 1055 global WPI source-only cohort not committed")
    if source["GBIF_new808_species_records_NOT_OBSERVED_in_stage0"] is not True:
        raise ValueError("New 808 source was selected using actual allflora outcomes")
    if source["original_previous_247_pilot_GBIF_all_Plantae_facet_summary_SHA256"]!=sha(root/OLD):
        raise ValueError("Prior 247 original allPlantae source not identical")
    original=source["source_all_1055_public_country_port_site_metadata_no_new_GBIF_response"]
    new=sorted([x for x in original if not x["previously_extracted_247_AllPlantae_pilot_source_port"]],
               key=lambda x:x["original_port_id"])
    if len(new)!=808 or len({x["original_port_id"] for x in new})!=808:
        raise ValueError("Fixed original 808 NEW WPI port IDs not available")
    mine=[row for i,row in enumerate(new) if i%SHARDS==shard]
    if len(mine)!=202:raise ValueError("Each of four deterministic source shards must have 202 ports")
    cache.parent.mkdir(parents=True,exist_ok=True)
    existing={}
    if cache.exists():
        for line in cache.read_text(encoding="utf-8").splitlines():
            if not line.strip():continue
            row=json.loads(line)
            k=row["original_port_id"]
            if k in existing:raise ValueError("Repeated raw facet port")
            existing[k]=row
    allowed={x["original_port_id"] for x in mine}
    if not set(existing).issubset(allowed):
        raise ValueError("Non-source GBIF plant facet present in shard cache")
    last_request=[0.]
    for i,row in enumerate(mine,1):
        port=row["original_port_id"]
        if port in existing:continue
        try:
            result,species=one_port(row,last_request)
            stored={**result,"original_speciesKey_count_pairs":sorted([[k,v] for k,v in species.items()]),
                    "GBIF_status":"COMPLETE","GBIF_error":None}
        except Exception as e:
            stored={
                "original_port_id":port,"source_country":row["original_source_country"],
                "source_quality_min50_faceted_records_min10_species_coverage70pct":False,
                "GBIF_status":"QUERY_FAILED_NOT_ZERO","GBIF_error":str(e)[:250],
                "original_speciesKey_count_pairs":[],
                "original_GBIF_query_and_response_SHA_pages":[]
            }
        with cache.open("a",encoding="utf-8") as f:
            f.write(json.dumps(stored,sort_keys=True,separators=(",",":"))+"\n")
        existing[port]=stored
        if i%30==0 or i==len(mine):
            print("ALL1055_GBIF_NEW808_SHARD",shard,i,"/",len(mine),
                  "query_success",sum(x["GBIF_status"]=="COMPLETE" for x in existing.values()),
                  "quality",sum(x["source_quality_min50_faceted_records_min10_species_coverage70pct"] for x in existing.values()),
                  "errors",sum(x["GBIF_status"]!="COMPLETE" for x in existing.values()),flush=True)
    if set(existing)!=allowed:raise ValueError("Original four shard site coverage incomplete")
    summary=[]
    for row in mine:
        o=existing[row["original_port_id"]]
        if o["source_country"]!=row["original_source_country"]:
            raise ValueError("Original WPI port source country changed")
        safe={k:v for k,v in o.items() if k!="original_speciesKey_count_pairs"}
        summary.append(safe)
    return {
        "schema":"gesui_ALL1055_GLOBAL_NEW808_DWCA_GBIF2011_24_Plantae_speciesKey_four_shards_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "source_original_1055_port_manifest_SHA256":sha(root/FROZEN),
        "old_247_GBIF_all_Plantae_source_SHA256":sha(root/OLD),
        "shard":shard,
        "number_of_source_shards":SHARDS,
        "expected_NEW808_site_total":808,
        "source_original_new_port_shard_site_count":len(mine),
        "source_original_2011_24_GBIF_plant_species_facets_cache_SHA256":sha(cache),
        "source_taxonKey":6,
        "new_gbif_query_failures_are_missing_not_zero":True,
        "original_public_2011_2024_Plantae_speciesKey_observations":summary,
        "total_success":sum(x["GBIF_status"]=="COMPLETE" for x in summary),
        "total_plant50record_supported":sum(x["source_quality_min50_faceted_records_min10_species_coverage70pct"] for x in summary),
        "no_original_observer_identifiers_or_photos_saved":True,
        "late_2026_LIVE_GBIF_index_MAY_DRIFT_from_old247":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--shard",type=int,required=True)
    p.add_argument("--cache",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if not 0<=a.shard<SHARDS:raise ValueError("Expected shard 0..3")
    root=Path(__file__).resolve().parents[1]
    result=run(root,a.shard,a.cache)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:result[k] for k in ("shard","source_original_new_port_shard_site_count","total_success","total_plant50record_supported")},indent=2))
if __name__=="__main__":
    main()
