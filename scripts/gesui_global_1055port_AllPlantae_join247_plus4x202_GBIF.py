#!/usr/bin/env python3
"""Merge 247 ORIGINAL prior pilot GBIF allflora sites + four NEW808 source shards.

Zero source-site reselection; original 1055 WPI sites and original 42 WWp
countries frozen BEFORE new GBIF records. No downstream ecological model
allowed if source-availability gate fails.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

BASE=Path("results/global_all_Plantae_WWp2015_20261110")
SOURCE=BASE/"ORIGINAL_ALL1055_42COUNTRY_WWp2015_BEFORE_NEW808_GBIF_SOURCE.json"
OLD=BASE/"ORIGINAL_247PORT_42COUNTRY_GBIF_ALL_PLANTAE_2011_24_RAREFACTION_SOURCE.json"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def build(root,paths):
    source=json.loads((root/SOURCE).read_text())
    old=json.loads((root/OLD).read_text())
    if source["schema"]!="gesui_original_1055_WPI_ports42country_808NEW_2011_24_ALLPlantae_PREOUTCOME_frame_v01":
        raise ValueError("Original 1055 source-first manifest changed")
    if source["original_previous_247_pilot_GBIF_all_Plantae_facet_summary_SHA256"]!=sha(root/OLD):
        raise ValueError("Pilot source GBIF response changed after original 1055 site selection")
    original=source["source_all_1055_public_country_port_site_metadata_no_new_GBIF_response"]
    idsource={x["original_port_id"]:x for x in original}
    if len(idsource)!=1055 or len({x["original_source_country"] for x in original})!=42:
        raise ValueError("Original 1055/42 source IDs not present")
    reuse={x["original_port_id"]:x for x in old["original_site_GBIF_facet_summary"]}
    if len(reuse)!=247 or any(not idsource[k]["previously_extracted_247_AllPlantae_pilot_source_port"] for k in reuse):
        raise ValueError("Original 247 allflora pilot subset cannot be reused safely")
    fresh={}
    shardhash={}
    for path in paths:
        p=Path(path)
        d=json.loads(p.read_text())
        i=d["shard"]
        if i in shardhash or d["number_of_source_shards"]!=4:
            raise ValueError("Duplicate or wrong original four-source shard receipt")
        if d["source_original_1055_port_manifest_SHA256"]!=sha(root/SOURCE) or d["old_247_GBIF_all_Plantae_source_SHA256"]!=sha(root/OLD):
            raise ValueError("Four GBIF shards used different source country/port frame")
        rows=d["original_public_2011_2024_Plantae_speciesKey_observations"]
        if len(rows)!=202:raise ValueError("Original 808 sites not split exactly into 202 each")
        for r in rows:
            pid=r["original_port_id"]
            if pid in fresh or pid in reuse or pid not in idsource:
                raise ValueError("Source 808 shard includes old pilot or duplicated/new foreign port")
            if idsource[pid]["previously_extracted_247_AllPlantae_pilot_source_port"]:
                raise ValueError("New source query re-fetched a prior pilot")
            if r["source_country"]!=idsource[pid]["original_source_country"]:
                raise ValueError("Changed original source country in new GBIF facet query")
            fresh[pid]=r
        shardhash[str(i)]={"source_receipt_SHA256":sha(p),
                         "source_original_facet_cache_SHA256":d["source_original_2011_24_GBIF_plant_species_facets_cache_SHA256"],
                         "success":d["total_success"],"quality_pass":d["total_plant50record_supported"]}
    if len(fresh)!=808 or set(shardhash)!={"0","1","2","3"}:
        raise ValueError("Four shard joins did not produce exact 808 previous-unqueried original WPI ports")
    results={**reuse,**fresh}
    if len(results)!=1055 or set(results)!=set(idsource):
        raise ValueError("Full 1055 original source port population incomplete or duplicated")
    rows=[results[k] for k in sorted(idsource)]
    bycountry=defaultdict(list)
    complete=Counter()
    oldnew={"pilot247":Counter(),"new808":Counter()}
    for r in rows:
        pid=r["original_port_id"]
        isold=idsource[pid]["previously_extracted_247_AllPlantae_pilot_source_port"]
        src="pilot247" if isold else "new808"
        if r["GBIF_status"]=="COMPLETE":
            complete["query_success"]+=1
            oldnew[src]["query_success"]+=1
        else:
            complete["query_error"]+=1
            oldnew[src]["query_error"]+=1
        if r["source_quality_min50_faceted_records_min10_species_coverage70pct"]:
            complete["predeclared_quality_success"]+=1
            oldnew[src]["quality_success"]+=1
            bycountry[idsource[pid]["original_source_country"]].append(pid)
        else:
            complete["quality_ineligible_or_API_error"]+=1
    countries_variation=sum(len(v)>=2 and len({idsource[k]["source_WWp_log1p"] for k in v})>=2 for c,v in bycountry.items())
    source_gate=(complete["predeclared_quality_success"]>=550 and
                 len(bycountry)>=35 and countries_variation>=30)
    return {
        "schema":"gesui_original_1055_WPI42country_2011_2024_AllPlantae_50records_247old_808new_GBIF_SOURCE_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_all1055_WWp_WPI_country_source_frame_SHA256":sha(root/SOURCE),
        "old_original_247_WPI_GBIF_all_Plantae_source_result_SHA256":sha(root/OLD),
        "four_new_2026_11_10_GBIF_allPlantae_source_shards_SHA256":shardhash,
        "previous_pilot_247_AllPlantae_GBIF_queries_NOT_REPEATED":True,
        "new_source_port_ids_expected":808,
        "new_source_GBIF_species_records_2011_24":"PUBLIC GBIF taxonKey=6 geoDistance=5km year=2011,2024 facet=speciesKey",
        "original_full_42country_all1055_site_count":1055,
        "original_GBIF_allPlantae_query_and_coverage":{"all_sites":len(rows),
                 **dict(complete),
                 "source_qualifying_countries":len(bycountry),
                 "qualifying_countries_with_2plus_WWp_distinct_sites":countries_variation,
                 "pilot247_vs_new808":{k:dict(v) for k,v in oldnew.items()}},
        "original_1055_port_all_Plantae_source_GBIF_facet_summaries":rows,
        "precommitted_full_1055_AllPlantae_response_support_GATE":{
            "min_quality_valid_ports":550,
            "min_quality_valid_countries":35,
            "min_countries_with_2plus_site_WWp_variation":30,
            "ALL1055_AllPlantae_WWp2015_macro_regression_ALLOWED":source_gate},
        "country_42_WWp_high_low_source_sample_not_selected_by_new808_GBIF_plant_sign":True,
        "all_2015_WWp_is_modelled_wastewater_generated_not_real_local_nutrient_discharge":True,
        "all_2011_24_Plantae_occurrences_are_GBIF_speciesKey_record_counts_not_census_abundance":True,
        "full_allflora_source_GBIF_v1_late2026_index_may_drift_with_old247_dates":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--shard",action="append",required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if len(a.shard)!=4:raise ValueError("Expected four precommitted source shard aggregate receipts")
    root=Path(__file__).resolve().parents[1]
    result=build(root,a.shard)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"source":result["original_GBIF_allPlantae_query_and_coverage"],
                      "model_allowed":result["precommitted_full_1055_AllPlantae_response_support_GATE"]},indent=2))
if __name__=="__main__":
    main()
