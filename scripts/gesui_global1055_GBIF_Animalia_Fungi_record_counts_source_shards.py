#!/usr/bin/env python3
"""GBIF original Animalia and Fungi indexed record COUNT comparators, source only.

Same 1055 original WPI sites, 2011–2024 occurrence years and 5km buffers as
prior all-Plantae. Same public API filters. Not species richness, not confirmed
negative-control biological organisms; counts proxy cross-kingdom RECORDING.
"""
from __future__ import annotations
import argparse,hashlib,json,math,time,urllib.error,urllib.parse,urllib.request
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

SOURCE=Path("results/global_all_Plantae_WWp2015_20261110/ORIGINAL_ALL1055_42COUNTRY_WWp2015_BEFORE_NEW808_GBIF_SOURCE.json")
PLANT=Path("results/global_all_Plantae_WWp2015_20261110/ORIGINAL_ALL1055_42COUNTRY_GBIF_PLANTAE_2011_24_50RECORD_RAREFACTION_SOURCE.json")
API="https://api.gbif.org/v1/occurrence/search"
PERIOD="2011,2024"
TAXA={"Animalia":1,"Fungi":5}
SHARDS=4
MIN_PAUSE=0.8
HEADERS={"Accept":"application/json","User-Agent":"gesui-global1055-2026-11-crosskingdom-recording-control/0.1"}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def blob_sha(b):return hashlib.sha256(b).hexdigest()

def request_count(lat,lon,taxon,last):
    param={
      "taxonKey":str(taxon),"year":PERIOD,
      "geoDistance":f"{float(lat):.7f},{float(lon):.7f},5km",
      "hasCoordinate":"true","hasGeospatialIssue":"false",
      "occurrenceStatus":"PRESENT","limit":"0"
    }
    url=API+"?"+urllib.parse.urlencode(param)
    for i in range(5):
        gap=max(0,MIN_PAUSE-(time.monotonic()-last[0]))
        if gap:time.sleep(gap)
        last[0]=time.monotonic()
        try:
            req=urllib.request.Request(url,headers=HEADERS)
            with urllib.request.urlopen(req,timeout=75) as response:
                data=response.read()
                if response.status!=200:raise ValueError("GBIF non200 response")
            obj=json.loads(data)
            if not isinstance(obj,dict) or type(obj.get("count")) is not int or obj["count"]<0 or obj.get("error"):
                raise ValueError("GBIF indexed count response malformed")
            return {"GBIF_status":"COMPLETE","original_indexed_record_count":obj["count"],
                    "request_url":url,"response_sha256":blob_sha(data),
                    "API_response_UTC":datetime.now(timezone.utc).isoformat(),
                    "zero_records_only_if_original_count_was_actually_zero":obj["count"]==0}
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,ValueError) as err:
            if isinstance(err,urllib.error.HTTPError) and err.code not in (408,425,429,500,502,503,504):raise
            if i==4:
                return {"GBIF_status":"API_ERROR_NOT_ZERO","original_indexed_record_count":None,
                        "request_url":url,"response_sha256":None,
                        "GBIF_error":(type(err).__name__+": "+str(err))[:240]}
            time.sleep(2.5*(i+1))
    raise RuntimeError("GBIF count source unreachable")

def run(root,shard,cache):
    source=json.loads((root/SOURCE).read_text(encoding="utf-8"))
    plant=json.loads((root/PLANT).read_text(encoding="utf-8"))
    if source["schema"]!="gesui_original_1055_WPI_ports42country_808NEW_2011_24_ALLPlantae_PREOUTCOME_frame_v01":
        raise ValueError("Original 1055 WPI WWp source no longer locked")
    if plant["original_all1055_WWp_WPI_country_source_frame_SHA256"]!=sha(root/SOURCE):
        raise ValueError("Original Plantae and WPI source identities mismatched")
    if source["all_original_country_WWp_source_sites_n"]!=1055 or source["original_country_count"]!=42:
        raise ValueError("Source all1055/42 frame does not match original positive hypothesis")
    sites=sorted(source["source_all_1055_public_country_port_site_metadata_no_new_GBIF_response"],key=lambda x:x["original_port_id"])
    taxa={x["original_port_id"]:x for x in plant["original_1055_port_all_Plantae_source_GBIF_facet_summaries"]}
    if len(sites)!=1055 or len(taxa)!=1055:raise ValueError("All 1055 source points not available")
    chosen=[s for i,s in enumerate(sites) if i%SHARDS==shard]
    if len(chosen)!=(264 if shard<3 else 263):
        raise ValueError("Unexpected original 1055 fixed port-ID 4-shard sizes")
    accepted={sid for sid,v in taxa.items() if v["source_quality_min50_faceted_records_min10_species_coverage70pct"]}
    if len(accepted)!=752:raise ValueError("Existing fixed full1055 Plantae source quality set changed")
    cache.parent.mkdir(parents=True,exist_ok=True)
    previous={}
    if cache.exists():
        for line in cache.read_text(encoding="utf-8").splitlines():
            if line.strip():
                v=json.loads(line)
                key=int(v["original_port_id"])
                if key in previous:raise ValueError("Duplicate source WPI port count query")
                previous[key]=v
    expected={s["original_port_id"] for s in chosen}
    if not set(previous).issubset(expected):raise ValueError("Extra non-source ports in GBIF cache")
    last=[0.]
    for i,s in enumerate(chosen,1):
        pid=s["original_port_id"]
        if pid in previous:continue
        responses={}
        for label,key in TAXA.items():
            responses[label]=request_count(s["original_port_lat"],s["original_port_lon"],key,last)
        item={
            "original_port_id":pid,
            "original_source_country":s["original_source_country"],
            "source_original_all_Plantae_quality_752_flag":pid in accepted,
            "source_2015_modelled_25km_WWp_log":float(s["source_WWp_log1p"]),
            "Animalia":responses["Animalia"],
            "Fungi":responses["Fungi"],
            "never_replace_API_error_as_true_absence":True
        }
        with cache.open("a",encoding="utf-8") as f:
            f.write(json.dumps(item,sort_keys=True,separators=(",",":"))+"\n")
        previous[pid]=item
        if i%30==0 or i==len(chosen):
            n=sum(all(z[label]["GBIF_status"]=="COMPLETE" for label in TAXA) for z in [{lab:v[lab] for lab in TAXA} for v in previous.values()])
            print("CROSSKINGDOM_COUNT_SHARD",shard,i,"/",len(chosen),"both_taxa_API_ok",n,"API_errors",len(previous)-n,flush=True)
    if set(previous)!=expected:raise ValueError("1055 original port 4-shard count source incomplete")
    total=Counter()
    for value in previous.values():
        for label in TAXA:
            total[(label,value[label]["GBIF_status"])]+=1
    result={
        "schema":"gesui_original_global1055_WWp2015_Animalia_Fungi_GBIF_count_only_4shard_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_WPI2019_1055_country2015_WWp_source_SHA256":sha(root/SOURCE),
        "original_2011_24_allPlantae_1055_response_source_SHA256":sha(root/PLANT),
        "original_query_taxonomy_ids":TAXA,
        "same_original_2011_24_year_and_5km_radius_as_original_Plantae":True,
        "original_fixed_country_source_752_plantae_quality_sites_not_reselected":True,
        "shard_id":shard,"shards_total":SHARDS,
        "original_1055_shard_expected_sites":len(chosen),
        "original_GBIF_crosskingdom_cache_SHA256":sha(cache),
        "per_original_source_port_taxon_indexed_record_counts":sorted(previous.values(),key=lambda x:x["original_port_id"]),
        "original_per_taxon_API_status":{label:dict(Counter(x[label]["GBIF_status"] for x in previous.values())) for label in TAXA},
        "GBIF_counts_are_recorded_occurrences_not_species_richness":True,
        "Animalia_and_Fungi_kinship_and_ecological_response_NOT_UNAFFECTED_negative_controls":True,
        "crosskingdom_2025_request_index_v1_not_time_frozen_archived_GBIF_download":True,
        "individual_GBIF_occurrence_records_not_acquired":True
    }
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--shard",type=int,required=True)
    p.add_argument("--cache",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if a.shard not in range(SHARDS):raise ValueError("Original GBIF recording-shard must be 0..3")
    root=Path(__file__).resolve().parents[1]
    out=run(root,a.shard,a.cache)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"shard":a.shard,"n":out["original_1055_shard_expected_sites"],
                      "status":out["original_per_taxon_API_status"]},indent=2))

if __name__=="__main__":
    main()
