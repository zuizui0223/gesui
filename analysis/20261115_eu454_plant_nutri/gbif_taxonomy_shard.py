#!/usr/bin/env python3
"""Original-source-locked EU H7 plant taxonomy acquisition (sharded; no effects).

Only look up source-original speciesKey taxon identities. Reuse the independently
obtained 1006-key PR12 mapping; do NOT query or fit W, HNS or site-level N scores.
No made-up trait values or re-queried GBIF ecological records.
"""
from __future__ import annotations
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import random
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

FACETS_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
SOURCE_FILE="GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json"
PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
PANEL_FILE="results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"
PR12_CROSS="GBIF_SPECIESKEY_EIV_EXACT_NAME_CROSSWALK.json"
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
TOTAL_PLANT_KEYS=11511
API="https://api.gbif.org/v1/species/"

def original_one(base,name):
    matches=list(base.rglob(name))
    if len(matches)!=1:raise RuntimeError("Original file missing or multiply sourced: "+str(matches))
    return matches[0]

def sha256(blob):
    return hashlib.sha256(blob).hexdigest()

def get_source(source_dir,cache_dir):
    path=original_one(source_dir,SOURCE_FILE)
    if sha256(path.read_bytes())!=FACETS_SHA:raise RuntimeError("Original H7 GBIF facets SHA256 changed")
    j=json.loads(path.read_text())
    if len(j)!=2808 or any(r.get("source_status")!="ok" for r in j):
        raise RuntimeError("Original source facet count/status changed")
    plant=[r for r in j if r.get("group")=="Plantae"]
    if len(plant)!=468:raise RuntimeError("Original 234x2 Plantae facets changed")
    # Original 468 source facets include 14 H7-period rows outside the fixed
    # true 454 analysis panel. Do NOT leak the extra 102 distinct taxa.
    panel_path=Path(PANEL_FILE)
    if sha256(panel_path.read_bytes())!=PANEL_SHA:
        raise RuntimeError("Genuine original 454 H7-period panel SHA mismatch")
    with panel_path.open(newline="",encoding="utf-8-sig") as handle:
        panel=list(csv.DictReader(handle))
    target={(str(x["h"]),str(x["period"])) for x in panel}
    if len(panel)!=454 or len(target)!=454:
        raise RuntimeError("Original analysis 454 rows missing/duplicated")
    selected=[f for f in plant if (str(f["HYBAS_ID_L07"]),str(f["period"])) in target]
    if len(selected)!=454 or len({(str(f["HYBAS_ID_L07"]),str(f["period"])) for f in selected})!=454:
        raise RuntimeError("Cannot reconcile authentic 454 cohort to 468 source original facets")
    keys=sorted({str(k) for f in selected for k,v in f["source_speciesKey_counts"].items() if int(v)>0},key=int)
    if len(keys)!=TOTAL_PLANT_KEYS:
        raise RuntimeError(f"454-cohort Plantae source-key union differs from original 11511: {len(keys)}")
    path=original_one(cache_dir,PR12_CROSS)
    previous=json.loads(path.read_text())
    if previous.get("source_SHA256")!=EIV_SHA or len(previous.get("species_taxonomy_and_matched_trait",[]))!=1006:
        raise RuntimeError("Prior independent EIV taxonomy source changed")
    cache={str(r["gbif_speciesKey"]):r for r in previous["species_taxonomy_and_matched_trait"]}
    if len(cache)!=1006:raise RuntimeError("Prior 1006-key mapping duplication")
    return keys,cache

def as_cached(k,a):
    if a.get("taxonomy_source_status")!="OK":
        return {"speciesKey":k,"status":"CACHED_PRIOR_NON_OK",
                "cached_original_status":a.get("taxonomy_source_status"),
                "from":"original_PR12","canonicalName":a.get("gbif_canonical_name"),
                "rank":a.get("gbif_rank"),"taxonomicStatus":a.get("gbif_taxon_status"),
                "phylum":a.get("gbif_phylum"),"kingdom":None}
    return {"speciesKey":k,"status":"OK","from":"original_PR12",
            "canonicalName":a.get("gbif_canonical_name"),
            "rank":a.get("gbif_rank"),"taxonomicStatus":a.get("gbif_taxon_status"),
            "phylum":a.get("gbif_phylum"),"kingdom":"Plantae",
            "legacy_pr12_exact":bool(a.get("exact_accepted_vascular_species")),
            "legacy_pr12_trait_NUTRI":a.get("trait_NUTRI")}

def fetch_taxon(k):
    err=None
    for attempt in range(5):
        try:
            req=Request(API+str(k),headers={
                "User-Agent":"gesui-academic-original-species-id-audit/20261115 (plant-trait source-only)",
                "Accept":"application/json"})
            with urlopen(req,timeout=22) as resp:
                body=resp.read(120000)
                status_code=getattr(resp,"status",200)
            a=json.loads(body)
            if not isinstance(a,dict) or not a.get("canonicalName") or not a.get("rank"):
                return {"speciesKey":k,"status":"INVALID_API_SCHEMA",
                        "http_status":status_code,"from":"new_GBIF_api",
                        "reason":"canonicalName or taxonomic rank missing"}
            if a.get("key") is not None and str(a["key"])!=str(k):
                return {"speciesKey":k,"status":"INVALID_API_KEY",
                        "from":"new_GBIF_api","returned_key":a.get("key")}
            return {"speciesKey":k,"status":"OK","from":"new_GBIF_api",
                    "canonicalName":a.get("canonicalName"),"scientificName":a.get("scientificName"),
                    "rank":a.get("rank"),"taxonomicStatus":a.get("taxonomicStatus"),
                    "phylum":a.get("phylum"),"kingdom":a.get("kingdom"),
                    "http_status":status_code}
        except HTTPError as e:
            err=f"HTTP {e.code}"
            if e.code in (400,404,410):
                return {"speciesKey":k,"status":"GBIF_KEY_NOT_FOUND","from":"new_GBIF_api",
                        "http_status":e.code,"reason":err}
            if e.code not in (403,408,429,500,502,503,504):
                break
        except (URLError,TimeoutError,OSError,json.JSONDecodeError) as e:
            err=type(e).__name__+": "+str(e)[:110]
        if attempt<4:
            time.sleep(min(24,2**attempt*2)+((int(k)%13)/13))
    return {"speciesKey":k,"status":"API_ERROR_AFTER_RETRIES","from":"new_GBIF_api","reason":err}

def main(args):
    if not 0<=args.shard<args.shards:raise RuntimeError("Invalid shard")
    keys,cache=get_source(args.source_dir,args.cache_dir)
    chosen=[k for k in keys if int(k)%args.shards==args.shard]
    if args.max_new_keys>0:
        cached=[k for k in chosen if k in cache][:8]
        new=[k for k in chosen if k not in cache][:args.max_new_keys]
        chosen=sorted(cached+new,key=int)
    n_cached=sum(k in cache for k in chosen)
    out=[as_cached(k,cache[k]) for k in chosen if k in cache]
    new_keys=[k for k in chosen if k not in cache]
    with ThreadPoolExecutor(max_workers=args.workers) as exe:
        futs={exe.submit(fetch_taxon,k):k for k in new_keys}
        for fut in as_completed(futs):
            out.append(fut.result())
    out.sort(key=lambda x:int(x["speciesKey"]))
    if len(out)!=len(chosen) or len({x["speciesKey"] for x in out})!=len(chosen):
        raise RuntimeError("Original source keys missing in shard output")
    status_count={}
    for x in out:status_count[x["status"]]=status_count.get(x["status"],0)+1
    manifest={
        "schema":"GESUI_EU454_GBIF_SPECIESKEY_TAXONOMY_SHARD_V1",
        "input_facets_SHA256":FACETS_SHA,"external_pr12_EIV_SHA256":EIV_SHA,
        "original_unique_specieskeys":len(keys),
        "original_all_source_468_facets_include_14_outside_cohort":True,
        "only_original_454_H7period_facets_selected":True,
        "shard":args.shard,"n_shards":args.shards,
        "pilot_max_new_keys":args.max_new_keys,
        "original_keys_assigned_to_shard":len([k for k in keys if int(k)%args.shards==args.shard]),
        "source_rows_returned":len(out),
        "reused_previous_taxonomy_keys":n_cached,
        "newly_queried_taxonomy_keys":len(new_keys),
        "source_status_count":status_count,
        "unresolved_api_keys_are_not_ecological_absences":True,
        "not_wastewater_or_biological_outcome_model":True}
    args.output.mkdir(parents=True,exist_ok=True)
    fname=f"TAXONOMY_KEYS_SHARD_{args.shard:02d}.json"
    (args.output/fname).write_text(json.dumps(out,ensure_ascii=False,separators=(",",":"),allow_nan=False)+"\n")
    manifest["output_file"]=fname
    manifest["output_SHA256"]=sha256((args.output/fname).read_bytes())
    (args.output/f"TAXONOMY_SHARD_MANIFEST_{args.shard:02d}.json").write_text(
        json.dumps(manifest,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_EU454_GBIF_TAXONOMY_SHARD="+json.dumps(manifest,sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-dir",type=Path,required=True)
    ap.add_argument("--cache-dir",type=Path,required=True)
    ap.add_argument("--shard",type=int,required=True)
    ap.add_argument("--shards",type=int,default=8)
    ap.add_argument("--workers",type=int,default=6)
    ap.add_argument("--max-new-keys",type=int,default=0)
    ap.add_argument("--output",type=Path,required=True)
    main(ap.parse_args())

