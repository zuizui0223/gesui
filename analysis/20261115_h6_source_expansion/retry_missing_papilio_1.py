#!/usr/bin/env python3
"""Selective retry ONLY one previously missing original GBIF source facet
from 30 new H7×two periods×six groups, preserve 359 SHA-verified successes.

This is strict API-source missingness recovery, never a biological zero or
new data-driven cohort choice. Outputs source-provenance patched 3shard files
and a literal 360-facet merged archive (not a causal/ecological model).
"""
from __future__ import annotations
import argparse, csv, hashlib, json, sys, time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from gesui_H7_original_fivefamily_plants_speciesKey_retention import task,FAMILIES,PERIODS

MISSING=("2070055610","2017-2019","Papilionidae")
WKT_SHA="47f597e0a798cd6155ee9fb69ce8916ed5629951d046d114feb9f667963310e4"
TAXA_SHA="afc524ebfa353e4d1c6ea51d208ba529a85fedb57c4f0d28d83545010bb76189"
QA_SCHEMA="GESUI_SOURCE_ONLY_30_NEW_H7_2ERA_SIX_ORIGINAL_FACET_ACQUISITION_V1"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def one(base,name):
    paths=list(base.rglob(name))
    if len(paths)!=1:raise RuntimeError("Original source file absent/duplicated "+name)
    return paths[0]
def key(x):return(str(x["HYBAS_ID_L07"]),str(x["period"]),str(x["group"]))
def dump(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,ensure_ascii=False,sort_keys=True,separators=(",",":"))+"\n")

def main(args):
    allrecords={}
    manifest={}
    originalsha={}
    by_shard={}
    for i in range(3):
        d=args.original/str(i)
        f=one(d,f"EXTRA_EEA_H7_REAL_SPECIESKEY_full_SHARD_{i:02d}.json")
        q=one(d,f"EXTRA_EEA_H7_REAL_SPECIESKEY_full_SHARD_{i:02d}_QA.json")
        audit=json.loads(q.read_text())
        if audit.get("schema")!=QA_SCHEMA or audit.get("source_shard")!=i or audit.get("source_only_mode")!="full" or audit.get("original_source_1219WKT_SHA256")!=WKT_SHA or audit.get("original_234_fixed_family_taxa_source_SHA256")!=TAXA_SHA:
            raise RuntimeError("Original full extra H7 shard source gate wrong")
        if sha(f)!=audit["original_output_source_JSON_SHA256"]:
            raise RuntimeError("Source archived shard original data SHA mismatch")
        rows=json.loads(f.read_text())
        if len(rows)!=audit["source_6_group_2period_facets_requested"] or len(rows)!=audit["source_H7_in_shard"]*12:
            raise RuntimeError("Original source shard lost facet")
        originalsha[str(i)]=sha(f)
        for r in rows:
            k=key(r)
            if k in allrecords:raise RuntimeError("Duplicated original source facet")
            allrecords[k]=r
        by_shard[i]=(f,q,rows,audit)
    if len(allrecords)!=360:
        raise RuntimeError("Original 360 source facet keys not present")
    miss={k for k,r in allrecords.items() if r["source_status"]!="ok"}
    if miss!={MISSING}:
        raise RuntimeError("Expected ONLY 1 previously missing original source key: "+str(miss))
    for k,r in allrecords.items():
        if k!=MISSING and r["source_status"]!="ok":
            raise RuntimeError("Previously success changed to missing")
    w=one(args.wkt,"HYBAS7_1219_ORIGINAL_NATIVE_POLYGONS_WKT_BEFORE_NEW_GBIF_COUNTS.json")
    if sha(w)!=WKT_SHA:raise RuntimeError("1219 original WKT SHA drift")
    native=json.loads(w.read_text())
    sites={str(x["HYBAS_ID_L07"]):x for x in native["basins"]}
    original_site=sites[MISSING[0]]
    if original_site["status"]!="all_original_polygons_WKT_queryable":
        raise RuntimeError("Original exact missing H7 source polygon not queryable")
    fam=one(args.taxonomy,"HYBAS7_234_REPEATED_EU_EQR_ORIGINAL_WKT_TAXON_FAMILY_PRE_NEW_OUTCOME.json")
    if sha(fam)!=TAXA_SHA:raise RuntimeError("Original frozen family taxonomy SHA drift")
    fixed=json.loads(fam.read_text())["verified_butterfly_family_taxa"]
    if set(fixed)!=set(FAMILIES):raise RuntimeError("Original butterfly source family keys changed")
    time.sleep(5)
    ret=task(original_site,MISSING[1],MISSING[2],fixed)
    if key(ret)!=MISSING:
        raise RuntimeError("Source-retried API returned wrong original key")
    report={
       "schema":"GESUI_SOURCE_ONLY_30_NEW_H7_SELECTIVE_SINGLE_MISSING_FACET_API_RETRY_V1",
       "source_original_three_shard_SHA256":originalsha,
       "original_missing_H7_period_group":list(MISSING),
       "original_source_original_facets":360,
       "original_source_success_facets":359,
       "original_source_missing_facets":1,
       "attempted_only_exact_missing_original_facet":True,
       "original_successful_359_facet_records_unmodified":True,
       "retry_return_status":ret["source_status"],
       "retry_original_query_count":ret["source_query_occurrence_count"],
       "retry_speciesKey_n":ret["source_speciesKey_count"],
       "api_missing_not_true_absence":True,
       "no_EIV_HNS_or_ecological_W_regression_fitted":True
    }
    if ret["source_status"]=="ok":
        allrecords[MISSING]=ret
    args.output.mkdir(parents=True,exist_ok=True)
    source_after={k:v for k,v in allrecords.items() if k!=MISSING}
    source_before={k:v for k,v in allrecords.items() if k!=MISSING}
    if source_after!=source_before or len(source_after)!=359:
        raise RuntimeError("Successful original 359 source facets were not unchanged")
    for i in range(3):
        f,q,rows,audit=by_shard[i]
        outrows=rows.copy()
        if any(key(x)==MISSING for x in rows):
            outrows=[ret if key(x)==MISSING and ret["source_status"]=="ok" else x
                     for x in rows]
            audit["original_selective_recovery_from_SHA256"]=originalsha[str(i)]
            audit["original_selective_retry_key"]=list(MISSING)
            audit["source_facet_ok"]=sum(r["source_status"]=="ok" for r in outrows)
            audit["source_facet_missing"]=len(outrows)-audit["source_facet_ok"]
            audit["source_missing_which"]=[list(key(r)) for r in outrows if r["source_status"]!="ok"]
        audit["prior_source_archive_validated_original_SHA256"]=originalsha[str(i)]
        dataout=args.output/"patched_three_shards"/str(i)/f.name
        dump(dataout,outrows)
        audit["original_output_source_JSON_SHA256"]=sha(dataout)
        qout=args.output/"patched_three_shards"/str(i)/q.name
        qout.parent.mkdir(parents=True,exist_ok=True)
        qout.write_text(json.dumps(audit,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    sorted_records=[allrecords[k] for k in sorted(allrecords)]
    merged=args.output/"MERGED_ORIGINAL_30H7_TWO_PERIOD_SIX_FAMILY_360_SOURCE_FACETS.json"
    dump(merged,sorted_records)
    report["final_original_source_complete"]=all(x["source_status"]=="ok" for x in sorted_records)
    report["final_source_success_count"]=sum(x["source_status"]=="ok" for x in sorted_records)
    report["merged_original_source_SHA256"]=sha(merged)
    (args.output/"SELECTIVE_ONE_FACET_RECOVERY_QA.json").write_text(
        json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_SOURCE_30_NEW_H7_SINGLE_FACET_RETRY="+json.dumps(report,ensure_ascii=False,sort_keys=True),flush=True)
    print("NEXT_TRAIT_GATE="+("SOURCE_ALL_360_COMPLETE" if report["final_original_source_complete"] else "ONE_SOURCE_STILL_MISSING"),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--original",required=True,type=Path)
    p.add_argument("--wkt",required=True,type=Path)
    p.add_argument("--taxonomy",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
