#!/usr/bin/env python3
"""Assemble exact original H7 plant speciesKey taxonomy and harmonized EIV NUTRI.

First report source-only taxonomy success and trait-score coverage. A
H7-period NUTRI score table is generated ONLY if hard preread source support
passes. No W, butterfly HNS, EQR, regressions or causal mediation used.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import requests

FACETS_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
EIV_URL="https://www.sci.muni.cz/botany/juice/NET/sources/EIV_Europe.csv"
TOTAL_UNION=11511
SHARDS=8
MIN_GBIF_SOURCE_SUCCESS=.99
MIN_SCORABLE=10
MIN_COV=.40
MIN_H7PERIOD=100
MIN_H6=40
MIN_COUNTRY=5

def sha_bytes(b):return hashlib.sha256(b).hexdigest()
def sha(p):return sha_bytes(p.read_bytes())
def locate(folder,name):
    arr=list(folder.rglob(name))
    if len(arr)!=1:raise RuntimeError("Source files missing / duplicate "+name+": "+repr([str(p) for p in arr]))
    return arr[0]
def simple_name(x):
    if not isinstance(x,str):return ""
    return re.sub(r"\s+"," ",x.strip()).casefold()
def pure_binomial(s):
    parts=s.split()
    return len(parts)==2 and bool(re.fullmatch(r"[a-z][a-z-]+",parts[0])) and bool(re.fullmatch(r"[a-z][a-z-]+",parts[1]))

def load_eiv():
    last=None
    for rep in range(4):
        try:
            r=requests.get(EIV_URL,timeout=(15,90),headers={"User-Agent":"gesui-academic-plant-habitat-trait-source-only/0.1"})
            r.raise_for_status()
            b=r.content
            if sha_bytes(b)!=EIV_SHA:raise RuntimeError("Published EIV_Europe original SHA256 mismatch "+sha_bytes(b))
            eiv=pd.read_csv(io.BytesIO(b),encoding="utf-8-sig",low_memory=False)
            if len(eiv)!=8908 or "NUTRI" not in eiv or "Unnamed: 0" not in eiv:
                raise RuntimeError("Original EIV European schema unexpectedly changed")
            eiv["name_norm"]=eiv["Unnamed: 0"].map(simple_name)
            eiv["num"]=pd.to_numeric(eiv["NUTRI"],errors="coerce")
            good=eiv["num"].between(1,12) & eiv.name_norm.map(pure_binomial)
            t=eiv.loc[good,["name_norm","num"]].copy()
            duplicate=t.groupby("name_norm")["num"].agg(["min","max"])
            ambiguous=set(duplicate.index[(duplicate["max"]-duplicate["min"])>1e-10])
            t=t.loc[~t.name_norm.isin(ambiguous)].drop_duplicates("name_norm")
            dictionary=dict(zip(t.name_norm,t.num.astype(float)))
            return dictionary,{"original_source_SHA256":EIV_SHA,
              "original_source_url":EIV_URL,"original_rows":len(eiv),
              "numeric_NUTRI_values_in_source":int(eiv["num"].notna().sum()),
              "unambiguous_exact_valid_binomials_with_NUTRI":len(dictionary),
              "ambiguous_duplicate_binomials_excluded":len(ambiguous),
              "source_verified_sha":True}
        except (requests.RequestException,OSError) as e:
            last=str(e)
    raise RuntimeError("Published EIV source request failed after retries: "+str(last))

def source_inputs(args):
    path=locate(args.source_dir,"GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json")
    if sha(path)!=FACETS_SHA:raise RuntimeError("Original 2808 source facets SHA mismatch")
    facets=json.loads(path.read_text())
    if len(facets)!=2808 or any(x.get("source_status")!="ok" for x in facets):
        raise RuntimeError("Original H7 source facet status/count changed")
    plants={(str(x["HYBAS_ID_L07"]),str(x["period"])):x for x in facets if x.get("group")=="Plantae"}
    if len(plants)!=468:raise RuntimeError("Not 234x2 source Plantae facets")
    allkeys={str(k) for r in plants.values() for k,v in r["source_speciesKey_counts"].items() if int(v)>0}
    if len(allkeys)!=TOTAL_UNION:raise RuntimeError("Union 11511 source plant speciesKey drift")
    if sha(args.panel)!=PANEL_SHA:raise RuntimeError("Original 454-panel SHA changed")
    with args.panel.open(encoding="utf-8-sig",newline="") as f:
        fullpanel=list(csv.DictReader(f))
    if len(fullpanel)!=454 or len({(x["h"],x["period"]) for x in fullpanel})!=454:
        raise RuntimeError("Original 454 H7-period sources missing/duplicated")
    if len({x["h"] for x in fullpanel})!=232 or len({x["h6"] for x in fullpanel})!=90:
        raise RuntimeError("Original H7/H6 mapping changed")
    return plants, allkeys, fullpanel

def load_shards(dir):
    all_tax={}
    shard=[]
    for n in range(SHARDS):
        fp=locate(dir,f"TAXONOMY_KEYS_SHARD_{n:02d}.json")
        manifest=locate(dir,f"TAXONOMY_SHARD_MANIFEST_{n:02d}.json")
        m=json.loads(manifest.read_text())
        if m["shard"]!=n or m["n_shards"]!=SHARDS or m["pilot_max_new_keys"]!=0:
            raise RuntimeError("Shards did not complete full original source mode")
        if m["input_facets_SHA256"]!=FACETS_SHA or m["external_pr12_EIV_SHA256"]!=EIV_SHA:
            raise RuntimeError("Shard original-source identity mismatch")
        if sha(fp)!=m["output_SHA256"]:raise RuntimeError("Shard content SHA mismatch")
        rows=json.loads(fp.read_text())
        if len(rows)!=m["original_keys_assigned_to_shard"] or len(rows)!=m["source_rows_returned"]:
            raise RuntimeError("Partial shard source rows")
        if any(int(x["speciesKey"])%SHARDS!=n for x in rows):
            raise RuntimeError("Taxon shard key misallocation")
        for x in rows:
            k=str(x["speciesKey"])
            if k in all_tax:raise RuntimeError("Taxon duplicated across shards "+k)
            all_tax[k]=x
        shard.append(m)
    if len(all_tax)!=TOTAL_UNION:raise RuntimeError("Original 11511 source GBIF IDs not all retained")
    return all_tax,shard

def taxon_score(a,dic):
    if a.get("status")!="OK":
        return None, "SOURCE_UNRESOLVED"
    if str(a.get("kingdom","")).casefold()!="plantae":
        return None, "NOT_PLANTAE_OR_UNKNOWN_KINGDOM"
    if str(a.get("rank","")).upper()!="SPECIES":
        return None,"NOT_GBIF_SPECIES_RANK"
    if str(a.get("taxonomicStatus","")).upper()!="ACCEPTED":
        return None,"NOT_ACCEPTED"
    if "tracheophyta" not in str(a.get("phylum","")).casefold():
        return None,"NOT_ACCEPTED_VASCULAR"
    n=simple_name(a.get("canonicalName"))
    if not pure_binomial(n):return None,"NOT_UNAMBIGUOUS_BINOMIAL"
    if n not in dic:return None,"EIV_NAME_OR_SCORE_NOT_MATCHED"
    return float(dic[n]),"EXACT_MATCH"

def make(args):
    plants,allkeys,panel=source_inputs(args)
    tax,shards=load_shards(args.shard_dir)
    if set(tax)!=allkeys:raise RuntimeError("Taxonomy keys not exactly original union")
    dic,eiv_source=load_eiv()
    for key,a in tax.items():
        if a.get("from")=="original_PR12" and a.get("status")=="OK":
            old=a.get("legacy_pr12_trait_NUTRI")
            new,_=taxon_score(a,dic)
            if old is None and new is None:continue
            if old is None or new is None or abs(float(old)-new)>1e-10:
                raise RuntimeError("Prior independently fixed PR12 exact NUTRI value mismatches fresh publisher dictionary for key "+key)
    cross={}
    for key in sorted(allkeys,key=int):
        a=tax[key];score,why=taxon_score(a,dic)
        cross[key]={"speciesKey":key,
                    "GBIF_source_status":a.get("status"),
                    "GBIF_source_origin":a.get("from"),
                    "gbif_canonical_name":a.get("canonicalName"),
                    "gbif_rank":a.get("rank"),"gbif_taxon_status":a.get("taxonomicStatus"),
                    "gbif_phylum":a.get("phylum"),"gbif_kingdom":a.get("kingdom"),
                    "EIV_NUTRI_exact":score,"trait_status":why}
    n_ok=sum(x["GBIF_source_status"]=="OK" for x in cross.values())
    source_ratio=n_ok/len(cross)
    rows=[]
    for x in panel:
        k=(str(x["h"]),str(x["period"]))
        if k not in plants:raise RuntimeError("H7-period missing original Plantae facet "+repr(k))
        facet=plants[k]
        spp={str(key):int(v) for key,v in facet["source_speciesKey_counts"].items() if int(v)>0}
        n=len(spp)
        num=int(facet["source_query_occurrence_count"])
        counts=np.array([spp[key] for key in spp],dtype=int)
        if sum(counts)>num:raise RuntimeError("Facet records unexpectedly greater than query count")
        scored=[(key,v,cross[key]["EIV_NUTRI_exact"]) for key,v in spp.items()
                if cross[key]["EIV_NUTRI_exact"] is not None]
        nscore=len(scored)
        sratio=nscore/n if n else None
        unresolved=sum(cross[key]["GBIF_source_status"]!="OK" for key in spp)
        unsupported=unresolved/n if n else 0
        if not all(int(v)>0 for _,v,_ in scored):
            raise RuntimeError("Zero original source occurrence count among included taxa")
        source_pass=num>=50 and nscore>=MIN_SCORABLE and sratio is not None and sratio>=MIN_COV
        rows.append({"h":k[0],"period":k[1],"h6":str(x["h6"]),
                     "country_holdout":str(x["country_holdout"]),
                     "all_plantae_source_records":num,
                     "source_indexed_speciesKey_records":int(counts.sum()),
                     "source_plantae_unique_speciesKeys":n,
                     "exact_NUTRI_scored_vascular_accepted_species":nscore,
                     "NUTRI_coverage_all_Plantae_speciesKeys":sratio,
                     "GBIF_taxonomy_unresolved_speciesKeys":unresolved,
                     "GBIF_taxonomy_unresolved_fraction":unsupported,
                     "full_source_trait_eligible":bool(source_pass)})
    passing=[x for x in rows if x["full_source_trait_eligible"]]
    h6=len({x["h6"] for x in passing})
    country=len({x["country_holdout"] for x in passing})
    pass_all=(source_ratio>=MIN_GBIF_SOURCE_SUCCESS and len(passing)>=MIN_H7PERIOD
              and h6>=MIN_H6 and country>=MIN_COUNTRY)
    status={}
    for a in cross.values():status[a["trait_status"]]=status.get(a["trait_status"],0)+1
    errors={}
    for a in cross.values():errors[a["GBIF_source_status"]]=errors.get(a["GBIF_source_status"],0)+1
    out={"schema":"GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V1",
         "source_facets_SHA256":FACETS_SHA,"source_H7_panel_SHA256":PANEL_SHA,
         "source_EIV":eiv_source,
         "original_Plantae_union_speciesKeys":len(cross),
         "original_cohort_H7period":454,"original_H7":232,"original_H6":90,
         "successful_species_taxon_sources":n_ok,"taxon_source_success_fraction":source_ratio,
         "source_status_count":errors,"exact_trait_join_status_count":status,
         "total_numeric_EIV_species_in_original_H7_union":
              sum(x["EIV_NUTRI_exact"] is not None for x in cross.values()),
         "quality_thresholds":{"raw_plantae_records_min":50,
               "numeric_EIV_species_min":MIN_SCORABLE,"all_Plantae_speciesKey_score_fraction_min":MIN_COV,
               "H7_period_min":MIN_H7PERIOD,"independent_H6_min":MIN_H6,
               "country_holdout_group_min":MIN_COUNTRY,"overall_API_source_success_min":MIN_GBIF_SOURCE_SUCCESS},
         "eligible_H7_period":len(passing),"eligible_H7":len({x["h"] for x in passing}),
         "eligible_true_H6":h6,"eligible_majority_country_groups":country,
         "hard_source_support_pass":bool(pass_all),
         "trait_coverage_all_plantae_q05_q50_q95":
               np.quantile([x["NUTRI_coverage_all_Plantae_speciesKeys"] for x in rows],[.05,.5,.95]).tolist(),
         "taxonomy_status_is_not_species_absence":True,
         "analysis_noncausal_and_no_W_or_HNS_opened":True,
         "GBIF_taxonomy_source_shards":shards}
    args.output.mkdir(parents=True,exist_ok=True)
    crosspath=args.output/"EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv"
    with crosspath.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(next(iter(cross.values()))))
        writer.writeheader();writer.writerows(cross.values())
    with (args.output/"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    (args.output/"EU454_FULL_EIV_SOURCE_GATE.json").write_text(json.dumps(
        out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    if pass_all:
        source_trait=[]
        for x in panel:
            key=(str(x["h"]),str(x["period"]))
            cnt=plants[key]["source_speciesKey_counts"]
            sc=[(float(cross[str(k)]["EIV_NUTRI_exact"]),int(v)) for k,v in cnt.items()
                if cross[str(k)]["EIV_NUTRI_exact"] is not None and int(v)>0]
            source_trait.append({"h":key[0],"period":key[1],
               "Plant_N_incidence":float(np.mean([v[0] for v in sc])),
               "Plant_N_record":float(np.average([v[0] for v in sc],weights=[v[1] for v in sc])),
               "n_scored_species":len(sc)})
        with (args.output/"EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv").open(
              "w",encoding="utf-8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(source_trait[0]))
            writer.writeheader();writer.writerows(source_trait)
    print("GESUI_EU454_PLANT_NUTRI_FULL_SOURCE_GATE="+json.dumps(
        {k:v for k,v in out.items() if k!="GBIF_taxonomy_source_shards"},
        ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)
    print("MODEL_GATE="+("PASS" if pass_all else "STOP_TRAIT_DATA_SUPPORT"),flush=True)
    return pass_all

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-dir",type=Path,required=True)
    ap.add_argument("--shard-dir",type=Path,required=True)
    ap.add_argument("--panel",type=Path,default=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"))
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    make(args)

