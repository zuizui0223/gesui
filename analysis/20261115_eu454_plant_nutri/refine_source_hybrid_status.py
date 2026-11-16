#!/usr/bin/env python3
"""Reclassify official 200-OK hybrid taxon names without assigning any trait.

Source-v1 audit STOP is retained unaltered. This source-only refinement
distinguishes GBIF retrieval success but non-binomial hybrid canonicalName
missing from a real taxonomy acquisition error. 2023 EIV original exact-match
criterion and 99% API-resolution threshold are NOT relaxed. No W or HNS read.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

FACET_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
SOURCE_KEY_N=11511

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def one(folder,name):
    a=list(folder.rglob(name))
    if len(a)!=1:raise RuntimeError("Missing/duplicate source "+name+": "+repr([str(x) for x in a]))
    return a[0]
def readcsv(p):
    with p.open(newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))
def writecsv(path,rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
def known_nonbinomial_hybrid(x):
    if x.get("source_http")!=200:
        return False
    if x.get("canonicalName") is not None:
        return False
    if str(x.get("rank","")).upper()!="SPECIES" or str(x.get("taxonomicStatus","")).upper()!="ACCEPTED":
        return False
    if str(x.get("kingdom","")).casefold()!="plantae" or "tracheophyta" not in str(x.get("phylum","")).casefold():
        return False
    name=str(x.get("scientificName") or "").strip()
    # This only classifies explicitly source-reported nothospecies/crosses.
    # No fuzzy or surrogate accepted-species name is invented.
    return bool(re.search(r"×|(?:^|\s)[xX](?:\s|$)",name))

def main(args):
    base=args.base
    v1=json.loads(one(base,"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (v1.get("schema")!="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V1" or
        v1.get("source_facets_SHA256")!=FACET_SHA or
        v1.get("source_H7_panel_SHA256")!=PANEL_SHA or
        v1.get("source_EIV",{}).get("original_source_SHA256")!=EIV_SHA or
        v1.get("original_Plantae_union_speciesKeys")!=SOURCE_KEY_N):
        raise RuntimeError("Original v1 provenance was not authentic")
    if (v1.get("hard_source_support_pass") is not False or
        v1.get("source_status_count",{}).get("INVALID_API_SCHEMA")!=181):
        raise RuntimeError("Original source STOP and 181 missing canonical trait statuses not retained")
    cross=readcsv(one(base,"EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv"))
    if len(cross)!=SOURCE_KEY_N or len({x["speciesKey"] for x in cross})!=SOURCE_KEY_N:
        raise RuntimeError("Crosswalk original 11511 source keys changed")
    invalid={x["speciesKey"] for x in cross if x["GBIF_source_status"]=="INVALID_API_SCHEMA"}
    if len(invalid)!=181:raise RuntimeError("Original invalid source set mismatch")
    qa_path=one(args.diagnostic,"GBIF_181_UNRESOLVED_SOURCE_DIAGNOSTIC.json")
    qa=json.loads(qa_path.read_text())
    if (qa.get("schema")!="GESUI_EU454_181_UNRESOLVED_TAXONOMY_SOURCE_ONLY_AUDIT_V1" or
        qa.get("original_invalid_n")!=181 or
        len(qa.get("source_taxon_raw_api_responses",[]))!=181):
        raise RuntimeError("Official 181-key source diagnostic incomplete")
    restored={};failed={}
    for x in qa["source_taxon_raw_api_responses"]:
        key=str(x["speciesKey"])
        if key not in invalid:raise RuntimeError("Source diagnostic unexpected taxon key "+key)
        if known_nonbinomial_hybrid(x):
            restored[key]=x
        else:
            failed[key]=x
    if len(restored)+len(failed)!=181 or len({str(x["speciesKey"]) for x in qa["source_taxon_raw_api_responses"]})!=181:
        raise RuntimeError("Duplicate or incomplete official source rechecks")
    # Inference only on source validity; leave ALL NUTRI numeric matches identical.
    for x in cross:
        key=x["speciesKey"]
        if key in restored:
            if x["EIV_NUTRI_exact"] not in ("",None):
                raise RuntimeError("Original hybrid had a numeric EIV trait (must not change)")
            x["GBIF_source_status"]="OK_NONBINOMIAL_HYBRID_NO_EIV_SCORE"
            x["trait_status"]="KNOWN_NONBINOMIAL_HYBRID_EXCLUDED"
            x["gbif_canonical_name"]=""
    n_ok=sum(x["GBIF_source_status"] in ("OK","OK_NONBINOMIAL_HYBRID_NO_EIV_SCORE") for x in cross)
    success_ratio=n_ok/SOURCE_KEY_N
    elig=readcsv(one(base,"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv"))
    if len(elig)!=454:raise RuntimeError("Original 454 source trait coverage rows changed")
    passrows=[x for x in elig if x["full_source_trait_eligible"]=="True"]
    n_h6=len({x["h6"] for x in passrows})
    n_country=len({x["country_holdout"] for x in passrows})
    hard_pass=(success_ratio>=.99 and len(passrows)>=100 and n_h6>=40 and n_country>=5)
    if len(passrows)!=451 or n_h6!=90 or n_country!=11:
        raise RuntimeError("Source trait eligible cohort altered on taxon source-status refinement")
    # Reopen real 454 panel for H7 keys ONLY, not environmental W/HNS data.
    panel=args.panel
    if sha(panel)!=PANEL_SHA:raise RuntimeError("Original H7 science panel SHA mismatch")
    pnl=readcsv(panel)
    if len(pnl)!=454 or len({(x["h"],x["period"]) for x in pnl})!=454:
        raise RuntimeError("Fixed original H7 science panel cohort changed")
    original_path=one(args.facets,"GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json")
    if sha(original_path)!=FACET_SHA:raise RuntimeError("Original GBIF source facet SHA mismatch")
    original=json.loads(original_path.read_text())
    if len(original)!=2808:raise RuntimeError("Original source original facets changed")
    plants={(str(x["HYBAS_ID_L07"]),x["period"]):x for x in original if x["group"]=="Plantae"}
    if len(plants)!=468:raise RuntimeError("Not original 468 source Plantae facets")
    by_key={x["speciesKey"]:x for x in cross}
    # If source gate passes, compute only 451 SOURCE-ELIGIBLE source-weight
    # NUTRI indices. Never assign N to 181 verified known hybrids.
    source_scores=[]
    if hard_pass:
        valid={(x["h"],x["period"]) for x in passrows}
        for row in pnl:
            key=(row["h"],row["period"])
            if key not in valid:continue
            facet=plants.get(key)
            if facet is None:raise RuntimeError("Trait source H7-period source facet missing")
            matched=[(float(by_key[str(k)]["EIV_NUTRI_exact"]),int(v))
                     for k,v in facet["source_speciesKey_counts"].items()
                     if int(v)>0 and by_key[str(k)]["EIV_NUTRI_exact"] not in ("",None)]
            if len(matched)<10:raise RuntimeError("An eligible source outcome lacks ten verified EIV taxa")
            score=float(sum(n for n,w in matched)/len(matched))
            weighted=float(sum(n*w for n,w in matched)/sum(w for n,w in matched))
            source_scores.append({"h":key[0],"period":key[1],
                "Plant_N_incidence":score,"Plant_N_record":weighted,
                "n_scored_species":len(matched)})
        if len(source_scores)!=451:
            raise RuntimeError("Scored source trait subset not identical to 451 eligible H7-periods")
    out={
        **v1,"schema":"GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA",
        "v1_original_hard_source_support_pass":False,
        "original_v1_source_OK_n":int(v1["successful_species_taxon_sources"]),
        "original_v1_source_success_ratio":float(v1["taxon_source_success_fraction"]),
        "official_source_recheck_SHA256":sha(qa_path),
        "official_source_recheck_181":181,
        "official_confirmed_valid_unscorable_nonbinomial_hybrids":len(restored),
        "source_unresolved_after_explicit_hybrid_reclassification":len(failed),
        "source_reclassified_incomplete_canonical_not_nutrient_matched":True,
        "source_hybrid_score_assigned":0,
        "successful_species_taxon_sources":n_ok,
        "taxon_source_success_fraction":success_ratio,
        "source_status_count":{"OK":int(v1["source_status_count"]["OK"]),
             "OK_NONBINOMIAL_HYBRID_NO_EIV_SCORE":len(restored),
             "STILL_UNRESOLVED":len(failed)},
        "original_EIV_NUTRI_scored_keys_unchanged":sum(x["EIV_NUTRI_exact"] not in ("",None) for x in cross)==4011,
        "original_plant_ELIGIBILITY_H7period_unchanged":len(passrows)==451,
        "hard_source_support_pass":bool(hard_pass),
        "status_resolution_change_after_QA_not_new_trait_matching":True,
        "not_independent_advance_preregistration":True,
        "no_ecology_regression_yet":True}
    if not out["original_EIV_NUTRI_scored_keys_unchanged"] or not out["original_plant_ELIGIBILITY_H7period_unchanged"]:
        raise RuntimeError("Any unplanned EIV trait score or cohort manipulation prohibited")
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"EU454_FULL_EIV_SOURCE_GATE.json").write_text(json.dumps(
        out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    writecsv(args.output/"EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv",cross)
    writecsv(args.output/"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv",elig)
    if hard_pass:
        writecsv(args.output/"EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv",source_scores)
    print("GESUI_EU454_HYBRID_SOURCE_RECLASSIFICATION="+json.dumps({
        "original_v1_source_success_ratio":out["original_v1_source_success_ratio"],
        "verified_noncanonical_accepted_hybrids":len(restored),
        "still_unresolved":len(failed),
        "source_success_ratio_without_fabricated_N":success_ratio,
        "original_4011_EIV_trait_matches_unchanged":out["original_EIV_NUTRI_scored_keys_unchanged"],
        "original_451_H7period_quality_cohort_unchanged":out["original_plant_ELIGIBILITY_H7period_unchanged"],
        "hard_source_gate_after_factual_taxon_status_refinement":hard_pass
        },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)
    if failed:
        print("UNRESOLVED_KEYS_EXAMPLE="+json.dumps(list(failed)[:20]),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--base",required=True,type=Path)
    p.add_argument("--diagnostic",required=True,type=Path)
    p.add_argument("--facets",required=True,type=Path)
    p.add_argument("--panel",type=Path,default=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"))
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())

