#!/usr/bin/env python3
"""Frozen-source EIV 1006-key *lower-bound* overlap with original 454 H7 Plantae.

Input EIV crosswalk was constructed independently for 36 earlier EEA source
sites, not fitted on 454 H7. Only 624 validated, EXACT, scored keys within
that restricted set are used here. No whole-community NUTRI means or W fits
are estimated, even if this lower-bound source support passes.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path

SOURCE_SHA = "251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
EIV_SHA = "ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"

def unique(root,name):
    a=list(root.rglob(name))
    if len(a)!=1:
        raise RuntimeError("Expected exactly one original source file "+name+": "+repr([str(z) for z in a]))
    return a[0]

def q(v,p):
    a=sorted(v);k=(len(a)-1)*p;lo=int(k);hi=min(lo+1,len(a)-1)
    return a[lo]+(a[hi]-a[lo])*(k-lo)

def main(args):
    original=unique(args.source_dir,"GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json")
    if hashlib.sha256(original.read_bytes()).hexdigest()!=SOURCE_SHA:
        raise RuntimeError("Source SHA mismatch")
    records=json.loads(original.read_text())
    if len(records)!=2808 or any(x.get("source_status")!="ok" for x in records):
        raise RuntimeError("Not the genuine 2808 source facets")
    plants={(str(x["HYBAS_ID_L07"]),str(x["period"])):x for x in records if x.get("group")=="Plantae"}
    if len(plants)!=468:raise RuntimeError("Original 468 Plantae facets changed")
    cross=unique(args.eiv_dir,"GBIF_SPECIESKEY_EIV_EXACT_NAME_CROSSWALK.json")
    eiv=json.loads(cross.read_text())
    if eiv.get("source_SHA256")!=EIV_SHA:raise RuntimeError("Independent original EIV checksum mismatch")
    taxon=eiv.get("species_taxonomy_and_matched_trait")
    if not isinstance(taxon,list) or len(taxon)!=1006:
        raise RuntimeError("Original PR12 1006 taxonomy metadata source invalid")
    roster={str(x["gbif_speciesKey"]):x for x in taxon}
    if len(roster)!=1006:raise RuntimeError("Duplicate original PR12 plant keys")
    known={}
    for key,rec in roster.items():
        n=rec.get("trait_NUTRI")
        if n is not None:
            n=float(n)
            if not 1<=n<=12 or not rec.get("exact_accepted_vascular_species") or not rec.get("EIV_name_exact_match"):
                raise RuntimeError("An unverified key was passed off as exact EIV match")
            known[key]=n
    if len(known)!=624:
        raise RuntimeError("Expected exactly 624 independent scored keys, got "+str(len(known)))
    with args.panel.open(encoding="utf-8-sig",newline="") as f:
        panel=list(csv.DictReader(f))
    if len(panel)!=454:raise RuntimeError("Original H7 panel count drift")
    if len({(str(x["h"]),str(x["period"])) for x in panel})!=454:
        raise RuntimeError("Duplicate original H7-period")
    rows=[];union_all=set()
    for x in panel:
        key=(str(x["h"]),str(x["period"]))
        if key not in plants:raise RuntimeError("Missing original plant facet "+repr(key))
        spp={str(k):int(v) for k,v in plants[key]["source_speciesKey_counts"].items() if int(v)>0}
        union_all.update(spp)
        n=len(spp)
        intersect=set(spp)&set(roster)
        score=set(spp)&set(known)
        c=len(score)/n if n else 0.
        rows.append({
           "h":key[0],"period":key[1],"h6":str(x["h6"]),
           "country_holdout":str(x["country_holdout"]),
           "original_plantae_speciesKey_n":n,
           "original_source_species_index_record_n":sum(spp.values()),
           "PR12_1006_catalogue_overlap_species_n":len(intersect),
           "PR12_624_verified_EIV_numeric_species_n":len(score),
           "verified_PRIOR_SUBSET_EIV_species_coverage_lowerbound":c,
           "lowerbound_passing_original_10_EIV_40pct_condition":len(score)>=10 and c>=.4,
           "not_full_EIV_join":True
        })
    passing=[x for x in rows if x["lowerbound_passing_original_10_EIV_40pct_condition"]]
    fr=[x["verified_PRIOR_SUBSET_EIV_species_coverage_lowerbound"] for x in rows]
    out={
      "schema":"GESUI_H7_EIV_PRIOR1006_LOWERBOUND_COVERAGE_V1",
      "source_H7_original_SHA256":SOURCE_SHA,
      "source_EIV_original_SHA256":EIV_SHA,
      "external_independent_PR12_matched_keys":len(roster),
      "external_independent_PR12_numeric_NUTRI_keys":len(known),
      "original_H7_period_count":len(rows),
      "original_H7":len({x["h"] for x in rows}),
      "original_true_H6":len({x["h6"] for x in rows}),
      "prior_1006_crosswalk_H7_union_species_overlap":len(union_all&set(roster)),
      "prior_624_scored_H7_union_species_overlap":len(union_all&set(known)),
      "original_all_H7_union_species":len(union_all),
      "species_based_subset_EIV_coverage_fraction_median":q(fr,.5),
      "species_based_subset_EIV_coverage_fraction_q05_q95":[q(fr,.05),q(fr,.95)],
      "lowerbound_EIV_complete_40pct_10species_H7periods":len(passing),
      "lowerbound_EIV_complete_40pct_10species_H6":len({x["h6"] for x in passing}),
      "lowerbound_EIV_complete_40pct_10species_countries":len({x["country_holdout"] for x in passing}),
      "prior_1006_subset_alone_exceeds_min_required_100rows_40H6_5countries":
        len(passing)>=100 and len({x["h6"] for x in passing})>=40 and len({x["country_holdout"] for x in passing})>=5,
      "caution":"Lower bound using only previously independently classified 1006 original GBIF species keys, not a complete EU454 source taxonomy crosswalk. Passing sites have at least this many valid scores; failing sites are UNKNOWN until all original species are taxonomically matched. No coefficient, plant N community mean or mediation was fitted.",
      "no_model_fitted":True
    }
    args.output.mkdir(parents=True,exist_ok=True)
    with (args.output/"EU454_PR12_EIV_LOWERBOUND_PER_H7.csv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (args.output/"EU454_EIV_PRIOR1006_LOWERBOUND_SUMMARY.json").write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_EU454_PRIOR_EIV_LOWERBOUND="+json.dumps(out,sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--panel",type=Path,default=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"))
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--eiv-dir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())

