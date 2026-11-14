#!/usr/bin/env python3
"""Reproduce genus-resampled host FAMILY breadth from Ferrer-Paris 2014.

This is a *non-geographic* larval host catalogue. Not measured host use at a
sewage outfall or in a HydroBASINS cell. Only distinct Latin-like binomials
from six original-family categories are used as source taxon candidates.
Taxonomic acceptance, synonymy, life stage, and original report validity
must be assessed separately before any ecological/geographic inference.

Usage:
  python scripts/gesui_larval_host_rarefaction_sixfamilies_20261114.py \
    --csv 20140912_CompiledButterflyHostPlantRecords_JRFP.csv --out /tmp/hosts
No network calls. Outputs per-butterfly summary and group-bootstrap diagnostics.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
import pathlib
import random
import re

FAMILIES=("Lycaenidae","Nymphalidae","Pieridae","Papilionidae","Hesperiidae","Riodinidae")
CORE=FAMILIES[:5]
REQUIRED=("Lepidoptera Family","Lepidoptera Name","Hostplant Family",
          "Hostplant Genus","Hostplant Name","Source","URL")
LATIN=re.compile(r"^[A-Z][A-Za-z-]+\s+[a-z][A-Za-z-]+(?:\s|$)")
PLACEHOLDER=re.compile(r"\b(?:sp|spp|cf|aff)\.?($|\s)",re.I)

def binomial_like(s):
    return bool(LATIN.search(s)) and not bool(PLACEHOLDER.search(s))

def rare_family(host_to_family,m):
    """Expected host family count among m DISTINCT reported host taxon names."""
    n=len(host_to_family)
    if n<m:
        return None
    freq=collections.Counter(host_to_family.values())
    ret=0.0
    for nf in freq.values():
        if n-nf<m:
            miss=0.0
        else:
            miss=1.0
            for j in range(m):
                miss*=(n-nf-j)/(n-j)
        ret+=1.0-miss
    return ret

def load(csv_path,drop_wikipedia=False):
    raw=pathlib.Path(csv_path).read_bytes()
    species={}
    invalid=collections.Counter()
    host_conflicts=0
    with pathlib.Path(csv_path).open(encoding="utf-8-sig",newline="") as fh:
        reader=csv.DictReader(fh)
        if tuple(reader.fieldnames or ())!=REQUIRED:
            raise RuntimeError("Unexpected original Ferrer-Paris column schema")
        rows=0
        for r in reader:
            rows+=1
            fam=r["Lepidoptera Family"].strip()
            fly=r["Lepidoptera Name"].strip()
            hfamily=r["Hostplant Family"].strip()
            host=r["Hostplant Name"].strip()
            source=r["Source"].strip()
            if fam not in FAMILIES:
                invalid["out_of_scope_or_unclassified_butterfly_family"]+=1
                continue
            if not binomial_like(fly):
                invalid["butterfly_name_not_latin_binomial_like"]+=1
                continue
            if not hfamily or hfamily.upper()=="NA":
                invalid["unknown_plant_family"]+=1
                continue
            if drop_wikipedia and source=="Wikipedia":
                invalid["wikipedia_sensitivity_excluded"]+=1
                continue
            key=(fam,fly)
            if key not in species:
                species[key]={"family":fam,"butterfly":fly,
                              "genus":fly.split()[0],
                              "all_hosts":{},"binomial_hosts":{},
                              "sources":set()}
            t=species[key]
            t["sources"].add(source)
            for target in [t["all_hosts"]]+([t["binomial_hosts"]] if binomial_like(host) else []):
                if host in target and target[host]!=hfamily:
                    host_conflicts+=1
                    continue
                target[host]=hfamily
    if rows!=43645:
        raise RuntimeError(f"Original source row count mismatch: {rows} != 43645")
    if host_conflicts:
        raise RuntimeError(f"Conflicting host-plant family labels: {host_conflicts}")
    return species,{"n_original_rows":rows,
                    "original_csv_sha256":hashlib.sha256(raw).hexdigest(),
                    "original_Git_blob_sha1":"ea93a5c0fdb17e361a9b6a9637a497a12e90493c",
                    "excluded_source_rows":dict(invalid),
                    "taxonomy_is_not_accepted_name_matched":True}

def percentile(x,p):
    vals=sorted(x)
    return vals[min(len(vals)-1,int(p*(len(vals)-1)))]

def genus_resample_mean(groups,rng):
    keys=list(groups)
    n=len(keys)
    total=0.0;size=0
    for _ in range(n):
        vals=groups[keys[rng.randrange(n)]]
        total+=sum(vals)
        size+=len(vals)
    return total/size

def summarize_species(species,m,seed):
    groups={}
    for family in CORE:
        g=collections.defaultdict(list)
        for t in species.values():
            if t["family"]!=family:
                continue
            val=rare_family(t["binomial_hosts"],m)
            if val is not None:
                g[t["genus"]].append(val)
        groups[family]=dict(g)
    rng=random.Random(seed)
    B=4000
    summary={};draws={}
    for family in CORE:
        g=groups[family]
        flat=[v for vals in g.values() for v in vals]
        if len(flat)<20 or len(g)<5:
            raise RuntimeError("Comparative taxon genus/species eligibility FAIL "+family)
        obs=sum(flat)/len(flat)
        draw=[genus_resample_mean(g,rng) for _ in range(B)]
        draws[family]=draw
        summary[family]={"n_butterfly_species":len(flat),
                         "n_genera":len(g),
                         "mean_expected_host_families":obs,
                         "descriptive95_genus_bootstrap":[percentile(draw,.025),percentile(draw,.975)]}
    contrasts={}
    for other in CORE[1:]:
        z=[a-b for a,b in zip(draws["Lycaenidae"],draws[other])]
        contrasts["Lycaenidae_minus_"+other]={
            "observed":summary["Lycaenidae"]["mean_expected_host_families"]-
                       summary[other]["mean_expected_host_families"],
            "descriptive95_genus_bootstrap":[percentile(z,.025),percentile(z,.975)],
            "bootstrap_fraction_positive":sum(q>0 for q in z)/B}
    return {"m_host_taxa_rarefaction":m,"families":summary,"contrasts":contrasts}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--csv",required=True)
    p.add_argument("--out",required=True)
    args=p.parse_args()
    out=pathlib.Path(args.out);out.mkdir(parents=True,exist_ok=True)
    full,qa=load(args.csv)
    without_wikipedia,qa_nowiki=load(args.csv,True)
    if len(full)!=3785:
        raise RuntimeError("Original candidate butterfly taxon count changed; re-audit taxonomy")
    records=[]
    for t in sorted(full.values(),key=lambda z:(z["family"],z["butterfly"])):
        bb=t["binomial_hosts"]
        records.append({
            "family":t["family"],"butterfly_name_as_reported":t["butterfly"],
            "butterfly_genus":t["genus"],"source_types_n":len(t["sources"]),
            "host_taxon_names_all_unique_n":len(t["all_hosts"]),
            "host_binomial_names_unique_n":len(bb),
            "host_families_distinct_binomial_n":len(set(bb.values())),
            "E_host_family_in_5_distinct_host_names":rare_family(bb,5),
            "E_host_family_in_10_distinct_host_names":rare_family(bb,10)})
    with (out/"FERRER_PARIS_2014_BUTTERFLY_SPECIES_HOST_FAMILY_RAREFACTION_SOURCE_DERIVED.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
    summary={
        "source":"Ferrer-Paris et al. (2014) compiled butterfly larval host data, doi:10.6084/m9.figshare.1168861.v1",
        "source_GloBI_repository":"https://github.com/globalbioticinteractions/ferrer-paris",
        "source_QA":qa,
        "without_wikipedia_source_QA":qa_nowiki,
        "source_observations_not_geographic_and_not_experimental":True,
        "not_measured_butterfly_or_host_occurrence_in_wastewater_basins":True,
        "not_source_taxon_accepted_name_verified":True,
        "model":"uniform hypergeometric rarefaction among distinct plant binomial names; genus-resampled descriptive CI of mean, 4000 replicates; NOT phylogenetic independent contrasts",
        "full_catalog":[summarize_species(full,5,20261114),summarize_species(full,10,20261115)],
        "without_Wikipedia":[summarize_species(without_wikipedia,5,20261116),summarize_species(without_wikipedia,10,20261117)]}
    (out/"FERRER_PARIS_HOST_FAMILY_RAREFACTION_SIXFAMILY_RESULTS.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True))
    print(json.dumps({"QA":qa,"full_catalog":summary["full_catalog"],
                      "without_Wikipedia":summary["without_Wikipedia"]},
                     ensure_ascii=False,sort_keys=True))
if __name__=="__main__":
    main()
