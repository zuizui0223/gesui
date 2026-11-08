#!/usr/bin/env python3
"""Read-only, SHA-locked transfer of existing chocho Lycaenidae resource results.

Not a new wastewater ecology test. Never relabel potential larval host resource
geography as butterfly occurrence/abundance, host quality, adult puddling, or
WWTP effect. The 25/22 family summary is pre-existing in chocho.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,urllib.request
from pathlib import Path

REF="e8b8f86c7a889b6a3224a4ef69ffe5e0c7738660"
REPO="zuizui0223/chocho"
SOURCES={
 "species":{
  "path":"data/frozen/figure_sources/anthropogenic_species_metrics.csv",
  "git_blob_sha1":"9f83d2d9ac8e46f1cfd4ded572c7a023273503dd"
 },
 "family":{
  "path":"provenance/reviewer_defenses/results/butterfly_taxonomic_family_expansion_v0.1.json",
  "git_blob_sha1":"95443e82e0c2ecb8e8a835e2c07ab8b62e7af27f"
 }}
FOCAL=["Lycaena phlaeas","Polyommatus icarus",
       "Brephidium exilis","Pseudozizeeria maha"]

def get_src(info):
    url=f'https://raw.githubusercontent.com/{REPO}/{REF}/{info["path"]}'
    request=urllib.request.Request(url,headers={
        "User-Agent":"gesui-read-only-chocho-frozen-resource-audit/1"})
    with urllib.request.urlopen(request,timeout=45) as r:
        raw=r.read()
    sha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if sha!=info["git_blob_sha1"]:
        raise ValueError(f'Original chocho Git blob mismatch: {info["path"]}')
    return raw,{"original_git_commit":REF,"original_github_blob_sha1":sha,
                "source_bytes_sha256":hashlib.sha256(raw).hexdigest(),"source_url":url}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    s,sinfo=get_src(SOURCES["species"])
    f,finfo=get_src(SOURCES["family"])
    species=list(csv.DictReader(io.StringIO(s.decode("utf-8-sig"))))
    family=json.loads(f.decode("utf-8"))
    if len(species)!=239 or family["resource_species"]!=239:
        raise ValueError("Mismatch in frozen species universe")
    if sum(float(x["introduced_added_units"])>0 for x in species)!=206:
        raise ValueError("Main positive-resource denominator changed")
    lyc=next(x for x in family["all_families"] if x["family"]=="Lycaenidae")
    if (lyc["species"],lyc["expanded_species"])!=(25,22):
        raise ValueError("Lycaenidae original family summary changed")
    by_name={x["species"]:x for x in species}
    out={}
    for name in FOCAL:
        r=by_name[name]
        b=int(r["native_resource_units"]);n=int(r["contemporary_resource_units"])
        d=int(r["introduced_added_units"])
        if n-b!=d:raise ValueError("Contradictory resource coverage: "+name)
        out[name]={
          "native_only_resource_region_units":b,
          "including_introduced_host_resource_region_units":n,
          "additional_region_units":d,
          "ratio":n/b if b else None,
          "percentage_added_vs_native_only":100*d/b if b else None,
          "resolved_host_species":int(r["resolved_host_species"]),
          "original_host_taxonomy_lower_bound_adequate":
              r["host_taxonomy_lower_bound_adequate"].lower()=="true",
          "not_butterfly_population_or_wwtp_exposure":True,
          "illustrative_species_chosen_posthoc":True,
        }
    result={
      "schema":"gesui_lycaenid_preexisting_chocho_resource_readonly_audit_v01",
      "publication_status":"EXISTING_CHOCHO_DATA__NOT_NEW_WASTEWATER_ANALYSIS",
      "source_original_commit":REF,
      "sources":{"species":sinfo,"family":finfo},
      "all_butterfly_species_n":239,
      "all_butterfly_resource_expanding_species_n":206,
      "Lycaenidae":{"species_n":25,"expanded_n":22,"expanded_pct":88,
        "median_added_region_units":lyc["median_added_units"],
        "median_host_family_count":lyc["median_host_family_count"],
        "comparison_not_significant_tested":True},
      "examples":out,
      "interpretation_constraints":[
        "Potential larval host resource geography is not larval performance, abundance or colonization",
        "Lycaenidae family 22/25 expansion is not distinct from other major butterfly families",
        "Pseudozizeeria maha has host_taxonomy_lower_bound_adequate=False and must not anchor primary quantitative inference",
        "None of these results measure WWTP location, sewage effluent, nitrogen concentration, salt, albumin, or butterfly puddling",
        "Species illustrations selected after observing results; not a predeclared confirmatory case study",
      ],
    }
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"CHOCHO_LYCAENIDAE_FROZEN_RESOURCE_AUDIT.json").write_text(
       json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    md=[
     "# Source-locked legacy chocho Lycaenidae findings (NOT wastewater results)","",
     "This audit uses only 239-species frozen resource geography and its independent family summary at exact chocho commit "+REF+".",
     "The two original Git blobs are checked byte-for-byte using Git SHA-1. No new GBIF or EEA ecology outcome is queried.","",
     "| Summary | Existing evidence |",
     "|---|---:|",
     "| Lycaenidae with introduced-host expanded potential resources | **22/25 (88.0%)** |",
     "| Family median additional resource regions | **73** |",
     "| All butterfly species with expanded potential resource opportunity | **206/239** |","",
     "| Illustrative species (posthoc) | Native-only regions | Contemporary including introduced hosts | Added | Original host-taxonomy quality |",
     "|---|---:|---:|---:|---|",
    ]
    for name,q in out.items():
        md.append(f'| *{name}* | {q["native_only_resource_region_units"]} | {q["including_introduced_host_resource_region_units"]} | +{q["additional_region_units"]} | {"adequate" if q["original_host_taxonomy_lower_bound_adequate"] else "NOT adequate"} |')
    md += ["","**Interpretation:** These are potential plant-host resource regions, not actual butterfly population change or wastewater responses. The extreme *Pseudozizeeria maha* ratio (16→275) has a non-adequate host taxonomy support flag and should **not** anchor the main claim.","",
      "Unresolved: whether local nutrient availability, adult organic nitrogen puddling, ant mutualism, host quality and effluent contamination interact. Every such connection is an **unmeasured hypothesis** at this stage.","",
      "Original sources:","",
      "- https://github.com/zuizui0223/chocho/blob/"+REF+"/"+SOURCES["species"]["path"],
      "- https://github.com/zuizui0223/chocho/blob/"+REF+"/"+SOURCES["family"]["path"],
      ""]
    (a.out/"README.md").write_text("\n".join(md),encoding="utf-8")
    print("GESUI_CHOCHO_LYCAENIDAE_RESOURCE_READONLY_AUDIT_SUCCESS",
       json.dumps({"species":239,"all_expand":206,"lycaenidae":25,
         "lycaenidae_expand":22,"pseudo_host_taxonomy_adequate":out["Pseudozizeeria maha"]["original_host_taxonomy_lower_bound_adequate"]}),flush=True)

if __name__=="__main__":main()
