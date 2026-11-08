#!/usr/bin/env python3
"""Independent offline EEA2022 matched-river/GBIF/Ellenberg NUTRI audit.

Never use this posthoc QC diagnostic as a preregistered nutrient-effect test.
Do not distribute entire third-party European EIV numerical species table.
"""
from __future__ import annotations
import argparse,hashlib,json,math,statistics
from pathlib import Path

SOURCE_P="54798e82f91585bdca26af9b86e7bd0d862deb1515a141d36628c1c273a66249"
GBIF_SHA="90ecb8d2a10936103c86698f8f863aae16603724c10fd1332f6a4e6b6b111bb9"
TRAIT_SHA="9382d876ef4919294e1913f1eac16fa6e983a5450dd34080a88b0880e2f1f3a8"
MAPPING_SHA="9c052ab23382fee4697e81630f01cd04e8ad4ad2aae2fb4eb2495b95f0024d80"
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads(path.read_text())
def exact_sign(pos,neg):
    n=pos+neg
    return min(1.,2*sum(math.comb(n,i) for i in range(min(pos,neg)+1))/2**n) if n else None
def same(a,b):
    assert math.isclose(a,b,rel_tol=0,abs_tol=1e-9),(a,b)

def audit(src,raw,trait,taxonomy,output):
    for p,h in ((src,SOURCE_P),(raw,GBIF_SHA),(trait,TRAIT_SHA),(taxonomy,MAPPING_SHA)):
        if digest(p)!=h:raise ValueError("Original source artifact changed: "+str(p))
    source=load(src);plants=load(raw);result=load(trait);mapped=load(taxonomy)
    assert len(source["pairs"])==18 and len(plants["plants"])==36
    assert result["fixed_source_pairs_before_GBIF"]==42377235104
    assert result["original_new_GBIF_run"]==42377454765
    assert result["source_EIV_qc"]["original_sha256"]==EIV_SHA
    mapping={int(x["gbif_speciesKey"]):x for x in mapped["species_taxonomy_and_matched_trait"]}
    assert len(mapping)==1006
    assert sum(x["exact_accepted_vascular_species"] for x in mapping.values())==911
    assert sum(x["EIV_name_exact_match"] for x in mapping.values())==624
    for obj in mapping.values():
        assert (obj["trait_NUTRI"] is not None)==bool(obj["EIV_name_exact_match"])
    sourceby={x["station_id"]:x for x in plants["plants"]}
    scored={x["station_id"]:x for x in result["site_trait_coverage_and_means"]}
    assert len(sourceby)==len(scored)==36
    for station,site in sourceby.items():
        scores=[(mapping[int(k)]["trait_NUTRI"],int(n))
           for k,n in site["specieskey_counts"].items()
           if mapping[int(k)]["EIV_name_exact_match"]]
        o=scored[station];num=len(scores)
        assert num==o["matched_numeric_EIV_N_speciesKeys"]
        fraction=num/site["distinct_species_keys"] if site["distinct_species_keys"] else None
        if fraction is not None:same(fraction,o["trait_coverage_fraction_all_Plantae_speciesKeys"])
        valid=(num>=10 and fraction is not None and fraction>=.4)
        assert valid==o["primary_EIV_site_quality"]
        if scores:
            same(sum(v for v,n in scores)/num,o["unweighted_mean_trait_N"])
            same(sum(v*n for v,n in scores)/sum(n for v,n in scores),o["GBIF_occurrence_weighted_mean_trait_N_sensitivity"])
    prior={x["pair_id"]:x for x in result["paired"]}
    rows=[]
    for p in source["pairs"]:
        name=p["pair_id"]
        a=scored[p["source_exposed"]["station_id"]]
        b=scored[p["source_control"]["station_id"]]
        old=prior[name]
        good=a["primary_EIV_site_quality"] and b["primary_EIV_site_quality"]
        assert good==old["eligible_primary"]
        d=a["unweighted_mean_trait_N"]-b["unweighted_mean_trait_N"] if good else None
        if good:same(d,old["primary_trait_source_positive_minus_zero"])
        rows.append({
          "pair_id":name,"block":p["spatial_block"],
          "source_EIV_rated_species":a["matched_numeric_EIV_N_speciesKeys"],
          "control_EIV_rated_species":b["matched_numeric_EIV_N_speciesKeys"],
          "source_EIV_coverage":a["trait_coverage_fraction_all_Plantae_speciesKeys"],
          "control_EIV_coverage":b["trait_coverage_fraction_all_Plantae_speciesKeys"],
          "source_GBIF_plant_records":a["all_GBIF_speciesKey_resolved_occurrences"],
          "control_GBIF_plant_records":b["all_GBIF_speciesKey_resolved_occurrences"],
          "source_quality_pass":a["primary_EIV_site_quality"],
          "control_quality_pass":b["primary_EIV_site_quality"],
          "pair_quality_pass":good,"EIV_mean_difference_if_primary_qualified":d
        })
    qual=[x for x in rows if x["pair_quality_pass"]]
    vals=[x["EIV_mean_difference_if_primary_qualified"] for x in qual]
    assert len(vals)==9
    same(statistics.mean(vals),result["primary_outcome"]["trait_source_positive_minus_zero_mean"])
    lower=sum(x["source_EIV_coverage"]<x["control_EIV_coverage"] for x in rows)
    higher=sum(x["source_EIV_coverage"]>x["control_EIV_coverage"] for x in rows)
    more_records=sum(x["source_GBIF_plant_records"]>x["control_GBIF_plant_records"] for x in rows)
    less_records=sum(x["source_GBIF_plant_records"]<x["control_GBIF_plant_records"] for x in rows)
    pos=sum(v>0 for v in vals);neg=sum(v<0 for v in vals)
    report={
      "status":"PASS_SOURCE_HASH_AND_INDEPENDENT_EIV_36_SITE_RECOMPUTATION",
      "source_original_GBIF_run":42377454765,"source_original_EIV_run":42378520227,
      "frozen_source_sha256":SOURCE_P,"GBIF_raw_sha256":GBIF_SHA,
      "EIV_result_SHA256":TRAIT_SHA,"speciesKey_trait_crosswalk_SHA256":MAPPING_SHA,
      "original_EIV_csv_sha256":EIV_SHA,
      "source_EIV_numeric_data_republication":False,
      "distinct_GBIF_speciesKeys":1006,"accepted_vascular_speciesKeys":911,
      "exact_EIV_N_rated_speciesKeys":624,
      "plant_sites":36,"matched_pairs":18,
      "source_positive_quality_sites":sum(x["source_quality_pass"] for x in rows),
      "source_zero_quality_sites":sum(x["control_quality_pass"] for x in rows),
      "predeclared_quality_pairs":len(qual),"predeclared_minimum_pairs":10,
      "predeclared_effect_inference_supported":False,
      "original_NUTRI_difference_exploratory":{
         "n":len(vals),"mean_source_minus_control":statistics.mean(vals),
         "median_source_minus_control":statistics.median(vals),
         "positive":pos,"negative":neg,
         "posthoc_exact_sign_p_descriptive_only":exact_sign(pos,neg)},
      "posthoc_trait_coverage_quality_diagnostic":{
         "source_coverage_lower_pairs":lower,"source_coverage_higher_pairs":higher,
         "pairwise_median_source_minus_control_fraction":statistics.median(
            x["source_EIV_coverage"]-x["control_EIV_coverage"] for x in rows),
         "posthoc_exact_sign_p":exact_sign(higher,lower)},
      "posthoc_GBIF_recorder_effort_diagnostic":{
         "source_more_record_pairs":more_records,"source_fewer_record_pairs":less_records,
         "posthoc_exact_sign_p":exact_sign(more_records,less_records)},
      "quality_pairwise_source_details":rows,
      "limitations":["Posthoc QC significance is not a prespecified ecology result",
         "Trait EIV NUTRI is species ecological preference, not river P or measured plant exposure",
         "GBIF indexed occurrences are not abundance or standardized effort",
         "No 2022 sewage-P causal attribution, no imputation of unknown botanical species",
         "Original n=9 does not meet required 10 qualifying-pair inference gate"]
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_INDEPENDENT_EIV_SOURCE_QC",json.dumps({k:v for k,v in report.items()
      if k not in ("quality_pairwise_source_details","limitations")},sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("source","raw","trait","taxonomy","out"):
        p.add_argument("--"+k,required=True,type=Path)
    a=p.parse_args()
    audit(a.source,a.raw,a.trait,a.taxonomy,a.out)
