#!/usr/bin/env python3
"""Verify 2018 EEA upstream WWTP reported-P source support in FR263→178→122.

Never equate 'no positive report' with an unpolluted control. This is an
identifiability/source-coverage test only; no ecological effect is fitted.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

SHAS = {
 "rz":"ff60c807be99468ae95e0c8cd6fabbb60a70a5e55aa3bf51f12dc5f67cdb9f9b",
 "fr263":"1824b0d9b84d35a460b4d1d0513659f93a16c59d7f5666a2d8707be1c6da9942",
 "fr178":"9a5d11a6c50ba49160da46097873c3d5eedda881e0a173f1aa90968c94950b11",
 "model":"d2a66c80527187ea6fcd7dbb6ec1cad8464740644ba17f4fbdd44bb6c1d44b2e"
}
SHARD_HASHES=[
 "d2a3e058f5922ec8abc877d308f7deeaa3ebeeaa272de611263d1d627da6830c",
 "be2da865e0b8d404e9cc197464c7aa1b37bd22850e0e8f8197cb2a0e0b40535d",
 "ac45cd0a1a01987081eb29f2e4080e08bc5848cb39f0ffaa37d8d4ca8a02d126",
 "c16d38fe55e35e6ad3ca50b2c3cb76f74530c5e5372d7e941eca265c01e92d9c"
]
SOURCE_SHA="5a074f5182711c4ab74ef9bdacfa3918e9fcb434dd0b75c67072a0f06459c3e1"

def checked_bytes(p,checksum):
    b=p.read_bytes()
    if hashlib.sha256(b).hexdigest()!=checksum:
        raise RuntimeError("Original source SHA mismatch "+str(p))
    return b
def j(p,sha):return json.loads(checked_bytes(p,sha))
def c(p,sha):return list(csv.DictReader(io.StringIO(checked_bytes(p,sha).decode("utf-8-sig"))))

def run(a):
    rz=j(a.rz,SHAS["rz"])
    if (rz["schema"]!="GESUI_2018_ORIGINAL_EEA_UWWTD_BENTHIC_RZ2018_250M_500M_SOURCE_ONLY_V1"
        or rz["status"]!="PASS_2018_RZ2018_ALL_ORIGINAL_BENTHIC_SITES_SOURCE_ONLY_NO_ECOLOGY"
        or rz["source_original_Github_2018_fixed_data_sha256"]!=SOURCE_SHA):
        raise RuntimeError("Original 2018 EEA UWWTD source lineage invalid")
    pop=rz["original_source_pilot_true_station_classification"]
    if (len(pop),sum(bool(q["source_2018_upstream_P_positive_report"]) for q in pop))!=(2639,204):
        raise RuntimeError("Original EEA2018 2639/204 population drift")
    bykey={q["site_key"]:q for q in pop}
    if len(bykey)!=2639:raise RuntimeError("Duplicate original 2018 source station")
    fr263=c(a.fr263,SHAS["fr263"]);fr178=c(a.fr178,SHAS["fr178"])
    model=j(a.model,SHAS["model"])
    if (model["schema"]!="GESUI_FR178_2018_2022_TWO_SCENE_MOSAIC_REAL_TP_H6_EXPLORATORY_V1"
        or model["status"]!="POSTHOC_EXPLORATORY_MOSAIC_RIVER_TP_NDVI_FITTED_NOT_ORIGINAL_MAIN_MODEL"):
        raise RuntimeError("The original genuine 122 numeric artifact changed")
    keys={
      "FR263":[q["EEA_identical_native_site_key_2018_2022"] for q in fr263],
      "FR178":[q["site_key"] for q in fr178],
      "FR122":[q["site_key"] for q in model["original_native_site_month_matched_values"]]}
    for k,n in (("FR263",263),("FR178",178),("FR122",122)):
        ids=keys[k]
        if len(ids)!=n or len(set(ids))!=n or not set(ids).issubset(bykey):
            raise RuntimeError("Source key coverage not exact at "+k)
    if not set(keys["FR122"]).issubset(keys["FR178"]) or not set(keys["FR178"]).issubset(keys["FR263"]):
        raise RuntimeError("The three source cohorts are not nested by true native EEA ID")
    positives={k:sorted(x for x in ids if bykey[x]["source_2018_upstream_P_positive_report"])
               for k,ids in keys.items()}
    if [len(positives[k]) for k in ("FR263","FR178","FR122")]!=[1,1,0]:
        raise RuntimeError("Report-positive counts differ from independently checked 1/1/0")
    only=positives["FR263"][0]
    if positives["FR178"]!=[only]:
        raise RuntimeError("One original reported-positive site lost before S2 opportunity cohort")
    optical={}
    for i,sha in enumerate(SHARD_HASHES):
        b=j(a.shards/f"FR178_2SCENE_MOSAIC_SOURCE_SHARD{i}.json",sha)
        if b["schema"]!="GESUI_FR178_S2_FROZEN_2SCENE_SAME_MONTH_CLOUD_MOSAIC_SOURCE_QA_POSTHOC_V1" or b["original_source_shard"]!=i:
            raise RuntimeError("Original no-aquatic-outcome source shard mutated")
        for x in b["observations"]:
            if x["site_key"] in optical:raise RuntimeError("Optical duplicate source ID")
            optical[x["site_key"]]=x
    if set(optical)!=set(keys["FR178"]):
        raise RuntimeError("Source-only 178 vs optical 178 ID mismatch")
    if {k for k,v in optical.items() if v["source_two_years_real_mosaic_pixel_eligible"]}!=set(keys["FR122"]):
        raise RuntimeError("Source-only NDVI filtering does not exactly match numeric 122")
    affected=optical[only]
    if affected["source_2018_2022_real_mosaic_months_n"]!=1 or affected["source_two_years_real_mosaic_pixel_eligible"]:
        raise RuntimeError("Unique 2018 reported-positive station optical source qualification differs")
    missing=[x for x in affected["original_monthwise_source_QA"] if int(x["year"])==2022 and int(x["month"]) in (6,8)]
    missed=sorted((int(x["month"]),int(x["band_counts"]["20_100m"])) for x in missing)
    if missed!=[(6,79),(8,78)]:
        raise RuntimeError("Observed source missingness changed for the single eligible reported-P river")
    rows=[{"stage":k,"n_exact_river_site":len(keys[k]),
           "n_original_2018_upstream_P_report_positive":len(positives[k]),
           "n_no_qualifying_positive_report":len(keys[k])-len(positives[k])}
          for k in ("FR263","FR178","FR122")]
    result={
      "schema":"GESUI_FR263_FR178_FR122_SOURCE_ONLY_2018_UWWTD_POSITIVE_P_IDENTIFICATION_V1",
      "status":"PASS_SOURCE_AUDIT_0_OF_122_ELIGIBLE_REPORT_POSITIVE",
      "source_2018_EEA_all_benthic_sites":2639,"source_2018_positive_upstream_P_report_sites":204,
      "stages":rows,"the_one_original_report_positive_site":only,
      "original_positive_report_site_dropped_after_2022_satellite_quality":True,
      "pixel_count_near_original_riverline_20_100m_2022_Jun_Aug":{"June":79,"August":78,"min_each":100},
      "only_one_valid_matching_summer_month_for_unique_report_positive_station":True,
      "reported_positive_false_does_NOT_imply_true_zero_sewage":True,
      "2018_upstream_annual_report_is_NOT_2022_local_delivery_or_contaminant_flux":True,
      "no_source_specific_exposure_contrast_in_FR122":True,
      "no_TP_NDVI_OR_causal_effect_estimated_by_this_source_audit":True,
      "original_input_byte_SHA256":SHAS,
      "original_source_4_mosaic_shard_sha256":SHARD_HASHES}
    a.out.mkdir(exist_ok=True,parents=True)
    (a.out/"FR122_SOURCE2018_UPSTREAM_UWWTD_P_POSITIVE_AUDIT.json").write_text(
      json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    with (a.out/"FR263_178_122_2018_UPSTREAM_UWWTD_REPORT_COVERAGE.csv").open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    print("FR122_WWTP_EXPOSURE_IDENTIFICATION="+json.dumps({
      "status":result["status"],"n_by_cohort":{r["stage"]:r["n_original_2018_upstream_P_report_positive"] for r in rows},
      "positive_site_excluded_Jun_Aug_pixel_counts":[79,78]},ensure_ascii=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    for arg in ("rz","fr263","fr178","model","shards","out"):
        p.add_argument("--"+arg,type=Path,required=True)
    run(p.parse_args())
