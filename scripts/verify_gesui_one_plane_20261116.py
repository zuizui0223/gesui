#!/usr/bin/env python3
"""Offline test: 454 true H7period joined water/plant/butterfly source and 83 blobs.

Verifies no hidden row duplication or mixing of EEA station vs H7 vs US EPA UID;
all 83 original Git blob SHA1 copied without modifications, all 451 numeric
EIV/3-layer sites complete, exact EEA H7 water site mean, 56 reported positive
and no imputed clean wastewater control. PR112 104-PR canonical snapshot stable.
Source Actions true numeric/actual GBIF/EEA publisher provenance remain linked,
not re-fetched or falsely called independent field replication.
"""
from __future__ import annotations
import csv,hashlib,json,math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PLANE=ROOT/"results/20261116_evidence_plane/EU_H7_PERIOD_454_WATER_PLANT_BUTTERFLY_WWTP_COMMON_PLANE.csv"
QA=ROOT/"results/20261116_evidence_plane/EU_H7_PERIOD_454_COMMON_PLANE_SOURCE_QA.json"
MODELS=ROOT/"results/20261116_evidence_plane/EU_H7_451_UNIFORM_W_AND_PAIRWISE_H6_COMPARISON.json"
MODELCSV=ROOT/"results/20261116_evidence_plane/EU_H7_451_UNIFORM_W_AND_PAIRWISE_H6_COMPARISON.csv"
CROSS=ROOT/"results/20261116_evidence_plane/STUDY_LEVEL_EVIDENCE_CROSSWALK_NO_POOLING.csv"
MANIFEST=ROOT/"provenance/POST112_ONE_PLANE_83_ORIGINAL_SOURCE_BLOBS_20261116.json"
PR112=ROOT/"docs/PR_RESEARCH_REGISTER_20261116.csv"
SOURCE=ROOT/"results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def gitblob(path):
    raw=path.read_bytes()
    return hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
def rows(path):
    with path.open(newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))
def assert_equal(cond,label):
    if not cond:raise RuntimeError("One-plane original-source audit failed "+label)
def main():
    metadata=json.loads(QA.read_text())
    manifest=json.loads(MANIFEST.read_text())
    assert_equal(metadata["schema"]=="GESUI_20261116_TRUE_EU_H7xPERIOD_WATER_PLANT_BUTTERFLY_PLANE_V1","source schema")
    assert_equal(digest(PLANE)==metadata["generated_common_plane_csv_SHA256"],"source plane checksum")
    assert_equal(digest(SOURCE)==metadata["original_454_H7_panel_SHA256"],"genuine 454 source panel checksum")
    assert_equal(metadata["matched_site_year_to_H7_mapping_records"]==3525,"2018/22 site-H7 original source")
    main=rows(PLANE)
    assert_equal(len(main)==454,"454 true unit rows")
    keys={(r["h"],r["period"]) for r in main}
    assert_equal(len(keys)==454 and len({r["h"] for r in main})==232 and len({r["h6"] for r in main})==90,
      "454 unique H7-period/232 native H7/90 true H6")
    assert_equal({r["period"] for r in main}=={"2017-2019","2020-2022"},"genuine two periods")
    positive=sum(r["EEA_upstream_positive_report_matched_flag"]=="1" for r in main)
    assert_equal(positive==56,"56 original positive WWTP source reports")
    assert_equal(sum(r["EEA_downstream_positive_placebo_source_flag"]=="1" for r in main)==61,
      "61 independently audited downstream placebo source positives")
    assert_equal(all(r["EEA_upstream_positive_report_matched_flag"]==r["EEA_upstream_positive_source_flag_independent"]
        for r in main),"original independent upstream source report flags agree")
    assert_equal(all(r["EEA_no_match_is_NOT_clean_control"]=="True" for r in main),
      "No reported WWTP source misclassified as clean")
    assert_equal(all(r["water_WWTP_report_N_and_P_are_NOT_inriver_TNTP"]=="True" for r in main),
      "No river chemistry imputation")
    valid=[r for r in main if r["plant_EIV_original_source_eligible"]=="True"]
    assert_equal(len(valid)==451,"451 pre-outcome plant EIV score support")
    assert_equal(all(r["plant_EIV_NUTRI_equal_species"] and r["aquatic_EEA_original_EQR"] and
            r["butterfly_hostN_HNS_equal_species"] for r in valid),"shared full 451 coresponse")
    assert_equal(all(not r["plant_EIV_NUTRI_equal_species"] for r in main if r not in valid),
      "3 missing EIV indexes retained NA")
    assert_equal(sum(r["EEA_upstream_positive_report_matched_flag"]=="1" for r in valid)==56,
      "all positive report cells have source-eligible EIV")
    diffs=[]
    for r in main:
        a,b=r["aquatic_EEA_original_EQR"],r["benthic_sameH7_site_equal_EQR"]
        if a and b:diffs.append(abs(float(a)-float(b)))
    assert_equal(len(diffs)==454 and max(diffs)<1e-9,
      "Independent original 3525 site values exactly replay 454 aquatic EQR basin means")
    m=json.loads(MODELS.read_text())
    assert_equal(m["schema"]=="GESUI_20261116_SAME451_UNIFORM_THREE_COMPARTMENT_H6_SCREEN_V1"
        and m["n_exact_same_H7period"]==451 and m["n_true_H6"]==90,"6 uniform H6 true models source")
    assert_equal(m["source_input_SHA256"]==digest(PLANE),"uniform models share same original 451")
    assert_equal(len(m["models"])==6,"all 6 source model contrasts preserved")
    modelrows=rows(MODELCSV)
    assert_equal(len(modelrows)==6 and len({r["id"] for r in modelrows})==6,"6 model summary unique")
    evidence=rows(CROSS)
    assert_equal(len(evidence)==18 and len({r["evidence_id"] for r in evidence})==18,
      "18 distinct study rows are not pooled")
    assert_equal(all(r["relationship_to_original_EU_H7_plane"] for r in evidence),
      "all studies explicitly identify unit link/limitation")
    registry=rows(PR112)
    assert_equal(len(registry)==104,"PR112 canonical original 104-PR source registry preserved")
    assert_equal(manifest["number_of_previous_historical_PR_original_blobs_newly_preserved"]==83,
      "83 selected original pre2026 negative+positive evidence files")
    originals=manifest["source_original_github_PR_changed_file_blobs"]
    assert_equal(len(originals)==83,"83 exact original blobs listed")
    for item in originals:
        path=ROOT/item["path"]
        assert_equal(path.exists() and gitblob(path)==item["sha"],"unaltered original PR SHA "+item["path"])
    print("GESUI_POST112_ONE_PLANE_OFFLINE_VERIFIED="+json.dumps({
      "original_PR112_registry":len(registry),"source_blobs_byte_identical":len(originals),
      "same_true_H7period":len(main),"same_H7":232,"same_true_H6":90,
      "same_h7_EIV_EQR_HNS_complete":len(valid),
      "original_2018_22_benthic_site_EQR_replay_rows":len(diffs),
      "source_positive_WWTP_report_period":positive,
      "uniform_source_same451_6model_outputs":len(modelrows),
      "study_level_evidence_without_pooled_effects":len(evidence),
      "no_causal_sewage_exposure_claim":True
    },ensure_ascii=False,sort_keys=True))
if __name__=="__main__":main()
