#!/usr/bin/env python3
"""Offline byte preservation + provenance/estimate-boundary smoke checks.

This is an integration integrity audit, NOT a new ecological re-analysis.
Original GBIF, EEA, EPA, HydroBASINS and Tichy 2023 complete raw source
artifacts are kept in referenced successful GitHub Actions and are not
implicitly available as part of the main repository.
"""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/"provenance/GESUI_VERIFIED_20261116_RESEARCH_BLOB_MANIFEST.json"
ORIGINAL_H7_SHA256="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
ORIGINAL_HNS303_SHA256="00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb"
SPECIES_CONTRIB_GIT_SHA1="02b0ade550c7cb2a7689757519a7fcf5c1f3703a"

def j(path):
    return json.loads((ROOT/path).read_text(encoding="utf-8"))
def close(actual,expected,tol=1e-8):
    if abs(float(actual)-expected)>tol:
        raise AssertionError("Expected "+repr(expected)+", observed "+repr(actual))
def within_zero(ci):
    if not (float(ci[0])<=0<=float(ci[1])):
        raise AssertionError("Negative/control interval incorrectly recast as positive: "+repr(ci))
def positive(ci):
    if not float(ci[0])>0:
        raise AssertionError("Positive result was not actually zero-excluding: "+repr(ci))

def verify():
    manifest=j("provenance/GESUI_VERIFIED_20261116_RESEARCH_BLOB_MANIFEST.json")
    if manifest["schema"]!="GESUI_MAIN_INTEGRATION_VERIFIED_20261116_GIT_BLOB_MANIFEST_V1":
        raise AssertionError("Original source-bounded main preservation schema changed")
    files=manifest["files"]
    if manifest["file_count"]!=98 or len(files)!=98:
        raise AssertionError("Exactly the 98 prespecified source blob copies required")
    if len({x["path"] for x in files})!=98 or sum(x["size_bytes"] for x in files)!=manifest["total_bytes"]:
        raise AssertionError("Duplicate path, unexpected size or source manifest drift")
    for x in files:
        path=ROOT/x["path"]
        data=path.read_bytes()
        actual=hashlib.sha1(f"blob {len(data)}\0".encode()+data).hexdigest()
        if actual!=x["git_blob_sha1"] or len(data)!=x["size_bytes"]:
            raise AssertionError("A source branch Git blob changed: "+x["path"])
        if not x["original_source_ref"].startswith("analysis/gesui-"):
            raise AssertionError("Unrecognized verified origin branch "+x["path"])
    h7=(ROOT/"results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")
    hns=(ROOT/"results/20261115_sync/GESUI_303_HOST_NUTRIENT_SCORES.csv")
    if hashlib.sha256(h7.read_bytes()).hexdigest()!=ORIGINAL_H7_SHA256:
        raise AssertionError("Original 454 H7-period numeric source SHA mismatch")
    if hashlib.sha256(hns.read_bytes()).hexdigest()!=ORIGINAL_HNS303_SHA256:
        raise AssertionError("Original 303 butterfly HNS taxon source SHA mismatch")
    with h7.open(newline="",encoding="utf-8-sig") as stream:
        rows=list(csv.DictReader(stream))
    if len(rows)!=454 or len({x["h"] for x in rows})!=232 or len({x["h6"] for x in rows})!=90:
        raise AssertionError("Original H7-period native 454/232/90 grain changed")
    with hns.open(newline="",encoding="utf-8-sig") as stream:
        if len(list(csv.DictReader(stream)))!=303:
            raise AssertionError("The source 303 butterfly host-score roster changed")
    q=ROOT/"results/20261115_hns_family_taxon/SPECIES_W_SLOPE_ATTRIBUTION.csv"
    b=q.read_bytes()
    if hashlib.sha1(f"blob {len(b)}\0".encode()+b).hexdigest()!=SPECIES_CONTRIB_GIT_SHA1:
        raise AssertionError("Original 186 scored butterfly source decomposition changed")
    shuffle=j("results/20261115_hns_trait_shuffle/PERMUTATION_RESULTS.json")
    if shuffle["source_species_n"]!=186 or shuffle["iterations"]!=9999:
        raise AssertionError("Original butterfly fixed within-family shuffles changed")
    close(shuffle["primary"]["p_two_sided"],0.0013)
    incidence=j("results/20261115_hns_trait_shuffle/INCIDENCE_VS_RECORD_SOURCE_VERIFIED_RESULTS.json")
    close(incidence["results"]["incidence_perm"]["permutation_p_two_sided"],0.0039)
    if not incidence["results"]["record_observed"]>0:
        raise AssertionError("Original HNS record co-composition sign was lost")
    eiv=j("results/20261115_eiv_specificity/EIV_NUTRI_FIVE_PAIRWISE_H6_COUNTRY_VERIFIED_SUMMARY.json")
    if sorted(eiv["results"]["comparison"])!=sorted(["LIGHT","TEMPE","MOIST","REACT","SALIN"]):
        raise AssertionError("Do not drop any of the five original opposing EIV controls")
    for item in eiv["results"]["comparison"].values():
        positive(item["NUTRI_after_other_H6_gain95"])
        if item["NUTRI_after_other_BH5_q"]>=0.05:
            raise AssertionError("Reported exploratory five-pair coefficient changed")
    stop=j("results/20261115_eiv_specificity/SOURCE_ONLY_SIX_PAIRWISE_WITHINH6_GATES.json")
    if stop["strict_six_trait_complete_case_gate"]["all_six_preoutcome_gate_pass"] is not False:
        raise AssertionError("Six-trait 40-H6 source STOP must not be silently converted to PASS")
    if stop["within_true_H6_geographic_replication_gate"]["within_H6_spatial_association_sample_gate_pass"] is not False:
        raise AssertionError("Original within-H6 30/40 spatial STOP must remain visible")
    fe=j("results/20261116_main_integration/183H7_40H6_H6FE_MAIN_SOURCE_RUN.json")
    if fe["source_job_conclusion"]!="success":
        raise AssertionError("Original true-H6 source execution not successful")
    close(fe["numerical_result"]["primary"]["beta"],0.11807615189599471)
    positive(fe["numerical_result"]["primary"]["beta95"])
    wn=j("results/20261116_main_integration/183H7_40H6_W_NULL_SOURCE_RUN.json")
    within_zero(wn["numerical_result"]["modeled_W_to_plant_N_within_true_H6"]["H6_cluster95"])
    within_zero(wn["numerical_result"]["modeled_W_to_butterfly_HNS_after_plantN_within_H6"]["H6_cluster95"])
    tempo=j("results/20261116_main_integration/EU219_H7_TEMPORAL_NULL_SOURCE_RUN.json")
    if tempo["numerical_result"]["n_exact_twoperiod_H7"]!=219:
        raise AssertionError("Temporal source repeated-H7 cohort changed")
    genus=j("results/20261115_h6_source_expansion/PLANT_EIV_GENUS_COMPONENTS_VERIFIED_SUMMARY.json")
    t=genus["derived_from_original_success_log"]["targets"]["hns_incidence"]["both_component_model_H6_cluster"]["source_components"]
    positive(t["between_genus"]["H6_cluster_description95"])
    within_zero(t["within_genus"]["H6_cluster_description95"])
    null=j("results/20261115_h6_source_expansion/PLANT_EIV_GENUS_GLOBAL_LABEL_SHUFFLE_VERIFIED_SUMMARY.json")["derived_from_original_success_log"]
    close(null["all_species_trait_shuffle"]["two_sided_permutation_p"],0.0012)
    close(null["within_genus_trait_shuffle"]["two_sided_permutation_p"],0.6864)
    oldnew=j("results/20261115_h6_source_expansion/OLD30_NEW10_EXTERNAL_TRANSFER_VERIFIED_SUMMARY.json")["derived_from_original_success_log"]
    within_zero(oldnew["new10_external"]["hns_incidence"]["gain95"])
    if not (ROOT/"docs/RESULTS_20261115_EPT_FOURTHCORNER_TRAIT_INCIDENCE.md").exists():
        raise AssertionError("Independent US aquatic EPT physical-habitat layer not preserved")
    report=ROOT/"docs/RESULTS_20261116_VERIFIED_MAIN_INTEGRATION.md"
    text=report.read_text(encoding="utf-8")
    for literal in ("下水", "未確認", "属内", "時間", "STOP", "EPT"):
        if literal not in text:raise AssertionError("Integrated full evidence boundaries missing "+literal)
    print("GESUI_VERIFIED_MAIN_PRESERVATION_PASS: 98/98 exact Git source blobs, original H7/303 source SHAs, five EIV controls, primary positive and null/STOP results.")
    print("This integrity PASS does not validate causal sewage effects or re-run original external input source archives.")

if __name__=="__main__":
    verify()
