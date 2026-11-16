#!/usr/bin/env python3
"""Offline original-blob and PR snapshot verification of gesui aquatic integration.

No external data access, no effect refit, and no original EEA/EPA ZIP needed.
Checks exact SHA1 Git-blob identities for every selected scientific source
and schema/coverage of the GitHub-sourced full PR register.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/"provenance/GESUI_20261116_AQUATIC_PR_SOURCE_BLOB_MANIFEST.json"
REGISTER=ROOT/"docs/PR_RESEARCH_REGISTER_20261116.csv"
VERIFIED_FILES=218
REGISTERED_PRS=104
REGISTER_RANGE=(2,110)
EXPECT_LABELS={"肯定":10,"否定":17,"STOP":13,"混合":43,"":21}
EXPECT_PRIORITY={14,16,17,19,56,57,70,71,72,73,74,75,76,77,78,79,80,81,82,83,97,99,100}
REQUIRED_PATHS={
 "results/eea2022_biology_bqe_20261108/ORIGINAL_EQR_BIOLOGICAL_RESULTS.json",
 "results/eea2022_bqe_falsification_20261108/BQE_2022_ORIGINAL_MISSINGNESS_AUDIT.json",
 "results/eea2018_independent_temporal_bqe_20261108/2018_ORIGINAL_TEMPORAL_BQE_REPLICATION.json",
 "results/eea2018_independent_temporal_bqe_20261108/SITE_DISJOINT_TWO_YEAR_MODEL_TRANSPORT_POSTHOC.json",
 "results/eea2020_source_only_20261108/RESULTS_AND_PREOUTCOME_STOP.md",
 "results/eea2016_thirdyear_source_20261111/2016_EEA_RIVER_BENTHIC_EQR_JRC_GREEN_MODEL_POINTNP_NUTS2_HELDOUT_RESULTS.json",
 "docs/RESULTS_20261113_FRANCE_I2M2_MEASURED_2022_RIVER_TP_EQR_273SITES_74H7_THREE_LAYERS.md",
 "docs/RESULTS_20261113_FR_REAL_2018_2022_TP_I2M2_263_SAME_SITES_SPATIAL_TEMPORAL_DECOUPLING_3BIOTA.md",
 "docs/RESULTS_20261113_FR263_SAME_SITE_TP_2018_2022_COMMON_MONTHS_I2M2_THREE_LAYERS.md",
 "docs/RESULTS_20261114_NRSA_THREECYCLE_EPT_AMBIENTNP_PHAB_HUC8_EXTERNAL_TRANSPORT.md",
 "docs/RESULTS_20261114_NRSA_3ERA_RARE100_AND_EPT_FUNCTIONAL_TRAITS.md",
 "docs/RESULTS_20261115_NRSA1314_PIECEWISE_PATH_VS_PHAB.md",
 "scripts/gesui_model_2022_eea_aquatic_EQR.py",
 "scripts/gesui_2018_independent_bqe_EQR_outcomes.py",
 "scripts/gesui_NRSA_3cycle_EPT_NOHUC_external_transport_models.py",
 "scripts/gesui_NRSA3cycle_EPT_FFG_HABIT_PTV_EXTERNAL.py",
 "analysis/20261115_piecewise_path/observational_paths.py",
 "analysis/20261115_joint_trait_env/state_and_penalty_audit.py",
}

def gitsha(data:bytes)->str:
    return hashlib.sha1(("blob "+str(len(data))+"\0").encode()+data).hexdigest()

def validate():
    m=json.loads(MANIFEST.read_text(encoding="utf-8"))
    if m.get("schema")!="GESUI_20261116_PR_REGISTER_AQUATIC_CURATED_SOURCE_BLOB_MANIFEST_V2":
        raise AssertionError("Manifest source-of-truth version changed")
    files=m.get("files")
    if len(files)!=VERIFIED_FILES or m.get("copied_git_blob_files")!=VERIFIED_FILES:
        raise AssertionError("Not all original 218 research assets included")
    if len({r["path"] for r in files})!=len(files):
        raise AssertionError("Duplicate selected result/science path")
    base=m.get("initial_main_sha")
    if not isinstance(base,str) or len(base)!=40 or not re.fullmatch(r"[a-f0-9]{40}",base):
        raise AssertionError("Original prior main SHA must be exact")
    size_total=0
    for row in files:
        name=row["path"]
        if name.startswith("/") or ".." in Path(name).parts or row["source_branch"].strip()=="":
            raise AssertionError("Unsafe path or missing original source branch")
        path=ROOT/name
        if not path.is_file():
            raise AssertionError("Selected original result not in checkout: "+name)
        raw=path.read_bytes()
        actual=gitsha(raw)
        if actual!=row["git_blob_sha1"]:
            raise AssertionError(f"Original-source Git blob mismatch for {name}: {actual}")
        if len(raw)!=row["size_bytes"]:
            raise AssertionError("Source original byte size drift for "+name)
        size_total+=len(raw)
    if size_total!=m["copied_size_bytes"]:
        raise AssertionError("Sum of 218 selected bytes differs from original audit")
    missing=REQUIRED_PATHS-{r["path"] for r in files}
    if missing:raise AssertionError("Required EU/USA original evidence missing: "+repr(sorted(missing)))
    if any(row["path"].startswith(".github/workflows/") for row in files):
        raise AssertionError("Old unreviewed source workflows not allowed to become active on main")

    with REGISTER.open("r",encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=REGISTERED_PRS or len({int(x["pr_number"]) for x in rows})!=REGISTERED_PRS:
        raise AssertionError("Full original 104-PR factual registry missing/duplicated")
    n=[int(x["pr_number"]) for x in rows]
    if (min(n),max(n))!=REGISTER_RANGE or n!=sorted(n):
        raise AssertionError("PR snapshot number span/order changed")
    for row in rows:
        if not re.fullmatch(r"[a-f0-9]{40}",row["head_commit_sha"]):
            raise AssertionError("Invalid source real PR HEAD commit SHA "+row["pr_number"])
        if row["judgment"] not in ("肯定","否定","STOP","混合",""):
            raise AssertionError("Unreviewed or invented source judgment")
        if row["pr_state"] not in ("open","closed","merged"):
            raise AssertionError("GitHub PR state invalid")
        if not row["pr_url"].startswith("https://github.com/zuizui0223/gesui/pull/"):
            raise AssertionError("PR crosslink invalid")
        if row["selection_in_this_integration"]=="priority-aquatic-selected" and int(row["pr_number"]) not in EXPECT_PRIORITY:
            raise AssertionError("Unexpected PR included without priority source classification")
    by={int(x["pr_number"]):x for x in rows}
    if by[106]["pr_state"]!="merged":
        raise AssertionError("Already merged PR106 status lost")
    actual_priority={int(x["pr_number"]) for x in rows if x["selection_in_this_integration"]=="priority-aquatic-selected"}
    if actual_priority!=EXPECT_PRIORITY:
        raise AssertionError("Aquatic original PR selection not identical to locked 23 PR numbers")
    labels=Counter(r["judgment"] for r in rows)
    if dict(labels)!=EXPECT_LABELS:
        raise AssertionError("Study outcome-judgment validation counts drift "+repr(labels))
    scientific_readme=(ROOT/"docs/EVIDENCE_LEDGER.md").read_text(encoding="utf-8")
    for marker in ["2018","2022","HUC8","PR登録簿","欧州","STOP","Chla"]:
        if marker not in scientific_readme:raise AssertionError("Missing primary and contrary evidence ledger topic "+marker)
    return {
        "schema":"GESUI_20261116_AQUATIC_218_ORIGINAL_BLOBS_104_PR_REGISTER_VERIFIED_V1",
        "verified_files":len(files),
        "verified_total_source_bytes":size_total,
        "verified_PR_rows":len(rows),
        "judgment_classification_counts":dict(labels),
        "original_branch_HEAD_SHA_values_all_40hex":True,
        "original_2018_2022_reported_WWTP_load_EQR_sources_present":True,
        "2016_modelled_NP_not_measurement_source_present":True,
        "three_original_NRSA_cycles_and_piecewise_negative_present":True,
        "unrelated_unreviewed_Draft_workflows_not_reactivated":True,
        "test_does_not_reestimate_biological_effects":True,
    }

if __name__=="__main__":
    result=validate()
    print("GESUI_AQUATIC_PR_REGISTER_SOURCE_BLOBS_QA="+json.dumps(result,sort_keys=True,ensure_ascii=False))
