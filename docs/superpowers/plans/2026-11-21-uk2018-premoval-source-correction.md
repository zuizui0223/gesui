# UK2018 PRemoval source and model correction plan

> **For agentic workers:** Execute the bounded correction continuously, using the existing frozen source protocol and independent code review. The user's continued instruction authorizes updates to Draft PR #137; main and other PRs stay outside this change.

**Goal:** Reconcile the 930 geographical candidates with 360 eligible TP observations, restore the omitted historical outlet coordinates, and report median and log-mean contrasts as different estimands.

**Architecture:** Retain the original source retrieval, pairing, 2018 cohort and thresholds. Add an aggregate-only support audit; repair the existing model and source receipts. Preserve the previous artifact byte-for-byte and publish the corrected output separately.

**Tech Stack:** Python 3.11, numpy, pandas, statsmodels, unittest, existing GitHub Actions.

**Spec:** `docs/PROTOCOL_20261121_UK_UWWTD2018_PREM0_REMOVAL_TO_EA_TP_SOURCE.md` and the documented pending 930-denominator audit in the existing result.

## Global constraints

- Original 2018 source, 250 m mutual uniqueness, one original outlet, 930 candidate pairs; no widening or outcome-based source selection.
- No source raw rows, station IDs, coordinates, or exact individual TP in public results.
- No ecological outcome, causal weights, imputed zero for missing TP, annual mass conversion, or causal treatment claim.
- Preserve original source thresholds and report their limitations; extra overlap diagnostics do not redefine them.
- Three-character UK region prefix is NUTS1; retain the full original uwwNUTS separately for source QA.
- Primary descriptive median uses all eligible numeric-TP sites; model complete cases are separately reported.

## Review focus

1. Omitted coordinates must be caught by a synthetic model with known flag coefficient and spatial confounding.
2. Missing covariates must not silently alter the primary descriptive population.
3. No-returned-record, censored and invalid-result cases must be distinguishable and source-page receipts retained.
4. Candidate/eligible counts must reconcile across source support tables without reading TP outcome values.
5. Unexpected source errors must leave an auditable JSON and fail execution, rather than look successful.

## Task 1: Preserve evidence and reproduce defects

- [x] Retrieve run 44850406458 artifact 14685161737; verify ZIP SHA256 e16750780254e4661b38dcebf0151a78743626d63d320972474fc9c63524e3a2.
- [x] Preserve original JSON SHA256 8c0ff6b0ceecebd7660644a5ff5fbdaf795b978a0aa3bbba5cc1aa7be10a49ff.
- [x] Add and run failing synthetic unittest cases for missing coordinates and complete-case drift before changing production code.

## Task 2: Correct the existing source comparison

Files: existing `eea2018_PRemoval_to_EA_TP_effluent_source_test.py`; new `test_premoval_source_model.py` and independent `premoval_source_support.py` / `test_premoval_source_support.py` in the same directory.

- [x] Carry historical outlet BNG coordinates through the frame and full adjusted model. Preserve old no-coordinate estimate as a labelled diagnostic.
- [x] Name coarse geography NUTS1 and validate its UK prefix. Keep full original code for QA.
- [x] Retain source-page receipts and all candidate fates; distinguish censored values, malformed values, source duplicates and repeated same-day observations. No substitutions for censoring.
- [x] Add source-only flag × region × fixed capacity-band support tables. Capacity bins are a new descriptive choice fixed before their results, not pre-registered confirmatory thresholds.
- [x] Retain EEA and EA registry provenance. Return nonzero after unexpected exceptions while writing the STOP JSON.
- [x] Run local targeted tests and syntax checks, then independent review.

## Task 3: Verify against the original sources and report

- [ ] Update the existing workflow to run the targeted tests before source retrieval; push only to existing Draft PR #137 with an expected-head lease.
- [ ] Inspect workflow output and final JSON status, source hashes/cohort counts and corrected model estimates; download and verify artifact digest.
- [ ] Save aggregate evidence and update result/protocol errata and PR summary. Preserve original median contrast, same-scale unadjusted/adjusted table and unresolved identity/selection limitations.

## Execution decisions

- The existing frozen protocol supplies the intended model; restoring East/North is a bug fix, not a new search for a favourable effect.
- The previous log-mean model without coordinates remains as a historical diagnostic. Neither it nor the new full model is directly comparable to a ratio of raw medians.
- Current API calls continue requesting 250 rows per page. This conservative requested size is supported; no claim about the current maximum across all formats is needed.
- A local source snapshot is used because the repository is accessed through the GitHub connector. Remote updates preserve all untouched tree entries and use a compare-and-swap head check.

## Verification ledger

- Native model regression tests: observed failures for omitted coordinates, complete-case descriptive drift, discarded record-fate audit and successful exit on source failure; corrected and green.
- Independent source-support helper: 17 tests, observed RED before implementation and GREEN after.
- Fresh independent review found three important gaps: source-cohort count alone, serialization outside failure guard, and adjustment gates suppressing primary descriptives. Each received a failing regression and correction.
- Review detail: export-duplicate counting applies after SA/material/unit filters; rejected-row counts remain returned records, not claims of independent samples. This scope is explicit in output metadata.
- Baseline station-set SHA is enforced before requesting TP. Full pair-mapping SHA is newly recorded; the earlier run has no pair hash, so retrospective pair identity is not claimed.
