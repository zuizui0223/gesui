# UK2018 permit-to-sampling identity audit implementation plan

> **For agentic workers:** Use superpowers:executing-plans to implement this continuation. The user has instructed continuing the existing source-first project; no new design approval or merge is requested.

**Goal:** Count source-reported historical permit/sample relationships among the unchanged 930 spatial candidates, preserving every unresolved edge.

**Architecture:** A pure Python helper handles exact identifiers, multiplicity and date screens. A source runner reads frozen official bytes and current WQE metadata privately, reuses the established spatial selector, and emits only aggregate JSON. An independent final review checks semantics and disclosure boundaries.

**Tech Stack:** Existing Python, openpyxl for reading only, numpy/pyproj and existing spatial selector; mdbtools for private Access extraction; GitHub connector for the existing isolated Draft branch.

**Spec:** `docs/PROTOCOL_20261121_UK2018_PERMIT_SAMPLE_ID_SOURCE.md`.

## Global constraints

- No TP/BIO outcome retrieval or regression; no annual mass calculation.
- Exact 2018 rows; complete regionalized permit/version/outlet/effluent keys.
- Preserve the 930 station-set hash and existing complete pairing hash.
- Publish aggregates and hashes only; raw downloads and database exports remain outside the repository.
- Current records and annual summaries do not prove sample-date-specific historical identities.
- Use the existing PR137 Draft branch with an expected-head lease; main remains untouched.

## Review focus

- Duplicate and repeated-header rows must not create extra sites or suggest independent P-permit groups.
- Same local sample code in different regions, or changed Midlands prefix, must not create a guessed join.
- Distinct versions/outlets must remain separate; no first-row selection.
- Missing/invalid/revoked/future dates must not qualify as verified 2018 activity.
- Partial paging, source changes or exceptions must not publish raw source fields or a false completed census.

## Task 1: Identity and date semantics

**Files:** Create `analysis/20261121_uk_prtr/permit_sample_identity_support.py` and `test_permit_sample_identity_support.py`.

**Interfaces:** `permit_key(row, source)` returns a regionalized 5-tuple or None; `ofwat_pairs(rows)` returns unique source-reported tuples plus aggregate QA; `date_overlap_status(row, year)` returns a controlled status; `summarize_candidates(candidates, pairs, corroborated, active)` returns aggregate multiplicity and corroboration counts.

- [x] Write synthetic tests for duplicate collapse, missing keys, region collisions, MI/MD source-documented correspondence, preserved permit leading zeros, distinct version/outlet retention, date-boundary failures, candidate membership and disclosure-safe aggregate output.
- [x] Run RED; implement only the required behavior; run GREEN.

## Task 2: Reproducible source census

**Files:** Create `analysis/20261121_uk_prtr/ofwat2018_permit_sample_identity_audit.py`; add one narrowly scoped Actions workflow.

**Consumes:** Task 1 helpers; existing UWWTD/WQE spatial selector.

**Produces:** aggregate source provenance, national workbook coverage, fixed-930 permit support and exact EA corroboration. Code accepts private verified caches for this run; default execution downloads only from the declared official sources or their UK National Archives preservation.

- [x] Require frozen XLSX, Access and EEA table hashes; parse selected columns only; record headers and page receipts.
- [x] Reconstruct the original 930 candidates and verify both station and complete-pair fingerprints before joining.
- [x] Check workbook and EA all/active/Effluents key consistency, duplicates, dates and source type. Keep incompleteness visible.
- [x] Run the new suite and relevant existing 25 source regressions; commit code/protocol before the study-cohort join.
- [x] Execute original-source audit, inspect status and aggregate margins, and retain the verified aggregate JSON.

## Task 3: Review, interpretation and delivery

**Files:** New results and provenance under `results/uk2018_permit_identity_audit_20261121/`; new Japanese result report in `docs/`; update PR137 summary.

- [x] Obtain a fresh code/scientific review; address important findings with failing-then-passing regressions.
- [x] Reconcile national and cohort counts independently. Explain annual-summary selection, cross-sheet duplicates, historical/current differences, and remaining UWWTD native-ID gap.
- [x] Prepare the corrected aggregate, provenance, Japanese report and prior-report addendum for the same Draft branch. Final publication verification is recorded in the PR delivery summary.

## Execution ledger

- Start base: `6bba31654520b4db1ee622b57eda6ad8a755c549`.
- Existing local directory is a selected-file snapshot, not a git checkout. Remote Git tree/commit/update-ref operations preserve unrelated files and the existing Draft branch. No local worktree is needed.
- Source discovery found a stronger public historical Ofwat workbook bridge than the initially considered current-register prefix construction. The explicit workbook IDs are therefore primary; the EA database provides cross-publication corroboration with shared EA WIMS lineage.
- Workbook schema/national counts were inspected before this protocol. Fixed-930 identity join has not yet been run at protocol creation; this is exploratory source research, not a preregistered confirmatory analysis.

- Task 1 complete: 17 new identity/source-boundary tests observed failing before implementation and passing after; all 42 tests in this source-analysis directory pass (17 new +25 existing). Python compilation passes.
- Original code run `44907070860` passed both original cohort fingerprints. Its identity counts are superseded: independent XML parsing found three legitimate alphabetic outlet identifiers that the first numeric-key validator discarded. Preserve that run and its correction history.
- Fix: retain opaque nonempty version/outlet/effluent strings, normalizing integer representations only. Two synthetic cases failed then passed. Add separate lost-pair counts (failing test before implementation) and a direct consent/Effluents conflict regression for the same joining path.
- Review: historical and current publications share EA WIMS lineage; revise independence language and explicitly bound English-region coverage. No declined-to-judge items.
- Reproduction: the local 930 station set and all 23 WQE page hashes exactly match the original run, but the hash containing projected floating-point coordinates differs across math runtimes. Preserve the local STOP, retain the strict original hash gate, and use the successful baseline-compatible Actions run for the complete mapping check. Record an additional native-pair hash and runtime versions; do not substitute it for the existing gate.

- Correction verification: 46/46 tests pass (21 identity/source-boundary +25 existing). The source-native key fix restores the three workbook identities found independently; pair-loss and conflict-path regressions pass.

- Adopted source run: `44907990658`, code `eb119f7530d8d69ae9d4184df7ef36b49420274a`, SUCCESS. All 46 regressions pass in Actions. Both original cohort fingerprints pass; all 23 WQE page hashes match run `44904429288`. Artifact `14708502882` ZIP digest and the extracted JSON digest were separately verified.
- Final source counts: 607 workbook pairs /606 points; 522 pairs /521 points corroborated across EA publications; 33 pairs /33 points also in the active table. National all-only support is 489 pairs /488 points, keeping multiplicity explicit.
- Fixed-cohort partition independently reconciled: 16 active corroborations +279 all-only corroborations +51 Ofwat-only associations +584 without a complete Ofwat association =930. All346 source-covered candidate points have one annual Ofwat key;295 are cross-publication corroborated.
- Independent date-boundary check: the earlier auxiliary 6,189 eligible keys included two revoked exactly on2018-01-01. The frozen half-open rule correctly retains6,187; no duplicate or key-normalization discrepancy remains. The separate427 invalid intervals comprise422 empty and5 reversed intervals. Active-table32 vs33 is exactly one restored alphanumeric outlet identifier.
- Japanese report and provenance document the annual-summary limitation, shared WIMS lineage, limited English-region source coverage, unresolved EEA physical outlet identity, missing annual volume and unchanged scientific gates. No TP/BIO model or daily-flow calculation was added.
- Delivery gate: publish only the explicit report/provenance/aggregate/addendum/plan files, then compare local Git blob SHA values with the final remote tree, inspect integrity checks, and recheck Draft status, base branch and unchanged main before reporting completion.
