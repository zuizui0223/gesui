# UK2018 permit outlet location implementation plan

> **For agentic workers:** Use superpowers:executing-plans for this authorized continuation; independent agents investigate original-source semantics and one fresh reviewer checks the implementation. Continue without another approval prompt for the existing source-only Draft scope.

**Goal:** Audit version-linked outlet locations for the unchanged295 corroborated points without promoting spatial agreement into formal EEA identity.

**Architecture:** A small pure helper validates National Grid cells and computes conservative distance classes. A separate runner reuses the accepted identity and cohort gates, then reads only permitted location fields from the same Access rows. It emits safe aggregateJSON through a narrowly scoped Actions workflow.

**Tech Stack:** Existing Python3.11, openpyxl, mdbtools, numpy/pyproj3.7.2 and accepted source readers; GitHub Git tree/commit/ref connector for the existing Draft branch.

**Spec:** `docs/PROTOCOL_20261121_UK2018_PERMIT_OUTLET_LOCATION_SOURCE.md`.

## Global constraints

- Keep original930 candidates and exact station/full-pair hashes; require the prior346/295 support partition.
- Use original5keys and half-open2018 dates; preserve opaque alphanumeric identifiers.
- Only `OUTLET_GRID_REF` supplies an outlet location. Other location fields never fill its gaps.
- Preserve grid resolution, missingness and conflicting versioned outlet geometry. No raw identifiers, names, coordinates or observations in public output.
- No TP/BIO/flow retrieval, models, mass estimates, river rerouting, candidate additions or main/PR124 changes.

## Review focus

- Geographic letter tile conversion must exclude I and reject out-of-domain tiles rather than wrap them.
- Coarse grids straddling a distance threshold must stay uncertain; use whole-cell min/max.
- Duplicate/cross-effluent records must not silently supply an arbitrary outlet position.
- Current metadata dates must not be presented as proof that a coordinate was surveyed in2018.
- Source/cohort failure and partial execution must yield a safeSTOP, not partial successful counts.

## Task 1: Grid and location semantics

**Files:** `analysis/20261121_uk_prtr/permit_outlet_location_support.py`, `test_permit_outlet_location_support.py`.

**Interfaces:** `grid_cell(value)` returns `(easting,northing,width)` or raises a controlledValueError; `distance_bounds(point,cell)` returns min/max; `distance_class(point,cell,radius)` returns within/outside/boundary. `location_index(rows,eligible_keys)` returns per-key field statuses and safeaggregateQA; `summarize_locations(pairs,index,matching,registry)` returns aggregate location and distance counts.

- [x] Write and observe failing synthetic tests: knownliteral tile/cell values, spaced/lowercase strings, bad letters/odd digits/off-grid,100m-cell threshold boundary, missing/duplicate/conflicting per-key and cross-effluent outlet geometry, correct outlet-vs-site selection and aggregate-only summaries.
- [x] Implement only those behaviors and run the new tests plus all46 accepted source tests.

## Task 2: Original-source execution

**Files:** `analysis/20261121_uk_prtr/ea2018_permit_outlet_location_audit.py`; `.github/workflows/gesui-UK2018-permit-outlet-location-20261121.yml`.

**Consumes:** accepted identity runner and Task1 helpers; frozen sourcebytes and candidate mapping.

**Produces:** safe aggregateJSON with schema, source receipts, exact baseline checks, version-linked outlet availability and fixed100/250/500m cell classifications.

- [x] Reuse the accepted source download/hash/read/join functions and preserve all failures. Validate the baseline support before reading new cohort results.
- [x] Add one synthetic end-to-end aggregation/controlledSTOP test for the runner's pure audit boundary, exercising real joins and source field selection.
- [x] Commit protocol and tested code before the295-point coordinate census. Execute in the baseline-compatible Actions environment and inspect status, artifact digest and arithmetic margins.

## Task 3: Review and deliver

**Files:** Japanese results in `docs/`; source aggregate/provenance in `results/uk2018_permit_outlet_location_20261121/`; additive prior-report note and PR137 summary.

- [ ] Obtain one fresh code/scientific review and independently reconcile adopted aggregate counts; fix important findings with RED/GREEN tests.
- [ ] Report public-source search coverage, coordinate semantics and remaining formalID/temporal gaps; do not reinterpret geographic disagreement as a measured ecological effect.
- [ ] Save explicit code/results/docs on the same Draft branch with expected-head lease; verify remote file hashes and CI and unchangedmain.

## Execution ledger

- Starting head `c864a7dad459797954746234fa99a51e58f87804`, tree `be8481a68019d16f3fbdbe2c50a81eb7c80de6ce`; PR137 remains openDraft stacked onPR136. Main `6a0ebac1274f23cd202e2f515176b6d77d264b54`.
- The local directory remains a selected-file snapshot, not a gitcheckout. Use remote tree operations and expected-head lease; no worktree/main mutation is needed.
- Ruling: keep the existing source-only workflow and inline implementation, with independent source investigators and one final fresh review. The user's continuation request covers this original-data audit; a new approval round would not resolve a missing scientific choice.
- Ruling: compare whole National Grid cells, not inferred metre-precise points. This conservatively preserves ambiguous boundaries but cannot account for undocumented historical coordinate errors.
- Pre-flight: Task2 consumes Task1 grid/index/summarization functions and the accepted identity functions; Task3 consumes the completed aggregate. No conflicting outputs or shared-file writes between investigators.
- Private prior/tmp downloads were pruned by workspace maintenance. Reacquire originalofficialbytes and verify existing hashes; never substitute a newer unverified source.
- Task1 complete: all16 geometry/location tests were observed failing against unimplemented interfaces; the implemented helper then passed all62 tests (16 new plus46 accepted).
- Task2 boundary verification: all4 audit/STOP tests were observed failing before implementation, then all66 source tests passed. Synthetic fixtures exercise the unchanged930/346/295/16/279/584 baseline and missing outlet geometry without replacing it from site coordinates.
- Read-only Git detection confirms this directory is not a repository checkout. The existing remote isolated Draft branch and expected-head lease are retained; no local worktree is created.
- Scientific note before cohort results: nationwide source exploration finds many coordinates ending in zeros and multiple distinct outlet IDs at the same recorded grid. Report notation resolution only; neither metre accuracy nor identity follows from matching coordinates.
- Initial protocol/code commit `762bfaa43e8133e8bc0fbbbd39d88e51316cf70d`; run `44913481948` SUCCESS with66 tests. ZIP SHA256 `fb2c6a39f88d168699828466a0a197130d95cad20c3560ab53e06b8343809e5a` matches GitHub artifact metadata. JSON SHA256 `4bc372eb38c6854bc8cfc82ba48c01f04f4b708f621fdf8e33b8d58ea1dee87a`. This is a pre-review-guard result, not the adopted final run.
- One fresh whole-implementation review found3 Important gaps: WQE coordinate-bearing pages were not pinned,295 pair identity was enforced by counts alone, and malformed/incomplete location exports could pass. All were reproduced through8 added regressions (initially2 failures and8 error/subtest outcomes). Fixes are applied in one review pass with each issue's tests checked after its change.
- Independent reconstruction from the accepted implementation and frozen raw inputs recovered295 station/five-key pairs with SHA256 `1ce2b5be966498d282a5f3ae23bff896fa27c9b63712d67baba3ea7881f308ae`. This additional fingerprint is newly reconstructed, not a claim about an old recorded field. All23 WQE raw-page hashes, original930 membership/native EEA mapping and every previous support count match. Local fullfloat projection differs; production keeps the strict accepted gate.
- Independent geometry implementation (no new geometry-helper reuse) matches all six threshold tables from the initial CI run. Actual counts will be adopted only after the strengthened run passes.
- Two Minor items deferred with bounded current risk: (1) a row missing part of its five-key is skipped by cross-effluent conflict logic; the SHA-locked actual source has zero invalid five-key rows, so this cannot affect this census; (2) direct output writes could leave stale/truncated output in an interrupted repeated local destination, whereas this workflow uses a fresh run-specific temporary path. No declined-to-judge review areas. No further reviewer dispatch is required; verify the fixes and source aggregate directly.
- All74 source tests pass after the three guard fixes. The parent reran the separate reconstruction and reproduced the exact295-pair fingerprint and all accepted baseline counts. Next verification is the strengthened production run on the unchanged strict spatial runtime gate.
