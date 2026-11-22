# Butterfly within-genus recording test implementation plan

> **For agentic workers:** Use superpowers:executing-plans for native implementation with source acquisition in parallel and one independent final review.

**Goal:** Run the source-verified, frozen new within-butterfly-genus test and report its actual result.

**Architecture:** Restore immutable raw sources; construct discordant same-genus pair strata; fit an exact within-stratum conditional likelihood with composite-score inference; retain positive, negative and failed checks.

**Tech Stack:** Python, NumPy, pandas, SciPy, unittest; existing GitHub repository.

**Spec:** docs/superpowers/specs/2026-11-22-butterfly-within-genus.md

## Global constraints

Same183 H7/40 H6 and original species; two periods clustered together; no imputed HNS or fabricated absence; no new outcome-driven thresholds; all raw-source hashes checked; no manuscript change or main merge.

## Review focus

- Plant genus and butterfly genus must never be interchangeable.
- Repeated GBIF records must not change the primary presence response.
- Conditional likelihood must omit arbitrary stratum intercepts and require both discordant outcomes.
- Cluster scores use weights once before squaring; shared pair observations are not independent replicates.
- Trait permutations keep species identities consistent across space/time, including initially tied species.

## Task1: Source and pair assembly

- [x] Commit specification before fitting a new biological effect.
- [x] Add behavioral tests for duplicate-record invariance, same-genus pairing, missing scores, stratum conditioning and pinned-source replay; watch relevant new behavior fail.
- [x] Implement `analysis/20261122_butterfly_within_genus/analyse.py` source builder and source-only gate. Input is the archived cache, output source eligibility summary; never fit in source-only mode.
- [x] Verify original4 community indicators within1e-10 and publish actual informative counts.

## Task2: Conditional likelihood and inference

- [x] Add hand-derived small-stratum likelihood and gradient tests, orientation-invariance and permutation-consistency tests; run RED before numerical implementation.
- [x] Implement dynamic-programming conditional log likelihood/gradient, deterministic rank handling, optimizer and H6 cluster sandwich.
- [x] Validate against independent enumeration/finite differences and a known synthetic effect before source fitting.
- [x] Run the frozen primary,399 H6 bootstrap,999 trait-label reference, genus influence and two fixed sensitivities. Log every failure; do not replace a failed model.

## Task3: Report and review

- [x] Save source hashes, source eligibility, aggregate results and Japanese report; no raw third-party bulk rehosting.
- [x] Run repo unittest suite and focused numerical checks; compare summary arithmetic to saved outputs.
- [x] Have a fresh reviewer audit the complete code/spec/results, resolve material issues, and preserve any numerical changes as implementation corrections.
- [ ] Publish this isolated branch and create a draft PR: blocked by automatic approval review pending explicit authorization for external GitHub publication. All analysis and review work is complete locally; no alternative publication route was attempted.

## Completion record

- Specification commit: `89c349cd3c1f63c328d24b01aa43408d3a6449b3`; implementation commit before the effect fit: `1071fd0c703cb371c3f4200087ce0a8d1d54371d`.
- Independent code/numerical/report review passed, with one minor attribution wording correction applied.
- Main result and all contrary checks are reported in `docs/RESULTS_20261122_BUTTERFLY_WITHIN_GENUS.md`.
- 23/23 tests passed; all 1,427 frozen resampling/influence fits matched independent enumeration to within 9.1e-7.
