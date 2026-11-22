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

- [ ] Commit specification before fitting a new biological effect.
- [ ] Add behavioral tests for duplicate-record invariance, same-genus pairing, missing scores, stratum conditioning and pinned-source replay; watch relevant new behavior fail.
- [ ] Implement `analysis/20261122_butterfly_within_genus/analyse.py` source builder and source-only gate. Input is the archived cache, output source eligibility summary; never fit in source-only mode.
- [ ] Verify original4 community indicators within1e-10 and publish actual informative counts.

## Task2: Conditional likelihood and inference

- [ ] Add hand-derived small-stratum likelihood and gradient tests, orientation-invariance and permutation-consistency tests; run RED before numerical implementation.
- [ ] Implement dynamic-programming conditional log likelihood/gradient, deterministic rank handling, optimizer and H6 cluster sandwich.
- [ ] Validate against independent enumeration/finite differences and a known synthetic effect before source fitting.
- [ ] Run the frozen primary,399 H6 bootstrap,999 trait-label reference, genus influence and two fixed sensitivities. Log every failure; do not replace a failed model.

## Task3: Report and review

- [ ] Save source hashes, source eligibility, aggregate results and Japanese report; no raw third-party bulk rehosting.
- [ ] Run repo unittest suite and focused numerical checks; compare summary arithmetic to saved outputs.
- [ ] Have a fresh reviewer audit the complete code/spec/results, resolve material issues, and preserve any numerical changes as implementation corrections.
- [ ] Commit and push this isolated branch; create a draft PR without merge, then report verified results and limitations to the user.
