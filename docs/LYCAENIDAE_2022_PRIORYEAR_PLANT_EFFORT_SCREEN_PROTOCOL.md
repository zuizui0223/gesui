# 2022 river WWTP × Lycaenidae: independently effort-screened exploration v1

**Protocol frozen prior to any 2022 butterfly photo outcome queries in this NEW cohort.** New cohort is distinct from the initial 24-pair 1 km test that failed the predeclared ≥12 adequately observed pairs requirement. That failure remains an immutable result.

## Question and causal boundary

Can a landscape-level, source-associated butterfly *photographic reporting share* be estimated in habitats where independent, **prior-year non-butterfly** observational effort is adequate? This does NOT test adult puddling, larval performance, sodium/protein enrichment, local N load, or wastewater causality. The independent effort proxy could itself correlate with sewage, urban access, plant visibility, or observer interests; filtering on it changes the estimand and can bias cross-environment comparisons.

## Completely source-only phase A (NO butterfly photo data)

- Original EEA Biology2026 2022 **source-only** benthic river monitoring sites in [frozen PR #14 registry](https://github.com/zuizui0223/gesui/blob/625af95fc12e01e9c796230ca4709be72225ab50/results/eea2022_biology_bqe_20261108/SOURCE_ONLY_BQE_SITE_REGISTRY.csv) (Git blob SHA1 `16f6b531e29ed09f070f9fb86f7d8feb899108ed`).
- Only original hydrologically **upstream ≤50 km** reported 2022 WWTP **both-positive measured** P and N masses vs zero **qualifying registered** sources; no original numerical EQR. No claim that registry source-zero is untreated/no sewage.
- Match exposed vs registered-source-zero stations **same country + 3° regional block**, using log(1+riverQ), log(1+upstream basin area), paired geographic distance [10,120] km, non-overlapping **5 km radii** (≥10 km all selected centers), no repeated station IDs. This is not a verified same river reach, and station proximity may include other hydrological networks.
- Deterministic SHA-based source-only exposed order, at most **16 pairs per country, 96 pairs in total**; no butterfly outcomes or observation data used.
- Save full eligible source support, every chosen site, controls, distance and QC SHA before any observational effort query.

## Prior-year *non-butterfly* observation opportunity phase B (NO 2022 butterfly results)

- Query official iNaturalist public v1 `/observations` at **each fixed site** within **5 km**, observation year **2021**, taxon **Plantae** only, research-grade, noncaptive, georeferenced, public geoprivacy, photo required, stated coordinate uncertainty ≤500m.
- This measures **2021 public plant photo reporting effort**, not number of plants, plant cover or nutrient availability. Plantae are deliberately separate from focal butterflies and their 2022 outcomes.
- Every original source-only candidate pair is evaluated without replacing failed/zero-plant buffers. Keep both sides if **≥8 plant-photo records each**, and `min(n1,n0)/max(n1,n0)≥0.20`, and the ratio is defined (both positive). Pre-specified threshold. Any missing API response means pair excluded from target analysis, not an ecological absence.
- To avoid posthoc outcome selection, sort surviving pairs by **minimum plant effort across two sites descending**, then by fixed source-only pair ID; at most **32 pairs**, and cap retained pairs to **8 per country**. This prioritizes comparable effort rather than Lycaenidae signal and changes generalization scope.
- Freeze selected pairs and source-only plant-effort query hashes in GitHub **before** any 2022 Lycaenidae/Pieridae/Papilionidae count query. If fewer than **12 retained pairs across at least 3 countries**, STOP before butterfly analysis, report effort feasibility only. The 12-pair/3-country gate is separate from the earlier 24-pair 1 km gate; cannot replace its conclusion.

## Same-year butterfly photo stage C (only after fully frozen A+B)

- At every retained fixed site, same spatial radius **5 km** and 2022, official iNaturalist descendant family IDs: Lycaenidae 47923, Pieridae 48508, Papilionidae 47223. Use identical photo, geoprivacy, research-grade, noncaptive, and ≤500m accuracy filters.
- Predefine derived response as **Lycaenidae / (Lycaenidae + Pieridae + Papilionidae) public photo reporting share** conditional on targeted family photos, NOT real butterfly relative abundance. Require ≥5 focal-three-family photos at **both** paired sites.
- A descriptive share result may be shown **only if ≥12 support-qualified pairs and ≥3 countries** remain. Otherwise STOP without inferential effect calculation. Save all zero/NA records and compare source-positive vs registered-source-zero total photo effort even when gate fails.
- If gate passes, report signed paired share differences, country-balanced and station-weighted averages, country-stratified directions, and country-cluster bootstrap intervals as **exploratory** only, with physically downstream registered-source placebo and within-country urban land-cover/host plant covariates **required before causal attribution**. No p-hacking / optimizing threshold / radius after butterfly results.
- A positive photo association would be an association of reporting, not puddling, not larval fitness, not a wastewater nutrient effect. The preexisting `chocho` 22/25 Lycaenidae *potential resource geography* remains unchanged and separate.

## Specificity and failure cases

Original 1 km 24 pairs had only **1/24 eligible** with ≥5 photos both sides and 32/36 photos in Luxembourg. This is not overwritten by a new larger-radius, non-butterfly-effort-screened design. In the new 5 km study, failures on 2021 plant effort, API access, 2022 butterfly photo support, country replication or independent urban/plant covariate availability must be reported equally.

Source: [gesui PR #20](https://github.com/zuizui0223/gesui/pull/20) original 1 km feasibility, and [chocho frozen resources](https://github.com/zuizui0223/chocho/tree/e8b8f86c7a889b6a3224a4ef69ffe5e0c7738660). Keep this branch Draft, no merge to `main` or changes in `chocho`.
