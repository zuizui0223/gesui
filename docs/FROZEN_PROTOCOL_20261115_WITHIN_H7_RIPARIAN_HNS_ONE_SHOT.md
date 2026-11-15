# FROZEN v1.0 — One-shot within-H7 riparian HNS test

**Freeze date:** 2026-11-15 JST
**Status:** PROTOCOL ONLY. Biological outcome **NOT FIT**; verified existing sources fail source/observation opportunity requirements.
**Repository:** existing unmerged PR #95 analysis branch. No new branch or PR. No edits to main or chocho.

## Study question and epistemic status

Do source-positive municipal WWTP reaches show a larger upstream-to-downstream change in the near-river versus inland gradient of **Lycaenidae host N affinity** than matched registered-source-negative reaches from the same H7?

This is a one-shot, hypothesis-directed protocol for a **new point-level analysis**, developed after prior exploratory H7 signals were seen. It is **not** a retrospective preregistration and, if 2022 observations overlap previous aggregated GBIF occurrence records, **not an independent cohort confirmation**. Source-coordinate and observation coverage gates must be frozen without opening the NEW local HNS outcomes. Existing H7 aggregate data cannot create exact river-distance observations.

An H7 fixed effect removes uniform between-basin variables. It does not remove the correlation between near-river habitat/urban access and W across basins; simply fitting W×river-nearness inside basins is not a wastewater-specific test.

## One primary ecological outcome

The sole primary Y is the **species-incidence-equal mean larval-host nutrient-affinity HNS** of accepted Lycaenidae species with numerically scored larval-host HNS. Use the previously frozen 303-key roster / **47 numerically scored Lycaenidae**. The exact original score source SHA256 is 00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb. Each taxon contributes once per site cell regardless of photo frequency; missing score does not mean zero HNS. No substitute all-butterfly, Pieridae or record-weighted HNS can replace the primary.

The date window is **2022-05-01 to 2022-09-30 inclusive**, synchronized only nominally to authentic positive **2022 annual reported** municipal WWTP N/P sources and actual discharge coordinates. The annual report is not measured receiving-water chemistry or confirmed larval-host/soil nutrient delivery. Only genuinely dated species-keyed, point-georeferenced GBIF records are eligible; raw 454-row H7-period taxon counts are NOT point observations.

## Geography — four cells per reach, source versus paired reference

Register source reaches using qualifying actual reported 2022 freshwater outfalls, directed NEXT_DOWN topology, original source coords, real receiving stream and exact native H7 ID. Exclude unresolved same-reach source/station ordering and disconnected river links. Pick each matching no-qualifying-reported-source pseudo-outfall reach **within the same H7** using river order, discharge, elevation, land cover and river access, with no butterfly outcome inspection. Reference means **no matched qualifying registered source**, not pristine or no sewage.

Before opening any butterfly HNS, fix reference-source matching rules: river order difference at most 1, flow ratio 0.5–2, altitude difference within 150m, WorldCover2021 urban and agricultural fractions within 10 percentage points, and all eight cell areas free from overlap with a second selected source/reference unit. If a valid source-reference match does not exist, exclude the pair with source-only reason, never change tolerance after looking at Y.

For each source/reference reach, use:
- **Upstream or downstream:** 1–5 true **river-network** km from actual outlet/pseudo-outlet; not Euclidean distance, and no ambiguous same-reach ordering.
- **Near river:** terrestrial point fully within 0–200m of that reach centerline.
- **Inland:** fully 500–1500m from any mapped river centerline, with equivalent relevant non-river land-cover opportunity. Discard the transition band 200–500m.
- **Geolocation:** original coordinateUncertaintyInMeters <=100m; its entire uncertainty disk must lie within the chosen distance band. Remove centroids, duplicate IDs, impossible coordinates, records without valid source date or accepted species match. Small unmapped tributaries remain a limitation.

Thus source plus reference, each with upstream/downstream × near/inland, gives **8 local cells per independently matched reach pair**. Never count 8 cells as 8 independent outfalls.

## Exactly one fixed inferential contrast

For reach j define:
D(j) = (HNS_down,near - HNS_down,inland) - (HNS_up,near - HNS_up,inland).

For each same-H7 matched positive-source/reference reach pair k:
T(k) = D(k, source-positive) - D(k, registered-source-negative).

The **only primary estimate** is equal-pair mean tau = mean_k[T(k)] on raw 1–9 HNS points. Estimate a **two-sided 95% CI** from a full-refit **genuine H6-cluster bootstrap** (2,000 replicates, seed 20261115), resampling the full matched pairs of all H7 inside the sampled H6 cluster. Report tau, CI, two-sided alpha=0.05 result, independent H6/source counts and all 8-cell raw observation denominators. Do not run a family-by-family search, multiple widths, regions, seasons, alternative distance cutoffs or ad hoc models and promote a favorable one. Normal multiplicity caveat remains because the current question was suggested by previous positive macro research.

**Interpretation:** tau>0 could be consistent with a directional source-associated river-near composition contrast. It remains observational: the reference may contain unreported sewage; city access, host-plant occurrence, land-cover heterogeneity, climate, recorder selection, other effluent and actual dilution can differ. It cannot be described as a measured causal wastewater→plant→butterfly indirect effect.

## Hard source-only stop gate BEFORE any new local HNS means are calculated

Require **at least 24 nonoverlapping matched source/reference pairs**, in **at least 12 genuine H6 parent basins and 3 countries**, with source/reference matching and actual upstream/downstream topology verified. This is a minimum opportunity threshold and does not automatically guarantee 80% power.

In **every** one of the eight local cells of a selected pair: >=20 unique source-georeferenced five-family butterfly speciesKey occurrence records; >=2 accepted Lycaenidae with frozen numeric HNS; >=70% Lycaenidae records with valid score; >=5 distinct observer-dates/standardized sampling events. If sampling-event/recorder/date metadata are not sufficiently available to measure opportunity, the cell fails instead of being imputed as sampled. Enforce non-overlap, season, geolocation, provenance and event/license eligibility at the source-only gate; choose candidate pairs without reading local numeric HNS.

No pooled visual reinspection or redefining 200m, 500–1500m, 2022 summer, 47 species, effort gates, H6 minimum or reference matching after a failure. If any gate fails: **STOP; no HNS outcome regression or interpretation**. The source-only feasibility check is not the one inferential biological run. Freeze eligibility manifest, retained IDs, raw source receipts and source QA SHA before seeing the new local HNS.

**Known pre-existing gate failures** (must not be ignored):
- DK standardized eBMS 2025: 10 strictly downstream sites but only **2 distinct WWTP outfalls/reaches**, upstream 4 sites; PR #42.
- DE 2022 registered-source river pairing: 73 river pairs/100 monitoring points, but only **5 original photo-effort-supported 5km matched quartets**, below its predeclared 8-gate, PR #36.
- EU 454 H7×period butterfly HNS uses aggregated speciesKey facet counts and does NOT retain occurrence-level river-distance coordinates, PR #95.

These sources do not presently establish 24 independently matched source-reference pairs with eight qualifying point-level cells. Accordingly **the biological primary test has NOT been executed** and no local wastewater-specific organism effect, including null, can be asserted.

## Fixed result interpretation if a future support-qualified cohort exists

Hypothesis sign is tau>0, but report two-sided estimate and CI. Predefined smallest biologically meaningful positive contrast for this **index** is +0.25 raw HNS points on its 1–9 scale; this is an explicitly provisional ecological judgment, not inferred from a new outcome. If CI crosses zero: **unresolved**, regardless of p and estimate. If upper 95% CI lies below +0.25: at most rule out a >=0.25 **HNS composition** contrast for these sampled sites, year and geometry, conditional on biases and model, NOT all terrestrial transfer.

A null HNS test cannot falsify adult puddling on mineral/organic/carrion resources, leaf N transfer or larval feeding, aquatic adult emergence, predator diets, small river-edge subsidies or biological response in other years/radii. Do not write "sewage never reaches land."

## Linked evidence (existing draft branches, not new publications)

- Exploratory macro W→HNS and source-direction caveats: https://github.com/zuizui0223/gesui/pull/88 , https://github.com/zuizui0223/gesui/pull/93 , https://github.com/zuizui0223/gesui/pull/95
- Denmark eBMS original source-only STOP: https://github.com/zuizui0223/gesui/pull/42
- Germany original source-only STOP: https://github.com/zuizui0223/gesui/pull/36
- Official GBIF sampling-event requirement: https://www.gbif.org/data-quality-requirements-sampling-events
- Prior aquatic insect deposition at stream edge is not species HNS: https://pubmed.ncbi.nlm.nih.gov/35523331/

No eBird. No TinyFish. No new inferential Actions dispatch from this frozen document. Keep Draft; no merge to main or chocho.
