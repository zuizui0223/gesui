# Denmark 2025 eBMS × wastewater: real 2021 habitat removes any apparent site-level positive Lycaenidae association

**Date: 2026-11-10. Status: exploratory POSTHOC, not a wastewater treatment or causal test.** This specific 2021 WorldCover adjustment was designed AFTER the original Germany/Austria WWTP-source positive photo association, 2023–25 Danish eBMS sample, and Denmark NEW399 aggregate 3-family-positive/5-family-negative sign reversal were known. No preregistration, out-of-time habitat test, or biological source intervention is claimed.

## Primary question

Could the mixed source-proximity association for Lycaenidae be explained by ordinary landscape habitat composition instead of a direct wastewater nutrient subsidy? Does the positive THREE-family pooled composition statistic persist at the **individual site level within shared geographic regions** once independently observed **built/cropland/grass/tree** cover is considered?

## Completely inspectable source material and QC

- **Butterflies**: licensed CC-BY-NC 4.0 Denmark eBMS original public v1.4 DwC-A, 2025 research observations, accepted species-to-family taxonomic mapping fixed independently of WWTP effects. Raw original archive ZIP SHA256 `2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b`, original **15,606 Event rows / 21,617 Occurrence rows**. Original SOURCE-FROZEN NEW399 sites exclude ALL older 107 analyzed coordinates. After source-only reporting gate, **286** sites have at least two dates with any original published Occurrence (near registered FW WWTP<=5km **70**, remote>=10km **92**, intermediate 5–10km **124**). Published ANY-Occurrence dates are a reporting opportunity proxy, not confirmed complete monitoring or a count of non-detections.
- **Source WWTP**: 2022 EEA official DK2022 annual registered both-positive N/P source, unique reporting plant with one freshwater outfall. **203 qualified sources nationally**; near sampled sites associated with only about 17 non-overlapping source-proximity components in the original NEW399 source frame; no measured downstream butterfly-receptor exposure or 2025 treatment change.
- **NEW independently sampled landscape information**: original **ESA WorldCover 2021 v200 10 m categorical raster**, real pixel classes sampled in **500m circles** around every 286 pre-selected public eBMS monitoring point. The source code reuses exactly the existing archived 3° WorldCover tile/COG class extraction helpers. Raster pixel counts, boundary tiles, missing-image QC, 2021 static class metadata and 286 site-hash joins were retained. **286/286 source sites passed**: no tile failure, >=500 circle pixels and >=200 known land pixels; classes and land-only denominator were examined. 2021 WorldCover classes are **NOT species-level larval host plants**, field grass biomass, actual 2025 conditions, nor WWTP-released inorganic/organic nutrients.

### Quantified source-site habitat imbalance (500m buffer, 2021; proportion of terrestrial known pixels)

| Original 2021 surface cover | WWTP-near 70 public sites | remote 92 public sites | near minus remote (percentage points) |
|---|---:|---:|---:|
| Built-up | **9.58%** | 3.73% | +5.85 |
| Cropland | 13.50% | **32.38%** | −18.87 |
| Grass | 33.95% | 32.34% | +1.61 |
| Trees | **42.23%** | 30.66% | +11.56 |
| Grass + trees | **76.18%** | 63.00% | +13.17 |
| Pixel fraction coded surface water (as % of all circle pixels) | 7.36% | 7.78% | −0.42 |

These are geographic distributions, not randomized treatment balance. In the six individually comparable original country×1° regions, local built, cropland and forest differences are heterogeneous and their directions vary; a single country-wide habitat gradient cannot represent all six.

## Within-geography and habitat-adjusted actual 2025 licensed butterfly family results

Only **six original 1° regions** contain both WWTP-near and remote reported source groups, yielding **139 original 2025 site-coordinate rows** (near=61, remote=78). The source component sample, taxon identities, report-day counts and region memberships were never adjusted to obtain a preferred sign. The three-family per-site composition analysis excludes six sites without any focal three-family individualCount (n=133). Five-family and absolute family rates retain 139 sites.

**Model (site-equal, exploratory):** near-vs-remote indicator + fixed six original country×1° region dummies + log(number of original any-published-Occurrence survey dates). A habitat-augmented version further includes 2021 fractions of built, cropland, grass and trees (the other land-cover classes are omitted to avoid exact compositional collinearity). Outcome either log((family individualCount+0.5)/documented dates) or log((Lycaenidae+0.5)/(the other two or four accepted families+0.5)), the latter ONLY if the relevant family total>0. The 0.5 continuity correction is an estimation choice, not invented biological abundance. An exponentiated near coefficient is a **transformed relative recorded rate or composition odds** under these models, not an effect of wastewater dose.

| Outcome, near/remote exponentiated source coefficient | Six-region + reporting effort | Add actual 2021 500m land cover | Descriptive 6-region bootstrap 95% (habitat adjusted) |
|---|---:|---:|---:|
| Lycaenidae recorded per survey date (continuity corrected) | **0.573** | **0.522** | [0.224, 0.988] |
| Lycaenidae / (Pieridae + Papilionidae), **3-family composition**, n=133 | 0.670 | **0.511** | [0.225, 0.985] |
| Lycaenidae / (Pieridae + Papilionidae + Nymphalidae + Hesperiidae), **5-family composition**, n=139 | 0.446 | **0.377** | [0.283, 0.492] |
| Pieridae recorded per documented date | 0.836 | 0.967 | [0.592, 1.510] |
| Nymphalidae recorded per documented date | 1.372 | 1.469 | [0.675, 2.359] |
| Hesperiidae recorded per documented date | 1.590 | 1.658 | [1.073, 2.387] |

Descriptive intervals are from 3,500 deterministic **country×1° region-cluster-with-replacement resamples**. There are only SIX independent broad geography clusters; the original near source IDs can be hydrologically connected and these intervals are NOT confirmatory hypothesis tests. No p-values are asserted. Site-level OLS treats original coordinates equally; it is a **different estimand** from pooled individualCount composition, exposure-weighted survey-date standardization and raw source area totals.

### Additional independent source diagnostics and sign stability

- A second pure-Python pivoted normal-equation implementation (independent of the first numpy least-squares code, requiring no new remote observation downloads) reproduced exactly **all six original outcomes × both region-only and habitat-adjusted specifications**. Original source JSON SHA values are verified.
- Omitting any one of the SIX original 1° regions leaves the habitat-adjusted near/remote ratios in the ranges: **Lycaenidae recorded-date ratio 0.355–0.643**, **three-family OR 0.374–0.644**, **five-family OR 0.333–0.407**. Each is below one. This is **descriptive leverage**, not a valid independent treatment replication.
- Nearby 5km geodesic **source components** inside the six regions total 15, but dropping one component can remove all near records from a region. Among 14 identifiable component-drop refits, habitat-adjusted three-family OR ranged **0.399–0.562**, five-family OR **0.285–0.407**. These component labels are **not independently verified hydrological networks or WWTP treatments**.
- The original pooled 2025 NEW399 three-family near/remote OR **2.270** was positive because Pieridae were particularly concentrated at remote reporting locations. The original same-1° equal-region OR was **0.908**, original date/area-balanced 3-family OR **1.819**, and the new **site-equal regional + actual habitat OR 0.511**. These different source/pool/weight/continuity estimands CAN legitimately have different signs without contradiction. The central finding is **sensitivity to sample, denominator and spatial/observer weighting**, NOT a stable sewage-mediated effect.

## Ecological interpretation and stop rule

The genuinely measured 2021 vegetation/land-cover opportunity differs strongly between sites close to registered outfalls and more distant public eBMS sites. Nevertheless, adding built/crop/grass/tree to a site-equal model **does not restore** the apparent pooled three-family positive; in this sample it shifts the source-near coefficient further negative. Thus 2021 coarse land cover alone cannot be declared the single cause of the positive pooled composition phenomenon, and posthoc covariate adjustment cannot distinguish habitat confounding from other selection processes. Source proximity remains *nonrandom* with respect to site and habitat.

Do **NOT** claim: phosphorus from sewage suppressed or benefited Lycaenidae, WorldCover tree fraction measures larval host plants, adult mud-puddling/salt/protein subsidies have been established, or a significant European macrolepidopteran wastewater law has emerged. The source registered outfall is from 2022, butterflies observed in 2025, no true 2025 effluent recipient dose or 2021/2025 host species at each site is measured, and only six shared region clusters support the comparison.

**Next meaningful test should be explicitly independent source-to-receptor exposure**: hydrologic connectivity, measured receiving-water nutrients/flow and actual host/microhabitat in standardized butterfly monitoring near >=10 independent wastewater sources and physically appropriate source-zero controls. Do NOT repeat the existing small Denmark sample under more posthoc denominators, cutoff tuning, or opportunistic photo abundance surrogates. Keep source and causal units distinct.

## Exact GitHub source receipts and reproduction

- [Original NEW399 source-only site/recording opportunity](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2025-new399-source-coverage-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json).
- [Real 2021 WorldCover 500m land-class pixel counts and 286/286 source QA](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2025-worldcover-habitat-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DK2025_NEW399_286SITE_WORLD_COVER2021_500M_NO_FAMILY_SOURCE.json) — [original 10m COG source retrieval SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42774511484).
- [Exact within-six-region 2025 recorded family models, 3500 source-region bootstraps and original source SHA](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2025-worldcover-habitat-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DK2025_NEW399_WORLD_COVER2021_LYCAENIDAE_3VS5_WITHIN6_REGION_HABITAT_QC.json) — [original frozen-source model SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42775422795).
- [Independent original site-cohort/family/land-class OLS regression coefficients and leave-area/source-component-out check](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2025-worldcover-habitat-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DK2025_NEW399_WORLD_COVER2021_FAMILY_REGION_INDEPENDENT_OLS_QA.json) — [second arithmetic audit SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42775758937).

**Study and original eBMS raw data remain on Draft only; CC-BY-NC 4.0 attribution; raw 10m rasters, private observer IDs or licensed individual observation rows not redistributed. No merge to main or unrelated chocho results.**