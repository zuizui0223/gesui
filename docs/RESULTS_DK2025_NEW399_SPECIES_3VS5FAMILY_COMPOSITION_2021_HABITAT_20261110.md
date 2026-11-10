# Denmark 2025 NEW399 eBMS: which butterfly species generate the wastewater-source 3-versus-5-family sign reversal?

Date 2026-11-10. **POSTHOC** species-level attribution of a known observational composition pattern. The preceding positive German/Austrian photograph share and 2025 Danish 3-family OR >1 vs 5-family OR <1 had already been observed. This is an exact arithmetic/species source composition analysis, NOT causal wastewater effects or a randomized source design.

## Original fixed sample and source QA

- Danish official eBMS v1.4 CC-BY-NC 4.0 Darwin Core Archive, original SHA256 `2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b`, **15,606 Event rows + 21,617 Occurrence rows**. Original source-only NEW399 sample excludes all 107 previously studied coordinates. Exactly **286/399** original sites pass the prespecified two+ survey dates with any published Occurrence. Near <=5km from a qualifying 2022 EEA registered freshwater WWTP: **70 sites, 539 linked dates**. Remote >=10km: **92 sites, 757 dates**. Intermediate 124 sites, 1,286 dates. Original 2022 registered WWTP source-distance and original five-family taxonomic classification are held fixed.
- Original GBIF strict species→family taxonomy was frozen by name before any butterfly counts: **65 accepted scientific-name records mapping to 63 distinct accepted usageKeys** because two pairs of original names coalesce under the same accepted species (*Aricia agestis*, *Thymelicus lineola*). Distinct accepted identifiers: Lycaenidae 20, Pieridae 7, Papilionidae 1, Nymphalidae 28, Hesperiidae 7. Every accepted taxon including those with zero records in a source group contributes to numerical reconciliation, not cherry-picked after comparing effects.
- Re-read the original licensed Event and Occurrence tables with exact SHA and core-ID linkage. All **399 original site×five-family counts** and source-linked survey dates reproduced EXACTLY, plus published source-group family totals: near [Lycaenidae 382; Pieridae 680; Papilionidae 0; Nymphalidae 3,822; Hesperiidae 372], remote [644; 2,602; 0; 3,544; 192]. No raw licensed site-species row, photograph, observer identity or sensitive coordinate is committed.

## Main new biological-composition finding: a few species dominate denominator geometry

| 2025 reported species | Near registered WWTP <=5km | Remote >=10km | Source-site distribution |
|---|---:|---:|---|
| *Pieris brassicae* (Pieridae) | **145** | **1,439** | 28 near vs 54 remote site coordinates |
| *Maniola jurtina* (Nymphalidae) | **1,828** | **870** | 61 near vs 59 remote |
| *Polyommatus icarus* (Lycaenidae) | **142** | **363** | 25 near vs 41 remote |
| *Plebejus argus* (Lycaenidae) | 63 | 0 | **all 63 at ONE near site**; potential leverage, not a source-system effect |

These are actual recorded adult individualCounts, not standardized true abundances. Original survey-date-normalized near/remote ratios for the first three species: *P. brassicae* **0.142**, *M. jurtina* **2.951**, *P. icarus* **0.549**. The source-near **Pieridae** has fewer records (680 vs remote 2,602) and its **remote** group is 55.3% *P. brassicae* (near group only 21.3%). The near-source **Nymphalidae** is 47.8% *M. jurtina*, vs remote 24.5%. Thus a pooled three-family ratio can rise when the denominator family is especially abundant elsewhere even if recorded Lycaenidae per survey day falls.

### Fully audited algebraic species deletion from BOTH groups (not biological species removal)

The original observed pooled **Lycaenidae/(Pieridae+Papilionidae) OR near/remote is 2.270**, versus **Lycaenidae/(Pieridae+Papilionidae+Nymphalidae+Hesperiidae) OR 0.771**. Papilionidae original counts happen to be zero in near/remote, so three-family denominator effectively equals Pieridae. Recompute the two source-group odds after subtracting one species from BOTH near and remote group observed counts, holding all other original counts fixed; this is exactly verifiable counterfactual arithmetic, NOT proof the species was a cause.

| Excluded species from observed near+remote composition | 3-family pooled OR | 5-family pooled OR |
|---|---:|---:|
| None: original observed | **2.270** | **0.771** |
| Remove *P. brassicae* | **1.289** | 0.614 |
| Remove *M. jurtina* | **2.270** | **1.065** |
| Remove *P. icarus* | **3.268** | **1.111** |
| Remove *P. argus* (one-site cluster) | **1.895** | 0.644 |

**Interpretation:** *Pieris brassicae* is a major driver of the pooled three-family positive (OR 2.27 becomes 1.29 without it); *Maniola jurtina* in the added Nymphalidae denominator alone suffices to change the sign of the pooled five-family statistic (0.771→1.065 after hypothetical removal). *Polyommatus icarus* also has enough remote dominance that excluding it switches five-family odds across unity; this proves the summary statistic is strongly composition-dependent, not that any species biologically drives effluent effects. Original all 63 accepted taxon IDs are included in machine-readable remove-one-species sensitivity.

## Added independent 2021 habitat check for three dominant source taxa

Retained original **six country×1° regions with BOTH source groups and 139 public 2025 locations** (61 near, 78 remote), with genuine pixel-sampled **ESA WorldCover2021 10m categories in 500m circles** for all 139. Source site family sums and species name mapping were independently checked against original 399-site outcomes; no new 2025 site recruitment or outcome-based thresholds. The three species were fixed by the above mass contributions, not selected to optimize a particular positive association. Descriptive OLS near coefficient from log((original species individualCount+0.5)/original ANY published-Occurrence survey dates), with area fixed effects and log dates; habitat model additionally built, cropland, grass and tree fractions, using source-equal public sites (not count weighting). Bootstrap resamples whole original **six 1° regions 3,000 times**.

| 2025 single species site-level near/remote transformed recorded-rate ratio | Geography/recording days | + WorldCover 2021 landcover | Six-region bootstrap descriptive 95% of habitat fit |
|---|---:|---:|---:|
| *P. brassicae* | **0.383** | **0.411** | [0.216, 0.760] |
| *M. jurtina* | **2.099** | **1.815** | [0.793, 3.641] |
| *P. icarus* | **0.788** | **0.662** | [0.473, 0.980] |

These models are explicitly posthoc and only six geographical resampling units. The result is consistent with **species- and landscape-specific recording responses**, NOT a universal WWTP butterfly effect. The original eBMS source data are not randomized and do not measure species detectability, actual 2025 river nutrient concentration, adult salt/protein puddling or real larval host abundance.

### Tentative natural-history interpretation, NOT a measured mechanism

- *P. brassicae* caterpillars use Brassica-type crops, gardens and wild crucifers; the independent 2021 WorldCover data find remote survey sites have more mapped **cropland** (32.38%) than near-source survey sites (13.50%). This supports a **plausible agricultural/habitat distribution alternative**, but it does **not** measure field Brassica host plants, pesticides or migration; importantly, the posthoc multivariable within-region fitted crop coefficient is NOT independently positive after other land-cover covariates, so **do not claim measured Brassica-mediated causation**. Supporting natural-history reference: [UK RHS cabbage butterfly host evidence](https://www.rhs.org.uk/biodiversity/cabbage-caterpillars).
- *M. jurtina* larvae use many Poaceae grasses; original near-source observations and more tree+grass mapped cover could reflect semi-natural or unmanaged habitats, but WorldCover tree/grass are not verified host species or grassland management. [European Butterflies reference](https://butterflies-europe.linnaeus.naturalis.nl/linnaeus_ng/app/views/species/taxon.php?id=114880).
- *P. icarus* larvae use legumes such as *Lotus corniculatus* and frequently associate with ants. Its measured 2025 recorded/day rate is LOWER near despite an overall three-family ratio >1. This already rules out interpreting three-family composition-positive as evidence of a universal Lycaenidae nutritional subsidy. [Naturalis Europe species ecology](https://butterflies-europe.linnaeus.naturalis.nl/linnaeus_ng/app/views/species/taxon.php?id=114729).

**Conclusive scientific boundary:** source proximity is a fixed 2022 registered EEA FW outfall location, not a 2025 identified physical butterfly effluent dose. Independent original 2025 HydroRIVERS flow check found strictly downstream candidate site coordinates are **10 but collapse to only TWO source outfalls/reaches**—not enough for a directed biological effect and not evidence all other butterflies avoid water. Demographic benefit, source P/N effects, ants, larval hosts and adult mud-puddling remain unmeasured. Preserve earlier DE/AT positive photo fractions while identifying the *reporting denominator* as a key structural confounder.

## Source and independent reproduction

- [Original licensed 65 accepted names→63 canonical GBIF IDs and all 2025 recorded species groups, all 399 site family reconciliation](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2025-species-composition-decomposition-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DK2025_NEW399_ORIGINAL_65SPECIES_3VS5_FAMILY_DECOMPOSITION.json) — [original eBMS v1.4 re-read SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42786494070).
- [Separate original-source only arithmetic audit all 63 canonical species and leave-one-species-out 3 vs 5 OR exactly reproduced](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2025-species-composition-decomposition-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DK2025_NEW399_SPECIES_ORIGINAL_63CANONICAL_INDEPENDENT_QA.json) — [independent QA SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42786902576).
- [Original 2021 WorldCover + original 2025 top3 species six-region habitat model](https://github.com/zuizui0223/gesui/blob/analysis/dk-ebms-2025-species-composition-decomposition-20261110/results/ebms_denmark_2023_2024_WWTP_source_only/DK2025_NEW399_TOP3_SPECIES_WORLD_COVER2021_6REGION_POSTHOC_QC.json) — [licensed original 2025 source re-read and 3,000 region bootstrap SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42787364808).

**KEEP separate Draft PR; no merge to gesui main or chocho manuscript.** Original CC-BY-NC eBMS archive raw records, pixel raster, private original observer IDs and raw 2025 site×species individuals are not redistributed.