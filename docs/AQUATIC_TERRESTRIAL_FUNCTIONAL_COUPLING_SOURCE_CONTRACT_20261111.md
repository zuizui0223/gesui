# Aquatic–terrestrial ecosystem-function coupling: **pre-outcome source contract** (2026-11-11)

## Scientific question
When are the ecological integrity of river benthos, potential export of aquatic insects, riparian plant resources, and butterfly/predator occurrence geographically **coupled** versus **decoupled** across nutrient and wastewater-source gradients?

**Important:** Direct **ecosystem services** (pollination, crop yield, pest suppression, water purification), aquatic insect **emergence**, and trophic transfer are not automatically measured by EQR, larval benthic counts, GBIF records or co-occurrence. All claims must be named as **potential functional linkages/proxies** until direct service measurements are available.

**No eBird. No inference that rivers and wetlands on different continents or catchments form a measured trophic chain.**

## Existing independent evidence (do not modify)
- EEA Europe's 2018 and 2022 upstream actual WWTP reported N/P adds small out-of-region information on benthic EQR. Not a sign of positive biological effects; report/source missingness and hydrologic uncertainty remain. `gesui` PR #14/#17/#56/#57.
- U.S. NWCA2021: real **wetland** plant and TN/TP measurements but WWTP-to-wetland flow routing support is insufficient; do NOT join nationally separated wetland plots to arbitrary streams to imply functional transmission.
- WPI/GBIF global botanical and butterfly photo associations primarily support recording co-geography and candidate host geography, not verified pollination or discharge causation.

## New first independent source: EPA NRSA 2018–2019
Official [EPA National Aquatic Resource Survey individual CSVs and metadata](https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys) are the *same* nationwide river visits with
1. benthic macroinvertebrate **taxon count / order / family / functional feeding group FFG** (BERW and BETB protocols must never be naively pooled);
2. site and survey geography (site identifier, visit, state, AG_ECO9, where available);
3. real water **total nitrogen** `NTL_RESULT` and **total phosphorus** `PTL_RESULT` with units and QA flags;
4. optional landscape/physical habitat controls in subsequent analyses.

Original official URL set (not copied/interpolated from elsewhere):
- Benthos: `https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_benthic_macroinvertebrate_count_-_data.csv`
- Site: `https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_SiteInfo.csv`
- Chemistry: `https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_water_chemistry_chla_-_data.csv`

## Phase 0: before *any numerical nutrient association*
Only read source headers, file hashes, row counts, unique native site+visit keys, sample-type support, taxonomic guild coverage and whether TN/TP numeric fields are **present/nonblank** (not their actual magnitudes). Do not fill missing with 0. No joining site rows to taxon rows unless source grain is validated.

The *source* gate for a full macro analysis is deliberately modest but independently defined here:
- at least 300 **unique** primary-protocol river site-visits with count, site record and both positive finite TN and TP (gate checks only availability, not nutrient values);
- at least 25 U.S. states and six distinct NRSA `AG_ECO9` units among those source-eligible site-visits;
- at least 60% of positive individual-count mass of primary-protocol taxon rows has nonempty `FFG`, and at least 80% has a known `ORDER` or `FAMILY` classification;
- unique scientific observation sample key and non-overlapping source modes; disallow double-counting taxonomic aliases or same organism through two sample types.

**Primary protocol**: NRSA `BERW` benthic count with 300-individual subsample `TOTAL300` if actually available and nonzero at a site. If source data show this is not possible at 300 sites, STOP; do not silently switch to `BETB`, mix protocols or relax thresholds. Other sample types may be separately audited, never silently included.

## Phase 1: conditional prespecified *river-only* functional proxy model
Only after Phase 0 passes, archive the frozen site manifest before opening TN/TP values. Define biological site response `EPT_individual_fraction`, `Chironomidae_individual_fraction`, `FFG_shredder_fraction`, `FFG_collector_gatherer_fraction` and aquatic benthic richness at identical source protocol; count zeros genuinely only after confirmed census. EPT taxa and chironomids are *potential*, not actual emergence.

Main environmental predictor separate `log10(TN mg/L)` and `log10(TP mg/L)` with QA, units, laboratory detection handling and site duplicates; never sum N and P unitlessly. Prespecified adjustments: ecoregion and state, year and site visit, stream size, conductivity, landscape cropland and urban cover where actually measured, and available field effort. Region-disjoint 5-fold, state-block bootstrap and source-year validation if available; no artifact of sample size or observation effort. Also fit negative-control taxonomic groups/metrics. Need enough states before inferential fitting, and audit variance in nutrient within strata. If sources fail, no coefficients are claimed.

## Phase 2: *joint geography* aquatic–terrestrial resource and pollinator bridge
At the same preselected river site coordinates and within 2017–2021, use independent public original GBIF GBIF 5 km buffers (10 km sensitivity), separately query:
- Plantae accepted speciesKey, including separately fixed *candidate butterfly host* taxa from `chocho` and GloBI interaction archive; log raw plant recording opportunity;
- Lycaenidae and other butterfly families using the same original GBIF observation quality and 2017–2021 time period;
- riparian candidate aquatic-insect consumers (Odonata and e.g. spiders/bats as non-equivalent comparators), **without** calling observed occurrence proof of predation at a site.

Source eligibility before any coefficient: >=250 original benthic-chemistry sites with complete original GBIF source requests, >=100 with >=50 plant reports and >=20 butterfly reports, >=15 states and >=5 eco groups with support; retain zero vs API missing distinction. If the gate fails, report source coverage only. Further analyses must control plant/butterfly recording effort and unrelated taxon reporting, avoid inflated spatial pseudoreplication, and test 5 vs 10 km separation.

In all phases: spatial correlation ≠ measured aquatic-to-terrestrial **flux**; no mediation/SEM claim without temporal source alignment, actual intermediary measurement, directed water delivery and falsification tests.

## Published methodological context
- Schürings & Olden 2026, DOI `10.1038/s41559-026-03041-1`: 14,177 US river EPT surveys and riparian insectivorous bird responses. **Our work does not reuse eBird** or advertise their result as ours; functional feeding groups and wastewater pressure would be the new question.
- EPA NRSA sample manual and data download/metadata: official EPA source above.
- USGS 2025 DOI and freshwater/terrestrial contaminant transfer studies: need actual insect emergence/pollinator function for causal transfers.

## Execution and provenance
All network requests in ephemeral GitHub Actions runner; store original source URLs, response byte SHA256, native site/visit key counts, missingness and literal QA decisions in the Actions artifact. No licensed/PII raw records committed. The initial commit contains **this contract before any new EPA source is downloaded**. Report source-only STOP, or preserve a source-locked candidate cohort before reading numeric nutrient or taxonomic effects.

**Keep Draft, no merge to `main` / `chocho`; do not overwrite PR #14/#17/#56/#57.**
