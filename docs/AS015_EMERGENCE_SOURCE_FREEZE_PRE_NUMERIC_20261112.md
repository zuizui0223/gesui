# Aquatic emergence calibration: source-first and effect-blind specification (12 November 2026)

**Purpose:** Determine whether aquatic insect community composition (EPT versus Chironomidae) is a reliable surrogate for the **quantity of actual adult emergence** exported from freshwater. This is an independent **mechanistic calibration of measurement**, **NOT** an EPA NRSA 2018–19 → Oregon river station link, nor a measured wastewater effect, nor a regional/macro replication.

**Authoritative original source:** Latorre-B, T.; Arismendi, I. (2026), *Adult aquatic insect sampling in Lookout Creek tributaries at the HJ Andrews Experimental Forest, 2022* (AS015 v2), EDI `knb-lter-and.12548.2`, DOI `10.6073/pasta/6e9b83e2547e3a82d92dab068476f1d7`. Official https://andrewsforest.oregonstate.edu/data/datacatalog/AS015 ; EDI PASTA+ API `https://pasta.lternet.edu/package/name/eml/knb-lter-and/12548/2` for original entity IDs; then `/package/data/eml/knb-lter-and/12548/2/<entityId>`. Open licensed CC BY 4.0, cite LTER data.

### What is known from the public *metadata before source outcome access*

- 2022-06-25 to 2022-09-01, **17 distinct first/second-order headwater streams in one Oregon Lookout Creek basin**, comprising **9 intermittent** and **8 perennial** reaches.
- Four emergence traps per reach, each 0.19 m², sampled at approximately three-to-four-day intervals.
- **1564 original taxon × site × sampling period rows**, nine fields; `SITE, STREAM_REGIME, SAMPLE_DATE, TIME_PERIOD, TAXA_ORDER, EMERG_ABUNDANCE, EMERG_BIOMASS` plus 2 source codes.
- Four reported categories: Ephemeroptera, Plecoptera, Trichoptera, and Diptera (specifically chironomids). Taxon amounts are **emergent adults per m² and g dry mass per m² per sampling interval**, not benthic counts and not annual flux.
- Metadata state 23 sequential sampling periods. **Expected mathematical design:** 17 × 23 × 4 = 1564; **must verify exact original keys/coverage** before using zero in a biological statement.

### Stage 0: original file source and outcome-blind cohort freeze

Read original EDI package metadata/entity names, fetch *actual original table* and require precise published MD5 `f3e6e400990da883ae0d341731431749` (metadata AS01501_v1.csv). No mirror, artificial rows or inferred zeros. Source QA on native columns, site × period × order completeness, site-regime consistency, date windows, valid numeric availability (not magnitude) only.

Source gate fixed **before magnitude access**: at least 15 of 17 streams, both hydrological regimes ≥6 streams, ≥18 unique sampled time periods per stream, at least 300 complete site × period × 4-orders events, ≥95% finite nonnegative numeric availability in both EMERG_ABUNDANCE and EMERG_BIOMASS (counts and biomass may have true numeric zero but missing is not zero), no duplicate primary keys. If source fails, STOP and no effect signs reported.

Store original byte MD5+SHA256, source URL and field names; freeze a deterministic hash-locked allowed native SITE/TIME_PERIOD/date/regime key manifest. Source stage MUST NOT compute any relation, mean, difference, correlation, or regime effect from numeric amounts.

### Stage 1: only on Stage 0 pass and using exactly its immutable manifest

Unit of *independent spatial inference* = **17 streams**, not 1564 individual taxonomic rows or 391 time-period panels. Unit of within-stream temporal description = sampling interval with full 4-order observations, **not** independent annual realizations.

Outcomes fixed:
1. per stream, across the same valid sampling periods, mean (per-interval) `EPT abundance fraction` = E+P+T adults / (E+P+T+Chironomidae adults); if total adults=0, EPT fraction **missing** for that interval, not 0.
2. per stream, mean *total emergence biomass* (g/m² per reported interval) and total emergence abundance (individuals/m² per interval); compare abundance and biomass, not replace one with the other.
3. `EPT biomass fraction` and `Chironomidae biomass fraction`, both with 0/0 missing; taxonomic resolution is limited to these four orders/family and **not** complete emergence across all aquatic insect taxa.
4. `EPT-dominant but low-biomass` versus `midge-dominant but high-biomass` discordance **as a descriptive pair of observed relative orderings**: predeclare low/high split by site-equal empirical median only as a secondary visualization. Do not claim a deterministic ecological type.
5. Hydrologic regime (9 intermittent, 8 perennial) contrast of *site-mean biomass* and EPT share separately, by an exact `choose(17,9)=24,310` label permutation assuming site exchangeability. Show 95% interval from 5000 stream bootstrap resamples **descriptive only**. Recheck elevation, dates and trap area to avoid spurious stream effects.

Primary correspondence: Spearman `rho(EPT abundance fraction, mean total emerging biomass)` over **17 site-level points**, plus second `rho(EPT abundance fraction, total abundance)` and ratio of biomass to abundance. No significance based on 391 repeated time events. For 17 stream sites, report exact rank permutation via Monte Carlo 20,000 labels/10,000 bootstraps if justified, or no inferential p if site dependence severe.

Within-stream temporal correspondence: demean EPT abundance fraction and log(1+total biomass) per stream, report site-cluster bootstrap and per-site signed effects without calling 391 events independent. Time trends and differences between regimes are exploratory, with repeated measures not pseudoreplicated.

### Critical interpretation boundary

This source can test whether **EPT composition is an adequate proxy for the amount of emergence**, which was an unsupported bridge in `gesui` PR #59. It **cannot** demonstrate that rising nutrients or WWTP discharges reduce/increase emergence. These 17 Oregon streams are in a *different landscape/year* from EPA NRSA 2018–19 sites; do NOT fit a merged EPT–TP–biomass row-level model, infer mediation, or transfer ecological slopes without an actual same-site measured chemistry and benthic source.

No eBird, no unsupported ecological-service/pollination claims, no undocumented switch of source/threshold. Results staged under new Draft branch, preserving PR #59 and original source/response data unchanged.
