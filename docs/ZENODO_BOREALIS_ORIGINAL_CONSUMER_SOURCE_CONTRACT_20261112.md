# Source-first extension: original publicly released consumer field data and wastewater contaminant transfer

**2026-11-12 JST.** Proposed additional evidence panels after PR #61. **Posthoc** by design: previous published paper results and HJA source/association results are known, but *no values from these newly accessed author deposits have been used to choose a favorable sign or new model*.

## A. Rideout et al. 2025: emerge aquatic insects → actual riparian birds

Canonical published original Dryad study: **doi:10.5061/dryad.sbcc2frfk**, original public Dryad dataset metadata version 6, nine named files and *original SHA256 digests*, collected separately from API (original list [Actions #43143134955](https://github.com/zuizui0223/gesui/actions/runs/43143134955)). Direct per-file Dryad API returned **HTTP 401** to anonymous Github runner and was not circumvented.

**New official alternative**: that exact original public Dryad dataset landing page explicitly lists the 2025 authors' accompanying Zenodo accession **doi:10.5281/zenodo.11186734** and version **10.5281/zenodo.11186733** under *Works referencing this dataset*. Search publisher-registered Zenodo API `https://zenodo.org/api/records/11186734` and `...11186733`, plus **only** public file URL advertised by that record.

**Source identity criteria:** author/title/year and article DOI match; compare individual file SHA256 to Dryad public file metadata for each:
- `sites_full.csv`: `3b4b38303fc3e70dcd00411742a20da8d77147cccc1ac56239736d9249e8afef`
- `sites_bird.csv`: `b4dc6834eb0b2d6a1906adf2a74dd1d8180765378e11f3e967302f271e4a5849`
- `emergence.csv`: `e36990e83a8a5f79ad06ca097911e849e2e6f3414b030d20753cf4376ac80dd4`
- `emergence_seasonality.csv`: `719552e0c0967f7e60db939fccb68a62d9cfff2d3f7daf2cff2`
- `bird.csv`: `5486ebd00866a6295cf4041f7dbaa873406c9c8daa0cc03ca8fb882065e1679f`
- `environment.csv`: `17c21017c5e385ecbd9e8fef89325854da533b1fbca2b754b1434f1b361d4bee`.
Zenodo may contain a software-only record instead of the original data or modified files; **never silently substitute it**. If exact original original SHA256 is not available, keep each source status "separate version, not original-reproduction", do not force data merge.
- Source gate: **>=20 unique native sites** available in both real 2022 emergence and 2022 bird field observations, original site IDs consistent, both insectivore metrics and biomass nonmissing for >=20 sites. Preserve original site groups; at least 3 sinuosity/habitat groups. Full 29 original sites (the paper analyzed 24 for bird fields) are NOT automatically 29 independent bird measurement units.
- Subsequent conditional analysis: compare actual 2022 emergence biomass/seasonality with field bird insectivore fraction on same native site, adjust sinuosity and original hedge area, no same-year causal nutrient/benthic link (chemistry sampled **September 2023**). No eBird.

## B. Direct wastewater microplastics via adult insects → riparian spiders

Wardlaw et al. (2026), *Environmental Pollution*, doi **10.1016/j.envpol.2026.129116**, released online 2026-09-05. **Eight original municipal WWTPs in Ontario Grand River catchment, 2022 and 2023**, with both up- and downstream larval/winged aquatic insects and Tetragnathidae spiders. Original paper **already reports** downstream larval Heptageniidae ~2.4x, adult Chironomidae ~3.7x, spider Tetragnathidae ~2x microplastics, variable by site. This is a **contaminant transport/harm** pathway, not a direct pollination/service benefit. Trophic transfer factors generally <1 adult→spider: no generic biomagnification.

The journal data statement says public Borealis V1 Canadian Dataverse repository on publication. Query only legitimate *public Borealis* Dataverse metadata/search API for original author/title/DOI, then if available follow public direct dataset DOI and native version, file SHA, original WWTP pair/time/taxon units.

Related but **distinct** published Wardlaw study (2026) *Limited evidence of elevated microplastics in lower trophic levels downstream of wastewater treatment plants*, doi **10.1093/etojnl/vgag111**, explicitly provides Borealis public data **doi:10.5683/SP3/TLDY6W**. Do not pretend it is the adult-insect/spider data; it may provide upstream-downstream source context or a within-study independent exposure falsification panel only if matched original WWTP IDs and dates are explicitly preserved.

**Pre-numeric gate**: >=6 independent original WWTP upstream/downstream pairs with adult insect **and** spiders measured in contemporaneous native period, actual data and real taxon IDs, original microplastics in consistent units and documented detection limits. If unavailable, source-only STOP, cite paper findings as original authors only. Never invent results or reconstruct plotted values.

## Across-panels decision framework

1. River N/P ↔ benthic EPT (US NRSA 1,220 original UID visits) is independently measured **but not source-specific sewage**;
2. Adult emergence quantity/temporal pulses (Oregon HJA 17 streams) is measured, but different site/year than US NRSA;
3. Realized predators (Ontario 29 sites, Quebec 20 sites) have field-source metadata but no native raw input yet;
4. Actual WWTP microplastics → adult insects → spiders (Grand River 8 WWTP) exists as a published observation, potentially a *different contaminant axis* not N/P trophic subsidy benefit.

**No cross-study row-level SEM, no coastal island analogy as proof, no cohort-free interchange of positive/negative associations.** The universal scientific question is: when does an apparent aquatic biological benefit or degradation result in positive/negative transfer to terrestrial consumers, according to **quantity × quality × timing × accessibility × contaminant transport**?

Only new unmerged Draft branch; preserve PR #61/60/59, gesui/main and chocho unchanged.
