# Denmark standardized butterfly monitoring × EEA wastewater: source-only protocol (2026-11-10)

**Status:** independent source-availability pilot, frozen before opening new Denmark butterfly occurrence-family outcomes. Motivating WWTP × Lycaenidae photo associations were known in Germany/Austria, so later hypotheses remain exploratory. Failed original Germany 73-pair and Sweden 2022 source gates are unchanged.

## Source data frozen before butterfly counts

1. **Standardized sampling effort:** eBMS Denmark transects, version 1.4 published 2026-02-10, GBIF dataset caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c, CC-BY-NC 4.0. Public source: https://ipt.brc.ac.uk/archive.do?r=ebms_denmark_transects&v=1.4 . Prior original current archive SHA256 2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b. Inspect the source Event core and meta.xml ONLY. Do not open any Occurrence extension/species outcome at this stage.
2. **WWTP source:** official public EEA 2022 UWWTD TreatmentPlants_2022/MapServer/0 and DischargePoints_2022/MapServer/0 REST. Source count 5,492 original plants with both positive measured N/P, 24,416 discharge points. Join nonempty original reporter+plant key uniquely with exactly one original freshwater (FW) direct outfall, Danish registered reporter DK2022, coordinate box W7 E16 S54 N58. Previous EEA source-only report recorded 69 qualified Danish outfalls; reverify original API records and result rows before any butterfly measurements.
3. The 2022 EEA discharge values are annual reported source mass, NOT receiving-water chemistry or adult puddling protein/salt, and not a known change in 2023/2024 effluent exposure.

## Event-only sample definition

- Independently count original 2023 and 2024 eBMS **Event** rows and unique calendar survey dates per original public site centroid, using meta.xml Event fields. Keep 2025 only as descriptive source availability. Count original site-year visit days, not event fragments.
- Strict repeated-site sample: same public latitude/longitude rounded to 5 decimals, at least **2 distinct visit days in each of 2023 AND 2024**, coordinates in Danish country rectangle 54–58 N, 7–16 E, explicit original coordinateUncertaintyInMeters <=500. If uncertainty missing, record it and exclude from strict local site proximity; never invent accuracy. Record repeated 2023/2024 sites before source-proximity analysis. Locality and parentEventID can diagnose site key quality; do not use them to invent duplicates or observer names.
- All source control and quality results must be recorded before any original butterfly occurrence table is opened, and never tune thresholds on signs of butterfly relationships.

## Fixed source-distance and geographic support gates

- Near primary is distance <=2km to closest qualifying source, secondary <=5km. Remote >=10km from every original **qualifying registered** outfall. 5–10km intermediate is neither near nor remote. Remote is not pollution-free.
- Primary: >=8 near <=2km and >=8 remote repeated-visit sites in >=3 original 1-degree geographic cells each containing both statuses.
- Secondary: >=12 near <=5km and >=12 remote repeated-visit sites in >=3 original 1-degree cells containing both statuses.
- If either radius does not meet its full criteria, STOP butterfly outcomes for that radius. No retuning or widening of radius after results. If both fail, archive source gap with zero butterfly outcomes read.
- If a gate passes: commit site list, original raw archive hash, year-specific visit dates/effort, WWTP source hashes and exact source distance band before fetching any Occurrence extension. Then require >=6 site-year occurrence count visits per source band in both years and >=5 distinct sites per source band with Lycaenidae/Pieridae/Papilionidae recorded each year before family-response coefficient interpretation, else STOP with explicit sample gate failure.

## Biological meaning of later study (still untested)

If all gates pass, an independent 2023-2024 standardized butterfly survey family composition near WWTP study could test whether Lycaenidae counts/relative share have a source-proximity association beyond survey effort, land cover and geography. 2023-2024 are not different exposure treatments; repeated annual surveys allow stability but not effluent causal change. WWTP 2022 reported P/N mass is not directly measured butterfly water/substrate protein/Na, larval host fitness, ant effects, or source-to-butterfly downstream river pathway. Separate independent family controls and standardized survey effort are required.

Previous original SeBMS 2022 matched pilot found 14 eligible Swedish WWTP freshwater sites, no qualified standardized transect sites within 2km and two point-count sites within 5km; no 2022 butterfly outcome was read. Do not erase this failure.

Store public source material with provenance and original CC-BY-NC 4.0 attribution. Keep new data on isolated Draft branch, no main or chocho changes.


## SOURCE-ONLY CORRECTION 2026-11-10 (before any species outcomes; do not erase history)

Initial automated source-only run [42757019282](https://github.com/zuizui0223/gesui/actions/runs/42757019282) stopped before producing an outcome or opening any eBMS Occurrence extension: it checked for **69** 2022 Denmark-reporter WWTP direct-FW outfalls, but the original EEA REST join and complete fixed Danish national envelope W7 E16 S54 N58 actually yielded **203**.

The earlier 69 count was only the number of Danish DK2022-reporter outlets falling inside a **Sweden-oriented overlap bounding box** (reported in the prior SeBMS feasibility inventory), **NOT the national Danish number**. This was a mistaken translation of an existing source geographic crop, not a result-driven attempt to enlarge a radius or choose favorable butterfly observations. Original EEA filters, 2022 source reporter DK2022, one-unique-direct-FW-outfall rule, Danish full bounding box and all 2km/5km near/remote/visit support gates remain EXACTLY as declared above. Corrected expected source-only EEA total is 203. The original failed workflow and both values remain in the history. No original Danish species abundances have been inspected at this correction stage.
