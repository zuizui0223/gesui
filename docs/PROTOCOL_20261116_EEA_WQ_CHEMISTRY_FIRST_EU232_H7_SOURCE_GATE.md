# Source-first: EEA 2026 WQ→original European flora/butterfly H7 opportunity (freeze before chemistry)

2026-11-16, exploration based on immutable PR #107 original source and main PR #106, BEFORE any new chemical station-year numeric values are inspected.

## Motivation and a changed sampling design

Prior PR #107 found only 9 original H7×period /8 H6 with strict *identical EEA benthic biology station and year* ambient TN and TP plus confirmed registered upstream WWTP annual N/P reporting, below minimum 30/20/10/3. Original alternative ≤1km EEA Water Quality station design still has a mathematical upper bound 16 report-positive H7×period, below minimum 30; do not try to relax the original criteria.

Here we deliberately **drop EEA benthic site-first selection**, not the water measurement identity checks. Select the whole official Water Quality monitoring station population inside exactly the prior 232 original HydroBASINS L7 polygons used for GBIF plants/butterflies. The new outcome, if sufficient source support, is the spatial association between ambient water TN/TP and *plant and butterfly habitat-nutrient-affinity composition* in the same H7×period. The outcome is **NOT** water benthic EQR, not wastewater effluent N attribution and not a field local foodweb.

## Original source geography locked before chemistry or ecological outcomes

- Original 454 genuine H7×period /232 H7 /90 native parent H6 panel SHA256: 1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca.
- Official HydroBASINS Europe level7 ZIP SHA256 b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44; 232 genuine polygons source mapped to 90 original native H6. All polygons valid, point-to-polygon shapely covered_by (exclude ambiguous border assignments).
- Official July 2026 EEA Water Quality ICM SpatialObjects source ZIP 2,797,036 bytes; SHA256 aa505be6f85d06c675e0d0fcb148a9148f80162f8e9590740a887cc95a3e349d. Source has ~69,820 valid fully identified WGS84 site locations. Require the original country/site identifier/scheme triplet, unambiguous coordinates and location inside ONE original H7 polygon.
- Map water station country to original H7×period majority country; **primary requires same country**, report cross-border any-country only as source opportunity sensitivity. No EQR/benthic station or closest distance requirement, and NO selection by NUTRI/HNS or W coefficients. No new GBIF queries or trait assignments.
- Exact source periods: 2018 annual ambient water quality ↔ original GBIF 2017–2019; 2022 annual quality ↔ original GBIF 2020–2022. This is a **3-year GBIF record window** associated with a single water-monitoring year, NOT same-day or same-year biological responses.

## Chemistry source, before any biological coefficients

- Official 2026 EEA Waterbase ICM AggregatedData original Part2 170,422,309 byte ZIP SHA256 315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a, **6,450,192 original rows**, verified same frozen original PR #107. Its annual-site statistics are from reported data where disaggregated samples are unavailable or an annual view is supplied. This **Part2-only** design is not comprehensive EEA Part1+Part2 chemistry; absence is unknown, not measured zero.
- Scan source with full four-part key (country, monitoringSiteIdentifier, monitoringSiteIdentifierScheme, calendar year); avoid cross-country identifier collisions.
- Only year=2018 or 2022, parameterWaterBodyCategory=RW, procedureAnalysedMatrix=W, total N determinand EEA_31615-01-7 and total P CAS_7723-14-0. Preserve canonical unit mg{N}/L, mg{P}/L, source finite mean>=0, resultNumberOfSamples>=3 and resultQualityMeanBelowLOQ not true. For each station×year×nutrient, require **exactly one** valid annual row. No nitrate/phosphate substitution, no missing zero-imputation, no mixing TN and TP from different stations.
- Freeze **source station-year** as the source-grain; when >1 qualified WQ station shares same H7×period/country, set each nutrient to the *median over these SAME qualifying TNTP station-year pairs*, with equal station weight. No species abundance, macroinvertebrate EQR, plant/butterfly score or local WWTP label is read for station choice.
- Primary chemistry panel = WQ station country matches original original GBIF basin country. Also count any-country panel, but not switch to it for a favored effect.

## Source-only hard gate (decide before any plant or butterfly fit)

In the primary country-matched **both strict TN and TP** station-year source panel, require at least:
- 30 original H7×periods;
- 20 distinct original H7;
- 20 genuinely distinct original parent H6;
- 3 original WQ country groups.

If any threshold fails: label **STOP_SOURCE_CHEMISTRY_PANEL_SUPPORT**; do not fit any biological model or silently relax a threshold because a coefficient looks favorable. Keep all source denominators and separate chemistry Part2-only limits. Do not require a WWTP-positive report in this new panel; such a requirement reintroduces the same ascertainment-limited population and would not be a source-independent ambient-chemistry test.

If all pass, an independent next stage may test:
1. Original same H7×period GBIF plant EIV NUTRI species-equal versus same ambient TN/TP, adjusted for **original land cover/river flow/country-method/source opportunity & W**, 5-fold H6 blocked Ridge10 and H6-cluster descriptive WLS. Analyze N and P separately; no P-vs-N causal mechanistic partition.
2. Butterfly host-nutrient HNS species-equal with plant NUTRI and observed ambient TN/TP, original source missingness and effort controls; no larval feeding or mediation claim. At least 20 H6 already, but substantial H6 leaveout power is limited. First pre-check independent H6 N and explicit new country and year before fitting.
3. An *independent* official WWTP-source exposure analysis requires actual outfall–WQ-station hydraulic upstream reach/flow/year alignment, measured source flow and credible no-exposure controls. Original HydroWASTE static basin outlet W alone and the source-reporting indicator do not qualify.

## Science limitations and reproducibility

The original water-quality measurements are not actual delivered wastewater fraction or terrestrial soil/leaf nitrogen. Basin polygons are geographically overlapping environments, not proof of transfer from the aquatic monitoring station to plants/larvae. GBIF species-equal composition indices are not actual abundance/density or ecosystem services. Original source date window mismatch, multiple comparison and geography selection are material. New source QA is not an independent ecological confirmation. Do not automatically merge prior Draft PR #107 or this branch into main/chocho.

Files: analysis/20261116_chemistry_first_eu232/source_first_h7_wq_geography.py and source_first_part2_TNTP.py. Run source-only geometry then chemistry in separate steps with the signed original source ZIP/sha, save source-only station/year cohort and compact aggregate to an Actions artifact, do not commit full EEA tables or raw individual measurements.
