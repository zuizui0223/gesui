# Pre-biology extraction contract: Austrian WISA monitoring station-year as the common wastewater–water–benthos–riparian unit

Created 2026-11-16 after final full-source EEA2022-282 source audit. **A prospective source/data acquisition gate only**; no actual WISA station rows or biological coefficients have yet been retrieved. It is NOT an executed ecosystem result and does not retroactively rescue the original EEA2022-BQE shortage.

## Why Austria / WISA, not more same-EEA H7 correlations?

The true original EEA UWWTD measured 2022 one-FW-outlet N/P national DE/AT source has original 282 source/riverline valid outlets (AT145/DE137) and source verified up/down ±1,500m rivers. True 2018/22 same-month original Sentinel2 red/NIR/SCL source QA on preselected48 yields 26 clear two-year source systems (AT14/DE12). But published EU EEA2022 biology **DE0/AT68** source sites and genuine 2022 WISE6 TP 301 geocoded field stations match 2km proximity on only one of 282 same-outlet source systems, with zero confidently same-outfall/D1 native reaches. We therefore have *incomplete source-table intersections*, NOT evidence of no nutrient transfer.

The official Austrian ministry's WISA cross-database query openly advertises these three scientific components accessible in a **common thematic query**:
1. H2O original physicochemical/chemical water quality at monitoring stations (TN/TP dates).
2. Federal benthos macrozoobenthos, phytobenthos and aquatic macrophyte survey records, with river-site/site-ID references.
3. EmReg-OW loads **aggregated to flow-water monitoring stations** (not necessarily a single specific WWTP outfall).

Government official landing: https://www.bmluk.gv.at/themen/wasser/wisa/daten/daten.html
Public JS frontend: https://wasser.umweltbundesamt.at/wisa-datenabfrage/
Related H2O official app: https://wasser.umweltbundesamt.at/h2odb/

**Source accessibility caveat:** These are public-facing JS applications; no programmatically verified station-year export, exact 2022 source records or complete discharge registry join has yet passed. Do not invent a REST API or promise full public bulk access. If public data access or terms require interactive query/registration, record that blocker as source accessibility, not a biological null.

## Phase 0 — source access and provenance (before any numeric ecology)

- Confirm whether official WISA/H2O outputs can provide legal reusable tabular downloads of station metadata, exact original coordinates/CRS, observation dates/year, river reach/site IDs, TN, TP with units/LOQ and source quality fields, macrozoobenthos/phytobenthos EQR or taxon tables, and EmReg-OW reported loads with original station/year.
- Preserve direct official API/export URLs, explicit access restrictions, filename original SHA256, original schema, permissions and row counts; don't copy proprietary or personal identifiers to public repo.
- Only TRUE original station/year records may be used. A "water body ecological status 2021" map is not an original 2022 site sample. Missing target-year chemistry or biological observation is unknown, not zero or control.
- If no independent same-station biology/water/annual emission source is available, emit SOURCE ACCESS STOP; DO NOT fit cross-boundary SEM using separate site averages.

## Phase 1 — exact source support (source-only, BEFORE NDVI or EQR outcomes opened)

Fix the true 2022 EEA N/P outfall source population of **145** previously geometry-passed AT original facilities, with genuine direct FW and measured annual N/P (not modelled W), and original 2022 plus 2018 satellite geographic source (original response-blind pilot24AT, 14AT source two-year local pixel pass). Keep ALL145 as source denominator and source observer-selected matching eligibility, not just posthoc selected satellites.

One potentially eligible facility must have:
- A confirmed WISA 2022 river monitoring station, with observed N and P **actual receiving river concentrations** in same year, actual latitude/longitude, original stable station ID, year and waterbody/taxon ID;
- Original MZB or phytobenthos 2022 recorded for the SAME physical monitoring site or a separate clearly source-mapped site within the SAME verified receiving river path, with biological observation timing and comparable assessment method;
- Original source annual load mapping from EmReg-OW to that WISA monitoring station with explicit attribution/multiple discharge entries, or separately a proven hydrologic directed EEA UWWTD 2022 outlet SOURCE → monitor path (and an actual source-negative reference reach that is not merely missing a report);
- A true outfall-local ±1,500m river path both bank circles and enough Sentinel-2 valid land pixels in 2022, with original source QA; the 2018 image serves only as historic land cover/phenology context unless 2018 genuine outfall mass and receptor water chemistry are independently available.

Prespecify **source opportunity minimum** of at least **20 independent outfall receiving-river systems** spanning ≥3 genuine distinct hydrological parent basins, with one real water chemistry station and one real benthic site in source year 2022. In the independent matched no-outfall reference pool require ≥10 genuine different receiving reaches with comparable river size/land use, same monitoring year. Record all original passed/failed source IDs, coordinates-quality and multi-outfall attribution. These minimums are DESIGN TARGETS, not data claims; revise only in a new declared cohort before numeric effects.

## Phase 2 — only if full source support

Primary outcome: true river benthic group condition/counted traits vs source-specific measured receiving-water TP/TN change, with upstream and reference river controls, year/method/effort/flow/land use, independent parent-basin CV.

Secondary terrestrial outcome: same-season true upstream/downstream riparian **stable landcover-masked green canopy index** from Sentinel-2 2022, original red/NIR/SCL scale+offset and local cloud pixel source QC. Only the source-year 2022 N/P reported load and 2022 water chemistry may enter the 2022 exposure. NDVI ≠ local plant biomass, biodiversity or net nitrogen assimilation.

**Cross-boundary hypothesis** (NOT yet identified): whether the direction/magnitude of water nutrient/benthos responses predicts co-located riverbank greenness beyond background, and whether genuine riverbank connectivity/hydrology moderates that association. Causal water→land material transfer requires its own evidence (leaf/soil N isotope, flood connectivity, insects etc.). Do not combine unrelated French TP/I2M2 sources with Austrian 2022 NP/Sentinel into a synthetic causal chain.

Only after source-first access/row count gating and a prospective analysis specification may the response-level regression be run. If STOP, keep distinct water and terrestrial results and do not create an unsupported wastewater impact paper.
