# Independent TP–plant river pilot: source-frozen before examining any new response

2026-11-07. New study **distinct** from the earlier fixed 16-outfall GBIF cohort, which is retained unchanged. The parent source-only cohort was selected without inspecting either measured TP concentrations or new GBIF plant assemblages.

## Frozen source selection

Use the exact six source-locked outfalls in `results/independent_deu_tp_source_cohort_20261107/STRICT_STATION_INDEPENDENCE_COHORTS.json` as selected by [source run 42371319376](https://github.com/zuizui0223/gesui/actions/runs/42371319376). Source: HydroWASTE v1 Git blob `0506f873f6541c5f324478894a9731290d722bce`, third-party HydroRIVERS Europe v1-derived GeoParquet SHA256 `5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f`, and Peters 2026 Figshare v5 `10.25400/lincolnuninz.30546371.v5`.

All six are in Germany, **three high and three low** modelled wastewater/river-flow fraction. They use **12 different TP stations**, six different HydroBASINS L12 units, and have pairwise outfall separation ≥57.28 km. The 16 original outfalls are excluded. Four-station joint TN–TP comparison is source-infeasible; this pilot is **TP-only**, not a surrogate for TN.

## Prespecified measured TP assessment

- Verify the original `TP.2_-_Harmonised.zip` (Figshare file ID 70896288, original MD5 `b596dbe9833c49426147a99f10f3fe31`). Restrict to the 12 preselected upstream/downstream TP station IDs; no nonpublic records, WRTDS interpolation or modeled TP substitution.
- Match only **positive finite measured TP concentrations** (mg/L) at an upstream and downstream station in the **same calendar month** from **2011-01-01 through 2022-12-31**.
- Each station-month is the median of observed concentrations; each outfall's descriptive TP gradient is the **median log(down/up) across matched months**, with counts of raw measurements and paired months. Require **≥12 matched months** before treating a source outfall as time-supported. Observations in a month are not independent replicates.
- Report TP gradient for all six, including signs and missingness. High/low and water-fraction association are descriptive only; no causal claim, no response-based switching of site IDs or years.

## Prespecified *new* GBIF plant survey

- At each of the 6 fixed outfalls, preselected directed river segments U2, U1, D1. **18 sites** fixed by the parent source-only geometry; do not replace any sites after inspecting plant outcomes.
- GBIF 2011–2022, kingdom `Plantae` key 6 (includes non-vascular plants), speciesKey facets with presence, georeference and spatial-quality filters, **1 km radius**. Report every site as an occurrence-index opportunity, not vegetation abundance/absence.
- Primary: hypergeometric expected richness after standardizing to **30 species-keyed records** per site, require each of U2/U1/D1 ≥30 species-keyed records and ≥6 taxa, plus max/min species-keyed records ≤10; `Y = [log(1+S_D1)-log(1+S_U1)] - [log(1+S_U1)-log(1+S_U2)]`.
- Prespecified low-depth sensitivity: 10 records and ≥4 taxa per site, clearly labelled weak observation support. Retain detected species-set Jaccard and distance differences as supplementary descriptive diagnostics, never raw migration claims.
- Join TP gradient to plant `Y` **only** if ≥4 outfalls satisfy both quality gates; report descriptive rho and denominator, no significance or general inference from n≤6. No response-driven threshold tuning.
- Report urbanization, cropland N, source plant effort, stream distance, time-alignment, other upstream WWTPs, seasonality as remaining potential confounders; no measured WWTP effluent chemistry or hydrological dilution at biological sites.

## Interpretation boundary

A positive TP river gradient would not demonstrate local sewage-P enrichment; TP stations are often 7–25 km from the outfall and may integrate multiple other sources. The WWTP fraction is modeled water discharge, **not observed P load**. A positive/negative plant association is only a small German spatial pattern. The previously positive `chocho` wastewater-volume/port-host associations and original `gesui` butterfly-free 48-site null are separate and preserved.

Outputs will be source MD5s, original station dates/values filtered to 2011–22, paired-month table, actual all-Plantae raw facets/site-species CSV, site-level estimator, and explicitly tagged cross-trophic cross-source *descriptive* correlations. No inference if the source/quality gates fail.
