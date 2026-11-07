# TN/TP source coverage decisions locked before reading concentration values (2026-11-07)

## The fixed old cohort: no power for a cross-outfall effect
Study population is the previously selected 16 HydroWASTE outfalls and 48 all-Plantae river sites, selected without nutrient or botanical outcomes. Source-locked Peters et al. *Global Nutrient Yields (2025 update)* Figshare `10.25400/lincolnuninz.30546371.v5` (CC BY 4.0) provides 19,424 TN and 41,639 TP site metadata rows, with 5,604 and 12,779 respectively meeting geographic and start/end-window filters in the pilot. These numbers are monitoring station *metadata*, not 2011–24 measured concentrations.

The source-first 500m nearest-HydroRIVERS-snap and directed NEXT_DOWN comparisons produced **only 1/16** fixed outfalls with BOTH up and down candidate stations for EACH of TN and TP within ≤25km river-network path; none within ≤10km. This fails the predeclared generalization/support gate; **do not regress the existing 16-outfall plant community index against TN/TP or change the first-cohort site choices**. See [source-run](https://github.com/zuizui0223/gesui/actions/runs/42369909177).

## Focal case preselection based ONLY on monitoring metadata and river topology
The sole doubly-supported outfall is **HydroWASTE WASTE_ID 24249**, Austria (outfall lat 48.560, lon 13.577; estimated water fraction about 1.710e-6 of local river discharge). The source network identifies:
- TN upstream `WQ016702` (site snapped to reach 20412137; upstream path ≈23.54 km), TN downstream `WQ016695` (reach 20412242, downstream path ≈8.59km).
- TP upstream `WQ976406` (same source reach as TN upstream; ≈23.54 km), TP downstream `WQ977059` (same downstream reach as TN downstream; ≈8.59km).
- TP upstream alternative `WQ977062` (observations start in 2021) is a **source-blind sensitivity**, not permitted to replace the primary pair depending on concentration sign.
- Outfall 24249 is modeled a **low** sewage-flow fraction source; even an observed downstream TN/TP increase would not identify sewage-N forcing rather than river-basin land use, other discharges, and flow.

## Predeclared measured-value inspection
1. Download only `TN.2_-_Harmonised.zip` and `TP.2_-_Harmonised.zip` from public Figshare v5, checking original file IDs and supplied MD5 (no WRTDS interpolated series or random forest fields). The source archives contain measured sample records; only the specified upstream/downstream station IDs may be extracted for this focal case. If the archive schema or licensing flag is incompatible, report an acquisition failure.
2. Primary common window **2012-01-01 through 2013-12-31**, the overlapping window of source site metadata for the primary TN pair, matched with TP on the same dates. Keep actual measurement dates (not merely first/last metadata).
3. Report original station/source/time/units/QC fields; exclude missing/negative/nonquantifiable measurements from numeric summaries and count excluded records. Do not assume TP and TN concentrations are N/P effluent output.
4. Calculate per-station monthly median of *observed* TN and TP in mg/L, then within-month `log(C_down)-log(C_up)` only if both stations have valid positive measurements in the **same calendar month** (no ex-post time-window widening). If fewer than six matched months in 2012–2013, label insufficient and provide counts only. Require matched months for the two parameters before comparing TN and TP signs.
5. Sensitivity: TP full source-overlap 2008–2022 with same original TP pair, but do not broaden TN's 2012–2013 primary window. Alternative TP upstream WQ977062 (2021–2022) is descriptive only and cannot substitute into primary.
6. Never claim causality or significance from one outfall; station spacing is markedly asymmetric, other urban/agricultural sources exist, model WWTP contribution is ~0.00017% of river flow. This case estimates a **river nutrient gradient**, not a globally wastewater-driven nutrient effect.

## New independent cohort (only after coverage failure is frozen)
A new monitored-source sample can be generated from the 3,569 source-only, geometrically eligible original outfall chains. Candidate selection must use only independent TN/TP monitoring source reach coverage, stream flow, area, country and quality before *new* GBIF plant outcomes are queried. New cohort and old 16 may not be mixed as an unbiased prospective validation of preexisting biological data.
