# Independent and explicitly posthoc sensitivity to spatial environmental controls, year and country

12 Nov 2026. Original adjusted source result [GitHub Actions #43227058486](https://github.com/zuizui0223/gesui/actions/runs/43227058486) is **already known**; all new diagnostics below are posthoc and cannot count as new field evidence.

Use unchanged original n454 H7×period 3-response sample, genuine verified 1219 original HydroATLAS L7 crop/urban/flow attributes, original year-specific reported WWTP N/P source and original equal basin period empirical EQR, Plant S50 and Lycaenidae. No requery, no reselect favorable sites or re-fit labels.

Must independently calculate original baseline B0 and static HydroWASTE B1 heldout MSE and **B2** (upstream cropland/use, urban/use, natural flow) → **B3** (+original HydroWASTE static) → **B4** (+original EEA reported upstream N). Compare exactly to first script for ALL THREE responses, not only EQR; fail if any original CV differs.

After reproduction, fixed posthoc subsets:
- 2017–2019 only, 2020–2022 only;
- FR only (n299; 4 positive reported annual N/P basin-periods), without FR (n155; 52 annual positive periods);
- exclude year-period rows where an original basin's station methods or source country metadata are mixed.
For each include n original H7, source-positive count, n countries, H7 hash fivefold MSE gains B2→B3 (static after upstream land), B3→B4 (year-reported N after land). Validate source positivity before claims; these checks share sites with the initial model and their direction is *not* independent replication.

Original 2018/22 EEA EQR methods differ by country, **country-out generalization and actual source attribution remain unresolved**. B2 uses original **upstream catchment** cropland and urban percent, as opposed to B2local, the source-protocol secondary local-scale sensitivity. More than one spatial source scale being examined is a limitation; neither test alone proves causation.

Preserve independent original result JSON and signed nulls, retain all original branch history, do not merge to main or chocho. No eBird.
