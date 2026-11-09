# Source-frozen 2022 WISE6 *same-river upstream–downstream* TP contrasts: eight source versus eight registered-source-zero Germany quartets

**Developed 9 November 2026, POSTHOC to the earlier broad-scale photographic Lycaenidae P-signal and the 168-reach WWTP-P environmental prediction.** This is a different nationwide original 2022 river-monitoring cohort, **not** selected from measured TP concentrations. The 8:8 quartet allocation was committed separately *before* inspecting these 32 station TP numeric outcomes, using original source-only river graph and WWTP outfall metadata. This is not a causal pollution experiment and does not test Lycaenidae butterflies.

## Original independent source-only comparator

The original Germany source-only selection is:
`results/wise6_2022_national_DEAT_TP_sourceonly/GERMANY_2022_WISE6_TP_8SOURCE_8ZERO_SAME_RIVER_MONITOR_QUARTETS_SOURCE_ONLY_LOCK.json`
on SHA-pinned Draft PR #30 source branch `analysis/wise6-2022-national-DEAT-river-TP-pairs-20261109`. Exactly **8 paired quartets × 4 unique original WISE6 2022 station IDs = 32 distinct sites**, of which **8 pairs** have at least one qualifying EEA 2022 both-positive measured N/P original one-freshwater-outfall on an intervening strictly directed river reach (including any tributary that flows into the traced intervening path), and **8 matched no-qualifying-registered-intervening-source controls**. The source comparison treats each quartet as one matched unit: one WWTP-interposed upstream/downstream pair versus one original source-zero upstream/downstream pair. Of eight quartets, 4 are matched in the same country×2° upstream block, 4 elsewhere; source matching was selected solely by original river path distance and country/region metadata. Only German 2022 river monitors are used here (not a DE/AT general result).

**Pre-exposure source facts:** physical HydroRIVERS v1 converted third-party GeoParquet directed topology; upstream/downstream WISE monitoring reaches strictly different original reach IDs connected by original NEXT_DOWN, directed path length <=50km. WWTP outfall strictly intervening on a traced reach/tributary between both original sample reach IDs, not necessarily on the main channel and not necessarily proven to alter sampled TP. **No same-reach outlet order guessed.** Original report of downstream-mixing network is third-party converted mirror and not confirmed against official HydroRIVERS. Group allocation does not use TP concentration, Lycaenidae photos or observer effort.

## Numeric TP is opened after this fixed source selection

Original EEA WISE6 DisaggregatedData 2026 original ZIP SHA256 `2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2`, original **102,569,057** individual public sample rows, water matrix W, river category RW, original Total phosphorus CAS_7723-14-0, original unit **mg{P}/L**, original 2022 dates. Stream 2GB compressed/36GB uncompressed original CSV without local full extraction. Filter **only the original 32 hashed site identifiers preselected**, never add stations with favorable chemistry; 13 of 32 were not part of the old 187 source/photo sites, so a new original original-source numeric read is required for complete quartets.

For each station:
1. Original `resultObservedValue` must be finite >=0 in original mg P/L. BelowLOQ TRUE or 1 **excluded** (no 0/halfLOQ substitutions); unknown flags flagged but may retain an explicitly finite measured result, other flag codes dropped. Multiple valid assays on same day collapsed to **median day TP**.
2. Qualify original station only with ≥**4 valid sampling dates in ≥3 calendar months** during 2022 and belowLOQ record fraction ≤**0.50**, otherwise station ineligible. Site-year TP is median of daily medians. Keep original station SHA; never publish individual sample IDs or full source ZIP.
3. A source-positive or source-zero upstream/downstream pair qualifies only when **both** source-locked station TP medians pass. To compare seasonally similar source sampling, require at least **2 calendar months with original quantified TP observations in common** between upstream and downstream. A quartet qualifies only if BOTH matched pairs qualify.
4. **After opening numeric validity but before testing any difference magnitude**, require **≥6 complete quartets**, each involving 4 original distinct WISE site IDs, across **≥4 independent source-defined country×2° geographic blocks**. Otherwise STOP WITHOUT estimating a source-interposition effect. Do not replace missing sites/quartets or relax LOQ/month cutoffs.
5. Report all 8 originally source-selected quartets and original 32 station quality flags and reason for missingness. Source-registered-zero is not a sewage-free stream.

## Fixed observational river difference contrast conditional on the above

Original station endpoint is annual 2022 measured TP in mg{P}/L. For each retained quartet, define source pair gradient:

`delta_source = log(TP_source_downstream + 0.001) - log(TP_source_upstream + 0.001)`

and registered-source-zero pair gradient:

`delta_control = log(TP_zero_downstream + 0.001) - log(TP_zero_upstream + 0.001)`.

The matched source-vs-control quartet contrast is **delta_source minus delta_control**, so positive means a larger *relative* upstream-to-downstream measured TP gradient in the pair with a reported qualifying WWTP introduced between stations, NOT that effluent caused it. The small +0.001 mg P/L constant is fixed *before* opening these new 32 original station numeric values; use only that one log transform, no threshold tuning.

Primary summary: **one equal-weight value per valid quartet**, mean and median differences, n positive/negative, median gradients by source arm. Exact two-sided matched sign-flip enumeration on retained quartets (up to 2^8=256 sign configurations) around zero is shown as a **descriptive randomization-style symmetry reference**, **NOT** a randomized WWTP p-value (source allocation is not randomized). Fixed seeded 9,999 quartet bootstrap replications for a descriptive 95% interval of the mean; country is DE only and block support may be 4–7 regions. Also report source reported annual original P mass and river path distance for context, but do not infer actual P concentration delivered.

**Main confounds:** Tributary dilution, river geometry/mixing, upstream land-use/industry and seasonally unmatched daily samples, even though both site-years are same 2022; source registry missingness and annual P mass reported; original network mirror not official topology-verified; site locations/LOQ and original expert `resultObservationStatus` quality flags not fully audited, and no direct source-speciation/protein/Na measures. A positive contrast indicates a new source-labeled observational TP change candidate, not a proven point-source causal impact. A negative contrast or gate failure is equally reportable.

**No butterfly species or actual Lycaenidae population/behavior data accessed for this study.** Previous 168-reach geographic upstream P prediction −1.88% and old coarse 1° Lycaenidae photo-reporting association are separate estimands. Historical 2018/2022 river EQR and `chocho` 22/25 potential introduced larval host-resource findings stay unchanged.
