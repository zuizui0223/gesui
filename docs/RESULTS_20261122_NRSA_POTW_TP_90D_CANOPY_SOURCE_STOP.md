# NRSA 2023–24 × exact-upstream POTW TP × riparian canopy: 90-day source GO/NO-GO (2026-11-22)

**Decision: STOP the sewage-specific canopy × effluent phosphorus benthic interaction. This is a SOURCE-INFORMATION STOP, NOT an ecological null effect.** No biological outcome coefficients were estimated or accessed in this workflow. The original source-selected candidate population was **195 native NRSA visits** with **657 different registered upstream outfall features** (815 site-feature edges), not 195 independent field sites or verified delivered sewage exposures.

## Reproducible evidence and source receipts

- Original *full* official EPA ECHO ICIS-NPDES FY2023, FY2024 and FY2025 DMR ZIPs, exact permit + external feature + external effluent monitoring location, TP parameter `00665`, finite numeric non-NODI and no censor qualifiers, monitoring period END **strictly before** native NRSA collection day.
- Same source-frozen original USGS NLDI up-river source census from [PR #74](https://github.com/zuizui0223/gesui/pull/74), original signed source artifacts [route source](https://github.com/zuizui0223/gesui/actions/runs/43449175435) and [registered source geometry](https://github.com/zuizui0223/gesui/actions/runs/43444463485). Current 2026 registry coordinates are not proof of historical discharge delivery.
- [**Completed 4-job Actions**](https://github.com/zuizui0223/gesui/actions/runs/45013286448), aggregate signed artifact `nrsa-prior-window-crossFY-canopy-method-gate` ID **14751310520**. This run reads ~**26.288 million** + **27.018 million** + **27.011 million** original DMR rows, joins with publicly released EPA NRSA benthiccount and PhysicalHabitat original tables using exact UID/site/visit/date, and merges 3 fiscal extracts by native UID.
- Original DMR ZIP SHA256 receipts: FY2023 `6283375ce59e0694256c656adbaf0f79a455af5451ee0d70ccecfa32c8928212`, FY2024 `a59c56565f475305b7dd28807643d0b01b42f4358c902067f71daec9dfec647f`, FY2025 `b6308f4256b7aaa0df41a02dd2046f0af96035f1dadf46a7eb87cfa667da8877`. Source EPA biological table SHA: benthiccount `5290d679f7eb444e97c8d45096beb4531b450186cef190580f41f4b9430b0805`; PHab `8a6094bb87f8a57cf45a35a1272b137e6a9327cefd158f4bee88b025eaf9b2f4`.
- Original report dates are period END dates, NOT the date results were submitted to EPA and NOT evidence the river site received those nutrients. Source archives may be subsequently refreshed by ECHO.

## Unique site, year and sample-method gate

| Same original source cohort, before biological outcome access | Native visits | Independent SITE_ID | HUC8 | States |
|---|---:|---:|---:|---:|
| Same-visit bona fide canopy + benthic method (all, BERW or BETB) | 194 | 154 | 138 | 41 |
| Any unqualified prior-90d TP concentration, mixed units/bases | 99 | 78 | 70 | 30 |
| Prior-90d TP mg/L, ANY statistical base | 59 | 45 | 41 | 19 |
| **Prior-90d TP mg/L, concentration-2, monthly-average `C2/MK`** | **42** | **31** | **27** | **12** |
| Same strict TP plus **BERW + valid midchannel canopy** | **7** | **5** | **5** | **4** |
| Same strict TP plus **BETB + valid bank canopy** | **35** | **26** | **23** | **11** |

The original 195 registered-upstream population contains **155** independent physical SITE_ID; **40 sites have two visits** and **no SITE_ID repeats between 2023 and 2024**. One original survey lacks a valid benthic-method pairing, explaining 195→194. All 195 native records match original EPA benthic metadata and PHab by exact UID, date, site and visit, with zero unmatched UID, date mismatch or site mismatch; there were zero missing canopy measurements within valid protocol pairs. The new strict cohort includes all registered-upstream sites, NOT simply the old **110 visits** with positive TP reports somewhere in the same calendar year; these are different eligibility definitions and must never be compared as if nested.

**Fixed source GO conditions from** `docs/PROTOCOL_20261122_NRSA_CANOPY_POTW_TP_SOURCE_FIRST.md`: at least 40 physical sites, 25 HUC8, 5 states, 20 reporting outfalls, ≥10 physical sites in each of 2×2 exposure×canopy cells and genuine exposure contrast. **FAILED**: 31 combined sites, 26 boatable sites /23 HUC8, 5 wadeable sites /5 HUC8; separate survey methods must not be pooled. Even before a 4-cell comparison, no method stratum passes. No defensible low-exposure control group or delivered in-stream dose was established.

## Sensitivity/alternative windows

- **Prior 30 days strict C2/MK mg/L**: 38 visits /29 sites /25 HUC8; BERW 4/3, BETB 34/26, source STOP.
- **Prior 90 days strict C2/MK mg/L**: 42 visits /31 sites /27 HUC8; BERW 7/5, BETB 35/26, source STOP.
- **Prior 365 days**: 42 visits /31 sites /27 HUC8 in current 3-FY data but the 2023 survey year has **91 visits with an unobserved earliest FY2022 segment**; DO NOT interpret 365-day support as a comparable full window unless FY2022 is also retrieved. Original 195 visits sampled between 2023-04-27 and 2024-10-03.

## Additional causal/exposure limits

- Of 815 native site–upstream feature edges, **99/195 visits (50.8%)** have ≥2 upstream registered outfalls (median 2, maximum 79). The source TP concentrations cannot be summed into a physically meaningful receiving-reach TP dose; site-specific outfall flow, longitudinal river flow and retention/dilution are unknown. Exactly **657** registered upstream outfalls span the 195 visits but only a smaller reporting subset has comparable DMR.
- **Censored qualifiers and nonfinite/NODI** values are excluded from the strict monthly-mean opportunity; do not reinterpret unreported, unitless or censored values as no sewage.
- Original source profiles in matched external-outfall DMR records: FY2023 2,418 concentration-type values, 1,498 mg/L, 898 unit-blank and 22 ug/L; FY2024 2,468 concentration-type values, 1,550 mg/L, 898 blank and 20 ug/L. These are full fiscal-match profiles, NOT 90-day site counts.
- Existing 2018–19 algae×canopy analysis [PR #78](https://github.com/zuizui0223/gesui/pull/78) did **not** support canopy moderation; it is not a sewage-specific test.

## Scientific decision

Stop the **U.S. WWTP TP × canopy × benthos** biological model under this fixed protocol; reporting an interaction from 5/26 sites or selecting another outcome/window after seeing results would be pseudoreplication/selection risk. This does **not** show that canopy cannot buffer wastewater impacts. Next scientifically distinct route is a larger same-UID **actual river-water TP × field canopy × benthic community** test using EPA NRSA 2023–24 with 2018–19 regional validation, explicitly labeled *ambient nutrient–habitat moderation, not sewage-specific*, or a purpose-built local WWTP matched-control study with effluent mass, river flow, shade and before/after monitoring. Do not revive the US wastewater-specific label without a new exposure/control source.
