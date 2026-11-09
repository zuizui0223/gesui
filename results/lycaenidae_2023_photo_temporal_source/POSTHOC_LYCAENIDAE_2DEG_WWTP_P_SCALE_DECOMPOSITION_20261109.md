# Wastewater P × Lycaenidae: between-region vs within-region test (2026-11-09)

**Status: exploratory POSTHOC.** The original positive 2022 P-source vs butterfly photograph study was already known. This analysis uses immutable SHA-checked public iNaturalist photograph aggregates and annual 2022 EEA wastewater source records; no new butterfly field observations or water/substrate chemistry were measured. Registered 2022 effluent P is a **static geographic covariate** in **all three** photograph years.

## Cohort and estimands

Same **42 country×1° cells** (DE 27 / AT 15) and **22 country×2° geographic blocks**, with 2021, 2022 and 2023 photographed Lycaenidae, Pieridae and Papilionidae family counts. The response is Lycaenidae vs the other two families **photographic count log odds** with 0.5 continuity correction, not butterfly individual abundance. The predictor is log(1+2022 official reported annual P t/year) from original reported WWTP freshwater outfalls, not measured receiving river TP or biochemical exposure at an adult puddling site.

**Within-block comparison** subtracts 2° group means from outcome and predictors. Only **11/22 blocks with >1 eligible cell, 31/42 source/photo cells** identify this slope; 11 singleton blocks are not fictitious local replicates. Within-block adjusted model includes 2021 independent Plantae photo effort, 2021 WorldCover built/crop terrestrial land classes and original country-terrestrial-cell area. A lat/lon-augmented version is reported separately.

**Between-block comparison** uses **22 equal-weight 2° group means**, adjusting country, latitude/longitude, the same plant-photo effort and land variables. Within and between are different models, not directly comparable causal effects.

## New slopes per log(1+annual P t/year)

| Photo year | Between 22 2° group means | Within 11 multi-cell 2° groups (unadjusted) | Within groups adjusted |
|---|---:|---:|---:|
| 2021 | −0.0176 | −0.0232 | −0.2989 |
| **2022** | **+0.2036** | **+0.0212** | **−0.0430** |
| 2023 | +0.0828 | +0.0439 | +0.0200 |
| 2021–2023 cellwise mean | +0.0896 | +0.0140 | −0.1073 |

**2022 between-region P** = +0.2036, **4,000 country-stratified whole-region bootstrap** descriptive 95% interval **[−0.0603,+0.4598]**, leave-one-geographic-region-out signs **22/22 positive**. The bootstrap includes zero. The 2021 between-region coefficient is slightly negative, so we cannot claim a stable cross-year between-region positivity.

**2022 within-region P, land-adjusted** = −0.0430, **4,000 within-block group bootstrap** interval **[−0.2780,+0.3776]**. Source P shuffled within 2° blocks in **3,999** exploratory placements: absolute-slope add-one two-sided fraction **0.702**. This is **not a confirmatory randomized-treatment p-value**, because land/source conditional exchangeability is not established. Adding within-region lat/lon gives P **+0.0080** (95% descriptive [−0.2850,+0.4721]).

**Three-year mean within-region P, land-adjusted** = −0.1073, descriptive 95% **[−0.2847,+0.0958]**, within-block shuffled-source fraction **0.192**. Three-year mean between-region P +0.0896 (4,000 country-stratified region bootstrap 95% [−0.1554,+0.3699]). The three-year means do not supply independent treatment-time or true population replication.

**Alternative 2022 within-region source indicators**, same land controls: official N mass **−0.0799**; WWTP source count **−0.0394**. Original N and P sources are highly collinear, so no chemically distinct P-specific benefit can be identified.

## Interpretation and scientific decision

The previously reported [original 44-cell 2022 positive P-photo association PR #24](https://github.com/zuizui0223/gesui/pull/24) is retained, along with the [WorldCover geographic-prediction sensitivity PR #25](https://github.com/zuizui0223/gesui/pull/25) and [2021–2023 prior-photo source negative controls PR #27](https://github.com/zuizui0223/gesui/pull/27). This **new conditional spatial test does not recover a consistent positive same-region relationship**. The existing positive is most evident at coarse geographic photographic reporting composition; it may reflect settlement, habitat, host distributions, observation/reporting biases or pollution co-location. It is not evidence that wastewater benefits Lycaenidae populations.

The independent [real measured river-TP intermediary test PR #28](https://github.com/zuizui0223/gesui/pull/28) also failed to support the proposed source P → river TP → Lycaenidae photo-composition route at landscape-cell resolution; it does not disprove site-specific effluent effects.

**Next evidence threshold:** standardized independent local butterfly surveys, valid WWTP upstream/downstream receiving habitat contrasts, actual effluent/substrate chemistry, larval host vegetation, ants, and photo observer effort. Published 1/5-km butterfly observation support gates failed; more 1° photographic associations alone cannot identify the adult water subsidy / protein puddling / larval fitness mechanisms. Do not merge these drafts into chocho's independent manuscript.

## Reproducible source-locked receipts

- [All numeric within/between P/N/count, original source SHA256, 4,000 block bootstraps, 3,999 within-region shuffles and LOO diagnostics](https://github.com/zuizui0223/gesui/blob/analysis/lycaenidae-2023-photo-time-negativecontrol-20261108/results/lycaenidae_2023_photo_temporal_source/POSTHOC_2021_2023_LYCAENIDAE_SOURCE_P_WITHIN_BETWEEN_2DEG_AUDIT.json)
- [Independent-source replay script](https://github.com/zuizui0223/gesui/blob/analysis/lycaenidae-2023-photo-time-negativecontrol-20261108/scripts/gesui_lycaenidae_2021_2023_posthoc_within_between_p_audit.py)
- [Replayed GitHub Actions run](https://github.com/zuizui0223/gesui/actions/runs/42616805707) (SUCCESS once completed). Branch remains Draft.
