# Observed EPT × periphyton chlorophyll cross-biota interaction with streambed embeddedness

**2026-11-14 | Exploratory source-frozen original EPA macroecology, POSTHOC relative to known previous EPT–Chla association.** Not independently sampled data, wastewater-effluent attributable exposure, adult insect emergence, river algal biomass/productivity, or an ecosystem service effect.

## Original data and non-optional constraints

Same authentic original [EPA NRSA 2018–19 National Aquatic Resource Surveys](https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys) native `UID` visits as [gesui #77](https://github.com/zuizui0223/gesui/pull/77): genuine observed periphyton Chla µg/cm²; BERW `TOTAL300` counted EPT fraction; same original measured water TN/TP; PHab riverbed embeddedness, midchannel canopy density, riparian vegetation, canopy fraction and actual wetted channel width. **946** unique visit UIDs / **880** physical SITE_ID / **603** HUC8 catchments / **47** states, fixed **before** interaction signs. No new missing-value imputation or site selection on effect direction.

Five prespecified interaction moderators were saved to the branch in [pre-numeric protocol](https://github.com/zuizui0223/gesui/blob/analysis/nrsa1819-ept-chla-context-interactions-20261114/docs/NRSA1819_EPT_CHLA_PHAB_INTERACTION_PRENUMERIC_PROTOCOL_20261114.md): `XCDENMID` canopy density, `XEMBED` channel embeddedness, `XPCAN` canopy presence fraction, `XGVEG` riparian vegetation ground cover, and `log_width`. Every model controls for all original habitat *main* effects, logTN/TP, state, eco9, lat/lon, year, season, stream order, conductivity and an extra Chla quadratic main term. No moderator cutoffs, null controls, outcome transformations or family correction were chosen from coefficient signs.

## Main physically co-sampled interaction results

HUC8-clustered OLS descriptive interaction β per +1SD log10 periphyton Chla × +1SD physical moderator. Holm correction for five fixed tests. Genuine 5fold whole HUC8 geographic heldout Ridge α10; refits performed with train-only standardization. A positive interaction coefficient here means the original Chla–EPT negative slope **becomes less negative** where moderator is higher.

| Original real PHab moderator × log10 Chla | β interaction (HUC8 cluster descriptive95) | Holm adjusted exploratory p | Change in heldout MSE vs quadratic/no-interaction reference |
| --- | --- | ---: | ---: |
| **XEMBED streambed embeddedness** | **+0.01921 [ +0.00693,+0.03150 ]** | **0.01084** | **+0.7404%** (fixed OOF whole-HUC8 bootstrap95 [−0.695%,+2.133%], uncertain) |
| XCDENMID midchannel canopy density | −0.00043 [−0.01489,+0.01403] | 1.00 | −0.479% (worse) |
| XPCAN canopy presence fraction | +0.00109 [−0.01367,+0.01585] | 1.00 | −1.058% (worse) |
| XGVEG riparian vegetated ground index | −0.00241 [−0.01494,+0.01012] | 1.00 | −0.218% (worse) |
| log10 XWIDTH actual wetted width | +0.00292 [−0.01065,+0.01648] | 1.00 | −0.995% (worse) |

All **five** interactions together worsen the quadratic reference heldout MSE **−1.135%** with group resample95 **[−2.697%,+0.444%]**. Adding the quadratic Chla main term itself to the original linear PR77 model slightly worsens heldout MSE by −0.289%. Do **not** treat many posthoc modifiers as independently confirmed ecosystem rules or replace original null results.

### Interpretation scale and actual environment

For the specific **single XEMBED interaction model**, using the authentic source population embeddedness quartiles and evaluating the conditional Chla slope around its source mean, the estimated EPT fraction change per 1SD higher log10 Chla is:

- At **lower embeddedness, 25th percentile XEMBED ~46.15%**: **−0.03761** EPT fraction (−3.76 percentage points).
- At **higher embeddedness, 75th percentile ~89.64%**: **−0.00507** (−0.51 percentage points).

The fraction is from standardized counted specimens, not absolute insect density, mass, aquatic insect flux, terrestrial predators, water self-purification or realized pollination services. This is a difference in **conditional observational slopes**, NOT the response of the same stream to adding algae or reducing sediment. Source original XCDENMID one value >100 and source XGVEG 22 values >1 are not clipped in primary and have explicit sensitivity tests.

## Entirely independent matrix and sensitivity check

Original authentic original biological river source/OLS/HUC8 outcome workflow [SUCCESS #43537987013](https://github.com/zuizui0223/gesui/actions/runs/43537987013). [Separate NumPy normal equation, HUC8-sandwich covariance, all five moderator β/95% and manually recomputed all heldout MSEs, plus 7 posthoc fixed-cohort sensitivities SUCCESS #43538310290](https://github.com/zuizui0223/gesui/actions/runs/43538310290). The code agrees within specified numerical tolerances, but uses the **same original EPA field observations**, not independent data.

| Original-source fixed sensitivity | n visits /HUC8 | Embeddedness×Chla β (cluster descriptive95) | Embeddedness interaction extra MSE gain vs curved Chla reference |
| --- | ---: | --- | ---: |
| Full fixed primary | 946 /603 | +0.01921 [+0.00693,+0.03150] | +0.740% |
| Unique physical site first visit only | 880 /603 | +0.01997 [+0.00745,+0.03249] | +0.984% |
| Original `TOTAL300` sample counts ≥250 | 775 /509 | +0.02231 [+0.00836,+0.03625] | +0.740% |
| Exclude source XCDENMID >100 entry | 945 /602 | +0.01913 [+0.00686,+0.03139] | +0.727% |
| Exclude original XGVEG >1 records | 924 /596 | +0.01985 [+0.00743,+0.03226] | +0.820% |
| 2018 native visits only | 487 /337 | +0.02052 [+0.00507,+0.03596] | +0.853% |
| 2019 native visits only | 459 /333 | +0.01761 [−0.00283,+0.03806] (uncertain) | −0.288% (worse) |

Thus duplicate visits, a small counted sample or single abnormal canopy-density value are unlikely to explain the *full-sample* interaction. However 2019 alone is imprecise, predictive benefit disappears and **fixed OOF group-bootstrap interval for full data overlaps zero**. Different years are not independent randomized interventions, and this is hypothesis generation not confirmed context-dependent mechanism.

## Ecological discussion, no causal overreach

One candidate explanation is that more embedded riverbeds constrain EPT assemblage composition such that Chla standing pigment supplies less additional information about EPT share, whereas where substrate is more open the negative cross-biota co-gradient is clearer. This could also arise from low-EPT floor effects, selective replacement by non-EPT taxa (e.g. chironomids), substrate/flow changes, light, contaminants, seasonal sample design or unmeasured nutrient inputs. **Do not claim reduced sediment directly changes EPT or algae.** A straightforward follow-up is to decompose relative chironomid vs EPT counts and their original taxon-ID richness *on the same exact original source visits*, checking whether apparent interaction is compositional rather than an absolute abundance/biomass phenomenon.

True municipal sewage effluent→river N/P→food-web causality remains unestablished because [PR #74](https://github.com/zuizui0223/gesui/pull/74) reported actual historical WWTP discharge source N/P dual-data coverage below precommitted eligibility, and this EPA 2018–19 panel is ambient stream TN/TP. 2023–24 original benthic publisher CSV links remain inaccessible from [PR #75](https://github.com/zuizui0223/gesui/pull/75); official NARS download table still lists them, but WQP4 exact native site-ID probes did not prove direct mapping. Nothing here measures emergence biomass or pollination.
