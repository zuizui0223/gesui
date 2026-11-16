# NRSA 2013–14 independent-cycle test of EPT–Chla riverbed moderation and true AFDM standing stock

**2026-11-14 | genuine independent EPA survey-cycle observations; confirmed source data, independent second-code numerical QA. DRAFT, deliberately preserves POSITIVES AND NULLS. This does NOT establish wastewater-derived nutrient causation or a realized ecosystem service response.**

## Source and comparability

Official [EPA National Aquatic Resource Surveys data table](https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys) lists original 2013–14 SiteInfo, WaterChem, Benthic Macroinvertebrate Count, Periphyton Chlorophyll-a, Periphyton Biomass AFDM and Physical Habitat Metrics sources. Official six EPA source URLs (2019 originals) were directly acquired with SHA256 receipts. **Original 2013–14 six table common native UID = 2,254 of 2,261 SiteInfo source UIDs**; benthic raw count n97,114 rows, including 58,667 `BERW` wadeable records and 38,447 `BETB` boatable records. Do not pool BERW/BETB. [Original source-only SUCCESS #43560419988](https://github.com/zuizui0223/gesui/actions/runs/43560419988).

Before new biological outcome magnitudes were opened, [independent six original-source unit/semantics audit SUCCESS #43560769667](https://github.com/zuizui0223/gesui/actions/runs/43560769667) froze **1,203 probability-sampled BERW native UID visits** (1,106 physical sites /754 original HUC8 /48 states). Source field units: original 2013 TN mg/L, TP µg/L (divided by 1,000 ONLY to express mg/L), *periphyton* Chla ug/cm², AFDM 2,198 original records mg/cm² and 60 original records bulk mg (not directly areal and therefore excluded from original areal AFDM models). NARS flags Q/ND/H/etc: primary strict exclude all nonblank original lab NARS flags rather than treat them as zero. 2013 SiteInfo had no published WGT_TP_CORE in this CSV: unweighted PROB-site associations are NOT national stream-mile prevalence estimates.

Official `BERW` (wadeable reach-wide benthic counting) and `TOTAL300` standardized counted-individual fraction is comparable at operational level to source original 2018–19 BERW. It is still **not absolute in-stream insect density, dry mass or measured adult emergence**. Likewise periphyton Chla is measured pigment stock not primary production; AFDM mg/cm² is standing organic mass not energy productivity. Source 2013 `XGV` is not assumed identical to 2018 source `XGVEG`, so the harmonized model only uses **six truly common original PHab variables**: XEMBED, XCDENMID, XPCAN, XFC_OHV, W1_HALL, and log10 original wetted width XWIDTH.

Exact same native UID and SITE_ID were verified across 2013 WaterChem, Chla, AFDM, PHab and original SiteInfo. WaterChem/PChla/PBIO `DATE_COL` agree with original SiteInfo dates; 2013 PHab `DATE_COL` is missing/blank on original eligible sites so its date cannot be *independently* verified, but the official original unique visit UID and physical site match. Original source is preserved, not fabricated.

## Actual separate 2013–14 survey cycle: EPT×Chla riverbed embeddedness

The earlier **2018–19** exploratory PR [#78](https://github.com/zuizui0223/gesui/pull/78) indicated a weak positive `log10 areal Chla × XEMBED` interaction on EPT counted fraction conditional on actual measured ambient TN/TP, PHab, field chemistry, spatial geography and seasonal covariates. The exact original 2013 source-eligibility criterion before new effects was `>=350` measured visits, `>=200` original HUC8 groups, `>=20` states.

Source-verified **2013–14 independent actual numeric cohort n1,116 native visits /706 HUC8 /48 states /1,028 physical sites**. Whole HUC8-disjoint fivefold Ridge α10 and HUC8-clustered descriptive OLS, +1SD logChla x +1SD XEMBED, original Chla quadratic nuisance term, true in-stream TN and TP, state, eco9, lat/lon, field season, Strahler order, conductivity and genuine common six PHab variables:

| Original independent EPA survey cycle | Native visit n | HUC8 | Embeddedness×Chla EPT-fraction partial β [HUC8 descriptive95] | HUC8 OOS extra MSE gain from interaction |
|---|---:|---:|---|---:|
| **2013–14** | **1,116** | 706 | **+0.01128 [+0.00024,+0.02232]**, exploratory p=.04530 | **+0.07346%**, fixed OOF HUC8 2000 bootstrap95 **[−0.613%,+0.790%]** |
| **2018–19 recalculated with same SIX PHab variables** | **946** | 603 | **+0.01957 [+0.00718,+0.03197]** | **+0.73305%** (still modest) |

[Actual 2013 original direct source field model SUCCESS #43561068908](https://github.com/zuizui0223/gesui/actions/runs/43561068908); [separate source-reconstituted manual NumPy harmonized 2013 vs2018 model/field QA SUCCESS #43561313134](https://github.com/zuizui0223/gesui/actions/runs/43561313134). 2018 original previous model had XGVEG also, so its earlier β+0.01921 is not strictly the same coefficient; this new β+0.01957 is the true common-six comparator.

### Critical geographic independence falsification (not optional, keep in main result)

Two EPA survey cycles share **419 of their HUC8 groups**. Of n706 HUC8 in 2013, **287 are exclusive to 2013**; of n603 in 2018, **184 are exclusive to 2018**. These are independent survey *cycles*, not fully independent sampled river basins; unweighted `n706+603` **must not be called 1,309 independent basins**. Union = **890 distinct HUC8**.

[Independent original-source exclusive HUC8 rerun SUCCESS #43561647073](https://github.com/zuizui0223/gesui/actions/runs/43561647073) found:
- **2013–14 HUC8 entirely absent from 2018 sample:** n**384** actual survey visits /287 independent HUC8/41 states: `Chla×XEMBED` β **−0.00210** HUC8 95 **[−0.02146,+0.01727]**, p=.832 and independent HUC8 heldout MSE extra **−1.192%** (worsens). This **does NOT reproduce** the original conditional slope.
- **2018–19 HUC8 entirely absent from 2013:** n229/184 HUC8 below earlier fixed independent support gate n300 for stable subpanel analysis: no effect model claimed.
- **2013 first unique physical SITE_ID only:** n1028/706 HUC8, β+0.01220 [0.00073,0.02367], extra OOS +0.172%, weak.
- **2013 BERW >=250 counted specimens:** n941/623 HUC8, β+0.00901 [−0.00321,+0.02122] and extra OOS −0.151%, uncertain.
- **2013 year alone:** 2013 n563/397 HUC8 β+0.01675 [+0.00176,+0.03174] but OOS −1.066%; 2014 n553/413 HUC8 β+0.00833 [−0.00848,+0.02514] and OOS −0.379%.

**Interpretation:** The *direction* is marginally repeated in a different field survey cycle when shared catchments are retained, but **does not generalize to catchments absent from the 2018 cohort**, nor robustly improve out-of-region prediction. No universal sediment-mediated bottom-up relationship is identified. The old positive exploratory beta from PR78 is preserved, and its independent spatial transfer failure is retained equally prominently.

## New true site-matched AFDM and community outcome

The original 2013–14 6-way native same-visit data provide **n1,094 actual areal AFDM/Chla/water chemistry/EPT/PHab visits in 693 HUC8 across 47 states** (other source 22 cases excluded due AFDM lab flag/unit/positive values). Original chemistry and periphyton original raw source units confirmed.

[2013 AFDM×XEMBED source numeric test SUCCESS #43561068908](https://github.com/zuizui0223/gesui/actions/runs/43561068908): EPT counted fraction vs `log10 AFDM × XEMBED`, with observed ambient TN/TP and common habitat controls, β **+0.00621 [−0.00599,+0.01842]**, exploratory p=.319, independent HUC8 OOS gain **−0.0230%**. This is NOT a confirmed dry organic mass moderation.

The original same-visit partial negative AFDM–EPT simple Spearman **ρ−0.1967**, Chla–AFDM **ρ+0.5728**, but nesting the original environmental and ambient nutrient covariates tells a different story.

[Actual original 2013 same-visit four nested resource models SUCCESS #43561534371](https://github.com/zuizui0223/gesui/actions/runs/43561534371):

| Additional original resource measurements beyond same baseline geography/PHab/TN/TP | EPT fraction HUC8 OOS MSE gain | Fixed heldout HUC8 2000-resample descriptive95 |
|---|---:|---:|
| + measured areal AFDM alone | +0.2130% | [−0.534%,+0.956%] |
| + measured areal Chla alone | +0.1448% | [−0.366%,+0.691%] |
| + both Chla and AFDM | +0.0417% | [−0.687%,+0.794%] |

Full simultaneous adjusted EPT counted-fraction partial slopes per +1SD source log variable:
- **AFDM β−0.01636 [−0.03708,+0.00437]** (zero-spanning).
- **Chla β+0.00237 [−0.01790,+0.02264]** (zero-spanning).
- TN β−0.01803 [−0.03994,+0.00389], TP β−0.00726 [−0.02682,+0.01231] both uncertain.

Observed N/P addition beyond geography+PHab on **log areal AFDM** HUC8 heldout MSE improvement **−0.7845%** [−1.337%,−0.270%] (degrades), and AFDM response TN β+0.03430 [−0.03142,+0.10003], TP β−0.01419 [−0.07039,+0.04200] both uncertain. For **log areal Chla** measured ambient N/P: MSE +0.454% [−0.924%,+1.771%] uncertain; TN β+0.12562 [+0.05232,+0.19891], TP β−0.02846 [−0.08525,+0.02832] source-conditional coefficients. Nutrient–Chla association is observational and not an effluent subsidy or net primary production rate.

[Independent separate manual NumPy Ridge10/OLS HUC8 test of *all same sample nested AFDM–Chla–EPT source contrasts* SUCCESS #43561707451](https://github.com/zuizui0223/gesui/actions/runs/43561707451) exactly reproduced the original same original field MSE and coefficient/interval outputs within tolerance. Not newly sampled independent ecology; independent numerical verification only.

## Scientific conclusion and next genuine data requirement

**Supported:** a large genuinely same-visit freshwater aquatic multi-biota/standing-stock macro dataset spanning two EPA cycles and 890 union HUC8, clear actual source-reported field quality, and habitat-associated conditional EPT–Chla patterns within the original 2018 source panel. **Not supported:** stable *new-catchment* generalization of the XEMBED×Chla interaction, a consistent experimental nutrient→areal AFDM→sensitive benthic insects cascade, actual wastewater-derived TN/TP dose, change in emerging insect dry mass, downstream riparian predators, pollination/water treatment or measured ecosystem service.

The historical exact WWTP effluent dual N/P source eligibility already failed in [PR #74](https://github.com/zuizui0223/gesui/pull/74). Do not re-label ambient EPA river TN/TP as sewage-specific, pool 2013 and 2018 source panels as independent watersheds, clip original units, or suppress the geographically distinct failure. Keep draft unmerged, original n946 PR78 effect preserved, no eBird.
