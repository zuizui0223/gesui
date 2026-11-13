# Genuine same-visit 2018–19 NRSA: River nutrients, periphyton pigment, benthic EPT and taxon-ID richness

**2026-11-14 | Draft, original data source-signed, exploratory posthoc observational macroecology.**

This study directly addresses whether **measured water-quality, algae pigment standing stock and benthic invertebrate community composition can be aligned at exactly the same original EPA survey visits**. The outcomes do NOT supply demonstrated wastewater-derived N/P, actual aquatic adult emergence or ecosystem services. The earlier [PR #74](https://github.com/zuizui0223/gesui/pull/74) source-specific WWTP source eligibility remains failed.

## Source-only native-UID census before new biological effects

- Authenticated, SHA-pinned EPA NRSA **2018–19 original SiteInfo, Water Chemistry, Benthic Macroinvertebrate Count (BERW), Periphyton Chlorophyll-a and Physical Habitat Larger Metrics** CSV files, not a literature aggregation.
- Original source-fixed n**1,221 native survey UIDs**, 1,125 distinct physical survey sites, 729 HUC8 independent watershed groups and 48 states. These all have the original source *table opportunity*, not necessarily valid quantitative responses. [Source-only original EPA SHA and identifier gate SUCCESS #43528570410](https://github.com/zuizui0223/gesui/actions/runs/43528570410). Original case `XGveg` was mapped using case-insensitive header normalization before opening physical covariate values, not posthoc fit tuning.
- The prior original [PR #73](https://github.com/zuizui0223/gesui/pull/73) stricter biological source/effort cohort has exactly **947 UID** cases /603 HUC8 /47 states with valid areal periphyton Chla and actual TN/TP and benthic EPT counts. Adding same-visit PHab leaves **946 visits, 603 HUC8, 47 states, 880 distinct physical SITE_ID**. This population and its two prespecified PHab columns sets were finalized before the cross-layer numerical effects.
- Original source metadata: periphyton Chla (ug/cm²), real river TN (mg N/L), TP (original µg/L correctly converted to mg/L), benthic BERW original fixed `TOTAL300` counted-subsample fraction EPT, distinct accepted original `TAXA_ID` keys with positive `TOTAL300` as *taxon-ID richness* (NOT species richness), not absolute density or adult insect emergence.
- **PHab broad**: XPCAN riparian canopy fraction, XGVEG ground vegetation cover index, XFC_OHV overhanging vegetation areal proportion and W1_HALL human disturbance. **PHab dense**: these + measured channel XWIDTH, mid-channel XCDENMID canopy density and XEMBED streambed embeddedness. Both tiers remain n946/603 HUC8; no opportunistic missing-data selection. 22 nonmissing raw XGVEG rows exceeded numeric 1; because this overlaps vegetation layers, values were NOT artificially clipped to [0,1].
- All cohort IDs, original CSV SHA256 signatures and source-failure/success receipts are in GitHub Actions artifacts. Original per-visit biological data and coordinates are not committed to the repository.

## Physically same survey visit, 946 NRSA UIDs / 603 HUC8 / 47 states: exploratory biological associations

Each endpoint has its OWN model with baseline: original state, eco9, geography, stream order, year, sampling season, river conductivity. The 5fold HUC8 hash grouped Ridge α10 ensures no site from the same watershed appears in training and validation. The two nested predictor groups: measured in-stream logTN+logTP and actual field PHab, then the independently observed other biological variable as extra covariate. All response-specific effect uncertainty is descriptive HUC8 clustered rather than experimental causality.

### Broad vs dense physical controls

| Observed response | TN+TP MSE improvement AFTER broad PHab | TN+TP MSE improvement AFTER dense PHab | Dense PHab improvement AFTER measured nutrients |
|---|---:|---:|---:|
| EPT original benthic counted fraction | +1.5363% | **+1.0835%** | **+11.1375%** |
| Periphyton areal log10 chlorophyll-a | −0.7772% (worsens) | **−0.9389%** (worsens) | **+3.3136%** |
| ln(1+native TAXA_ID richness) | +1.4820% | +1.4471% | +3.2599% |

All are *predictive* differences relative to different response-specific reference models, NOT comparable biological magnitudes.

**Raw cross-biota Spearman** (same n946 visits):
- Benthic EPT counted fraction vs Chla log: **ρ−0.04857** — weak, no claim of strong universal cross-biota co-response.
- Benthic unique taxon-ID richness vs Chla log: **ρ+0.07128** — weak.
- EPT counted fraction vs taxon-ID richness: **ρ+0.28936** — correlated but not interchangeable.

**Dense PHab fully adjusted observational partial slopes** (all compared per +1 SD of the named original predictor):
- **EPT fraction vs other measured log10 Chla: β−0.019864** (HUC8 descriptive95 **[−0.033648,−0.006080]**) after physical and TN/TP. Adding Chla to the geographic EPT predictor after nutrients + PHab produces **+1.0857%** MSE gain, descriptive HUC8 out-of-fold bootstrap95 **[−0.0507%,+2.2629%]**. So *directional partial relationship but independently weak prediction*, not demonstrated biomass→EPT mechanism.
- **Log10 areal Chla vs measured EPT fraction:** β−0.05642 **[−0.09504,−0.01779]**. EPT extra heldout prediction +1.0990%, bootstrap95 **[−0.0443%,+2.3086%]**. This is a symmetric observational co-state check, NOT independent replication or directionality.
- **EPT fraction vs TN:** β−0.03402 **[−0.05844,−0.00960]** after dense PHab/TP/Chla. TP β−0.00054 [−0.02232,+0.02124] includes zero.
- **ln taxon-ID richness vs TP:** β−0.03562 **[−0.06536,−0.00588]** after dense PHab/TN/Chla and original counted-individual effort control. TN β−0.02721 [−0.06157,+0.00716] includes zero.
- Extra TN+TP AFTER dense PHab EPT fraction **+1.0835%** HUC8 descriptive fixed heldout bootstrap95 **[−0.4293%,+2.6530%]**. After PHab, the Chla extra nutrient prediction **−0.9389%** [−1.4357%,−0.4546%] worsens. These do not prove absence of nutrient effects, only poor extra prediction for this particular outcome/model.
- Even after exposure/physical stratification, ecological sites and organism detection processes are not randomized, and the measured habitat itself may be a mediator or collider. Do not interpret PHab adjustment as proving nutrients biologically inert.

## Numerical and source audits

- Original genuine source/UID availability [SUCCESS #43528570410](https://github.com/zuizui0223/gesui/actions/runs/43528570410).
- Main raw EPA same-visit nutrient×benthos×Chla×PHab 2-tier/3-response geographic OOS models [SUCCESS #43528870509](https://github.com/zuizui0223/gesui/actions/runs/43528870509).
- **Independent manual NumPy Ridge + OLS HUC8 sandwich, all 30 heldout models, all partial slopes/intervals and the fixed whole-HUC8 bootstrap source cohort exactly agree with the original code** [SUCCESS #43529245688](https://github.com/zuizui0223/gesui/actions/runs/43529245688). It is a code audit on the same field data, not independent ecological observation.

## What can and cannot be concluded for sewage and services

**Supported as exploratory macroecology:** genuine same-river visit TN/TP and multiple original freshwater biotic states; strong added information of physical habitat for EPT and Chla, weaker cross-biota nutrient effects; negative conditional EPT–Chla co-occurrence signal but not strong predictive gain. Native taxon-ID richness rather than photos makes the biodiversity endpoint closer to standardized field surveys.

**NOT demonstrated:** wastewater caused richness loss; wastewater raised or lowered biofilm biomass; P caused functional group changes; measured EPT fraction represents insect production/emergence or actual bird/bat/spider food subsidy; chlorophyll-a equals primary productivity; butterfly observation equals pollination service. There is no verified at-site year-specific WWTP exposure sufficient for multi-trophic effect identification in [PR #74](https://github.com/zuizui0223/gesui/pull/74). This is a 2018–19 source cohort distinct from the 2023–24 original areal AFDM [PR #75/#76](https://github.com/zuizui0223/gesui/pull/76) and must not be synthetically merged as a temporal causal chain.

## Next source work

The EPA 2023–24 original Count/MMI/Taxa files remained unreachable (HTTP 404) from publisher links in [PR #75](https://github.com/zuizui0223/gesui/pull/75). The original NRSA 2023/24 QAPP states biological data will eventually enter WQX / Water Quality Portal; availability there must be tested by true monitoring-location and visit keys, not assumed. The same 2018 n946 result stands whether or not the 2023 original website recovers.

Keep PR unmerged/DRAFT, no eBird or publication outside user-requested repository changes.
