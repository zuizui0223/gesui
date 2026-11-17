# 2026-11-17 | 原本ZIP SHA照合：厳格162種の国別一般化、陸上水分形質 vs 水域EQR

**Status: REAL original signed artifact verified + NEW posthoc arithmetic on the ORIGINAL six country-heldout MSE tables. NO new fitted ecological model, new field observations, independent replication, significance test, or causal conclusion.**

## Provenance
- Source: PR #119 [GitHub Actions #43901046421](https://github.com/zuizui0223/gesui/actions/runs/43901046421), original artifact ID **14240413667**, ZIP name `gesui-20261117-H7-original162-same-HOST-N-F-plant-moisture-real-aquaticEQ`; original ZIP **SHA256 a22ea230c8e6beb4533bc696504194879927e9fd6a553377fba5255ce1464bc2** was freshly rechecked from the downloaded real ZIP, *not* copied from PR prose.
- Internal member `GESUI_STRICT_SAME_HOST_NF_H7_PLANTF_EQR_CROSSDOMAIN_RESULTS.json`; source model schema `GESUI_20261117_H7_BUTTERFLY_SAME_ACTUAL_FERRERPARIS_LARVALHOST_EIV_NF_MATCHED_PLANT_AND_WATER_EQR_V1`.
- True complete strict source: **162** published butterfly species with numeric N/F scores derived on the exact same documented 2014 larval-host plant taxa; **447 genuine EU H7 × period**, **89 H6 parents**, 11 source country labels; original 2017–19/2020–22 source composition, not butterfly population abundance or measured feeding.
- Country-holdout analysis has an **eligible six-country subset: 424 H7×period / 81 H6** (FR 298, HR 16, IT 48, LT 30, NL 12, NO 20). Therefore these reported six-country gains do **not** cover all 447/89 source H7/H6 or all 11 source country labels. **FR accounts for 298/424=70.3%** of the six-country site-period observations.
- Original stage used H6-purged country withholding to fit its models; the present step only recomputes aggregations of their already-saved country MSEs. Country holdouts share some training sites, so they are NOT six independent studies or inferential units for a meta-analysis.

## Audited country-level genuine heldout MSE improvements

`gain = 100 × (1 − MSE_added / MSE_baseline)`, using **the same existing original heldout predictions** for each comparison.

| Fully withheld country | n H7×period | Plant-community MOIST added after plant N / butterfly host N | Benthic EQR added after both plant traits / host N |
|---|---:|---:|---:|
| France FR | 298 | +20.618% | −1.930% |
| Croatia HR | 16 | +26.867% | +0.805% |
| Italy IT | 48 | +45.872% | +5.364% |
| Lithuania LT | 30 | +39.579% | +1.492% |
| Netherlands NL | 12 | +61.625% | −3.504% |
| Norway NO | 20 | +13.493% | −11.435% |
| **Country-equal ratio of MSE** (original PR119 independently reproduced) | 424 | **+34.087%** | **−0.740%** |
| **H7×period-count weighted ratio of MSE** (NEW reaggregation) | 424 | **+28.796%** | **−0.985%** |
| **H7×period weighted, drop FR from already withheld country-error table** (NEW reaggregation) | 126 | **+36.783%** | **+0.174%** |

The new H7×period weighting is *not* a rerun of a 424-point pooled validation model: it pools the **previously withheld-country prediction errors** with site-period weights. The original country-equal ratio is `1 - mean(MSE_added) / mean(MSE_baseline)`, not the arithmetic average of the six country-specific percentage gains (which would be plant +34.676% and EQR −1.535%). The distinction is preserved explicitly.

### Single-country omission diagnostic (posthoc; **no retraining**)

When each original country set of prediction errors is removed from the six-country aggregation, the **H7×period-weighted** plant F gain ranges **+24.335%** (drop LT) to **+36.783%** (drop FR). Corresponding EQR gain ranges **−1.804%** (drop LT) to **+0.174%** (drop FR). This is a *sensitivity calculation only*, not a resampling CI. It shows that the existing cross-country predictive asymmetry is not a peculiarity of assigning equal weights to tiny NL vs large FR, yet makes no claim of independent replication.

### Temporal falsifier still fails

Original strictly matched 216 repeated H7 periods: `Spearman(Δbutterfly-host MOIST, Δplant MOIST) = +0.05016`; `Spearman(Δbutterfly-host MOIST, Δaquatic EQR) = −0.05938`. The cross-country positive plant-habitat correspondence is mainly spatial: **no synchronized change** was demonstrated between the two three-year windows.

## Ecological interpretation and next falsifiers

1. **Supported observational pattern:** terrestrial community moisture-habitat affinity tracks butterfly community's documented larval-host moisture habitat affinity across eligible countries, even after plant N and butterfly host N controls. The **plant** and **larval host** EIV scores draw from the SAME European habitat-indicator reference, so shared genus composition, spatial climate, sampling and functional habitat preferences remain alternative explanations to any local realized host use.
2. **Unsupported pathway:** benthic EQR and static wastewater geography add essentially no consistent extra **country-transfer** prediction. This is not evidence of *zero physical river-to-land nutrient transfer*, an ecological boundary, sewage deterrence, or a sewage causal effect.
3. **Urgent next independent falsifier:** genus-aware source-first trait randomization / between- vs within-genus decomposition for **plant moisture** (analogous to existing PR #105 **nutrient**-trait genus check), and common-H6 support check **before** any local model. No claim of within-genus F ecological filtering until scored host species and source GBIF site occurrence matrices are genuinely reconstructed. The PR #105 N result was largely between-genus and its conclusion must not be silently transplanted to F.
4. **Water-specific falsifier still uncollected:** same 2022 outfall-direction matched source N/P and local receiver TP/EQR plus same-site standardized **riparian** biota/soil or emergence/spider consumption. The Austria67 S2 natural-land 5/10 source STOP and 16-outfall local-GBIF 250m 0/16 STOP remain in place.
5. **Publication gate:** retain primary aquatic observational source/condition signal, use terrestrial trait matching as a separate compartment with exact domain/temporal sample declarations; **no** `WWTP→river TP→plants→butterflies` path model or ecological-service effect from these data.

### Reproduction
A pure-Python standard-library recomputation is in `scripts/verify_pr119_country_transfer.py`. It reads the original signed ZIP, requires its exact published SHA and expected schema, verifies every country's original MSE gain and the published equal-country aggregate, then reports the new H7-weighted and country-exclusion sensitivity. The original third-party source/data are **not rehosted**; download artifact ID 14240413667 from original Actions. A generated local JSON result can be treated as a derived audit receipt, not a new source.
