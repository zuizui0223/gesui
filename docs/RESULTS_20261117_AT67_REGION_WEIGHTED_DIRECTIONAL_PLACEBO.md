# 2026-11-17 | Austria67: genuine 2022 aquatic source signal under region weights and source-direction placebo stress

**Source-signed artifact verification + new posthoc reweighting of the ORIGINAL five region-heldout MSEs, not a new fitted biological model or independent validation dataset.** This corrects overstatement of how strongly the downstream-source placebo proves causal direction.

## Original source provenance, physical support and units

- [Original Austria67 official same river 2022 TP / EEA benthic EQR × upstream annual WWTP N/P vs P source #43890958636](https://github.com/zuizui0223/gesui/actions/runs/43890958636), artifact ID `14237141308`; actual ZIP SHA256 verified **f92c7eb204bf86b83a4da9a1ec2cc37dd199a1eb324cdb176f175a5a6d2d732f**.
- [Original 2021 WorldCover / source river landcover #43888461430](https://github.com/zuizui0223/gesui/actions/runs/43888461430), artifact ID `14236375421`; actual ZIP SHA256 verified **09246a4b11cd638fcc5382c00cf6ced9b97668709b9374fd87f827102074837d**.
- Real 2022 67 official Austrian **water monitoring sites** and genuine river annual TP / EQR. Original five 3° geographic withheld blocks: 15_3 n9, 15_4 n14, 15_5 n21, 16_4 n13, 16_5 n10. 2021 WorldCover 1km is a landscape **covariate**, not 2022 riverbank vegetation change.
- The EEA annual facility-reported upstream P and N tonnes/year / mean river flow (Q) are source/flow potential concentration-equivalent predictors, NOT measured source fractions actually received at water or banks. Same-source upstream N/Q and P/Q Spearman **+0.94527**. Annual reported loads, ambient river TP, modeled W and basin EIV remain different environmental variables.

## Same original region-purged validation errors, arithmetic only

Original five-region error-equal source comparison `100×(1−ΣMSE_added/ΣMSE_base)`; NEW source-site-weighted error aggregation `100×(1−Σ(n_region × MSE_added)/Σ(n_region × MSE_base))`. **No regional model refit, sample filtering, new region choice, biological outcome change or permutation p-values.**

| Predictor added after original WorldCover + actual river lnTP unless specified | Original five-region equal MSE gain | NEW 67-site-weighted MSE gain | Site-weighted gain range removing any one region's ORIGINAL heldout errors |
|---|---:|---:|---:|
| Eligible upstream reported **P/Q** | **+8.388%** | **+9.263%** | **+5.683% to +11.731%** |
| Eligible upstream reported **N/Q** | **+12.199%** | **+13.332%** | **+7.914% to +17.635%** |
| Eligible upstream **P/Q added after binary source-presence** | **+10.440%** | **+11.725%** | **+8.438% to +12.835%** |
| Binary eligible upstream-source reporting presence alone | **−3.418%** | **−3.886%** | **−5.000% to −0.107%** |
| **Downstream P/Q negative-control geography** | **−13.208%** | **−9.912%** | **−14.335% to −0.551%** |
| Downstream N/Q negative-control geography | **−2.551%** | **−2.021%** | **−3.034% to −1.110%** |

The original **P/Q gain is positive in all five region-specific holdouts** (+7.46, +1.77, +19.83, +1.04, +11.52%). **N/Q is 4/5 positive**: region `16_4` has **−20.20%** even though the region-equal and site-weighted source results remain positive. N and P source proxies are highly collinear, and adding P/Q **after N/Q** had original MSE gain **−1.317%**; N/Q after P/Q **+2.898%**. Thus one cannot identify nutrient-specific water treatment.

### Novel caution from source reweighting: the DOWNSTREAM P/Q placebo is fragile

Original downstream reported P/Q placebo worsens EQR prediction by **−13.208%** across five region-equal heldout MSEs, largely due to extreme deterioration **−47.350% in region 16_5**. Omitting **16_5** from the same heldout-errors table changes the *site-weighted* downstream deterioration to only **−0.551%** (region-equal **−1.654%**). Other regions show smaller or positive placebo improvements. This does **not** invalidate the upstream P/Q predictive signal, but it **weakens the simple narrative** that a dramatic downstream-negative placebo by itself proves true hydrological wastewater direction or sewage causality.

## Manuscript implications and gates

- Positive conservative **aquatic result**: original 2022 facility reported source/flow predictor provides spatially transportable incremental EQR prediction beyond *local* landcover, actual river annual TP and geographic background within these 67 stations. Sensitivity to original region weighting is limited for upstream P/Q, even when any original region is omitted from error aggregation.
- Important **source attribution weakness**: external N/Q and P/Q predictors covary; no independently measured receiving-site WWTP fraction, N/P separation, full upstream 50km land use, source intervention or pre/post change. Five region heldouts are too few to treat descriptive fixed-error intervals as robust independent causal inference. Original reverse-source placebo asymmetry is region-sensitive.
- No concurrent riverbank plant biomass, vegetation trait composition, leaf tissue N/P, aquatic emergence→spider uptake, or butterfly puddling/life-stage outcome. The separate near/outer S2 vegetation support 5/10<prelocked6/10 remains technical STOP; original 1km WorldCover **cannot** substitute.
- In a single manuscript, retain authentic aquatic evidence as a substantial within-water panel and terrestrial botanical–butterfly trait matching as a separate region/population panel. No `effluent→river→plants→butterflies` SEM or ecological-service positive/negative effect.

Pure standard-library independent aggregate arithmetic is provided in `scripts/verify_at67_region_transfer.py`. Download the two original ZIPs from the exact source runs above; the original government datasets and raster assets are **not** copied to the repo. New JSON receipt is in `results/20261117_at67_region_transfer/`.
