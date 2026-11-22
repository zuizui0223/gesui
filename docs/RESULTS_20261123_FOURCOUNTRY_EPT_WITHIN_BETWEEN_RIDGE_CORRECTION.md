# Four-country original fine sediment–EPT transfer: numerical re-evaluation and within-vs-between-country decomposition
**2026-11-23 JST | Signed source-model audit (POSTHOC), not preregistered confirmatory ecological hypothesis.**

## Main correction to prior PR #146

Our previous independent four-country *descriptive* cross-country holdout of original normalized EPT relative fraction against original visual <2 mm fine sediment cover reported only **+1.45%–+2.40%** MSE improvement at Ridge `alpha=10`. **This magnitude reflected extremely strong shrinkage**: each of the three training countries was assigned equal *total* sample weight 1, yielding train total weight **3**, while slope penalty was 10 on raw fines fraction `s ∈ [0,1]`. Therefore prior descriptions of only a **weak ecological relationship because improvement was small** were unwarranted. The original fourcountry source, source taxonomy harmonization, domain definition, response and heldout countries are all **unchanged**. Under an explicitly fixed POSTHOC sensitivity grid, unpenalized `alpha=0` gives much larger, still only **descriptive** heldout MSE improvement.

### Exact original 4-country pointwise holdout MSE gains (%)

| Whole heldout target country | Original samples | Unpenalized α=0 | α=1 | **Original α=10** | α=100 |
|---|---:|---:|---:|---:|---:|
| Australia | 3099 | **+25.354** | +12.233 | **+1.814** | +0.190 |
| Brazil | 232 | **+29.614** | +16.155 | **+2.396** | +0.251 |
| New Zealand | 342 | **+27.643** | +10.241 | **+1.448** | +0.151 |
| UK | 2818 | **+33.962** | +12.460 | **+1.785** | +0.186 |

The same four-country-holdout source, original normalized outcome and training-country-equal weights apply to every column. The alpha grid was frozen in [posthoc protocol](https://github.com/zuizui0223/gesui/blob/analysis/external-fourcountry-within-between-transfer-audit-20261123/docs/PROTOCOL_20261123_FOURCOUNTRY_WITHIN_BETWEEN_TRANSPORT.md) **before running the new sensitivities**, but this was designed after seeing the α=10 gains; **do not select α=0 afterward and call it a preregistered result or publish a new universal effect size**.

### Common measured fine-cover 0–50% range sensitivity

| Holdout country | Source valid n after 0–50% restriction | α=0 MSE gain | α=10 MSE gain |
|---|---:|---:|---:|
| Australia | 2300 | +15.207% | +0.271% |
| Brazil | 101 | +38.864% | +0.334% |
| New Zealand | 213 | +5.654% | +0.163% |
| UK | 2575 | +23.108% | +0.276% |

The common-domain restriction simultaneously **changes the source sample and available fine-cover variability**, so smaller original α=10 gains there cannot prove absence of an ecological fine-sediment response. Brazil has only 101 sampled analyses in the restricted domain and no independent site count, so its +38.864% **is not a reliable independently replicated national site effect**.

## Source ecology vs country composition — exact algebraic gain audit

The original `M_pool` model uses *pooled* covariance with one training-country vote. This covariance can be decomposed:
`Cov_train(s,y) = mean_country Cov_within_country(s,y) + Cov_across_3_country_means(s_bar_c,y_bar_c)`, and similarly for sediment variance.

For any heldout country and fixed training beta, the heldout gain versus the training grand-mean predictor is **exactly**
`MSE0 - MSE1 = 2 beta Cov_target(s,y) + 2 beta (mean_s_target - mean_s_train)(mean_y_target - mean_y_train) - beta^2 [Var_target(s) + (mean_s_target - mean_s_train)^2]`.

**`target` EPT means appear only in this AFTER-THE-FACT decomposition, never in training or predictions.** Each contribution below is normalized to the original heldout baseline MSE and is measured as percentage points of total **MSE improvement**, not as EPT changes or ecological regression slopes.

| Heldout target | **Within target-country covariance contribution** | **Country-mean-offset contribution** | Predictor penalty | Reproduced original α=10 gain |
|---|---:|---:|---:|---:|
| Australia | **+1.573** | +0.275 | −0.033 | +1.814% |
| Brazil | **+0.067** | **+2.376** | −0.048 | +2.396% |
| New Zealand | **+1.420** | +0.047 | −0.018 | +1.448% |
| UK | +0.780 | **+1.026** | −0.021 | +1.785% |

**Inference from exact predictive mechanics**: Australian and New Zealand full-domain gains are dominated by real **within-country observed covariance** (which was negative for all four target countries). Brazil's positive heldout gain under alpha10 is *overwhelmingly* due to **country mean offset**, and UK's is more than half driven by country means. Therefore **the original 4/4 positive target-MSE gains alone should NOT be cited as 4 equivalent demonstrations of the same individual-stream-level sediment response**.

### A fixed within-training-country-only predictor retains positive heldout gains

Define `M_within` using ONLY original within-training-country demeaned `Cov_j(s,y)` and `Var_j(s)` for the slope, with the *same* training-country-equal grand mean intercept. No target-country EPT values used in fitting.

| Heldout country | Original pooled Ridge α=10 gain | New `M_within` α=10 gain | M_within α=0 gain (sensitivity) |
|---|---:|---:|---:|
| Australia | +1.814% | **+1.170%** | +25.312% |
| Brazil | +2.396% | **+2.160%** | +29.902% |
| New Zealand | +1.448% | **+0.814%** | +23.793% |
| UK | +1.785% | **+1.390%** | +31.077% |

Original source country-balanced within-country covariance and actual heldout country within-country covariance were **negative in all four countries**; the fixed within-only trained beta was negative for every training set. This is coherent *directional ecological information* for the relative EPT share, **but not independent physical-site statistical significance**, not a causal effect, and it overlaps the previously published McKenzie et al. EPT/fine sediment conclusion. The `M_within` comparison still can benefit from mean-offset correction at the target; its positive heldout gain is not identical to the mathematically isolated within-target gain term above.

Visit-equal train weighting but preserving the same total weight 3 (to isolate sample weighting rather than regularization): α10 heldout MSE gains Australia +1.452%, Brazil +2.056%, NZ +1.424%, UK +1.206%; signs remain positive, magnitudes change. All values are **descriptive**. The taxonomy upper-bound in original UK `Insecta` uncertain source label changes the original MSE gain by <0.001 percentage point.

## Original source and scope guards

- **Original independent archival source**: [Loughborough Figshare doi:10.17028/rd.lboro.24615570.v1](https://doi.org/10.17028/rd.lboro.24615570.v1), source original `MathersEtAl-2023-PairedVisualFineSediment-Data.xlsx`, **CC BY-NC 4.0**; original MD5 `1e4d6e37b73c663406e5567dd79d9452` and source SHA256 `a02573bef0e11ea19c9465a78ebecc130fb9c39d76d35a79863c2efd5debcd48`. Precisely 6,491 original paired taxon/sediment `ANALYSIS_ID` observations (AU3099, BR232, NZ342, UK2818).
- **Original taxon response is normalized relative abundance**, not raw integer counts or biomass, and original 4country workbook does not supply a physical site identifier or within-site repeat grouping. Therefore **NO site-cluster independent bootstrap95, sample-wise p, or generalization to 6,491 independent waterways is supported**. Country-level independent domains only n=4; country-heldout MSE is a conditional descriptive property of these source samples.
- **Original source taxonomy harmonized before original EPT–sediment fit**; 302 names, only UK `Insecta` uncertain <0.02% normalized mass; original source-only fixed unknown mass gate passed, no relaxation. [Original PR #146](https://github.com/zuizui0223/gesui/pull/146) documents all strict source QA.
- **The original paper already reports EPT indices as strong indicators and regional variation**: McKenzie et al. (2024), *Global Change Biology* 30:e17084, [doi:10.1111/gcb.17084](https://doi.org/10.1111/gcb.17084). Do not claim first demonstration of deposited fine sediment/EPT association, global universality or newly discovered environmental response.
- External original visual channel-bed *deposited fines <2mm %* **is not** USA EPA field `XEMBED` (which mechanically assigns sand/fines as 100% embedded). Original EPA historical HUC8 source-status interaction cannot be independently replicated without same exposure and source-history definition. **No delivered wastewater TP/N or water-to-land nutrient transfer ever measured here.** Keep sewage-specific source GO failure unchanged [PR #139](https://github.com/zuizui0223/gesui/pull/139).

## Source-signed and computational independent audit

1. [Original signed 4-country α=10 fit from publisher original source #45032490453](https://github.com/zuizui0223/gesui/actions/runs/45032490453), artifact ID **14758661327**.
2. [All 16 original α=10 numeric results REPLAYED and within/between target exact covariance identity + fixed α grid + country-vs-visit weighting #45083417535 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45083417535), original output artifact **14778218065**. Source SQL/publisher raw files not committed.
3. **[Independent raw publisher XLSX re-read, 2×2 weighted normal equations (not original fitting code), within-country source centering, all alpha-grid heldout gains and exact target MSE identity #45083615144 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45083615144)**, independent numeric artifact ID **14778103720**. All combinations and source decompositions reproduced to numeric tolerance. The independent solver uses the SAME original observational data, so this checks math, not field-level replication.
4. Code: `scripts/gesui_fourcountry_EPT_within_between_country_gain_audit_20261123.py`, `scripts/audit_gesui_fourcountry_EPT_withinbetween_manual_WLS_20261123.py`, with immutable original source protocol [PROTOCOL_20261123_FOURCOUNTRY_WITHIN_BETWEEN_TRANSPORT.md](https://github.com/zuizui0223/gesui/blob/analysis/external-fourcountry-within-between-transfer-audit-20261123/docs/PROTOCOL_20261123_FOURCOUNTRY_WITHIN_BETWEEN_TRANSPORT.md). **Posthoc** scope and negative/limited results preserved.

## Publication decision

**Source/numeric GO** for a strictly *descriptive supplementary cross-source ecological sensitivity* that replays a known prior sediment/EPT relationship and audits its apparent transfer mechanics. **Scientific novelty and causal GO: NO** for standalone global sediment/EPT paper: original McKenzie paper has already addressed the geographic ecological phenomenon, site-grain identity is incomplete, and we cannot compare absolute XEMBED to visual fine-cover physical scales or test historical EPA source-status interactions. The contribution may be useful as an **audited cautionary case of sampling-domain composition and regularization in ecological transportability analysis**, but additional truly independent source/measurement coverage is required to elevate it to a new universal ecological theory.

**Do not add these four-country results to the primary *sewage* paper as independent WWTP pollution exposure data.** 
