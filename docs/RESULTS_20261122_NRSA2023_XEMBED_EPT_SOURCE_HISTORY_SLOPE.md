# 2023年NRSA EPT群集：河床埋没度への応答傾斜は元EPA解析標本の出現履歴で異なるか

**2026-11-22 / POSTHOC / observational; no causal sewage source inference.** The single exploratory interaction estimand was fixed **before the new numerical test** in [PROTOCOL_20261122_NRSA2023_XEMBED_EPT_SOURCE_STATUS_INTERACTION.md](https://github.com/zuizui0223/gesui/blob/analysis/nrsa2023-ept-xembed-source-moderation-20261122/docs/PROTOCOL_20261122_NRSA2023_XEMBED_EPT_SOURCE_STATUS_INTERACTION.md), but the question was selected after earlier positive source-novelty transfer gap findings, so cannot be called preregistered or independent confirmatory testing.

## Ecological question and prior support
The original prior EPA NRSA 2008–09/2013–14/2018–19 source-trained EPT predictor was more useful when the whole 2023–24 target HUC8 had been *represented somewhere* in those older eligible source panels, although **every** 2023 target HUC8 was excluded from training. The previous 2023 source-locked comparison found the prediction gap survived matching on original source geography, stream order, water TN/TP, XEMBED, channel width, true particle granulometry, and NHDPlusHR modeled flow/urban category ([PR #143](https://github.com/zuizui0223/gesui/pull/143)). Another fixed analysis found no clear primary E/P/T **order composition** differences ([PR #144](https://github.com/zuizui0223/gesui/pull/144)).

**New direct test**: Does the observed within-2023 relationship of original TOTAL300-counted EPT fraction to original measured bed XEMBED have a *different slope* between the two source-history groups, after environmental source-only matching?

## Original data and independent sample units
- Original actual EPA 2023–24 source SHA on five published files from signed [PR #141](https://github.com/zuizui0223/gesui/pull/141), original BERW/WADEABLE 1,165 source-selected visits/1,046 SITE_ID/693 original HUC8.
- Earlier original 2008/13/18 three-era census rebuilt: 2,648 visits and 1,030 source HUC8, old training after removing **ALL 693 target HUC8** is 1,185 visits/551 HUC8. Target sample old-source-ever 479 HUC8 (879 visits) versus source-never 214 HUC8 (286 visits). **Historical-source-never is NOT previously never observed rivers; it only refers to analytical source inclusion.**
- Use strictly the *original response-blind source-matched HUC8 pairs* selected **before** new slope results: 214 `grain` pairs source SHA `8c78c9060e0d512666062760f0467baaa54f6cfabeeaf76867ad609e3f1ea8e2` across 9 eco9, and 209 `hydrology/urban` pairs source SHA `4d21361aa86707039d4afd62cfa11c70f629c50ea01e1a0663e4c00930f254b2`. These two populations heavily overlap and are not independent experiments.
- One HUC8 = one observed EPT fraction (average of its original visit ratios) and one observed environment mean; grain-match pairs do not reuse HUC8. No EPT outcomes used to select pairs. Original EPA TOTAL300 count and provenance from original source code, not group interpolations.

## Single primary matched-pair source-history × XEMBED contrast

Let `D_i = meanEPT(novel source HUC_i) - meanEPT(previous source HUC_i)`. Define `M_i = [mean XEMBED(novel)+mean XEMBED(previous)]/2`. Under a simple linear model for each source-history group, the slope of D on M measures the **conditional difference in the two groups' XEMBED–EPT slopes**. Control pair differences in original XEMBED, measured stream logTN/logTP, field width, latitude, longitude and stream order. All variables scaled by original 2023 all-693 HUC8 source-only standard deviation. This is a model-dependent association, NOT causal sediment effect or a randomized source status.

| EPA original 2023 HUC8 paired test | Source grain 214 matched HUC8 pairs | Hydrology/urban 209 matched HUC8 pairs (sensitivity) |
|---|---:|---:|
| Standard deviation of original XEMBED | 25.612 units | Same fixed original 693 HUC source SD |
| **NOVEL−EVER slope difference, EPT percentage points per 1SD XEMBED** | **+3.759 pp** | **+4.270 pp** |
| HC3 pair-robust description 95% | **[+0.225,+7.292] pp** | **[+0.709,+7.831] pp** |
| HC3 exploratory p-value, unadjusted | **0.03707** | **0.01875** |
| 2000 pair-resample model-refit 95% | **[+0.167,+7.215] pp** | **[+0.894,+7.864] pp** |
| Descriptive adj-R² of paired-difference model | **−0.0013** | **+0.0341** |

- Unadjusted original EPT fractional slopes in the SAME 214 matched HUC pairs: old-source-ever **−13.423 EPT percentage points per 1SD XEMBED** [−15.820,−11.025], source-novel **−8.454 pp** [−11.450,−5.459]. These are *unadjusted*, not the slopes decomposed by the paired regression. Both indicate lower EPT fraction in more embedded/granulometrically finer original source bed, while the source-novel association is shallower.
- Primary 214 pair data remain adequate after source gate: 9 eco9, matched-pair mean embeddedness 10–90% ≈32.93–98.03 in EPA units; range ≈2.54 original source SD.
- Prespecified diagnostic controlling the *leftover paired differences in PCT_FN/PCT_SAFN/LSUB_DMM* changes primary γ to **+4.334 pp per XEMBED SD** [HC3 +0.798,+7.870], p=0.01630, paired bootstrap [+0.725,+7.827]. Despite this, XEMBED source highly correlates with grain fraction, so controlled estimates may reflect multicollinearity/selection rather than an independent depositional process.

## Critical falsifier: heldout-pair prediction does NOT establish universal moderation

On the SAME 214 source-frozen pairs, entire HUC pairs are assigned to deterministic fivefold CV (heldout counts 41, 51, 43, 38, 41), train-only standardization, Ridge alpha 10, baseline paired difference controls versus baseline plus source-slope feature `mean(XEMBED)`.

- Pair-heldout baseline D MSE: **0.06221863**
- Adding the matched XEMBED-moderation term pair-heldout D MSE: **0.06055467**
- **Heldout pair-equal MSE improvement +2.674%**, but fixed OOF matched-pair 2000 bootstrap descriptive95 **[−0.976%,+5.866%]** includes zero.

**The local slope difference is statistically detectable in this posthoc matched model, but stable prediction of the group difference in previously unused HUC8 pairs is not yet supported.** The low adjusted R² means the model explains little of pair-to-pair differences in EPT relative fraction. This is essential to the ecological interpretation, not a nuisance to conceal.

## What this result explains, and what it cannot
- It is **consistent with** the prior transfer asymmetry: older-source-ever HUC8 have more strongly negative EPT–XEMBED slope than source-novel HUC8, so a historical single slope learned within source-selected older years may travel poorly to previously unrepresented source areas.
- It **does not identify why**: source status is historical data eligibility, not a treatment, and may proxy differences in river-type histories, habitat composition, genetic/taxonomic pool, catchment land use, substrate sampling frame and unmeasured interactions. A simple `grain PCT_FN/SAND` match does not control all sediment dynamic processes; original XEMBED is also a particle-composition proxy.
- Neither a significant HC3 coefficient nor a bootstrap interval on 214 posthoc pairs establishes independent time-out transfer of the **source-moderation model**. In particular source-novel 214 HUC8 have already informed the present descriptive test; the 5fold pair split is internal cross-validation, not an independent year or continent.
- **No direct measured wastewater delivered nutrient load, sewage mechanism, aquatic insect emergence or land consumer flux.** Do not mix this standalone source-domain benthic-model analysis with the EPA POTW effluent SOURCE-STOP [PR #139](https://github.com/zuizui0223/gesui/pull/139).

## Audit and repository state
- Original EPA genuine 2023 source and three-era original rebuilding + exactly paired matched original XEMBED response test: [Actions #45029392553 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45029392553), signed aggregated artifact ID **14756503777**, original code `scripts/gesui_NRSA2023_matched_XEMBED_source_status_slope_20261122.py`.
- **[Independent EPA source+numerical QA Actions #45029595902 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45029595902)**、独立計算のsigned aggregate artifact **14756852797**。別のスクリプト`scripts/audit_NRSA2023_XEMBED_EPT_source_status_independent_numpy_20261122.py`で、旧3期原本ソース、2023真正SOURCE UID、source-only HUC8-pair SHAを再監査し、**statsmodelsを使わない手計算のNumPy正規方程式・HC3サンドイッチ共分散・214/209ペアの2,000回bootstrap、scikit-learnを使わないRidge10の全5fold train-only標準化／pair-heldout予測誤差**を独立再計算。主なgamma **+3.758670121346144 pp**、HC3[+0.225429710121706,+7.291910532570584]、粒径追加条件・209水理感度・CV +2.674% [−0.9763,+5.8664]まで元の推定結果と数値一致。**計算の再現であって独立EPA現地データの再現ではない**。
- **Research decision:** retain as *exploratory moderation signal* with a negative/nonrobust predictive guardrail, not as universal biogeographic source-gating law or sediment causal mechanism. Stop hypothesis expansion until external independent survey method/era or mechanistic sediment data are available.
