# Empirical persistence diagnostic — 2026-11-25

Protocol fixed at d8fb30ddf73e8468fa60e79bc77f7f1f7262d6ea before calculation. This is a descriptive audit on previously opened observations, not another hypothesis screen or independent replication. Verified 20 ZIP SHA256 values (including prior archives, matched FR73 and third source gate). Reconstructed 2020–22 butterfly inclusion on the exact common PR158/160 roster using scipy hypergeom.sf; maximum difference from sealed 2020–22 level 2.75e-11. Original 263 same-site means independently reconstructed and matched FR73 TP/EQR table. Exact balanced variance decomposition and independent weighted pairwise identities checked for 558 variable/taxon/cohort rows. Both host catalogues use the previous fixed common 46 accepted taxa (263 original and 519 Clarke links). No model fitting or p-values.

## 1. Four-window persistent geographic structure
For each variable/species separately, H6-equal and within-H6 H7-equal weights, four equal windows. Variance about grand mean = persistent differences among H7 means (H6 means plus within-H6 H7 means) + common period mean variation + H7×period residual. Median across nonconstant species; component medians need not sum to one. Weighted Pearson correlation of observed first vs last spatial levels is separate from the variance decomposition.

| Cohort | Index | persistent H7 mean share | H6 mean share | H7×period residual share | first–last r |
|---|---|---:|---:|---:|---:|
| 103 H7/53 H6 quality | butterfly inclusion | .70328 | .62006 | .28890 | .55596 |
| same | Ferrer-Paris mean host inclusion | .62660 | .54173 | .32685 | .45644 |
| same | Clarke mean host inclusion | .66166 | .56285 | .29780 | .46227 |
| 176 H7/75 H6 all | butterfly inclusion | .67943 | .60129 | .31424 | .59428 |
| same | Ferrer-Paris mean host inclusion | .64667 | .54596 | .31269 | .46047 |
| same | Clarke mean host inclusion | .67560 | .58056 | .28512 | .51805 |

46 nonconstant species in every row. H6 geography accounts for substantial variation; first/last correlations remain positive at median but far from unity. This does not prove static niche constraints or that all species are persistent.

## 2. Exact same-catchment intersection, two windows
2017–19/2020–22 terrestrial vs 2018/2022 physical river measurements. In full FR73 (24 H6), persistent H7 mean share logTP .97046, EQR .90328; cross-year r .95187/.82811. In the 61 H7/23 H6 intersection of France73 and common inclusive 176, the same quantities:

| Variable | persistent H7 mean share | cross-year r | nonconstant taxa |
|---|---:|---:|---:|
| measured river logTP | .96902 | .95397 | 1 |
| river EQR | .90894 | .84552 | 1 |
| Ferrer-Paris mean host inclusion | .80690 | .71445 | 46 |
| Clarke mean host inclusion | .80961 | .73074 | 45 |
| Lycaenidae inclusion | .84073 | .75630 | 42 |

Quality intersection only 33 H7/18 H6: logTP .97132/.95789, EQR .85244/.77080, old host .83450/.73039, Clarke .80178/.71592, butterfly .85142/.84843 (nonconstant taxa 46/45/41 for old/Clarke/butterfly). Smaller support and descriptive only. These are matched H7 and windows, NOT colocated riparian plant/larva plots. Finite-window means absorb stochastic variation, especially with two windows; do not interpret the fraction as reliability or a true latent ICC, or compare two-window and four-window fractions as process effects.

## 3. Observation-process context
H6-equal median ratio of identified indexed records on quality 103 H7, 2020–22 to 2023–25: plants .460, pooled butterfly records .939; equal-H7 medians .263 and .461. All176 H6-equal median .448 and .757. Counts of GBIF records do not estimate abundance; collection/reporting/index changes could affect opportunity, rarefaction addresses count but not all observation differences. This does not demonstrate observation noise caused the temporal null.

## Manuscript interpretation
Persistent geographic structuring is a documented property of the observed variables, alongside substantially variable temporal residuals. This makes the observed spatial/temporal contrast plausible and provides a common descriptive bridge, not an identified mechanism. Candidate explanations include persistent source/habitat geography, species pools, observation process and differing response timescales. The variance shares and correlation are not tests of these mechanisms. Do not pool across variables, infer wastewater flux, or claim equivalence from the previous nulls. Place matched two-window diagnostic beside Figure 2 and four-window panel beside Figure 4 or supplement. Preserve the previous source-specific EQR-network nonreplication and all change/coherence failures.

## Reproduction
Local derived-only archive gesui_persistence_diagnostic_20261125.zip SHA256 e6b454db4beadce50550f47d59ce33d3bd5a4e15e8c588c08a432994d2741c6b; audit.py SHA256 7ebf993e07b9dcdba60f82362104bf7fae6d779925fc4ade02330b52c4215416; variance_summary.csv SHA256 2487fb1fab759ccd7316c08437b7aca5f3cb3029b15f59389392029253c40fba. Full code, variance rows, record ratios, receipt in conversation ZIP. Report committed only; no remote replay, no rehosted third-party archives, no main merge.
