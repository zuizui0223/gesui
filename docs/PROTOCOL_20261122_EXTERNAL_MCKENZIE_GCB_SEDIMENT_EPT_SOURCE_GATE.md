# Independent external sediment–EPT source opportunity: 4-country Loughborough dataset
Date 2026-11-22. **Source-first gate — no biological values accessed before schema & grain audit**.

## Why this source and strict non-identity to EPA
McKenzie et al. (2024) *Global Change Biology* 30:e17084, doi:10.1111/gcb.17084, openly reports **6,491 paired aquatic invertebrate and percent deposited fine sediment samples** in **Australia (3099), Brazil (232), UK (2818) and New Zealand (342)**, data doi: **10.17028/rd.lboro.24615570** (Loughborough Figshare). It established regional differences of fine sediment–invertebrate assemblage response. Those findings are already published, NOT a new discovery by our project.

A rigorous independent-sourcing question is whether the **measured EPT individual fraction** (not EPT richness) can be harmonized with the true original EPA `TOTAL300` EPT fraction and whether any old-data sediment slope / source domain prediction generalizes. EPA `XEMBED` is method-defined percent of grains embedded, automatically 100% for sand+fines; external variable is visually estimated **percent channel bed area covered by sediment <2mm**. **Never compare regression coefficients as if same physical variable.** The 4-country source has no EPA 2008/13/18 HUC8 historical source membership or same source UID, so it CANNOT replicate the history-status interaction or classify novel-EPA HUC8.

## Source GO before model fit
1. Obtain API original metadata from `https://api.figshare.com/v2/articles/24615570` (or institution-specific equivalent) and publisher-provided file listing, name, file size, SHA/MD5 if available, license, DOI, version. Never invent file names or data.
2. For any candidate original quantitative files, download the publisher's exact file via supplied file URL. Record SHA256, sizes, column names, number of records, and identify the **grain** (per site vs per sample vs per taxon) without opening biology numeric values if possible.
3. Establish real sample ID that links sediment percent to abundances, country labels, sample year/site identifiers and taxonomic E/P/T classification. Record whether individual **abundance** (not only taxon occurrence/richness or precomputed trait scores) is available and how E/P/T taxonomy is determined.
4. Source GO only for countries with **at least 100 independent original sites, at least 150 paired actual samples, >=25 distinct fine sediment % values and >=20 percentage points of observed sediment coverage**. Do not treat repeated visits as independent sites. For holdout-country prediction, require at least **2 countries GO**, identifiable response/measurement equivalence and feasible country-level categorical mapping. Do not pool raw counts or infer absent missing data.
5. If only precomputed ecological diversity indices are supplied without EPT individual ratio or source site ID, **STOP** original EPT fraction replication. We can still report a descriptive independent paper comparison, but NOT numeric transfer validation.
6. If 4-country DOI is inaccessible, record failed endpoints and HTTP status, do not create fake results or use paper's published summary numbers as if newly computed.

## Analysis ONLY if source GO
- Fixed endpoint original EPT **individual fraction** (E/P/T true abundance divided by all positive counted benthic animals; unclassified taxa in denominator).
- Predictor original visual % fine sediment; age/country/site controls original fields if provided; country-blocked cross-validation and whole-site disjoint split.
- Primary predictive question **source-domain robustness**, not a new description that sediment affects EPT. Comparison: baseline country-appropriate constant/available environmental predictors vs adding sed % (bounded), within-country leave-groups-out and country-heldout. Also report country slopes with site-cluster robust intervals. A source lacking original site IDs can only support per-sample descriptive correlation, not independent holdout.
- **Never** describe any result from this 4-country source as independent testing of HUC8 old-source-status×XEMBED; only a separately sourced conceptual comparison.
- Review licensed reuse, authors, DOI and original data provenance. No external raw site records committed; only aggregate manifest and code on branch.
