# Original four-country sediment–EPT sampling design audit, 2026-11-23

**SOURCE ONLY** before any new biological outcome analysis. This is posthoc to country-heldout results in [PR #147](https://github.com/zuizui0223/gesui/pull/147). The original GCB paper already reports regional sediment–EPT response differences (McKenzie et al., *GCB* 30:e17084), and the 6,491 Figshare XLSX rows are **normalized taxon relative abundance**, not independent physical sites.

## Strict new source provenance
Publisher original [Loughborough dataset v1](https://doi.org/10.17028/rd.lboro.24615570.v1), CC BY-NC 4.0:
- Original XLSX `MathersEtAl-2023-PairedVisualFineSediment-Data.xlsx`, file 43299387, source SHA256 `a02573bef0e11ea19c9465a78ebecc130fb9c39d76d35a79863c2efd5debcd48`, MD5 `1e4d6e37b73c663406e5567dd79d9452`.
- ORIGINAL accompanying source descriptions `MathersEtAl-2023-PairedVisualFineSediment-DataSourceSummaries.docx`, file **43299384**, publisher **MD5 `08d81e4d97af34053b453a1c1e4f8c7d`**, **279,236 bytes**. This is an original documentary source, not generated project data.

First verify the exact original DOCX MD5/bytes, extract **document metadata, original structure, limited source-method hints, source references and linked sites** without copying the original manuscript verbatim. Determine from it whether physical sampling site IDs, repeated visits, country-specific sampling design, sampling date, and taxonomy normalization are defined. Explicitly cross-check primary paper source methods. Original Excel `ANALYSIS_ID` merely joins taxonomy and fine-cover columns; numeric sequential sample IDs in AU/BR/UK are NOT independent sites by assertion. NZ dotted IDs (`Death.GRES1.42417` etc.) may contain source/dates; DO NOT classify their prefixes as true independent physical sites until the source document defines the naming system and repeated visits.

## Conditional sensitivity only if site grouping authenticated
If original data source documents explicitly define site grouping that can be joined to original exact 6,491 `ANALYSIS_ID`, perform a true grouped leave-site-out (or nested source-unit heldout) **only within countries with provable identity**. Source GO is at least 30 source sites and at least 100 sample rows within that country, with grouped IDs reproducibly and unambiguously extracted from the original documented code. Other countries remain SOURCE STOP for site-independent statements. Never impute site IDs, dates or negative/zero counts from arbitrary delimiters.

If site grouping cannot be authenticated from documents, output **SITE SOURCE STOP**. Preserve earlier published four-country country-heldout MSEs as **descriptive conditional predictions only**. Do not compute independent-site bootstrap CI/p, claim 6,491 independent streams, or infer EPA `XEMBED` historical source-domain × EPT slope replication. These datasets have different sediment concepts and no sewage treatment plant delivered nutrient dose.

## Audit of source paper novelty and attribution
McKenzie et al. original article explicitly analyzed EPT indicators and regional sediment response; any result repeating their within-country association is background, not independent novel ecology. The most substantive novel contribution at present may be algorithmic: separating shrinkage/within-country source covariance/country mean-offset in model transport. Even that is exploratory and source-specific. Reject broad claims of universal ecological prediction until truly site-independent validation exists.

No raw original DOCX/XLSX committed; only source hashes, source audit code, minimal documented provenance, conditional source GO/STOP verdict and aggregate signatures. All PRs remain drafts, `main` unchanged.
