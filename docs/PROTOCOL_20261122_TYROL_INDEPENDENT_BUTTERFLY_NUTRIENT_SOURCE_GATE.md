# Independent butterfly × plant NUTRI data: South Tyrol 93-site source gate

Date: 2026-11-22. **Exploratory prior-informed test**; the Guariento et al. 2023 paper already reports plant EIV N–butterfly feeding breadth via fourth-corner (p=.01). It is not an independently preregistered new discovery. First perform source audit, and do not fit any HNS-by-nutrient effect until the independent source support checks below pass.

## Source version frozen before HNS support check

- Original study: Guariento et al., *Biodiversity and Conservation* (2023), DOI `10.1007/s10531-022-02498-3`, [CC BY 4.0 publisher landing page](https://link.springer.com/article/10.1007/s10531-022-02498-3).
- Official Springer supplementary workbook: `10531_2022_2498_MOESM3_ESM.xlsx`, byte SHA256 **`e0774f313c5e008e7d37cef7d20d9f822f1cc02d7d0ba9f0b3301d5a9cf9202c`**, acquired in [Actions #45017500460](https://github.com/zuizui0223/gesui/actions/runs/45017500460), artifact `14752133756` contains **only a workbook inventory**; no original XLSX is stored in this repository or uploaded as artifact.
- Source workbook's original sheets: `Community` (93-site x 100-butterfly-species counts); `Variables` (the same 93 site_code and plant EIV nutrient `Nutrients (Elenberg plant indicator)`, light, site mean temperature, elevation, site year, land use); `Traits`; `Species list`; `Metadata`. Metadata `site_code` must match. The counts come from standardized 30-minute surveys conducted 3–4 times per site (not one repeat across years).
- Original accepted-numeric HNS roster: `results/20261115_sync/GESUI_303_HOST_NUTRIENT_SCORES.csv`, SHA256 **`00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb`**. 303 original keys, 186 accepted numeric non-alias-conflict scores. HNS is a different European host-plant affinity score, **not a measured leaf-N concentration or a direct field fertility measurement**.

## Fixed source checks

1. Original workbook SHA and HNS table SHA must match. Only exact accepted binomial butterfly species names matched after NFC/NBSP whitespace normalization; never impute, split, or arbitrarily synonym-map cryptic `Colias hyale/alfacariensis`, `Pyrgus malvae/malvoides`, `Leptidea sinapis/juvernica`, other `Genus sp.`/slash taxa. Original publisher's ambiguous species allocations are a provenance caveat.
2. `Community` exactly 93 sites and 100 butterfly species columns; `Variables` 93 unique site_code, exact 1:1 join. Plant N must be finite and recorded for at least 40 sites, and have at least 10 distinct numeric values; no land-cover/land-use proxy for N.
3. A qualifying genus must have two **different numeric HNS** species, each actually recorded on at least **two distinct survey sites**. Keep all original raw species counts; do not discard genera or site by fitted effect direction. Require **>=10 qualifying butterfly genera** to proceed to an exploratory combined-genus numeric test. This count gate is conservative and borrowed from prior gesui support practice, not a universal statistical power criterion for this new study design.
4. Report how many qualifying genera have two species co-recorded at any site. Original 2019–2021 sites are **cross-sectional**, not a within-site temporal change. Do not pool these 93 sites with 183 watershed/40 H6 GBIF data as independent macro-replicates.
5. Report source effort and missingness limitations including HNS host-score coverage >=80% as a source sensitivity; if <10 genera remain, do not silently change the threshold or fit and report a primary coefficient on that restricted subset.

## Model only if source gate passes

Infer a cross-sectional **within-genus relative individualCount composition** effect; NOT true species population change. Use observed site × species counts, preserving true sampling time and native site_code. Fit a conditional-multinomial species composition model within each genus and site:

`Pr(species=s | genus=g, site=i) ∝ exp(alpha_species + HNS_species × [beta*z(EIV_N_i) + gamma*z(EIV_L_i) + delta*z(T_mean_i) + eta*z(elevation_i) + sum_l rho_l*landuse_il])`.

Condition on genus total observed counts per site and use site-level geographic/land-use-block standard errors or site cluster bootstrap, *not* individual butterfly records as independent replicates. Rank check for N × HNS after controls; if near collinear, mark identification STOP. `HNS_s` is globally fixed to the accepted score table, no posthoc habitat-specific reassignment. Define the primary effect as `beta` before opening it; only claim generality if held-out site/landuse prediction improves. Fit the exact same dataset with and without `beta`, show direction, interval, and both favorable/unfavorable holdouts. Add a species-within-butterfly-genus HNS-label randomization only if exchangeability assumptions defensible; otherwise descriptive.

Any calibration or historical positive selection makes this a **new exploratory replication** even if numbers agree. The prior paper's significant positive broader larval feeding breadth association is already known, so neither broad feeding niche nor plant EIV N association alone can be claimed as newly discovered here.

## Known residual risks

- Standardized observations are still adult-flight snapshots; 6 sites had only three visits. Original species matrix may include adult individuals with species assigned proportionately from uncertain field identifications; check fractional counts.
- Observed butterfly HNS scores are derived from European-known host plants, which may differ from the actual local larval host resources; site EIV N is Pignatti et al. 2005 regional plant indicator affinity, not a measured wastewater effluent, field N deposition, leaf nitrogen or soil chemistry.
- Land-use, elevation, light and N can be strongly confounded. No credible wastewater exposure or aquatic organism measurements in this 93-site dataset.
- If no eligible support, report **SOURCE STOP**, not a negative beta or absence of an association.