# Aqua-terrestrial realized-consumer field bridge: pre-outcome source assessment

**2026-11-12 / Draft / source-first before re-reading numeric responses from two newly identified public field datasets.**

This branch extends earlier evidence, but does **not** merge unrelated datasets into a synthetic causal chain. Novel source products are read as separate **field panels** with explicit years and scales. This protocol is drafted *after* the prior 17 Oregon stream EPT-vs-emergence analysis and after checking the original papers' reported findings; any new models are therefore **posthoc / independent reanalyses of published field observations**, not prospectively independent confirmations.

## Scientific objective
Is the degree of coupling between (1) stream water/landscape stress, (2) larval benthic quality or aquatic emergence quantity/size, and (3) actual riparian consumer ecology determined by **prey quantity, prey quality, and terrestrial habitat jointly?**

Quantities and traits must be kept distinct. Realized functions require measured consumers (e.g., dietary aquatic prey fraction via isotopes), not inferred from GBIF records or generic presence. No eBird or other opportunistic bird count datasets.

## Independent field panel A — published 2025 Dryad, South Nation watershed, Ontario

Official publisher: Dryad [doi:10.5061/dryad.sbcc2frfk](https://doi.org/10.5061/dryad.sbcc2frfk). Original paper: Rideout et al., *Quality versus quantity: response of riparian bird communities to aquatic insect emergence in agro-ecosystems*, 2025, doi:10.3389/fsufs.2024.1484377. The authors collected **29 field sites** (agricultural ditches and forested/natural watercourses), not eBird. Reported data include:
- original site locations and watercourse type;
- benthic macroinvertebrate sample and metabarcoding indicators;
- measured emergence insect dry biomass by repeated trapping in **May–July 2022**;
- **field** riparian bird point-count community collected June **2022**;
- water chemistry and benthic site sampling in **September 2023** (NOTE **one-year mismatch**, so cannot assert simultaneous water-nutrient→insect→bird mediation);
- original landscape/vegetation measures.

**No ecological numerical data has been read through the connector before this source contract.** The paper's previous conclusions are known and must be cited, not relabeled our first discovery. Any reassessment must use the same 29 site IDs, not 29 independent watersheds, and acknowledge detection limits on single-visit bird counts.

## Independent field panel B — published 2026 Dryad, Québec, 2 rivers / 20 sites

Official publisher: Dryad [doi:10.5061/dryad.v6wwpzh8z](https://doi.org/10.5061/dryad.v6wwpzh8z), published 2026-03-10. Gagnon, Beauchamp, and Harvey, *Urban-driven homogenization of aquatic subsidy size structure cascades to riparian predator communities*, Oikos doi:10.1002/oik.12080. Original site design: **10 stations along each of just TWO streams**; spatial autocorrelation is strong. Field sampling May–June **2024**, source files contain:
- original individual emerging aquatic insects (taxon/family, body length, observed abundance);
- original spiders (family-level counts, e.g. Tetragnathidae), other riparian arthropod field data;
- field water quality including conductivity/dissolved O₂, and independently derived imperviousness;
- spider (Tetragnatha) and aquatic/terrestrial prey stable isotope δ¹³C and δ¹⁵N, and original MixSIAR diet fraction outputs.

The authors reported smaller emergent prey with impervious land cover but mostly stable overall biomass, with aquatic insect contribution to spider diet varying substantially. Those are **previously published results**, not new discoveries here.

## Stage 0: acquire original official dataset bytes and audit model feasibility, no reanalysis of outcome signs

Use original **Dryad public API v2** GET and immutable version + file IDs (documentation: `https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.XXXX` -> version API -> per-file metadata). If API blocked, use **official Dryad** public download URLs discovered on its own dataset page; no guessed mirrors, no generated data.

- Pin dataset DOI, dataset ID, published version, individual original file names, SHA256 of each file's original bytes, raw ZIP member names and each member's SHA256; license is Dryad CC0 unless file-specific explicit restriction says otherwise. Do not assume a license just from domain.
- Find original README(s), capture raw source schema (header names), row counts, source-native candidate SITE / DATE / TIME and sample IDs (NO model coefficients).
- Require for panel A **≥20 exact distinct original site keys** jointly available across an **emergence biomass site summary and the original field bird consumer record**, including at least 3 habitat/sinuosity groups or ≥15 sites with original numeric landscape and field-source metadata. Original water chemistry / benthic data coverage and temporal mismatch reported independently, NOT forced to pass all layers for a consumer-only test. Stop before consumer models if source-to-site joins fail.
- Require for panel B **≥15 distinct same-native SITE IDs** with emergence measurement and original spider/diet-related measurements; presence of *actual* spider diet mixing-model source field or reliable isotope records must be demonstrated independently, not inferred from original paper. **Only 2 original rivers**; source gate cannot grant catchment-level generalization.
- If site-matched support fails or outcomes are unavailable, stop **without** relaxing denominators, grouping site rows differently, switching to unrelated species, or inventing values.

Output source-only QC JSONs with fixed original site IDs and file hashes in temporary GitHub artifacts; do not commit raw external research data unless explicitly licensed and necessary.

## Stage 1: conditional, limited sample-size analyses (never merge A and B site rows)

**Panel A, if support PASS:** site-level mean/sum emergence biomass as per README *versus* actual riparian bird **insectivore proportion** and species richness from field surveys; compare a priori environment/sinuosity or forest/hedge controls with+without emergence quantity (and if the original benthic dataset supports **matched** 2023 quality, test EPT richness fraction as a **time-misaligned ecological descriptor**, not contemporaneous mediator). Use leave-one-site-out CV, original ecological groups, Spearman plus 10k site-level rank permutations; no p-values treating individual birds as independent field replicates. Original insectivore assignment must follow taxon/guild defined in original data, not a guessed trait. Avoid auto-stepwise model selection on n29.

**Panel B, if support PASS:** site-level emergent insect mean size, size range, abundance/biomass vs actual spider family counts or modelled isotope aquatic reliance from same site, with imperviousness and longitudinal position as covariates. Do not fit 20 independent-watershed random effect / claim replicated catchment-level inference (only two streams). Show descriptive within-each-river effect contrasts and leave-one-river-out transfer if supported; the site-level evidence is within-river spatial association. Stable isotope consumer samples, mixing-model output and prey source samples may be unpaired temporally; read native dates.

**Cross-panel synthesis:** compare *definition and direction* of quantity-vs-quality associations in the two genuine **independent field programs** and the 2022 HJA n17 emergence study, plus EPA NRSA2018–19 n1,220 measured river nutrients and benthic EPT. Distinguish source stressors: agriculture, urban imperviousness, and sewage are **different** and there is no measured WWTP dose in these panels. Do not combine slopes or create a fabricated globally pooled full-chain SEM.

## Publication claim boundary
This branch seeks a genuine empirically supported link **emerging aquatic insects ↔ riparian consumer abundance/diet** in observational field studies. It does **not** prove N/P from municipal WWTP caused predator/flower/butterfly pollination outcomes. Before claiming a global water-land *ecosystem service* generalization, need original independently replicated catchments, same-year source routing, measured emergence flux and realized terrestrial responses.

Keep Draft stacked on `analysis/aquatic-emergence-calibration-as015-20261112`. Do NOT merge `gesui/main`, alter older positive PRs or `chocho` manuscript.
