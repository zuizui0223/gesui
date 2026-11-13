# Spatial transfer robustness for same-basin wastewater–benthos–plants–butterflies macro analysis

**Frozen on 2026-11-13 before obtaining the new HydroBASINS Level 6 source crosswalk or new geographic-block model results.** This is **POSTHOC** relative to PR #65/#66, because prior H7-hash 5-fold improvements and their signs have been seen.

## Scientific question

Does the observed additional out-of-basin predictability of EEA river benthic EQR from (1) the **static modeled HydroWASTE outlet-pressure** and (2) the **sparsely reported 2018/22 upstream N-load index** survive validation where **geographically neighboring H7 subbasins must be held out together**, or when source countries are omitted completely? Equally test the original **Plantae 50-record rarefied indexed richness** and **Lycaenidae five-family record composition** in exactly the same source cohort.

### Same original sample and outcomes

Use the immutable original prior [232 native HydroBASINS L7 river basins ×454 basin-periods](https://github.com/zuizui0223/gesui/actions/runs/43215754868) from EU2018 and EU2022 EEA actual field EQR, accepted original GBIF S50 & five-family Lycaenidae. Preserved originals:
- full matched original source per-basin/period values with original country, BQE method, year, source GBIF record opportunities, local area, native river station modelled flow and prior HydroWASTE static outlet intensity;
- identical native [HydroATLAS level7 cropland, urban, natural outlet discharge](https://github.com/zuizui0223/gesui/actions/runs/43225986249) original SHA256 `ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686`;
- original 2018/2022 upstream reported N/P station annual mass converted to hypothetical concentration equivalents, [source-only original](https://github.com/zuizui0223/gesui/actions/runs/43222864099) SHA256 `012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26`;
- original land-use-adjusted B2/B3/B4 three-response numerical saved results [source-success #43227058486](https://github.com/zuizui0223/gesui/actions/runs/43227058486) from PR #66's successor branch. **The original H7-hash predictive gains remain unchanged and must first reproduce exactly**.

No choice of new outcomes, year, local/upper land-use definition, region, source reporting threshold, source butterfly family denominator or original 454-row sample based on results.

### Source Stage 0: authentic Level 6 geography, before new spatial-block results

Acquire official original HydroBASINS v1c **EU Level6** and **EU Level7** shapefile archives from `https://data.hydrosheds.org/file/hydrobasins/standard/hybas_eu_lev06_v1c.zip` and `...lev07...`; require original Level7 ZIP original source SHA256 `b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44`; verify Level6 official source original SHA256 after download. Read all 232 original H7 polygons from official Level7 ZIP, and all authentic Europe Level6 polygons. Match each source H7 polygon via geometrical intersection against Level6 using highest original *area intersection*, require fraction ≥95% assigned to one unique Level6 parent, no nearest-centroid or guessed parent ID substitution. Record source polygons with invalid original geometry and source coverage. Under 95% no forced geometry repair; if H7 valid/geo overlap fraction deficient, keep source unmatched.

H6 SOURCE gate: **≥220 of 232 original H7 basins** get a unique ≥95%-area-cover parent, and ≥30 independent Level6 parents among included H7, and at least five Level6 parents per held-out 5-fold group. If FAIL, STOP official Level6 blocked ecological fitting; report exact source and the original H7 results only. If PASS but <232 H7 assigned, preserve the original 454 baseline sample for strict comparison if possible; do NOT silently drop unmatched H7 from baseline. Source limitation must be explicitly documented before any exploratory sensitivity with a different common-case cohort; prefer require **232/232** for fair unchanged sample, else stop the *strict* spatial test.

Topological H6 grouping is the geographic spatial validation unit; neighboring Level6 polygons can still be adjacent and share a country; this is **stronger than H7 hash** but is not fully independent ecoregion transfer.

### Models and evaluation, fixed after source confirmation

On exactly the unchanged 454 source basin×period rows, build same original features and three Y families:
- **B2**: country FE + native assessment-method FE + period + GBIF original Animalia/Plantae/five-family reporting count opportunities + original EQR site count + H7 area + original site DIS_AV_CMS + HydroATLAS ORIGINAL **UPSTREAM** cropland fraction, upstream urban extent fraction, `log1p` natural L7 outlet flow. Preserve country/method categories and weights from prior model. Run train-only scaling and Ridge alpha 10.
- **B3**: B2 + `W_static`: modelled cumulative HydroWASTE discharge to **basin OUTLET** divided by upstream area.
- **B4**: B3 + original report-source annual upstream `log1p(N input mg/L-equivalent)`, where source positive only 56/454 basin-periods. Negative source index is absence of positive matched reporting not zero WWTP.

Primary robustness: **5fold source-H6-parent-heldout Ridge10** with all H7×periods sharing an H6 parent assigned to the same fold, deterministic SHA256(H6_ID) mod 5, fold-count / basin-count exact. Compare OOF MSE B2, B3, B4 separately for each of three responses. Report B2→B3, B3→B4 signed MSE improvement and **original H6-parent-bootstrap 2,000 descriptive intervals**, effective # unique H6 parents. Never interpret source predictions as causality or inferred services.

Country transfer: **leave-one-original-country-out (LOCO)**. Assign exactly one stable primary country to each original H7 basin across BOTH time windows via majority native source country, deterministic alphabetical ties; source H7 with inconsistent country labels must be grouped together. Omit country fixed effects for LOCO because heldout countries are unseen; use identical non-country baseline for all 3 Y, retain source EQR method labels with fixed categories (report unseen method columns) and original recording/physical variables. Report OOF MSE on country blocks, for each source country, aggregate **row-weighted** and **country-equal-weighted** gains, show performance in every country with ≥10 basin-periods. Countries with 1–5 rows are descriptive only. A source effect cannot be called pan-European if country-heldout prediction fails in major countries.

Important: Source N/P positive reports are heavily uneven (only 4 positive H7-period rows France vs 52 outside); no source-positive→source-negative causal comparisons, no simulated nonreports.

Keep independently reported original parent-PR successes, null results, and no `main` / `chocho` merge. No eBird. No new source effect coefficients should be advertised before H6/geographic block and country validation output.
