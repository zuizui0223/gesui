# Final population-density sensitivity and manuscript package

Fixed before ecological fitting. Existing outcomes are already known; this is an exploratory sensitivity on the same observations, not independent confirmation or preregistration on unseen data. One additional population covariate; no alternative epochs, resolutions, footprints, outcomes or subsets selected after results.

## Exposure audit and source

Original W is `log1p(upstream modelled HydroWASTE discharge / (upstream km2 + 1))`, as reconstructed in `analysis/20261115_h6_source_expansion/verify_original_1219_covariate_replay.py`. HydroWASTE uses 2015 WorldPop to fill some served-population/discharge information (https://essd.copernicus.org/articles/14/559/2022/). Thus do not invert W to create population or claim statistical independence of all underlying inputs.

Use separately produced JRC GHS-POP R2023A, epoch 2020, Mollweide ESRI:54009, 1 km population counts. GHSL derives residential population from GPWv4.11 and built-up mapping (https://human-settlement.emergency.copernicus.eu/ghs_pop2023.php). It is a distinct product from the WorldPop layer used in HydroWASTE, but common census and settlement information can remain. Static 2020 represents the vicinity of both 2017–19 and 2020–22 windows; no annually observed population change is claimed.

Use the already cached official HydroBASINS Europe level 7 v1c polygons and NEXT_DOWN connectivity. For each of the 232 original H7 catchments, collect all upstream H7 polygons including itself. Compute population in these polygons by exact fractional-cell overlap in the equal-area raster CRS. Divide total population by summed polygon area (km2); use `log1p(upstream_population_density_2020)` as the sole new term. This matches upstream rather than local-only human concentration. It does not fully control local human activity.

## Source checks before fitting

Verify archive/raster/polygon and specification hashes, source CRS, pixel size, population units and nodata. Check unique H7 IDs and acyclic upstream graph; do not double-count an ancestor within a target. Validate summed native SUB_AREA against target UP_AREA within max(2 km2, 1% UP_AREA); report any failures. Polygon validity repairs, if needed, must be recorded and cannot change projected area by more than 0.1%. Require at least 95% valid raster-area coverage of the complete target upstream area, finite nonnegative population, and positive area; do not replace nodata with zero. Do not extrapolate population over missing area.

If coverage removes targets, retain only eligible original H7 and compare both models on exactly that same subset. Require >=186 original H7, >=72 original H6, >=8 original country-holdout groups and a full-rank model; otherwise report source insufficiency without relaxing gates. Eligibility is determined without HNS values. No special treatment of favourable countries or taxa.

## Fixed models and decision

Reproduce the EU451 controls from specification commit20b2cb6: crop, urban, riverflow, naturalout, area, nEEA, nAnimal, nPlant, nButterfly, primary_record_coverage_raw, primary_host_completeness, NUTRI_coverage_all_Plantae_speciesKeys, log1p_eiv_scored_taxa; country_method and period categories. Response remains primary_incidence_mean across all five butterfly families. No Plant_N_incidence/EQR/family-share adjustment or alternative outcome.

On the common eligible sample, estimate W coefficient and heldout W increment (a) with original controls and (b) with the same controls plus log1p upstream population density. Reuse original H6-equal WLS, H6 cluster-robust 95% interval, H6 hash 5-fold Ridge alpha=10 and 2,000 fixed-error H6 bootstrap draws, seed20261115. Training-only preprocessing. Original-data W SD is fixed across the paired models. Report change in coefficient descriptively, not as a causal percentage explained.

Population-adjusted support requires both positive coefficient-interval lower bound and positive heldout MSE-gain-interval lower bound. One criterion alone yields that limited support. Neither criterion yields uncertain positive independent W information, not proof of no sewage effect or proven agricultural causation. Report a supported negative direction if present. Inspect W/population correlation and residual W variability for identifiability, without changing the model based on those diagnostics.

## Finish the paper package

After this sensitivity, stop adding ecological models. Prepare a Japanese abstract and three review-ready figures: (1) landscape-only versus population-adjusted W–butterfly association and heldout gain; (2) French water–land spatial correspondence and uncertain temporal changes; (3) aquatic chemical/source prediction evidence with samples and uncertainty distinguished. Aquatic evidence is contextual validation within water, not proof that butterfly exposure is measured correctly. Do not call nonsignificant river paths absence of river transport. Keep Virginia biosolids outside this paper's ecological analysis and sludge/adult nutrition as unmeasured mechanisms.
