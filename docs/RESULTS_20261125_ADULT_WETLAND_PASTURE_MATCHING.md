# Wetland and pasture proxies for adult butterfly resource traits (2026-11-25)

Protocol fixed 3ac942288c70ee613af0ce50c3735aee6e7c4c2f. Original source acquisition [45568408341 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/45568408341), signed artifact 15001585693 SHA256 1243704d7dea61297e4a2a8043bf139540d9f32952a567ad7f35027c99e7e042. Native BasinATLAS Level07 DBF SHA256 1c73db5c07cb02ee6d10b963a0fc9fa4321fab55ff43721e8d6cd5c1054fb113. Full 57,646-row native field extraction; exact independent local replay of 1219×8=9752 legacy values. All 232 original H7 have new fields. Wetland sg2 excludes lakes/reservoirs/rivers; only 15/232 nonzero integer-percent values; pasture nonzero 223/232. Pasture is circa-2000 EarthStat, wetland legacy GLWD; not contemporaneous animal dung, carrion, mineral concentration or wet microhabitat. Original publisher redirect host changed; initial source job stopped, corrected exact delegated Cloudflare host and retained native hash. No GLW4 livestock or CORINE raw data acquired.

Same prior 357 rows/202 H7/78 H6 broad168, 357/203/78 common44, prior source support and all prior adult/host indices reproduced (max error 2.78e-16). Nonzero wetland 13/202 and 14/203 H7; pasture 195/202 and 194/203. Fixed log1p local percent, prior background/host/effort/method/country/period/W plus open-water difference, forest, elevation, annual temperature and precipitation, population density. Ridge10 H6-equal, train-only scaling, fixed-loss 2000 H6 bootstrap; paired adult An/Mi labels shuffled within source butterfly family/genus with support gate.

## Primary matching (increment over background + hosts + other resource)

| Endpoint | Ferrer-Paris or broad | Clarke |
|---|---|---|
| Broad168 Mi within-family wetland | -2.5753% [-5.9018,-.4124] | not applicable |
| Broad168 An within-family pasture | +.1114 [-2.3483,2.8324] | not applicable |
| Common44 Mi total wetland | -1.3843 [-4.0048,.4290] | -.7808 [-2.8060,.7199] |
| Common44 An total pasture | +.0104 [-3.3116,3.2648] | -.3040 [-3.7968,2.7142] |

All primary unsupported. Broad within-family Mi wetland permutation p=.9334, An pasture p=.4684; common44 total Mi wetland p=.943/.6956. An genus-shuffle STOP: three variable genera, seven movable taxa (fixed minimum ten); no p-value assigned. No wetland biological absence inferred.

## Exploratory cross-country secondary signal, not a primary pass

Common44 An within-genus pasture increment with 3 country holdouts (318 rows, all shared H6 excluded from training) +17.8086% [3.9633,27.4837] old, +12.8202 [2.1629,22.7009] Clarke. Individual country gains: FR +3.9633/+2.1629%, IT +5.0046/+3.0508%, NO +27.4837/+22.7009%. Three-country fixed-loss intervals do not establish generality; Norway dominates gain. Same component H6 holdout +1.3755 [-3.2725,5.6838], +.1983 [-5.6453,4.3775]; signed standardized-exposure slopes +.021959 [-.002252,.046171], +.013571 [-.011158,.038300]. No supported primary or taxon-specific permutation. Other nominal positives include wetland vs An between genus (nonmatching), low support, no multiplicity adjustment. Ten of 198 comparisons nominal positive. Do not rescue via secondary results. No causal relationship or actual feeding, dung, carrion, sodium measurement. Preserve previous host identity and ant range findings and all nulls.

## Reproduction

Local complete archive gesui_adult_landscape_20261125.zip SHA256 421aff1dfcbd90d2437136c8b7e2a98eedf81264ba2322642ae5c5a0a35b8384; audit.py 099ad6d6d194da7422e57b388904b3f2ad3678fa11b49f05db5161c82a8d2dd3; summary.csv 08adfda84da3bc072147efbc5b9479546db773cc5b3b49c9ad66e4fb547c0f61; receipt.json cefc86d5272b0df9c8a498f6309fc5fc000d3030e84396b78c23e7d1ffb419ca. 198 comparisons, 108 slopes, 12744 loss rows, 39992 null draw rows, 25704 saved primary prediction rows. Independent ridge solver discrepancy 5.58e-15, loss replay 1.25e-13. Source acquisition code and protocol and report in GitHub; numerical code and complete tables in local ZIP; no remote numerical replay or main merge. Next distinct independent source: direct FAO GLW4 ruminant density and finer wet habitat, not further threshold or endpoint screening on this panel.