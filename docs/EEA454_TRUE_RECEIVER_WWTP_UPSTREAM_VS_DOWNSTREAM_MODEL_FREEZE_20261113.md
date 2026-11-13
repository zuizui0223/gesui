# Same-receiver vs basin-outlet wastewater macro test: post-source exploratory contract (13 Nov 2026)

**Defined after confirming HydroRIVERS original source gate and previously known WWTP–EQR associations, before estimating this new comparison. This is exploratory/posthoc, not an unobserved prospective test.**

## Evidence and fixed population

Source original published Europe HydroRIVERS v1 SHP ZIP SHA256 `500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5`, DBF SHA256 `87b80b198b8fc512591f234e92b7fa5067fa87271246ff3a0334a7df1d0a972c`. Real **938,544** regional `HYRIV_ID` reaches and acyclic `NEXT_DOWN` network. Original HydroWASTE published 58,502 sources, of which **27,143** outfalls have original reach IDs on Europe network. No raw HydroRIVERS SHP or DBF is redistributed; derivative source features only.

Source original EEA **3,525/3,525 actual site-year reaches** matched 2018/22 original study, **2382 distinct monitoring `HYRIV_ID`** and the **same 454 basin-periods/232 H7** as PR #65/#66/#67. Source-only untouched receptor exposure [run #43340832143](https://github.com/zuizui0223/gesui/actions/runs/43340832143), derived H7-period JSON SHA256 `d50eba0f934278cb28eb403af7b84d4da5a19af209a6d05be312a2ee44ad6d6d`. Genuine **301/454** basin-periods with ≥1 registered **strictly upstream on other reach**, 300/454 within **50 km conservative FULL reach travel distance**. **168/232** H7 basins have strict-upstream registered sources. **438/454** H7 periods contain original registered facility on distinct downstream ≤50km reaches as geographical direction comparator. Same-reach ambiguous 471 original site-year records, **explicitly excluded from upstream**.

Status `Operational` is confirmed for only 3,767/58,502 original global facilities; rest mostly `Not Reported`. Source direct flow observation `QUAL_WASTE=1` is not facility operation proof. HydroWASTE WASTE_DIS is **m³/day** modelled/mixed quality, not true 2018/22 measured N/P nor pollutant-specific emission. The stream network and WASTE_DIS are STATIC original source products.

## EXACT prior observed-response population and covariates, no selection allowed

Same **454** authentic 2018/2022 H7×period rows and original **90 H6 true hydrological parent groups**, created from official HydroBASINS EU Level6 (16/232 native H7 required outcome-blind `make_valid` of 207 originally invalid parent polygons), matching source SHA `372b800b85abe2f395241987cea2000ff901cf5aff31bf1d1fb825304f2495fc`. Read previous original [90-H6 source result #43338763077](https://github.com/zuizui0223/gesui/actions/runs/43338763077), and **numerically reproduce prior B2,B3,B4 for ALL three outcomes before any new receptor analysis**.

B2 baseline = country and original EEA BQE method/period, source H7 area, EEA site count, original EEA site modeled discharge, original GBIF Plantae/Animalia/five butterfly-family indexed observation opportunities, HydroATLAS actual original **UPSTREAM cropland/urban percentages and natural modeled flow**. Use train-only standardization, Ridge α=10, 5 source-H6-parent-heldout folds. Both years of same H7 grouped together.

Three separate fixed responses in same rows:
- real original EEA monitored benthic EQR;
- original GBIF Plantae speciesKey expected rarefied richness S50 from 50 raw indexed specimens, not real flora;
- original Lycaenidae log odds among five accepted butterfly family indexed records, not abundance or pollination.

## New predictor contrasts and directional negative control

From genuine source-only per-EAA-river station receptor manifest, aggregate once per native H7×period with equal physical-site weights (already frozen, no cherry-picked removal):

- `W_old_basin_outlet`: original HydroWASTE **upstream-to-BASIN OUTLET**, log1p(original modelled discharge / H7 upstream area+1), retained unchanged for exact replication.
- `W_receiver_all_strict`: `log1p(original strict-upstream OTHER-reach registered WASTE_DIS m³/day / (EEA site's modeled mean DIS_AV_CMS × 86400))`, then averaged per native H7×period at prior source stage. This represents source opportunity for discharge that is geometrically upstream of the EEA receiving *reach*, not actually measured nutrient concentration at the interior station.
- `W_receiver_50km_conservative`: same but only distinct upstream reach sources whose **full source-reach and intervening-reach lengths** sum ≤50km; original method may exclude some true upstream sources farther within their reaches, so it is a conservative exposure proxy.
- `W_DOWNSTREAM_50km_control`: original distinct **DOWNSTREAM** ≤50km WASTE_DIS divided by original station modeled river flow, a geographic association comparator, **not an equally sampled causal placebo**.
- Additional source-quality sensitivities: strict 50km operational-only and `QUAL_WASTE=1` original-flow m³/day divided by upstream total H7 area (not same exact station fraction; metric marked separately). No undeclared new exposure category.

Original reported annual N/P station-source mass is separate. Use initial source-all vs downstream-specific results BEFORE considering positive or negative nutrient signs.

Compare nested models on same 454:
1. **B2** alone;
2. B2 + old outlet W;
3. B2 + strict receiver all upstream;
4. B2 + receiver upstream 50km;
5. B2 + downstream 50km index;
6. B2 + old outlet W + receiver upstream all;
7. B2 + receiver all + downstream 50km;
8. B2 + receiver all + original reported annual N for orthogonal source metadata (only 56 H7 periods positive);
9. B2 + receiver upstream 50km + original operational-only source quality as source-reporting sensitivity.

Report all signed H6-parent heldout 5fold MSE deltas for each of **three** outcomes; 2000 H6-parent bootstrap descriptive 95% intervals; original country-heldout with source country dummy removed as a *different model estimand*. No tuning years, sample, taxon families, original method, 50km window or penalty according to effect signs. No biological models if source 454 H7/3525 station reach origin QA had failed (it passed unchanged).

## Key safety / scientific boundaries

No biological significance is inferred from the origin of static source points; the original 27,143 registered regional WWTP include mostly status-not-reported and estimated discharge and station-reach positional uncertainty. Even source upstream-specific predictive gain cannot prove sewage **caused** EQR, butterfly, plants or food-web effects. If downstream control is stronger, retain as an observation process / urban geography confound. Country transfer failure from PR #67 must not be erased. This is exploratory physical-scale refinement, NOT an independently replicated causal field design.

Keep this as Draft on top of PR #67 and never merge `gesui/main`, `chocho`; no eBird.
