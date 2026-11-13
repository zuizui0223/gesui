# Macro watershed source-direction test: actual upstream to the monitored river reach versus basin-outlet wastewater geography

**13 November 2026 · New authentic hydrographic-source and three biota response analysis · Strictly observational POSTHOC · DRAFT / UNMERGED.**

## Main question

Does the earlier positive predictive increment from modeled registered WWTP flows calculated **to a subbasin outlet** actually track estimated discharge that is topologically **upstream of the benthic monitoring reach**, or is it largely information about the wider river network and settlement/registration geography? This is answered on the **same 232 European HydroBASINS Level7 subbasins and same 454 H7×2017–19 /2020–22 rows** as preceding PR #65/#66/#67, not from literature example pooling or a reselected sample. All original outcome signs and land-use source methods are preserved.

## Authentic new physical original source (before biological outcomes)

Original official [HydroRIVERS v1 Europe & Middle East SHP ZIP](https://www.hydrosheds.org/products/hydrorivers), source ZIP SHA256 `500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5`, original native reach DBF SHA256 `87b80b198b8fc512591f234e92b7fa5067fa87271246ff3a0334a7df1d0a972c`. **938,544 unique original river reaches** with full acyclic directed `HYRIV_ID → NEXT_DOWN` graph. Original HydroWASTE v1 exact publisher 58,502 registered treatment plants with modeled/mixed-source `WASTE_DIS` **m³/day**, 27,143 have official `HYRIV_ID` within this original Europe network; most plant operation statuses originally `Not Reported`, not verified active.

**Every one of 3,525 original EEA2018/2022 monitoring station-year source entries in the unchanged n454 H7-period macro frame** maps to original HydroRIVERS network (2,382 distinct station river reaches); original station source reach identifiers were retained in both original EEA Git blobs, including same method/year/time/source key. This is a real river network, not proximity-only Euclidean buffer or shared H7 polygon assumption.

**Prior outcome-blind source routing gate** [SUCCESS #43340832143](https://github.com/zuizui0223/gesui/actions/runs/43340832143) used only upstream topology and discharge-source original rows, without opening EQR, plant, or butterfly numerical outcomes. Derived source manifest SHA256 `d50eba0f934278cb28eb403af7b84d4da5a19af209a6d05be312a2ee44ad6d6d`. Exactly original **454/454 period rows** pass; **301/454** have ≥1 registered treatment facility at a distinct reach upstream of at least one monitored station; **300/454** within 50km on the original network counting FULL upstream reach lengths. Thus 168 distinct original H7 have a verified original *registered source upstream network position*, not a verified operating source/pollution dose. Registered sources on the **same** receiver reach are excluded due to unknown within-reach upstream/downstream ordering (471 of original 3,525 station-years have an unresolved same-reach registered source).

**Directional contrast source:** 438/454 H7×periods have a registered downstream facility within 50km, but only 434 have **positive modeled downstream WASTE_DIS**. This downstream observation is **NOT a balanced, randomized or source-matched placebo**: stations have far more facilities downstream because river and settlement geography differ. Still, it offers an important sensitivity check for geographic source enrichment.

## Main new outcome analysis — SAME original 454 H7×period rows, SAME physical/environmental controls and 90 original H6 parent holdouts

Baseline B2 = original EEA year/method/country/mean station river flow/site count, original GBIF Animalia/Plantae/five-family record opportunities, native H7 area, authentic HydroATLAS **UPSTREAM cropland fraction, urban percentage and natural modeled outlet discharge**. Group entire original H6 parent catchments (90 source-native groups, 5 folds), not just random H7 IDs. The original stricter 232 H7→H6 geometry source gate failed 216/232 without repair, then a **separate outcome-blind explicit `make_valid` sensitivity corrected only original invalid parent geometries** and restored 232/232; do not call those original data intrinsically valid or the modified source a new prospective trial.

Add separately:
- old static modeled **W_to_H7_BASIN_OUTLET** per upstream area;
- new genuine strict **W_to_EEA_MONITOR_RIVER_REACH** (source reach topologically distinct and upstream) / native receiver DIS_AV_CMS;
- new conservative fully-traversed **50km upstream** receiver index;
- **50km downstream** direction-control modeled wastewater facility index, not an equal physical placebo.

**Observed signed heldout MSE gains**: same original 454 H7×periods, n232 H7 and 90 H6 parents; separate outcome models, with original type-specific units:

| Predictor added to identical B2 original baseline | Real EEA benthic EQR | GBIF botanical speciesKey-rarefied 50-record S50 | GBIF Lycaenidae five-family indexed log odds |
|---|---:|---:|---:|
| **Old modeled WWTP source to H7 basin OUTLET** | **+11.129%** | −0.196% | −0.762% |
| **NEW modeled WWTP strictly UPSTREAM of monitored reach** | **+1.434%** | −1.297% | +0.449% |
| NEW strictly UPSTREAM of monitored reach within full-reach 50km | +1.275% | −2.414% | +0.252% |
| NEW distinct source reach **DOWNSTREAM** within 50km | **+5.782%** | +0.614% | +1.977% |

EQR original H6-parent-resample descriptive 95% MSE gain for old H7 basin-outlet model **[+3.944%,+18.267%]**; NEW strictly-upstream receptor model **[−0.580%,+3.577%]**, new upstream 50km **[−0.414%,+3.273%]**, and downstream direction-control **[−1.841%,+12.865%]**. Only outlet-model 95% estimate excludes 0, and that model is *not receiver-level delivered exposure*. This is not proof of zero true nutrient effect.

**Additional/independent source information**:
- Adding strict receiver-upstream source **on top of the old basin-outlet model** improves original EQR H6 heldout MSE by only **+0.464%**, original H6 bootstrap95 [−0.738%,+1.732%].
- Adding **downstream source** on top of original basin-outlet index improves original EQR H6 heldout MSE **+3.456%**, bootstrap95 [−2.341%,+9.560%].
- Adding original EEA actually-reported annual N source index on top of strict receiver-upstream HydroWASTE index **worsens** EQR H6 MSE **−0.491%** (original H6 bootstrap95 [−1.157%,+0.108%]), whereas adding reported N on top of OLD basin-outlet source shows +0.687% with 95 crossing zero. Original annual reported N source is positive only in 56/454 H7-year rows and is strongly selected by country and operational reporting.

**Country-heldout differences** (stable original H7 majority source country, n11 groups; country FE omitted because country was unseen): Original EQR static H7-outlet gain **+0.953%** row-weighted; new TRUE strict-upstream EQR **+1.544%**, strict upstream 50km **+1.116%**; **downstream 50km +11.432%** row-weighted. This cannot establish pollutant directionality: the downstream point inventory is much denser, heterogeneous land/river/country methods and facility statuses are uncontrolled. Original PR #67 country-equal-weight original outlet effect **−7.989%** remains a major negative transfer result.

The three response-layer results indicate that registered wastewater geography contains more stable extra predictive information about EEA biological **quality index** than the opportunistic GBIF botanical and butterfly **recording-based** proxies; this is a difference in prediction and observation processes, not an observed ecological service dissociation.

[Main true original receiver-specific upstream/downstream 3-layer macro CV, SAME 454 observations, H6 and country blocked — SUCCESS #43341517792](https://github.com/zuizui0223/gesui/actions/runs/43341517792). The full machine-readable result JSON is in the Actions artifact; no raw licensed HydroRIVERS data redistributed.

## Separate independent numerical and source-quality audit

A **different sklearn Ridge/SVD** numerical implementation verified the full original data source and **all 3 response families × 12 source models × 2 geographic evaluation designs =72 original MSE values** to numerical tolerance. This is numerical reproducibility, **not a second ecological dataset**.

Critical original source-quality audit:
- Original receiver study 3,525 native station-year observations, 301 original H7×periods with ≥1 distinct upstream registered facility, **299 period rows with a positive-modelled WASTE_DIS amount inside 50km**.
- Original downstream 50km registered presence is 438 H7-period rows, **positive modelled wastewater flow in 434**.
- Original **QUAL_WASTE=1** (=original actually reported WWTP flow quality code) source pressure within 50km was positive in only **68/454** original H7×periods.
- Original **STATUS exactly `Operational`** upstream within 50km is positive in **0/454** H7×periods, not because no sewage plant works, but because source operation-status coding is mostly `Not Reported`. Thus a sensitivity restricted to status-confirmed operating plants has no variance and cannot support a negative ecological conclusion.
- Observational registered-source upstream vs downstream log-pressure correlation in source H7-period rows ≈ **−0.0088**, upstream pressure vs old outlet source log-pressure ≈ **+0.299**: important evidence that the old registered outlet index is NOT interchangeable with topologically sourced receiver-level exposure.

[Independent sklearn 72-model numerical agreement and original reach-status source quality audit — SUCCESS #43342095959](https://github.com/zuizui0223/gesui/actions/runs/43342095959).

## Interpretation, implications for the user's global macro programme

**Do not state that wastewater N/P causes a consistent ~11% deterioration of benthic diversity.** The +11.1% is an incremental **predictive MSE gain** from an original modeled **basin-outlet geographic index**; it is not magnitude or sign of biological damage, and the corresponding genuine receiver-upstream model gives only a small, uncertain +1.4% prediction gain. Downstream facility geography even predicts 5.8% better than baseline in H6 blocked folds (despite being physically unable to be a straightforward nutrient input to an upstream reach). This pattern is consistent with reporting/urbanization/river-network structure contributing substantially to the original positive geographic predictive result. It does NOT rule out true causal sewage effects or directly quantify their direction.

**Same spatial/temporal unit achieved**: EEA benthic real EQR, GBIF rarefied indexed Plantae, Lycaenidae relative family source counts, original upstream land/urban/flow and HydroWASTE point inventory all align to 454 H7-period rows. But a mechanistic water→insects→plants→butterflies or pollination *ecosystem service cascade* remains **unobserved** (no same-unit actual emergent insects, plant hosts, predator diets, pollen transfer or seed set). The earlier n169 exact-site cross-year three-layer co-response was near zero and is retained.

**Necessary future physical-quality data**: (i) georeference same-reach outfalls and stations along the line with validated downstream orientation, site positional uncertainty and same-date plant status; (ii) source 2018/2022 WWTP effluent N/P mass for ALL source-receivers and routing with retention/other sources, not registered presence; (iii) measured original receiving-river water TN/TP and standardized habitat/survey source; (iv) country/external basin validation. Reject binary "no source reported means no effluent." Do not substitute an index response for a realized service.

**Preserve** current PR #67/66/65, prior unmerged source branches, source SHAs and all positive and null results. This report is POSTHOC, unmerged Draft only; no `gesui/main` or `chocho` changes and no eBird.
