# Real source-locked river Plantae dataset (16 outfalls, 48 locations)

**Original analysis:** https://github.com/zuizui0223/gesui/actions/runs/42365305227 (SUCCESS). Source geometry: https://github.com/zuizui0223/gesui/actions/runs/42364866919, fixed 16 original selected outfalls.

Actual GBIF 2011–2024 Plantae (not butterfly-host-only) speciesKey observations at 16 directed U2/U1/D1 river triplets; 48/48 sites had plant records; 1,402 unique indexed speciesKeys across all sites. Each source file has a content SHA256 in SOURCE_RECEIPT.json. No additional GBIF request is made by archive preservation.

**Prespecified strict response quality:** 12/16 triplets contain >=30 species-keyed records and >=6 indexed plant species at all three sites, with max/min site species-keyed record imbalance <=10. Exact record-rarefaction to 30, then richness curvature Y = [log1p(S_D1)−log1p(S_U1)] − [log1p(S_U1)−log1p(S_U2)].

Modeled wastewater water fraction vs Y: Spearman **rho +0.125874**, n=12, **exploratory p=0.69668**. Sensitivity with rarefaction to 10 records: rho +0.173626, n=14, p=0.55277. Neither shows a clear nutrient-mediated vegetation response. The primary high-minus-low mean Y is -0.059678 (7 high vs 5 low qualified), not positive.

Companion Jaccard raw species-set turnover is strongly affected by geographic segment lengths and GBIF detection. It must not be interpreted as a source-specific effluent effect or species dispersal without matching distance and observation intensity.

**Caveats:** the third-party HydroRIVERS v1 Europe conversion has not been independently verified against original official bytes; hydrogeomorphic source similarity does not measure N/P load. U1 and D1 are on different river reaches, but different distances, settlements and preexisting upstream discharges remain uncontrolled. GBIF records are not vegetation cover/abundance, and Plantae includes nonvascular taxa. This is the first independently acquired all-plant directional river observation table, **not an N causal test**. The older chocho port/host positive findings remain unchanged.
