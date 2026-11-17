# 2026-11-17: First actual 2018–2022 Sentinel-2 10m river-station near–far NDVI (8 native H6, 8/8 QA passed)

**New source-authentic measured satellite pixels, no river chemistry / benthic / sewage biological effects fitted.** This is distinct from original HR-VPP TPROD, the earlier 2018-only RZ land-cover study, and the prior Austria S2 10-site source STOP.

## Original unchanged source and pre-outcome eligibility
- Original source universe = 263 identical French physical river WQ/I2M2 2018+2022 field sites, 73 H7 /24 H6.
- 24 H6 were each assigned a single source field station by SHA256 original EEA site_key, before optical NDVI/river TP/EEA EQR inspection. 2018/2022 June-July-August actual Element84 Earth Search Sentinel-2 L2A red, NIR, SCL source objects [Actions **#44008867555**](https://github.com/zuizui0223/gesui/actions/runs/44008867555) 24/24 qualified 2+ months; original source JSON SHA256 **22414dec230bc9bf9c308d823f166303cc843143c6ade922af6cb44b1f2c3f3c**.
- Next pixel smoke = exactly 8 of fixed 24 H6 selected by SHA256(original H6 ID). Scene for each year/month preselected by original tile cloud fraction (<=40%), at most 2 scenes/month/year; original red/NIR asset `raster:bands.scale` and `offset` **read from the actual item, do not hardcode 2018 offset=0**. SCL 4/5 real vegetation land, ESA WorldCover2021 class 10 tree,20 shrub,30 grass,90 herbaceous wetlands; fixed **station-centered** inner 0–250m and outer 250–500m rings. Minimum 100/200 clear land pixels, 40/100 natural clear pixels, >=10 each common natural class per ring/year; at least two **matched calendar months** 2018 &2022 for each site; smoke GO >=6/8 source H6. These are source-quality thresholds, **not** ecology coefficients or a model selected on NDVI values.

## Genuine real original optical pixel results
[Original Actions #44009758637](https://github.com/zuizui0223/gesui/actions/runs/44009758637) completed COG window extraction of all eight fixed H6; original data JSON was written **before** a final print-only Python error (`json.dumps(...,flush=True)`). The workflow UI concluded failure because the final log statement raised `TypeError`, **not** because real source pixel extraction failed. All eight source station-year pairs had >=2 actual matched months **8/8 PASS**; seven passed 3 paired Jun-Aug months, one passed 2. Original artifact #14295254324 exists and includes complete per-scene class-wise original genuine reflectance/clear-land pixels, JSON SHA256 **f6b2d24e26e6519ab6f2d84b0aa9ec842e5c9147fbcba2d10090ec7be7acdb0e**. 50 attempted source scenes, 3 insufficient clear/natural pixel attempts; no apparent original COG access failures in those scene attempts.

Observed exploratory 8 original per-site NDVI **delta of (near−outer) from 2018 to 2022**, source matched same original calendar months:
| source H6 | paired original months | 2022−2018 near-minus-outer NDVI |
|---|---:|---:|
| 2060022970 | 3 | −0.028817 |
| 2060552460 | 3 | −0.005973 |
| 2060499090 | 3 | +0.019190 |
| 2060552470 | 3 | −0.081304 |
| 2060564740 | 3 | +0.134770 |
| 2060021030 | 3 | +0.019895 |
| 2060441280 | 3 | −0.035093 |
| 2060020730 | 2 | +0.060142 |

Raw unweighted 8-site mean **+0.01035**, median **+0.00661**, 4 positive and 4 negative. This cannot be generalized as an ecological effect, a water nutrient subsidy or regional habitat response; no TP/I2M2 model was fit for this pilot. A fixed 2021 class mask may also select on landcover occurring between 2018 and 2022.

## Next full24 genuine source
A separate exactly source-fixed full24 H6 original same-calendar-month real pixel QA was created with print bug corrected, [workflow on Draft PR #120](https://github.com/zuizui0223/gesui/blob/analysis/gesui-20261117-onepaper-pathway-identification-decision/.github/workflows/gesui-20261117-FR263-full24H6-S2-real-pixels.yml). Pre-outcome full24 threshold >=18/24 real cloud-free pixel-supported H6. Before completion, **do not assert** the full24 result or permit an all-263 land ecology coefficient.

## Boundary and causal qualification
The current 8 pilot NDVI bands are **field-station-centered circles**, not actual lateral distance to the river. True HydroRIVERS source geometry QA on a source-pinned thirdparty v1 conversion passes only 192/263 authentic riverline station-reach geometry; official original ZIP retrieval returned 403 and stays failed. Additional source geometry confirmation and original same-source years/quality masks are required before claiming a physically realized cross-compartment mechanism. Annual satellite optical index is not biomass in mass units, not foliar nitrogen, and not biological water→land resource flux.
