# 2026-11-16 Austria GZÜV official water monitoring station/source audit, 2022 common data access

**Verdict: authentic public national station coordinates, exact IDs and crosswalk PASS. Same-year 2022 measured river TN+TP / macrozoobenthos / wastewater receiving load / riparian vegetation source chain remains UNVERIFIED. No new ecology coefficient fitted.**

## Why this matters for the single-paper goal

Before this work, original gesui PR116 showed a genuine source-year 2022 reported N+P facility outfall / river line/remote sensing study in AT+DE: 1,156 official source facilities, 282 physically validated river paths (AT145, DE137); 48 geographic independent satellite pilot, 26 source systems with seasonally comparable 2018/2022 local sensor pixels. Yet only one of 282 outfalls had a nearby real 2022 EEA benthic site and WISE6 2022 TP site within 2km (zero verified same receiving source reach). Missing same-location aquatic data is a source support problem, not an ecological null. Separate French actual TP→I2M2 and broad-European plant nutrient NUTRI→butterfly HNS results cannot be concatenated into a causal wastewater food chain.

## Official BMLUK and national Umweltbundesamt source

- Government WISA (https://www.bmluk.gv.at/themen/wasser/wisa/daten/daten.html) explicitly advertises public crossdomain data querying H2O physicochemical river-quality data, macrozoobenthos/phytobenthos and EmReg-OW emission load aggregates by monitoring stations. H2O five-step public quality-query app (https://wasser.umweltbundesamt.at/h2odb/fivestep/abfrageQdPublic.xhtml) and WISA frontend require JavaScript; accessible pages do not demonstrate a documented automated export of original 2022 observations. Aggregated EmReg source at stations is not a single source-specific wastewater delivered mass.
- The *GZÜV Messnetz Fließgewässer* public dataset was independently discovered through documented EU public data catalog search, original BMLUK Austrian government dataset ID b054ce75-f0b1-449e-b11c-9603a629c9c5. Public original distribution is Austrian Environment Agency share https://docs.umweltbundesamt.at/s/SNJqq7TEprxJtpJ exposing genuine **FWMST.gpkg**. The public site legitimately returns a single file with SQLite magic rather than a ZIP even though source public menu URL mentions accept=zip; no unpublished API was accessed.
- Original authentic GeoPackage SHA256 **4b87a495cd71bf2c413494b80158e8cfc91bb0de3e6e2559c186298388986fcf**, 675,840B, original CRS EPSG:31287, official FWMST point layer **3,788 national water-course monitoring station records**. Each original station has unique government PK source ID and unique internationally compatible PKINT ID, e.g., PK FW10000027 / PKINT ATFW10000027. The original field BEGIN_DAT is identically **2025-12-09** across source rows; it does NOT confirm 2022 monitoring or prove a station start date (likely dataset export metadata). Raw third-party geo file was not copied to Git or Actions artifacts. [Original GZÜV archive and location verification success #43863117767](https://github.com/zuizui0223/gesui/actions/runs/43863117767). Previous public metadata search and real share access [#43861131035](https://github.com/zuizui0223/gesui/actions/runs/43861131035), [#43861576087](https://github.com/zuizui0223/gesui/actions/runs/43861576087) and [#43861925144](https://github.com/zuizui0223/gesui/actions/runs/43861925144) SUCCESS.

## Original 2022 annual-report-positive outfalls to genuine historical GZÜV station coordinates

Use exactly original **145 AT** facilities with both actual 2022 annual reported measured N and P (original signed #43856642187), outfall physical river path already source validated and not chosen based on response strength:

| Great-circle outfall radius | Distinct source-positive 2022 facilities with any historical GZÜV station | Distinct official stations within all buffers |
|---|---:|---:|
| 500 m | **41/145** | 46 |
| 1 km | **77/145** | 103 |
| 2 km | **117/145** | 225 |
| 5 km | **144/145** | 710 |

These are point-distance source opportunities only, NOT directed same-reach 2022 measured TN+TP, actual benthos, treated water delivery or controls.

## Exact station identity against original EEA2022 actual benthic site IDs

Original signed in-main source real 2022 EEA benthic site registry Git blob SHA1 325f199e985d9862ffee5b0edb205c61a7be3833 has genuine Austria **68** benthic station IDs. Joining EEA code ATFW... to the original BMLUK GZÜV PKINT verifies **68/68 exact station-ID matches**, independently of numeric ecological scores. Government GPkg-to-EEA site coordinates have median **0.643 m** difference; **67/68 within 5m**. One physically important coordinate discrepancy remains: PKINT **ATFW54110117** is offset **954.614 m** between EEA2022 and national station geodata. Resolve original source coordinate correctness before considering that station co-located; do not silently promote same ID to the same 250m river-bank patch.

When counting AT145 source-positive outfalls within2km of the national *GZÜV stations known to overlap these EEA68 2022 benthic sites*, 2 AT facilities meet geometric proximity; original EEA2022 site geographic coordinates yield only 1 AT facility near source. Neither demonstrates being on the true physical downstream river reach, and at least one location source differs by ~955m. This is not a new robust wastewater field sample. [Original ID and coordinate QA run #43864306874 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43864306874).

## GO / STOP criteria and next exact source acquisition

**GO as an original national monitoring station registry + EEA interoperability dataset:** actual station spatial mapping and stable code linkage are proven. This justifies a targeted original 2022 chemistry/MZB table source investigation in WISA/H2O using GZÜV PKINT.

**NO-GO as an empirical 2022 wastewater-to-river-to-terrestrial causal outcome paper:** National source map 3,788 stations are historical locations, not actual 2022 sampled stations; national WISA water/biota per-PKINT original 2022 dates/values are not yet acquired. Annual EEA outlet 2022 tonnes/year is not measured receiver-water TN/TP, EmReg station average is not single WWTP attribution, satellite NDVI not plant biomass or biodiversity. Full 2022 receptor/source, landcover, no-interposed-discharge reference river and true independent hydrologic blocks must precede causal model.

Next source test: use only publicly documented H2O five-step official export or a permitted WISA cross-domain data download to retrieve a **source-only** original station-PKINT×2022 table (presence/metadata of TN+TP samples, macrozoobenthos sample dates, source sample methods) for original **225 national station IDs** found within2km of original AT145 positive outfalls. No numeric ecological response fitting until exact 2022 source counts and downstream connectivity pass. Source can be requested from Austrian Environment Agency at the official published national surface-water contact if export is interactive or restricted; remain UNSENT without user authorization. Missing original data is a source STOP, not an ecological null or proof data do not exist.

The original 4 source action scripts and exact GeoPackage SHA checks are on this review-only branch; main and chocho untouched.


## Source-to-request actionable 117×225 official PKINT crosswalk (NEW)

[Original official river GPkg/public download, 145 source annual positive N/P outlets, and compact source-only 2km PKINT pairs verified in Actions #43865663130 SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43865663130). This generated an inspectable **225-row CSV** mapping **117 original reported-positive EEA2022 outfall facility IDs** to **225 distinct government GZÜV internationally prefixed PKINT station IDs** (exactly one pair per station ID in this original geometry). Deterministic derived-source CSV SHA256 **040d610b6172d165787775406ea50c0f73f34db1db3fcb6c8576898736146dde**. Fields include original outfall ID, genuine PK, PKINT, geodesic outfall-to-station metres, actual 2022 facility reported annual N and P tonnes, and explicit two false-to-infer flags. The original full external GeoPackage is not rehosted. Source-list counts are identical to the previously verified 41/77/117/144 facilities at 0.5/1/2/5km.

The 225 official **PKINT** station IDs form a bounded original source request or legitimate public WISA/H2O query set for genuine **calendar 2022 TN+TP sample-day and macrozoobenthos source availability**. The CSV itself is NOT a measured 2022 biology/pollution time series, not a sampling absence table, and not a cross-ecosystem causal result. A statement that 117 WWTPs actually affected 225 monitoring stations would be false because only geographic proximity, not exact directed receiver transport, is verified.
