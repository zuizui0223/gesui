# Source-locked HJA AS015 seasonal subsidy continuity (2026-11-12)

**Scope and timing:** Exploratory functional decomposition chosen AFTER earlier AS015 17-stream EPT fraction/mean biomass correlations, but before looking at newly constructed **within-stream EPT vs Chironomidae temporal covariance/continuity**. Same original source v2 has MD5 `f3e6e400990da883ae0d341731431749`, SHA256 `c8ae00e91e7b210fbc86f5b08c6047ceaeade30bb2d410c3748f495fea0de4b2`, and pre-outcome sample/period manifest SHA256 `51354ce5cc3ab3acc63aecc28dd93534438ab3ffaf09b86ab2b0994b51874923`. Original source acquisition GH Actions #43125631204, original 17-stream adult emergence summary #43126089434, independent validation #43126473521. No new sites, no changing source/year, no multiple end-point rescue.

**Question:** Even when adult EPT *relative* abundance does not predict total emergence biomass across streams, is the **time stability of total aquatic→terrestrial adult emergence biomass** buffered by offsetting seasonal peaks between EPT and Chironomidae (Diptera) within the same stream?

Only 17 streams, all within one HJ Andrews Oregon basin, summer 2022, each 23 interval-specific collections of four orders (E+P+T+D). This is *potential prey subsidy through adult emergence*, not actual prey capture by spiders/birds, measured pollination, sewage discharge, nutrient exposure, or global wetland service.

## Predefined calculations
Per complete native `SITE × TIME_PERIOD × SAMPLE_DATE` event:
- `EPT_t` = sum observed EMERG_BIOMASS for Ephemeroptera, Plecoptera, Trichoptera, original g dry mass/m² per collection interval.
- `D_t` = Diptera/Chironomidae EMERG_BIOMASS.
- `T_t = EPT_t + D_t`, valid measured zero retained.
- Interval emergence biomass can be zero; log(0) not used.
For each 17 site (23 periods), calculate:
1. `Cov(EPT_t,D_t)` across 23 periods with **sample covariance**; compute `C = -2 Cov(EPT,D)/(Var(EPT)+Var(D))` (positive = complementary/asynchronous fluctuations reduce total observed temporal variance relative to sum of component variances; zero = uncorrelated; negative = covariance amplifies total variance). If both component variances zero, C undefined and site is not used.
2. `φ = Var(T) / (SD(EPT)+SD(D))²` (Loreau–de Mazancourt type synchrony index, 0 to 1 when defined); values near 1 correspond to more synchronous variation; note `φ` is low also when different variation amplitudes, and does not by itself establish beneficial stabilization.
3. `CV_T=SD(T)/mean(T)`, fraction of the total 23-event emergence biomass in the **largest 5 of 23 events** (fixed ~22% time intervals), and the maximum one-event fraction; reflect pulsed availability, not adult flight duration, consumer feeding.
4. adult `EPT biomass share` summed across 23 intervals and `EPT adult abundance share` summed across periods, and their difference to quantify different size/biomass weighting.
5. Within-stream chronology check, original 23 dates exactly, natural ordering from original TIME_PERIOD (not a priori lag with other taxa). In an **exploratory deterministic circular shift** null, hold each 23-element EPT biomass series and each 23-element Diptera biomass series fixed and rotate Diptera through 23 circular time shifts, recompute `Var(EPT + shifted D)`. Report actual shift0 percentile within shifts (not an inferential test, only assesses alignment sensitivity while approximately preserving each group's autocorrelation).
6. Regime I=9 vs P=8 stream site means of C, φ, CV and pulse fraction, with no claim of independent replicate catchments; descriptive differences only. For uncertainty bootstrap streams (n=17) 10,000 times or report wide intervals, not n391 significance.
7. Across streams, correlate the mean total biomass with CV and EPT share to detect whether high biomass is concentrated in pulses; this is exploratory and can reflect common denominator.
Do not interpret `EPT` or `Diptera` as complete aquatic guilds or actual terrestrial consumer assimilation.

**Main test:** count streams with `C>0` and their median C; summary of pulse fraction; plotability from stored aggregate; report **all** signs (including if positive fraction is small or regimes differ). No ecological effect significance threshold or data dredging by selecting favorable streams.

**What needs additional original data to become an actual linked ecosystem service:** same river & date measured nutrient/source loads, larval-to-adult correspondence, riparian predator or pollinator diet, and original independent multiple catchments; DO NOT translate this seasonal decomposition into a wastewater benefit or damage.
