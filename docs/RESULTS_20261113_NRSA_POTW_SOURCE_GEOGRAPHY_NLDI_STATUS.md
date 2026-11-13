# EPA NRSA 2018–2024 × POTW: source-only geography and flow-network audit (2026-11-13)

**Status: Source audit; hydrological six-shard job running at time of this entry. Not an observed sewage nutrient, biodiversity, biomass or ecosystem-service effect.**

## Source-authentic original-data results completed

Original source universe from [PR #73](https://github.com/zuizui0223/gesui/pull/73): exactly 2018/19 **1,221 original pre-outcome-frozen BERW native UID visits** and 2023/24 **2,129 original source-only matched nutrient/periphyton native UID visits**. EPA NARS SiteInfo byte-pinned originals. ECHO current weekly NPDES *individual External Outfall* coordinates verified and original national ZIP SHA retained in [Actions source artifact 43444463485](https://github.com/zuizui0223/gesui/actions/runs/43444463485).

| Fixed geography/opportunity, NOT effluent delivered | 2018–19 | 2023–24 |
| --- | ---: | ---: |
| Frozen unique native NRSA UID | 1,221 | 2,129 |
| Registered municipal external outfall straight-line within 5km | 65 | 372 |
| Within 10km | 163 | 878 |
| Within 25km | 434 | 1,626 |
| HUC8 with site within 25km | 264 | 810 |
| US states with site within 25km | 41 | 48 |
| Site metadata with original native COMID (any source distance) | **1,221** | **858** |
| Contemporary-date-*plausible* 2026 ECHO registered POTW original individual external outfalls nationally | **5,491** | **26,048** |

**Historical source selection bias**: 2026 ECHO register is not the full 2018 nor 2023 historical inventory. Permit-effective/expiration/termination dates imply *plausibility* not proof of operation or complete historical reporting. The 2018–2023 increase in geographic candidate rate **must never be interpreted as increased wastewater discharge or causal biodiversity change**. Current register missing or no nearby registered facility is neither a historical zero nor a negative control.

**Local position is not hydrologic exposure:** NHDPlus native COMID is not verified actual site coordinate on channel; 2023/24 retains only 858 native linked COMIDs. Never silently replace missing COMID with arbitrary closest modeled reach. Source-based routing only for matched COMID subset, all exclusions counted.

## API route source QA and exact identity rules

[USGS NLDI compatibility probe success #43444057539](https://github.com/zuizui0223/gesui/actions/runs/43444057539) found NLDI feature sources `npdes`, `epa_nrsa` and true COMID upstream-tributary/downstream-mainstem navigation. Example NLDI `npdes` identifiers are **permit ID + original outfall code**; EPA original facility license (9 chars) and original numbered outfall (3 chars) are not interchangeable with registered facility location. The exact byte-level permit/outfall match is mandatory, with padding of original ICIS numeric 1–3-character feature codes only according to EPA data dictionary; ambiguous or unrelated IDs omitted.

[6-shard directed topology source stage](https://github.com/zuizui0223/gesui/actions/runs/43444789943): source-only 25km geo candidate screen then 50km river UT and DM from native COMID, mark independent-upstream, downstream, same-reach ordering unresolved, NLDI API errors independently. The stage **does not** read DMR concentrations, mass, NRSA ecological values, or ecosystem services, and should not be declared successful until every shard and independent merger pass.

## Expected next source gate (not yet established)

Independent 6-shard audit verifies unique UID, unique HUC8, unique administrative permit and **exact permit/outfall feature**. To justify actual DMR nutrient assessment require >=40 distinct source-upstream NRSA visits, >=25 independent outfall features, >=10 HUC8, >=5 states; then report exact historical native DMR 00600(TN)/00665(TP) outfall/year source and true physical units, never `NODI` as zero. This is still only **reporting opportunity**, not verified downstream dose, until actual measured stream flow/hydrolocation are connected.

Original existing [PR #73](https://github.com/zuizui0223/gesui/pull/73) EPA instream TN/TP observed effects (EPT2018 relative-count fraction negative TN, 2023 areal AFDM weak/null, Chla2023 negative TP) are never relabeled as wastewater-specific outcomes. No eBird; no merged main or other independent branch touched.


## 2026-11-14 independent six-shard QA — COMPLETE (replaces previous 'running' status)

[All six official USGS UT/DM jobs SUCCESS #43444789943](https://github.com/zuizui0223/gesui/actions/runs/43444789943); [separate independent complete-population/exact ID merge SUCCESS #43449175435](https://github.com/zuizui0223/gesui/actions/runs/43449175435). The merger compared all original native NRSA UIDs and every exact permit+feature against the pre-outcome 25km manifest, and refused absent or duplicated shards.

| Pre-DMR true directed NLDI registered source | 2018–19 | 2023–24 |
| --- | ---: | ---: |
| Original source-blind 25km + native COMID candidates queried | 434 | 647 |
| USGS both UT and DM success | 432 | 645 |
| API unavailable (not zero) | 2 | 2 |
| **Distinct NRSA visit UIDs with at least one upstream DIFFERENT-reach registered municipal EXO NPDES source** | **36** | **195** |
| Verified upstream site×outfall source pairs | 61 | 815 |
| Unique upstream EXO feature pairs (permit+feature) | 55 | 657 |
| Unique upstream administrative permits | 50 | 402 |
| HUC8 containing source-associated NRSA sites | 30 | 139 |
| States with source-associated NRSA sites | 12 | 41 |
| Site-visits with genuine upstream outfall within 10km straight distance | 18 | 89 |
| Downstream different-reach source-site candidate pairs | 302 | 642 |
| Ambiguous overlapping UT+DM source pairs (excluded from upstream) | 12 | 23 |
| Frozen **minimum 40 site UIDs /25 outfalls /10 HUC8 /5 states before historic DMR source scan** | **FAIL** | **PASS** |

Scientific decision:
- **2018/19 STOP:** n36 is below the original minimum 40; do not relax this after seeing sign or claim true WWTP-specific stream nutrients. Retain existing 2018 EPT/TN spatial observational result separately.
- **2023/24 PROCEED TO HISTORIC DMR SOURCE-ONLY:** n195 at 657 potential individual outfalls geographically and topologically linked to actual NRSA station COMIDs. Importantly this is still a 2026 registry, and NLDI NPDES reach is a modelled source locator, not verified discharge during 2023 or 2024. All 815 *site×source* candidates must NOT be counted as independent sites, and 657 features are NOT 657 reporting in-period.
- **Not yet established:** actual year-specific 2023/24 monitoring period N/P mass or water discharge, effective effluent dose, 40 truly measured source-positive sites, 40 actual comparable confirmed no-discharge sites, or any source-induced EPT/biomass/pollination outcome. The true DMR report may fail after the source-only gate.

The full original fixed site-to-upstream-feature candidate identifiers live solely in the Actions artifact, not copied into publication-ready biological claims. Same-reach ambiguity and two-direction duplicated routes are excluded, and 2+2 missing USGS results remain NA.

## 2026-11-14 authentic 2023–24 historical DMR: SOURCE GATE FAIL for same-feature N+P, no ecological WWTP coefficient

The source-only historical effluent monitoring question was investigated exactly on the original source-frozen 195 NRSA visit UIDs, 657 NLDI/ECHO-verified individual upstream EXO registered outfalls, with original corresponding 2023 or 2024 visit year and original original 2026 ECHO registry caveat. Original EPA national ICIS-DMR FY2023, FY2024, FY2025 full ZIPs were downloaded and streamed; their original zip SHA256 and entire provenance are in [full fiscal source SUCCESS #43449939212](https://github.com/zuizui0223/gesui/actions/runs/43449939212). Calendar periods are obtained from `MONITORING_PERIOD_END_DATE`, not wrongly equating a fiscal-year archive with a calendar year. Original EPA **permit + external-outfall feature** and gross/net effluent monitoring location are jointly required, and NODI/non-numeric never becomes reported zero.

[Exact native source + 3-fiscal independent join and unit-type separation SUCCESS #43450461565](https://github.com/zuizui0223/gesui/actions/runs/43450461565):

| Original NRSA visit-year / source-side outcome (NOT ecological response) | 2023 | 2024 |
| --- | ---: | ---: |
| Source-frozen real NRSA 2023/24 survey UIDs with potential true upstream *registered* municipal EXO | **91** | **104** |
| Distinct candidate source EXO features | 283 | 374 |
| Individual source features reporting numeric, strictly positive **effluent N** in **concentration type** | 36 | 68 |
| Individual source features reporting numeric, strictly positive **effluent P** in concentration type | 85 | 110 |
| **Same exact outfall** with positive reported **N AND P concentration-type** | 35 | 66 |
| Actual NRSA visit UIDs connected to a source with both positive reported N/P concentration-type | **25** | **38** |
| Same exact outfall with positive N AND P **quantity-type** report | 20 | 30 |
| NRSA visit UIDs connected to a source with both positive reported N/P quantity-type | **13** | **17** |
| Original predeclared minimum 40 NRSA visit UIDs for verified source-specific N/P | **FAIL** | **FAIL** |

**Quantity-type DMR value is not automatically discharge mass per year:** underlying `VALUE_TYPE_CODE`, `DMR_UNIT_DESC`, reporting frequency, monitoring location and statistical base must be inspected before mass standardization. Concen­tration reports remain concentration. Original zip data sources completed with 26,287,632 rows FY2023; 27,018,548 FY2024; 27,011,856 FY2025, respectively, not selected by downstream biological response. The report covers only original 657 registered outfall candidates with verified up-river NLDI route and exact-year NRSA visit, not every available historical WWTP.

### Scientific decision

- Original **2018–19** NLDI source candidate UID n36 fails the *pre-DMR* minimum40. Stop original N/P source effect test.
- Original **2023–24** NLDI source candidate UID n195 passes the *pre-DMR* gate, but same-source *both N and P reported in same site-year* fails in both original survey years (2023 n13 quantity / n25 concentration, 2024 n17 / n38). Thus **full original WWTP measured N/P→in-stream nutrients→AFDM/Chla/EPT biodiversity causal or mediation models are NOT FIT**.
- A missing/reporting-zero source is **not** a verified no-wastewater site. The previously defined minimum40 properly monitored non-exposed comparator sites is also **not established**.
- The earlier 2018/19 TN↔EPT and 2023/24 TN/TP↔AFDM/Chla results are preserved as **measured ambient river nutrients**, not attributable WWTP effluent effects. No increased-biomass or ecosystem-service benefit/harm conclusion follows from the new DMR coverage.
- **Prospective route:** source-complete discharge records and exact recipients, properly reported quantity units, actual stream discharge/dilution, and measured in-stream N/P+biological endpoints on enough common visits. Do not tune buffers, quality rules, environmental endpoints or falsely pool 2018 with 2024 merely to pass positivity.
