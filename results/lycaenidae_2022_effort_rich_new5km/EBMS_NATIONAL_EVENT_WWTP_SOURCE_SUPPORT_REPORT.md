# eBMS national public 2022 standardized visits × official EEA WWTP source compatibility: independent audit

**2026-11-08, source-only results.** Public licensed SeBMS and Denmark eBMS DwC-A sources were downloaded and hashed. **No butterfly species-level occurrence/count records were read.** Do not confuse source feasibility with no ecological wastewater effect.

## 1. Swedish/Danish national published sampling-event data exist

| Public archive | Original source ZIP SHA256 | Original full Event rows | 2022 Event rows | Unique 2022 public coordinate points |
|---|---|---:|---:|---:|
| Sweden SeBMS CC0 | `384a2e5145df29268fe2de79fde9e196fc113e68adfc922a81acd58069437348` | 39,466 | **2,677** | **448** |
| Denmark eBMS CC-BY-NC | `2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b` | 15,606 | **125** | **22** |

Both official originals contain a Darwin Core `Event` core and original event dates/coordinates, without needing to open the occurrence/species-count extension. Sweden has both `Fixed route Pollard walk transect` and `Point site count` protocols; Denmark's `samplingProtocol` and explicit `samplingEffort` fields are absent in the published event core. Central eBMS data access requires coordinator approval; no restricted data retrieved.

Sources: [Swedish published source](https://www.gbif.org/dataset/be77e203-486c-4651-91b9-8347968b728c), [Danish source](https://www.gbif.org/dataset/caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c), [verified public Event source-only audit](https://github.com/zuizui0223/gesui/actions/runs/42454579954).

## 2. Original EEA 2022 WWTP freshwater source country code repaired before species outcomes

Original public EEA **2022** source REST reported 5,492 WWTP plants with both positive measured N/P discharge annual masses and 24,416 discharge points; source-only combination and filtering retained unique direct-freshwater outfall codes.

Original Swedish report code was **`SE369`**, not `SE_` or `SE-...`; this was discovered from EEA original reporter-code distribution at source-only stage, not after opening SeBMS species counts. The initial source filter therefore returned 0 but was corrected to exact `SE369`. Corrected original 2022 public report layer has **14** qualifying unambiguous Swedish FW outfalls in the study extent. The identifier itself has **no 2022 suffix**: its individual laboratory sampling/reference year is not independently verified, although drawn from the EEA 2022 REST release.

[source-only EEA×SeBMS corrected SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/42455098633); [full source-only corrected QC JSON](https://github.com/zuizui0223/gesui/blob/analysis/lycaenidae-effort-rich-river-cohort-20261108/results/lycaenidae_2022_effort_rich_new5km/SEBMS_2022_SITE_VS_EEA_WWTP_SOURCE_ONLY_QC.json).

## 3. The prespecified 2km/5km versus 10km source proximity gates fail — STOP before butterfly counts

Strict Swedish public site-quality gate: 2022 in fixed Europe bounding box, **≥2 independent dated visits per original source location**, a single consistent source coordinate and **reported uncertainty ≤500m**.

- Out of 2,677 original Swedish 2022 event rows, original protocols were **864 fixed Pollard-walk transect visits** and **1,813 point-count visits**. Under the specified 500m uncertainty and repeated-visit gate, **no fixed Pollard-walk site** qualifies; 2022 source-only position uncertainty is too coarse for a 2–5km proximity comparison on those routes.
- **Point-site surveys:** 226 distinct source-quality sites / 1,720 visits survive the same quality gate.
- From 14 original qualifying Swedish 2022 EEA measured-N/P FW outlets, nearest-point comparison at those 226 sites yields **0 within 2km**, **2 within 5km**, **215 at ≥10km**. The intermediate 5–10km band comprises 9 sites.
- Predeclared source-only support gate: **≥20 near ≤2km and ≥20 far ≥10km**, across ≥2 common counties. Primary **FAILS**. Predeclared ≤5km secondary gate (≥20 near) also **FAILS (2 sites only)**.

No species-specific SeBMS occurrence extension, Lycaenidae counts, abundance/per-visit test, puddling or sewage biological effect was calculated. The population of Swedish public transects close to available *reported measured N/P freshwater outlets* is too thin for this source definition and temporal support; this is not evidence that wastewater has no effect.

## 4. What is preserved and recommended

- The historical `chocho` 22/25 Lycaenidae expansion in **potential introduced-host resource geography** remains supported and separate. Its values are not "sewage-driven butterfly growth."
- The `gesui` 2018/2022 bottom-fauna biological-EQR prediction gains are preserved and refer to entirely different monitoring units.
- The previous 1km photo opportunity support failure (1/24) and the new 5km 2021-plant-effort-screened 2022 photo support failure (4/17) both remain in the evidence ledger.
- The next defensible way to seek a new wastewater–Lycaenidae *biological* result is **independent transect effort and count records in countries with sufficiently reported verified WWTP discharges and matching source-year**, with independent urban/land cover and host plant availability. Central eBMS raw count access requires a scientific data request; alternatively locate openly licensed national event datasets with enough same-year visits and source overlap. Do not loosen the 2022 Swedish proximity cutoffs after this failed support gate.
- Total N/P reported effluent mass is not organic protein, adult sodium drinking, plant nitrogen uptake, ant partner abundance, or local nutrient concentration. None is identified here.

All results stored on stacked Draft PR #21; no main merge or modification of chocho.
