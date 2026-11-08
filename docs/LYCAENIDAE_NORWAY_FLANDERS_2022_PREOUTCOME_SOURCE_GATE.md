# Open Norway and Flanders standardized butterfly sampling events × 2022 WWTP: pre-outcome source gate (2026-11-08)

**Pre-outcome protocol; no butterfly occurrence, count, family response, puddling image or EQR numerical score read.** New independent country source-coverage audit, distinct from 2022 Sweden and opportunistic iNaturalist pilots. Original `chocho` 22/25 Lycaenidae potential introduced-host geographic resource expansion and `gesui` 2018/2022 EQR predictive gains are held fixed and not modified.

## Source files and licences (publisher original only)

- **Norway:** Norwegian Institute for Nature Research, *Bumblebees and butterflies in Norway*, GBIF dataset `aea17af8-5578-4b04-b5d3-7adf0c5a1e60`, CC BY 4.0. Official original Darwin Core Archive `https://ipt.nina.no/archive.do?r=butterflies_bumblebees2020`. Published release dated 2026-05-04, stated temporal coverage 2009–2025; area-representative Lucas-grid network, 89 grid sites total, each with 20 sections of 50m, ~three visits per summer. Original entire ZIP checksum to be recorded; do not assume all published events occur in 2022 or all field positions are precise.
- **Flanders:** INBO *Meetnetten.be - Transects for butterflies in Flanders, Belgium*, GBIF dataset `63938753-1fec-4c08-ae39-e9f8a6576521`, CC0. Original DwC-A `https://ipt.inbo.be/archive.do?r=meetnetten-butterflies-occurrences`. Published 2026-02-27, records until 2023-05-20. Survey of priority species plus other observed butterflies. **Sensitive original sites are generalized to 1, 5 or 10 km UTM cells**, so point-distance exposure must use published coordinate uncertainty. Never try to de-obfuscate.
- **Flanders historical secondary:** INBO *Butterfly monitoring scheme in Flanders, Belgium*, GBIF dataset `81b909c7-5612-4148-8e70-99fea2cdb572`, CC BY-NC 4.0. Original `https://ipt.inbo.be/archive.do?r=bms-occurrences`, ends 2020; **2022 source-year match unavailable**, so metadata-only historical comparator, not 2022 replication.

## Stage A (before any 2022 WWTP matching)

1. Fetch only publisher-original public ZIP; lock downloaded SHA256, `meta.xml` SHA256, declared licence, original Event `rowType`, `core` source, event date/site/protocol/geolocation/uncertainty fields and archive member inventory. **Never read or extract occurrence extensions or individual butterfly counts.**
2. Treat `eventID` as unique event identifier and nested `parentEventID` as potential grid/route hierarchy, **not** automatically independent sites. Keep full exact hierarchy class metadata (one-way hashed if necessary) only in derived QC. A visit must have a specific YYYY-MM-DD date; repeated survey events at same grid/route on different dates establish repeat visits. Do not count a root Event/year period or overlapping sections as independent replicates.
3. Source-only independent monitoring units: original nonblank stable grid/location ID when supplied; otherwise one coordinate location rounded to 4 decimal places + protocol and parent identifier, marked **provisional**. Count raw 2022 events, distinct 2022 dates, spatial footprints, grid/route groups and metadata-defined coordinate precision. If no independent grid site ID exists, report uncertainty and do not promote 20 sections to 20 biological replicates.
4. Strict 2022 eligibility: at least **2 distinct dated visits** in the same monitoring unit; published public coordinate uncertainty **≤500 metres**, non-obscured; group original site coordinates consistent to 500m. Missing source accuracy is **unknown**, not zero. Stated coarse Belgian 1–10km public generalization fails 500m whenever documented uncertainty exceeds 500m. Do not attempt to reconstruct withheld locations.
5. Lock the source-only station inventory and site-key/visit hashes and record counts *by* national scheme and country. No automatic long-term species abundance inference.

## Stage B (only after stage-A site information is archived)

1. Use EEA official **2022** source REST with `uwwNDischargeMeasured > 0 AND uwwPDischargeMeasured > 0` and unambiguous **one** published discharge point, `dcpWaterBodyType=FW`; original reporter+plant key exactly joined. Source 5,492 eligible reporting records and 24,416 discharge points before joins (locked historic result). **Do not assume national reporter ID prefix conventions**; report code frequency and geographically plausible country-coded rows to determine official valid prefixes before looking at butterfly species outcomes. Discharge observed mass is reported annual output, not an independently verified same-date chemical assay.
2. Compute geodesic distance from each independent 2022 repeat-visit site to its country's **2022 source-only qualified freshwater discharge points**, never to a WWTP building. No river connectivity assumed for terrestrial butterflies, and air/soil nutrient exposure is **not** proven by Euclidean distance.
3. **Fixed primary near ≤2km versus registered-source remote ≥10km**; prospectively fixed second near ≤5km sensitivity. Apply full position uncertainty ≤500m to both, plus a conservative boundary margin: near if center_distance + 0.5km ≤2/5km, far if center_distance −0.5km ≥10km; leave uncertain boundary sites **unclassified**. Do not select nearest from reported outlets of a different country or include coastal/nonfreshwater or multiple-outlet facilities.
4. The unit of replication is an **independent Lucas grid / transect monitoring site**, not 50m sections or 3 visits. To permit an outcome model require **≥10 independent strict 2022 near sites and ≥10 far sites**, each with ≥2 dated visits, across ≥2 matched subregions/counties with both near and far available. For Norway 89 grid squares total, report available 2022 independent grid count; for Flanders report coarse-location exclusions explicitly. The same fixed rule applies to both countries.
5. **STOP** if source accuracy, 2022 year coverage, reporting-source availability, or site matching gates fail. A 2020-only Flanders release cannot substitute for a 2022 validation, and cannot be silently paired with 2022 WWTP mass.

## Stage C conditional, not authorized by source discovery alone

Only after BOTH A and B succeed and results are committed BEFORE opening any `occurrence.txt` butterfly identities/counts: prespecify per-visit effort, taxa and detectability, land cover/host plant background, spatial blocks, sampling-frame eligibility and negative controls. Compare actual Lycaenidae counts standardized by visits and total butterflies, with within-country/source-only controls. Even if positive, no adult protein/Na puddling, ant mutualism or wastewater causal nutrient claim without their direct independent measurements.

## Exact source references

- https://www.gbif.org/dataset/aea17af8-5578-4b04-b5d3-7adf0c5a1e60
- https://www.gbif.org/dataset/63938753-1fec-4c08-ae39-e9f8a6576521
- https://www.gbif.org/dataset/81b909c7-5612-4148-8e70-99fea2cdb572
- https://github.com/zuizui0223/gesui/pull/21

**Scientific status:** new independent source-only feasibility design, NOT a confirmatory test of wastewater ecological effects.
