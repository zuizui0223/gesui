# Third NRSA cycle (2008–09): original benthic, nutrient and riverbed source-first audit

**2026-11-14. Independent-third-cycle 2008–09 source contract locked BEFORE any new numeric ecology values.** This extends [2013–14 vs 2018–19 independent survey comparison](https://github.com/zuizui0223/gesui/pull/79); preserve its critical **new-HUC8 generalization failure** and previously positive conditional source signs. Do not merge existing drafts or infer sewage causal effects.

[Official EPA NARS data download table](https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys) lists original **NRSA0809** siteinfo, CHEMW water chemistry, `NRSA 0809 Benthic Macroinvertebrate Count`, `NRSA 0809 Benthic MMI/Condition Class`, `NRSA 0809 Physical Habitat Larger Set of Metrics`. **It does not list a NRSA 2008–09 original periphyton Chla/AFDM dataset** alongside 2013/2018 periphyton datasets; therefore do NOT invent 2008 algae or set it to 0. Original source period-specific ecological comparison is `XEMBED (bed embeddedness)` and observed nutrients TN/TP vs counted benthic EPT relative share/quality, **NOT `Chla×XEMBED` replication**. The third cycle can falsify a purported invariant substrate + nutrient benthic association.

## Source-only gate
Official publisher primary URLs from the EPA table are pinned by downloaded raw file SHA:
- Site `https://www.epa.gov/sites/default/files/2015-09/siteinfo_0.csv`
- CHEMW `https://www.epa.gov/sites/default/files/2015-09/chem.csv`
- Physical `https://www.epa.gov/sites/default/files/2015-09/phabmed.csv`
- Original benthic counts `https://www.epa.gov/sites/default/files/2016-11/nrsa0809bentcts.csv`
- Benthic MMI `https://www.epa.gov/sites/default/files/2015-09/bentcond.csv`.

Raw original headers and associated publisher metadata are essential. Do not combine `BERW`/boatable methods, `CHEMB` with `CHEMW`, sites marked HAND/NAWQA with probability site without documented rules, or duplicate native UID/visit without verifying identity. Separate original 0809 MMI from native counted EPT fraction (MMI is an integrity *score*, not species richness or fish production). Original CHEMW `NTL` and `PTL` are both originally **µg/L** per [EPA metadata](https://www.epa.gov/sites/default/files/2015-09/chem.txt); this differs from 2013/2018 N source units. Standardize only with source-aware conversion.
- Pre-numeric support gate: >=300 source-valid unique UID with count chemistry physical + source HUC8, >=150 distinct hydrological groups and >=20 states. If native HUC8 absent in source and cannot independently link frozen USGS HUC8, do NOT fabricate HUC8 from coordinates, nearby reach, or county. Stop biological model and produce source QA if gate fails.
- If gate passes, original measured EPT fraction vs `XEMBED`, TP and TN, and same-era original published benthic MMI (if native UID link), full PHab, co-season and geographic predictors; independent HUC8-heldout 5fold, HUC8-cluster explanatory intervals. Important: **no direct 2008–09 EPT–Chla comparison possible from original listed data**, no downstream animal/plant ecosystem service inference, and no verified wastewater-specific effluent load.
- Never treat unreported N/P, missing benthic count, inaccessible source, or 404 as biological zero. Preserve all direction and negative/null results.
