# 2026-11-11 source-only NRSA join-key clarification — before nutrient/FFG ecological estimates

The source-only first EPA official download succeeded: NARS NRSA2018/19 count, site and water CSV original source SHA256:
- benthic `4e9b3a0798b25f67c7a96d4800b94f86544bbd623abaf119358e067b158324a7`
- site `1b5b2119c881290191edeb47f1d80a7711f10440e9851574cfa86fd57be770e3`
- water `cb9a3349ae3708834eab533e71423d74c1778bf76aba0737f0ae1bcc8f50baa1`

In the **first precommitted** site+VISIT_NO join, **1,220** primary BERW + positive TN/TP matches / 48 states / 9 eco groups, 97.1% FFG **positive taxon-row** coverage; gate **STOP**, due to **two** duplicate site+VISIT_NO keys. The follow-up metadata-only audit established that each duplicate composite key corresponds to **two distinct native EPA UID values** in both official site and chemistry source, with identical state/ecoregion/geographic attributes, and the matching site/water UID itself is consistent. Therefore site+VISIT_NO is not the authoritative unique observation key in the source. It was wrong to interpret the four rows as duplicated identical observations, but also wrong to declare all duplicate keys biologically invalid.

**Source-only repair rule, fixed before all new nutrient-versus-taxon effects:** Rejoin by the official **UID**, validate its uniqueness independently in each site and chemistry CSV, check site/visit agreement, and ensure `SAMPLE_TYPE=BERW` + positive original `TOTAL300` support. Keep all **distinct UID visits**, do not pretend physically repeated `SITE_ID` is independent river. Check study-year, original georeferencing, 48-state/ecoregion source and nutrients' literal units. Do not choose a preferred numerical TN/TP observation, taxon, radius, or ecological result. Continue only if original thresholds (>=300, >=25 states, >=6 ecogroups, >=60% FFG **individual-mass**, >=80% order/family **individual-mass**) are met with one-to-one UID mapping; otherwise STOP.

Note the original chem TP is **micrograms P/L** while TN is **milligrams N/L**; future regression must convert units before any mechanistic comparison. Laboratory flag/LOD handling and sub-sample comparability are independent next-stage QA. No wastewater source/delivery has been established, and no actual emergence, prey transfer, host-plant establishment or pollination service has been measured.

Audit receipts: [original gate run SUCCESS but STOP](https://github.com/zuizui0223/gesui/actions/runs/43032890833), [native UID/individual count mass diagnostic SUCCESS](https://github.com/zuizui0223/gesui/actions/runs/43033418114). This is an explicit post-source **key correction** with no ecological outcomes seen; not represented as a preregistered confirmatory test.
