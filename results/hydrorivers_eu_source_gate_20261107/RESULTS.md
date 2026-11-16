# Source-only directed stream connectivity recovery

Original [GitHub Actions success](https://github.com/zuizui0223/gesui/actions/runs/42364866919) (2026-11-07). Official HydroRIVERS server denied access (403); an independently hosted v1 EU GeoParquet was SHA256-frozen and the original HydroWASTE 58,502 WWTP source was verified against its Git blob.

- 4,713 eligible European outfall river IDs; **4,713 matched to the third-party HydroRIVERS v1-derived mirror**.
- 4,497 valid four-reach directed chains, 3,569 sites meeting pre-outcome spacing, 16 outfall units selected source-only (4 each AUT,CZE,DEU,DNK; 8 nominal high, 8 low wastewater-water-fraction).
- Selected outfall-to-matching-line offset median **0.02503 km**, maximum **0.04674 km**, verifying geographic agreement of the copied network and discharge points at these sites.
- No GBIF, plant community, nitrogen concentration or river quality was used to select these sites. Coordinates and river reach IDs are in candidate_site_design.json.
- **Critical validation limitation:** the third-party GeoParquet conversion has been hash-verified as distributed but NOT independently compared to the official HydroRIVERS original bytes/topology. This is a source-linked candidate network, not proof of environmental N exposure or insect/plant dispersal.
