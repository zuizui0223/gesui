# Results — same-design state/change and fixed-consumer aquatic correspondence (2026-11-25)

Protocol frozen at ef008c25f330abb0ea9ddbbe623b3cc1d045e7a1 before fitting, after the previous T3 endpoint was opened. Actual local run completed; 17 original ZIPs checked against signed digests; same exact four-window roster. Historical source original 990 + 186 plant family keys and two graphs on 46 accepted common Lycaenidae; canonical extracted Clarke 519-edge text and extra-key family text matched connector-derived FNV32 (noncryptographic transfer check). Same original prediction design and H6 folds. 16,652 common direction predictions matched sealed PR158 baseline, change, null-change and now at max 1.67e-16. New prediction CSV SHA256 df4749cc0c77a117a8027788d5428eec92ce539bcf5525a3c31c568804be0eae sealed before later species counts read by this script. Independent Ridge solver max diff 3.21e-15. No new source observations, no tuning.

## A. Controlled state/change comparison (identical design and target)
All percentages are H6-equal reduction in MSE; 5000 fixed-loss H6 resamples, nominal 95%, no model uncertainty or multiplicity adjustment. Quality 103 H7/53 H6, all 176/75; 46 species. Both graphs on identical rows. The same historical design as PR158 is used for all models, not comparing different training periods.

| Cohort | Graph | State vs baseline [95%] | State vs null-state | Change vs baseline | Both vs state | Both vs baseline |
|---|---|---|---|---|---|---|
| quality | old | 0.1259 [-0.0466,0.2950] | 0.0862 [-0.0304,0.2087] | -0.0686 | -0.0881 [-0.1750,-0.0163] | 0.0378 |
| quality | Clarke | 0.1001 [-0.0963,0.2911] | 0.0368 [-0.0667,0.1450] | -0.0152 | -0.0050 [-0.0301,0.0221] | 0.0951 |
| all | old | 0.2213 [0.0317,0.4300] | 0.1536 [0.0736,0.2420] | -0.0063 | -0.0477 [-0.0975,-0.0084] | 0.1737 |
| all | Clarke | 0.1328 [-0.0599,0.3402] | 0.0599 [-0.0265,0.1415] | -0.0214 | -0.0243 [-0.0695,0.0106] | 0.1085 |

No robust both-catalogue increment under the stronger design. State point estimates more favorable, but not a formal test of a biological constraint versus driver distinction. Preserve PR156 original frozen state forecast positive and PR158 change/coherence negative; don't interpret A as independent replication.

## B. Earlier river EQR and late resource compatibility, consumer fixed at 2017–19
On the identical A score roster, normalize 2017–19 inclusion within common 46 taxa and hold constant. For each fold train-only early plant geography/frequency/family exact mean null, resource index = weighted observed minus null adjacency times rarefied plant inclusion. Predict mean of 2020–22 and 2023–25 or their difference. Controls 2014–16 & 2017–19 index, plant S50/host mass, early consumer mass, recording effort, static geography/land use/WWTP W, EEA site number and country×method. H6-disjoint spatial validation; later plant responses used in training regions, NOT forward heldout time validation. No later butterfly composition in target, but later source eligibility and common predictor coverage affect roster.

| Sample | Graph | EQR -> late mean [95%] | EQR -> late difference [95%] | W -> late mean given EQR | W -> late difference given EQR [95%] |
|---|---|---|---|---|---|
| quality | old | +4.9005 [-0.4529,10.2203] | -1.3595 [-4.5432,1.1489] | +0.6651 | +9.4585 [-1.9272,22.7830] |
| quality | Clarke | +1.0837 [-0.0853,2.4151] | -1.3575 [-3.0360,0.1530] | -3.3862 | +7.6987 [-1.0520,17.8993] |
| all | old | +16.0645 [6.3672,26.0732] | +0.1828 [-0.9722,1.7133] | -0.5432 | +1.0796 [-7.1252,7.3699] |
| all | Clarke | -0.2002 [-3.5585,2.3941] | -0.9238 [-5.7941,4.1829] | +9.3166 | +11.9439 [1.4341,21.1850] |

Old all late mean EQR coefficient negative all 5 folds (-.0852 to -.0203), quality negative all 5; Clarke all inconsistent. The isolated nominal W positive in Clarke all late difference is not robust across catalogue/quality, W is modeled geography not measured dose. The long-term EQR hypothesis fails cross-catalogue rule. Do not promote it as general aquatic-land coupling or wastewater effect; both state and change biological mechanism remain uncertain. No new 2016/2024 EQR, no independent field validation, no causal direction. Keep nulls.

## Archive
Local reproducibility ZIP SHA256 766c780dde0a85c30f5520fd599c851295b516a705a2a25357728e03c495c5c4. 48 output rows and full fixed H6 losses, water target and coefficient tables, historical predictions, exact local source texts, source receipt and code. audit.py SHA256 07a26a75971b3cbdabf59d87532a3a1ac8d7892e3ef53f1a25c0c93a909c1a2e; results.csv SHA256 8c88a972d0fd8d3db74e2a55aba30bbe665331bf13d3e40d6279ef19883f45a4. Full numerical files and code in conversation ZIP only; this report and protocol committed. No raw source redistribution. Remote replay not run; main unchanged.
