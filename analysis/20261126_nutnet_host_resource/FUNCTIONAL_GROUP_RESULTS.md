# NutNet functional group and host identity — retrospective extension
2026-11-26. Original source Chen et al. 2025, DOI 10.1038/s41467-025-59166-7; Figshare v4 DOI 10.6084/m9.figshare.26412295.v4 (CC BY 4.0). Original source artifact 15085929755 SHA256 23dfc003ad8fe09399d52e24cc3b3e27c6a3fe853b1ac87a1a847bcfdf6f64da. All nine member checks verified. Previous baseline DiD and EIV quality STOP unchanged.

Native groups: 72 sites, 70 exact coordinate pairs, 216 matched blocks, 444 plots, 4934 positive-cover rows. Europe18 subset: 16 coordinate pairs, 54 blocks, 114 plots. Each site has three source blocks (not necessarily labelled 1,2,3). Average plots within arm and blocks within site; sites equally weighted. Six blocks have two plots per arm. Native labels GRASS, LEGUME, FORB, GRAMINOID, WOODY, CACTUS, UNKNOWN (103 missing group rows retained), not exact plant families or biomass. Site bootstrap 20000, coordinate bootstrap 10000; 19999 conditional paired block label swaps with Holm adjustment across the three focal native cover endpoints. Posthoc relative to previously viewed aggregate outcomes.

| Outcome | Global72 contrast [nominal 95%] | Europe18 contrast [nominal 95%] |
|---|---|---|
| total cover | +10.05352 [5.03353,15.10625] | +10.83519 [.66102,21.09449] |
| GRASS | +8.54625 [2.84459,14.30302] | +19.70648 [10.19056,29.95188] |
| LEGUME | -1.17731 [-3.45903,1.22179] | -3.64074 [-8.81676,.39630] |
| FORB | +3.04542 [-1.61832,8.05220] | -.86296 [-9.97051,8.12252] |
| WOODY | -.14898 [-2.14082,1.74698] | -4.54630 [-10.62597,-.45370] |
| GRASS minus remaining cover | +7.03898 [-2.77966,16.86048] | +28.57778 [9.88509,48.23718] |
| total label richness | -2.10880 [-2.96759,-1.35185] | -2.61111 [-4.59259,-.83333] |

GRASS Holm p=.0012 global, .00015 Europe; 49/72 and 15/18 positive. ~85% of global net gain is GRASS, accounting ratio not a statistically established global dominance contrast. European GRASS share +.106408 [.020785,.183964], global share uncertain; global grass-minus-rest contrast also uncertain. Europe is nested, not independent validation.

Identity-only host endpoint distinct from EIV: fixed prior 94 consumers, all accepted documented host edges 653/1250 old/Clarke instead of scored-only 484/966. Median exact-key cover capture .97986577/.99084859 control/NPK. 47 matched blocks and 13 sites/7 countries meet new taxonomy-only 70% cover criterion; EIV quality remains zero sites. Full18 primary and quality13 descriptive post-treatment subset, selection risk. Host-union cover +7.59/+9.19 full18 (both uncertain). Mean host cover per listed host full18 -.01678 [-.12596,.07331]/+.01462 [-.05709,.07830]; quality13 -.04098 [-.18853,.08096]/-.00458 [-.09675,.07567]. Fraction listed hosts recorded full18 -.003407 [-.007737,.000556]/-.004002 [-.007448,-.000966]; quality13 -.005676 [-.010678,-.001323]/-.004888 [-.009247,-.001092]. Nominal joint breadth decline only in post-treatment selected subset, not established on full18. Unmapped unknown; outside union not verified nonhost. No local consumer census.

Exploratory fixed strict guilds: both catalogues all accepted hosts Poaceae-only 12 consumers (11 Nymphalidae), Fabaceae-only 14 (9 Lycaenidae). Full18 normalized host cover: grass-only +.104622 [-.106743,.333184]/+.128579 [-.169986,.470576]; legume-only -.136684 [-.494180,.087259]/-.071634 [-.234040,.029137]. Paired difference +.241307 [-.089603,.702494]/+.200213 [-.086257,.509591], both Holm p=.2119. Quality13 +.295347 [-.150436,.921314]/+.151442 [-.206299,.543759], both p=.3730. Direction consistent, differential response unsupported. Family and guild confounded, no actual feeding or butterfly population observations.

124 group/identity rows and 24 guild rows. Independent 56 point checks max absolute error 3.55e-15; prior scored endpoint replay <4e-15 and equivalent host weight sums <9e-16. Independent check initially assumed all block labels 1,2,3; corrected to actual source labels, no main estimate change. Local reproducibility ZIP gesui_20261126_nutnet_groups.zip SHA256 8afee81fe0280c2ddbd265524287da2da76cc598cfad573159ed4db45c20f6a4. Original archives not rehosted. No remote numerical replay. NPK with initial micronutrients is not sewage or separate nitrogen/phosphorus effects; cover is not biomass. No main merge.
