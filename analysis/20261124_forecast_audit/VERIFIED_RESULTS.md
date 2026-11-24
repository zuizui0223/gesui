# Verified third-period consumer forecast and posthoc influence audit (2026-11-24)

Original protocol and model unchanged. Historical-only train and later scoring completed [run 45394379127](https://github.com/zuizui0223/gesui/actions/runs/45394379127), prediction artifact 14920242609 (SHA256 a7745f667fd27f6185d7e1ce3fd537612967e669621f632fbc263707c309a7b9), score artifact 14920636600 (SHA256 59113276d8e53f3ab06550fd614bfaaef44ebf982059c8a727768499fac8fbfe). Prediction file SHA256 43f6a31e4e6855977997cf0d071d996958c1bcc9c81d0073f74f18c2b400c394, 63,910 prediction rows sealed before T2 scoring. Four source shard ZIPs, gate ZIP and original host ZIP also downloaded and all eight archive hashes verified. Independent local scipy hypergeometric survival computation reproduced all 1092 H6 loss rows, max absolute error 1.9831e-13. 217 eligible H7/88 H6 and quality 143/70; common Lycaenidae 195/84 and quality 133/67. No third-period refitting.

## Prespecified outcome, percentage reductions in MSE
| Cohort | Family/catalogue | baseline gain [95%] | structured null gain [95%] | persistence gain [95%] |
|---|---|---|---|---|
| quality | common Lycaenidae old | 0.3259 [0.1526,0.5040] | 0.2037 [0.0624,0.3628] | 7.6340 [0.2382,13.9905] |
| quality | common Lycaenidae Clarke | 0.2615 [0.0492,0.4643] | 0.1404 [0.0047,0.2902] | 7.5743 [0.1424,13.9276] |
| quality | Hesperiidae old | 0.2067 [0.0400,0.4025] | 0.1825 [-0.0168,0.4200] | 7.5882 [-0.7300,15.9083] |
| quality | Lycaenidae old | 0.2790 [0.1037,0.4559] | 0.2064 [0.0638,0.3557] | 6.4141 [-0.8900,12.4884] |
| quality | Nymphalidae old | 0.0735 [-0.1847,0.3522] | 0.0659 [-0.0790,0.2061] | 0.1917 [-7.4388,6.8332] |
| quality | Papilionidae old | 3.0400 [1.4696,4.7206] | 3.2884 [1.8021,4.8706] | -6.4818 [-20.3473,2.6335] |
| quality | Pieridae old | -0.3455 [-1.0194,0.3037] | 0.0200 [-0.3273,0.3435] | 7.9747 [1.2821,13.1435] |
| all | common Lycaenidae old | 0.4110 [0.1593,0.6722] | 0.2835 [0.1187,0.4518] | 1.5718 [-6.4795,8.5237] |
| all | common Lycaenidae Clarke | 0.3620 [0.0678,0.6848] | 0.1929 [0.0286,0.3604] | 1.5233 [-6.5603,8.4725] |
| all | Hesperiidae old | 0.1046 [-0.1163,0.3129] | 0.1117 [-0.1367,0.3548] | 5.1622 [-4.0923,13.2539] |
| all | Lycaenidae old | 0.3821 [0.1400,0.6269] | 0.2766 [0.1084,0.4508] | 1.4483 [-6.4981,8.1048] |
| all | Nymphalidae old | 0.0238 [-0.2564,0.3114] | 0.0208 [-0.1213,0.1596] | -0.5550 [-7.2403,5.3053] |
| all | Papilionidae old | 2.9394 [1.6524,4.4156] | 3.0516 [1.8704,4.4589] | -5.9902 [-18.1123,2.5784] |
| all | Pieridae old | 0.2344 [-0.5033,1.0510] | -0.1080 [-0.4195,0.1937] | 7.3814 [0.9111,12.7297] |

Primary intersection passes both common catalogues against baseline and structured null, though small gains and Clarke lower bound close to zero. New endpoint is future butterfly indexed species inclusion, not earlier EQR-network endpoint (Clarke non-replication retained). No independent wastewater effect or tested causal lag.

## Posthoc audit (fixed predictions, no refit)
Paired H6 bootstrap 5000 seed 20261124, Bonferroni-8 within each cohort (4 comparisons x2 references), same H6 support, family minus Lycaenidae. Papilionidae positive in both cohorts, others uncertain. Quality baseline +2.7610 pp [0.5649,5.2715], structured null +3.0820 [1.0226,5.3911]; all +2.5573 [0.7213,4.7202] and +2.7750 [1.0712,4.8435]. Targets differ by species and baseline MSE; gains not directly comparable ecological effects. 10 Papilionidae vs 52 Lycaenidae, 26 Hesperiidae, 97 Nymphalidae, 20 Pieridae. 

Species leverage: average error within H7, H6, then equal H6 by species; rank absolute loss reductions, recompute after deleting high contributors. Quality:
| Panel | comparator | original gain | remove top 1 species | remove top 3 species |
|---|---|---:|---:|---:|
| Papilionidae | baseline | 3.0400 | 0.1626 | -0.9090 |
| Papilionidae | structured null | 3.2884 | 0.4891 | -0.7932 |
| common Lycaenidae old | baseline | 0.3259 | 0.2825 | 0.2031 |
| common Lycaenidae old | structured null | 0.2037 | 0.1607 | 0.1008 |
| common Lycaenidae Clarke | baseline | 0.2615 | 0.2051 | 0.1039 |
| common Lycaenidae Clarke | structured null | 0.1404 | 0.0984 | 0.0238 |

Papilionidae Iphiclides podalirius accounts for 82.19% of gross positive loss reductions vs structured null (8 old catalogue plant keys). Other contributors Zerynthia polyxena 9.59%, Z. rumina 4.53%. Species exclusions are posthoc outcome-ranked fixed-loss calculations, not independent confirmation or refitted models. Negative values and rare species with unstable relative percentages retained. Geographic influence: quality common Lycaenidae Clarke null gain +0.1404, remove top 5 H6 +0.0032, remove 10 -0.0497; old +0.2037 -> +0.0642 -> +0.0088. Papilionidae null +3.2884 -> +2.1017 -> +1.2430. Therefore distinct species and geographic leverage, no general strength ranking.

Taxonomy key audit: 6 new family-key combinations, 9 new-key butterfly records among all eligible H7, no known key seen in a different family. Original matched family key shares all exceed 99% at every nonzero family/site (minimum Lycaenidae 99.32%, Nymphalidae 99.86%); this is continuity of old keys not a live taxonomic synonym re-resolution. New keys remain in denominators. All source and outcome checks local; no new CI rerun. Source taxonomy full re-resolution, direct EQR/WWTP attribution, independent third-period change-on-change, and actual feeding not established. main unchanged.
