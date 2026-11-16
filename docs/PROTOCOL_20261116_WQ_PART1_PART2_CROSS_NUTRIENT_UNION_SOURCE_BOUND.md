# Corrected precomputed upper bound for Part1/Part2 station-year TN/TP union (2026-11-16)

Correction to the rejected naive bound: adding original Part1 valid BOTH H7-period n37 to original Part2 valid BOTH H7-period n5 is NOT an upper bound of a per-nutrient union. Same station-year might have valid total N only in Part1 and valid total P only in Part2 (or vice versa). This cross-table complementarity must be allowed without counting it as two independent measurements.

An exact but cheap source-only *upper bound*, not exact union, is:

If the final station-year per-nutrient union has TN+TP both, then at least one of:
(A) original Part1 already has both at that station-year, or
(B) original Part2 has >=1 qualified TN or TP at the same station-year.

Therefore every qualifying H7×period in a correct Part1-first/Part2-fallback union must belong to the union of:
S1 = the genuine Part1 BOTH-qualified H7×period set (size exactly 37 primary, 24 genuine H6), and
S2 = the genuine Part2 ANY-valid-TN-or-TP H7×period set (unknown until original source scan).

Hence upper bound on distinct H7×period <= 37 + |S2|, upper bound genuine H6 <= 24 + distinct H6 of S2. This upper bound can be loose due to geographic overlap and station-year mismatch, but it must NEVER undercount a valid nutrient-wise union. If either bound cannot meet 100 H7×period /40 H6 frozen science gate, it is mathematically unnecessary to redownload and rescan 1.64GB/96.6M Part1 again solely to compute union, and biology effects remain STOP.

Independent Part2 source-only step:
- Reuse exact original WQ station/year/H7 polygon/source country primary frozen candidate SHA256 c5cde7bd8134b254c5e95b06849292fa260fd6f01832e6b6d2e1994a903c5910.
- Read official original 2026 Part2 170,422,309 byte ZIP source SHA256 315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a; 6,450,192 genuine rows, years 2018,2022, same station identifier+scheme+country and year, RW/W, actual total TN EEA_31615-01-7 or total TP CAS_7723-14-0, canonical elemental mg/L, original >=3 samples, valid finite LOQ-gated mean, exactly one valid annual report for the particular nutrient.
- Record all year-matched primary country WQ source site-year keys for which **either** TN or TP is qualified, regardless of whether both qualify. No original EQR, PlantN, HNS, effect sizes or signs read.
- Before reading any biological coefficients, compute strict source-only S2 counts distinct H7period and genuine H6, and apply loose *sufficient* impossibility bound to original science minimum 100 H7period/40 H6. If bound below threshold, stop. If bound could reach, one must perform the full original 96.6M Part1 per-nutrient source extraction and exact source-priority station-year union; cannot claim a source gate passed from an upper bound.

This correction leaves old original Part2-only STOP and Part1-alone source PASS/science STOP intact. Never treat absent source values as zero or environmental absence. Source documentation and next actions remain a Draft PR only.
