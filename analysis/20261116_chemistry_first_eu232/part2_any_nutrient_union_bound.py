#!/usr/bin/env python3
"""Rigorous source-only upper bound for Part1+Part2 nutrient-wise TNTP union.

A valid strict Part1-first, Part2-fallback SAME station-year BOTH TNTP pair
must either come from Part1 BOTH (old n37 H7period,24 H6), or involve >=1
strict Part2 nutrient at that exact station/year. Therefore a conservative
geographic upper bound is |Part1 BOTH H7period| + |Part2 ANY H7period|; H6
bound similarly. Do NOT replace "Part2 ANY" with "Part2 BOTH" (invalid).
No EQR / NUTRI / HNS values or wastewater causal effects opened.
"""
from __future__ import annotations
import argparse,csv,collections,json,sys,zipfile,io
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_chemistry_first_eu232"))
from source_first_part2_TNTP import load_frozen_geography
sys.path.insert(0,str(ROOT/"analysis/20261116_reported_np_source"))
from acquire_eea_2026_part2_exact_site import download_part2,parse_strict,CODES

SOURCE_ZIP_SHA="315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a"
SOURCE_CSV_N=6450192
P1_SCHEMA="GESUI_20261116_CHEMISTRY_FIRST_PART1_WQ_EU232_STRICT_SITE_YEAR_TNTP_V1"
P2_ANY_SCHEMA="GESUI_20261116_CHEMISTRY_FIRST_PART2_ANY_NUTRIENT_UNION_BOUND_V1"

def load_p1(p):
    q=json.loads((p/"CHEMISTRY_FIRST_EU232_PART1_SOURCE_GATE.json").read_text())
    if q.get("schema")!=P1_SCHEMA or not q.get("source_support_gate_pass"):
        raise RuntimeError("Original validated Part1 source gate not PASS")
    s=q.get("primary_original_GBIF_basin_country_matched_qualified")
    if not s or s["H7period"]!=37 or s["true_H6"]!=24 or s["strict_valid_original_WQ_siteyears"]!=621:
        raise RuntimeError("Original source Part1 37/24/621 documented source cohort changed")
    return q,s

def scan(zipfile_path,original):
    qualified=collections.defaultdict(lambda:{"TN":[],"TP":[]})
    stats=collections.Counter()
    with zipfile.ZipFile(zipfile_path) as z:
        paths=[x.filename for x in z.infolist() if x.filename.lower().endswith(".csv")]
        if len(paths)!=1:raise RuntimeError("Original Part2 source CSV not unique")
        with z.open(paths[0]) as stream:
            f=io.TextIOWrapper(stream,encoding="utf-8-sig",newline="",errors="strict")
            header=f.readline()
            delimiter=";" if header.count(";")>header.count(",") else ","
            fields=next(csv.reader([header],delimiter=delimiter))
            required=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
              "phenomenonTimeReferenceYear","parameterWaterBodyCategory",
              "procedureAnalysedMatrix","observedPropertyDeterminandCode",
              "resultUom","resultNumberOfSamples","resultMeanValue","resultQualityMeanBelowLOQ"]
            if not set(required).issubset(fields):
                raise RuntimeError("Official Part2 2026 source fields changed")
            rows=csv.DictReader(f,fieldnames=fields,delimiter=delimiter)
            n=0
            for q in rows:
                n+=1
                v=q["phenomenonTimeReferenceYear"]
                if v not in ("2018","2022"):continue
                key=(q["countryCode"],q["monitoringSiteIdentifier"],
                     q["monitoringSiteIdentifierScheme"],int(v))
                if key not in original:continue
                stats["source_any_parameter_at_selected_WQ_siteyear"]+=1
                if q["parameterWaterBodyCategory"]!="RW" or q["procedureAnalysedMatrix"]!="W":
                    continue
                nutrient=None
                if q["observedPropertyDeterminandCode"]==CODES["TN"]:nutrient="TN"
                elif q["observedPropertyDeterminandCode"]==CODES["TP"]:nutrient="TP"
                if nutrient is None:continue
                valid,reason=parse_strict(q,nutrient)
                if reason:
                    stats[nutrient+"_source_unqualified_"+reason]+=1
                else:
                    qualified[key][nutrient].append(valid)
                    stats[nutrient+"_valid_source_annual_rows"]+=1
            if n!=SOURCE_CSV_N:raise RuntimeError("Official original Part2 rows unexpected "+str(n))
    return qualified,stats,n

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    part1,p1counts=load_p1(a.part1)
    geom,candidates=load_frozen_geography(a.geography)
    if geom.get("source_candidate_CSV_sha256")!=part1.get("frozen_H7_WQ_station_year_CSV_SHA256"):
        raise RuntimeError("Both original Parts need identical pre-chemistry station selection")
    tmp=a.output/"OFFICIAL_2026_PART2_TEMP_SINGLE_NUTRIENT.zip"
    origin=download_part2(tmp)
    if origin["official_2026_original_part2_ZIP_SHA256"]!=SOURCE_ZIP_SHA:
        raise RuntimeError("Official Part2 2026 original source ZIP SHA changed")
    qualified,stats,total=scan(tmp,candidates)
    tmp.unlink()
    allrows=[]
    primary=[]
    fullN=0;fullP=0
    for key,site in sorted(candidates.items()):
        nutrients=qualified.get(key,{"TN":[],"TP":[]})
        TN=(len(nutrients["TN"])==1)
        TP=(len(nutrients["TP"])==1)
        if TN:fullN+=1
        if TP:fullP+=1
        if not (TN or TP):continue
        r={"h":site["official_original_HYBAS_ID_L07"],
           "h6":site["original_true_HYBAS_ID_L06"],
           "period":site["original_GBIF_source_H7_period"],
           "year":int(site["source_WQ_calendar_year"]),
           "country":site["source_WQ_country"],
           "original_majority_country_matches":site["WQ_station_country_equals_original_GBIF_basin_majority_country"].lower()=="true",
           "part2_TN_strict_valid":TN,
           "part2_TP_strict_valid":TP,
           "part2_any_strict_nutrient":True}
        allrows.append(r)
        if r["original_majority_country_matches"]:primary.append(r)
    if (fullN,fullP)!=(56,114):
        raise RuntimeError("Original verified Part2 56 TN /114 TP source site-year qualities not reproduced: "+
            str((fullN,fullP)))
    def ctr(rows):
        return {"station_years":len(rows),
                "H7period":len({(v["h"],v["period"]) for v in rows}),
                "distinct_H7":len({v["h"] for v in rows}),
                "true_H6":len({v["h6"] for v in rows}),
                "country_groups":len({v["country"] for v in rows}),
                "with_TN":sum(v["part2_TN_strict_valid"] for v in rows),
                "with_TP":sum(v["part2_TP_strict_valid"] for v in rows)}
    full=ctr(allrows)
    p2=ctr(primary)
    loose={"original_Part1_BOTH_H7period":p1counts["H7period"],
           "original_Part1_BOTH_H6":p1counts["true_H6"],
           "Part2_any_original_H7period":p2["H7period"],
           "Part2_any_original_H6":p2["true_H6"],
           "unconditional_valid_source_union_H7period_UPPER_BOUND":
               min(454,p1counts["H7period"]+p2["H7period"]),
           "unconditional_valid_source_union_H6_UPPER_BOUND":
               min(90,p1counts["true_H6"]+p2["true_H6"]),
           "exact_union_NOT_computed":True}
    notenough=(loose["unconditional_valid_source_union_H7period_UPPER_BOUND"]<100 or
               loose["unconditional_valid_source_union_H6_UPPER_BOUND"]<40)
    if primary:
        p=a.output/"CHEMISTRY_FIRST_PART2_ANY_NUTRIENT_H7_SOURCE_KEYS.csv"
        with p.open("w",encoding="utf-8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(primary[0]));w.writeheader();w.writerows(primary)
    out={"schema":P2_ANY_SCHEMA,
         "status":"STOP_ECOLOGY_EVEN_AFTER_POSSIBLE_PART1_PART2_COMPLEMENT" if notenough else
                  "CANNOT_STOP_FROM_LOOSE_UPPER_BOUND_EXACT_PART1_PART2_UNION_REQUIRED",
         "official_part2_original_SHA256":SOURCE_ZIP_SHA,
         "source_genuine_original_Part2_CSV_rows_scanned":total,
         "same_original_2026_spatial_station_year_source_sha":geom["source_candidate_CSV_sha256"],
         "Part1_independent_original_source_gate":p1counts,
         "Part2_ANY_TN_or_TP_full_source_WQ":full,
         "Part2_ANY_TN_or_TP_primary_WQ_countrymatched":p2,
         "independent_original_source_Part2_qualified_sample_rows":dict(stats),
         "valid_Part1_per_nutrient_then_Part2_fallback_union_geographic_bound":loose,
         "original_fixed_ecology_minimum":{"H7period":100,"true_H6":40,"countries":5},
         "strict_science_support_provably_unreachable_from_source_union_upper_bound":bool(notenough),
         "source_Part2_BOTH_qualified_H7period_n5_NOT_A_VALID_BOUND_FOR_CROSS_PART_UNION":True,
         "no_plant_or_butterfly_effect_fitted":True,
         "no_sewage_chemical_source_apportionment":True}
    (a.output/"CHEMISTRY_FIRST_EU232_PART2_ANY_NUTRIENT_UNION_UPPER_BOUND.json").write_text(
        json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_PART1_PART2_CORRECTED_UNION_UPPER_BOUND="+json.dumps({
       "status":out["status"],"Part2_any_primary":p2,
       "Part1_BOTH_primary":p1counts,"safe_upper_bound":loose,
       "source_finite_Part2_TN_and_TP":(fullN,fullP),
       "source_gate_disposition_without_ecology":bool(notenough)},sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--geography",type=Path,required=True)
    p.add_argument("--part1",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
