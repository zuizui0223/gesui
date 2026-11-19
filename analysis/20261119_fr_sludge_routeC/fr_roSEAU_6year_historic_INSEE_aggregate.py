#!/usr/bin/env python3
"""Six original ROSEAU source years: union of COG-verified candidate communes.

The raw per-station agricultural sludge mass is not allocated to the listed
communes. Only historical INSEE-verified code/name candidates are combined.
"""
import argparse,csv,hashlib,json
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict

YEARS=tuple(str(y) for y in range(2017,2023))
SCHEMA="gesui_France_sludge_2017to22_embedded_INSEE_source_year_specific_COG_v1"

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def summarize(base):
    per={}
    candidates={}
    receipts={}
    for year in YEARS:
        q=list(Path(base).rglob(f"FR_ROSEAU_{year}_HISTORIC_INSEE_COG_SOURCE_QA.json"))
        if len(q)!=1:raise ValueError("Original COG source receipt missing/duplicated for "+year)
        v=json.loads(q[0].read_text())
        if v.get("year")!=year or v.get("schema")!=SCHEMA or not v.get("no_plant_NUTRI_butterfly_HNS_EQR_used"):
            raise ValueError("Original year or source-outcome contract mismatched")
        c=list(Path(base).rglob(f"FR_ROSEAU_{year}_COG_VALIDATED_CANDIDATE_COMMUNE_CODES.csv"))
        if v.get("status")!="OFFICIAL_YEAR_COG_AND_ORIGINAL_SLUDGE_CODE_NAME_VALIDATED":
            per[year]={"status":v.get("status"),"source_error":v.get("error"),"ready":False}
            candidates[year]=set()
            continue
        if len(c)!=1:raise ValueError("Public code list missing for a successfully verified year")
        with c[0].open(newline="") as f:
            reader=csv.DictReader(f)
            rows=list(reader)
        if any(r.get("year")!=year or r.get("interpretation")!=
            "candidate_in_station_spreading_commune_list_and_station_positive_tMS_NOT_field_application" for r in rows):
            raise ValueError("Nonconforming independent source-only municipality code list")
        codes={r["INSEE_code"] for r in rows}
        if len(codes)!=len(rows):raise ValueError("Duplicated candidate INSEE within one year")
        n=v.get("source_counts",{})
        if len(codes)!=n.get("distinct_positive_codes_matching_both_COG_code_and_normalized_name"):
            raise ValueError("Original source year counts not equal to candidate INSEE list")
        per[year]={"status":v.get("status"),"ready":True,
                   "orig_lab_assay_rows":n.get("original_repeated_lab_rows"),
                   "distinct_STEU":n.get("distinct_STEU_source_stations"),
                   "positive_STEU":n.get("positive_agricultural_spreading_mass_stations"),
                   "original_station_commune_code_pairs":n.get("all_original_station_commune_code_pairs"),
                   "positive_station_commune_code_pairs":n.get("positive_station_commune_code_pairs"),
                   "distinct_positive_verified_code_name_pairs":len(codes),
                   "distinct_positive_departement_prefixes":len({x[:2] for x in codes}),
                   "COG_code_name_conflicts":n.get("COG_code_name_conflict_pairs_all_source"),
                   "source_hash":v.get("original_ROSEAU_CSV_SHA256")}
        receipts[year]={"original_source_receipt_sha":sha(q[0]),"COG_verified_list_sha":sha(c[0])}
        candidates[year]=codes
    early=set().union(*(candidates[y] for y in YEARS[:3]))
    late=set().union(*(candidates[y] for y in YEARS[3:]))
    union=early|late
    posyears=defaultdict(set)
    for y in YEARS:
        for code in candidates[y]:posyears[code].add(y)
    return per,receipts,posyears,{
        "full_6year_COG_verified":len(receipts)==6,
        "2017to19_positive_candidate_codes_across_years":len(early),
        "2020to22_positive_candidate_codes_across_years":len(late),
        "candidate_codes_in_both_periods":len(early&late),
        "2017to22_unique_candidate_codes_across_years":len(union),
        "distinct_source_application_departments_code_prefix":len({x[:2] for x in union}),
        "numeric_only_min_50_candidate_communes":len(union)>=50,
        "numeric_only_min_5_department_prefixes":len({x[:2] for x in union})>=5,
        "positive_in_both_predeclared_periods":bool(early) and bool(late),
        "confirmed_commune_year_APPLICATION_EVENTS":False,
        "annual_STEU_mass_not_allocated_to_commune_code":True,
        "H7_source_geo_crosswalk_NOT_completed":True,
        "causal_sludge_to_plant_butterfly_source_GO":False,
        "prior_NUTS2_point_diffuse_source_stop_remains":True
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    per,receipts,candidates,gate=summarize(a.input)
    a.out.mkdir(parents=True,exist_ok=True)
    summary={"schema":"gesui_FR_ROSEAU_2017to22_COG_validated_candidate_commune_SOURCE_only_aggregate_v1",
             "UTC":datetime.now(timezone.utc).isoformat(),
             "years":per,"source_sha_receipts":receipts,"gate":gate,
             "no_biological_values_read":True}
    (a.out/"FR_ROSEAU_2017TO22_INSEE_SOURCE_CANDIDATE_COMMUNES_GATE.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    with (a.out/"FR_ROSEAU_2017TO22_PUBLIC_CANDIDATE_COMMUNE_YEAR_PRESENCE.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["INSEE_code","candidate_source_years","source_interpretation"])
        for code,years in sorted(candidates.items()):
            w.writerow([code,";".join(sorted(years)),
                        "listed_in_STEU_spreading_communes_with_positive_STEU_agricultural_tMS_NOT_field_exposure"])
    print("GESUI_FR_SLUDGE_SIX_YEAR_INSEE_CODES",json.dumps(gate,sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__": main()
