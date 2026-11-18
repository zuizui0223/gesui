#!/usr/bin/env python3
"""Source-only chronological gate across six published Grand River WWTP data matrices.

Original exact 8-facility source: Wardlaw et al. 2026 two Borealis DOIs.
Year key for each organism site is independently present in isotope samples.
Effluent collection-year mapping comes from *the published paper's methods*
(2023 at seven WWTPs, 2022 only Caledonia), not inferred from time=1,2,3.
No response magnitude, contaminant concentration, mediation model, or p-value
is read or computed. Stop original full-chain >=6 facilities/2 visits gate.
"""
import argparse, collections, csv, hashlib, json
from pathlib import Path

SCHEMA="GESUI_20261118_GRAND_RIVER_8_WWTP_SIX_MATRICES_CHRONOLOGICAL_GATE_V1"
ORIGINAL_JSON_SHA="4b25863efa087fcb05090f14c66cb93b9b194f041a2796db2880e66f9a8e72c3"
SITES=("Ayr","Caledonia","Guelph","Hespeler","Kitchener","Paris","Preston","Waterloo")
YEAR_BIOTA={"Ayr":2023,"Caledonia":2022,"Guelph":2023,"Hespeler":2022,
            "Kitchener":2022,"Paris":2023,"Preston":2023,"Waterloo":2022}
YEAR_EFFLUENT={site:(2022 if site=="Caledonia" else 2023) for site in SITES}
BRANCH={"Ayr":"Nith","Guelph":"Speed","Hespeler":"Speed",
        "Caledonia":"Grand","Kitchener":"Grand","Paris":"Grand",
        "Preston":"Grand","Waterloo":"Grand"}
MIN_MATCHED_MONTHS=2
MIN_FACILITIES=6
TAXA=("Larval Hydropsychidae","Adult Hydropsychidae","Tetragnathidae")
def sha(raw):return hashlib.sha256(raw).hexdigest()

def run(args):
    raw=args.original.read_bytes()
    if sha(raw)!=ORIGINAL_JSON_SHA:raise RuntimeError("Original signed six-table metadata SHA changed")
    src=json.loads(raw)
    if src["schema"]!="GESUI_20261117_BOREALIS_GRAND8_TWO_DOI_SOURCE_LOCKED_WWTP_SITE_YEAR_STAGE_GATE_V1":
        raise RuntimeError("Original source schema mismatch")
    files=src["source_files"]
    if len(files)!=6 or files["effluent"]["original_real_record_n"]!=77 or files["cross_species"]["original_real_record_n"]!=906 or files["isotopes"]["original_real_record_n"]!=918:
        raise RuntimeError("Original publication file source sample sizes differ")
    if set(src["source_exact_8WWTP_site_names_common_to_six_original_files"])!=set(SITES):
        raise RuntimeError("Exact all-eight six-matrix WWTP site-source IDs do not match")
    isoyears=collections.defaultdict(set)
    for entry in files["isotopes"]["source_stage_x_year_x_WWTP_site_x_updown_counts"]:
        k=entry["keys"]
        isoyears[k["Site"]].add(int(k["Year"]))
    if {k:sorted(isoyears[k]) for k in SITES}!={k:[v] for k,v in YEAR_BIOTA.items()}:
        raise RuntimeError("The original published isotope site/year cohort differs")
    eff=collections.Counter()
    for entry in files["effluent"]["source_stage_x_year_x_WWTP_site_x_updown_counts"]:
        k=entry["keys"]
        eff[(k["Site"],str(k["time"]))]+=int(entry["original_sample_records"])
    stage=collections.Counter()
    for entry in files["cross_species"]["source_stage_x_year_x_WWTP_site_x_updown_counts"]:
        k=entry["keys"]
        stage[(k["Site"],str(k["Time"]),k["Location"],k["Taxa"])]+=int(entry["original_sample_records"])
    original_times=("1","2","3")
    rows=[]
    for site in SITES:
        match=[]
        for t in original_times:
            every_stage=all(stage[(site,t,side,taxon)]>=2
                            for side in ("upstream","downstream") for taxon in TAXA)
            if YEAR_BIOTA[site]==YEAR_EFFLUENT[site] and eff[(site,t)]>=1 and every_stage:
                match.append(t)
        efftimes=[t for t in original_times if eff[(site,t)]>=1]
        rows.append({
          "original_facility":site,"original_named_branch":BRANCH[site],
          "original_biological_field_year":YEAR_BIOTA[site],
          "effluent_field_year_from_published_2026_methods":YEAR_EFFLUENT[site],
          "original_effluent_field_month_timepoint_count":len(efftimes),
          "original_effluent_available_timepoints":",".join(efftimes),
          "same_year_three_stage_original_specimens_and_effluent_months_n":len(match),
          "same_year_three_stage_original_timepoint_ids":",".join(match),
          "passes_two_matched_source_months":len(match)>=MIN_MATCHED_MONTHS
        })
    anysite=sum(x["same_year_three_stage_original_specimens_and_effluent_months_n"]>=1 for x in rows)
    full=sum(x["passes_two_matched_source_months"] for x in rows)
    if anysite!=5 or full!=4:raise RuntimeError(f"Original year/source support changed: {anysite=}, {full=}")
    args.out.mkdir(parents=True,exist_ok=True)
    path=args.out/"GRAND8_ORIGINAL_SIX_MATRIX_EFFLUENT_STAGE_YEAR_SOURCE_ONLY_QA.json"
    result={
      "schema":SCHEMA,
      "status":"STOP_ORIGINAL_6_WWTP_MIN_YEAR_MATCHED_TWOVISIT_FULL_CHAIN_4_OUT_OF_8",
      "source_only_metadata_sha256":ORIGINAL_JSON_SHA,
      "original_effluent_source_doi":"10.5683/SP3/TLDY6W",
      "original_insect_adult_spider_source_doi":"10.5683/SP4/U0MN9S",
      "published_field_year_methods_doi":"10.1093/etojnl/vgag111",
      "one_connected_Grand_River_watershed_three_named_branches":True,
      "total_original_distinct_WWTP":8,
      "n_facilities_with_at_least_one_valid_sameyear_exposure_3stage":anysite,
      "n_facilities_with_at_least_two_valid_sameyear_exposure_3stage":full,
      "predefined_min_full_chain_original_facilities":MIN_FACILITIES,
      "original_full_chain_source_gate_pass":full>=MIN_FACILITIES,
      "effluent_replicate_500mL_subsamples_are_not_independent_timepoints":True,
      "effluent_time_code_1_2_3_not_a_calendar_year":True,
      "source_no_numeric_microplastic_concentrations_loaded":True,
      "source_not_a_failure_of_actual_material_transfer":True,
      "no_effluent_source_to_spider_causal_mediation_fit":True,
      "original_facility_year_stage_support":rows}
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    with (args.out/"GRAND8_FACILITY_2022_2023_SAME_YEAR_STAGE_QA.csv").open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    print("GESUI_GRAND8_ORIGINAL_EFFLUENT_STAGE_TIME_GATE="+json.dumps({
      "status":result["status"],"anyyear_n":anysite,"two_visit_n":full,
      "min_required":MIN_FACILITIES},ensure_ascii=False),flush=True)

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--original",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    run(a.parse_args())
