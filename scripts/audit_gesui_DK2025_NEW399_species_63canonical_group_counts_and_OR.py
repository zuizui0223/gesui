#!/usr/bin/env python3
"""Pure-Python independent arithmetic audit of original 2025 species/family OR.

Reads only immutable source-locked original 2025 species×source aggregate and
previous official 399×family source totals. Does NOT import the species parser,
query any taxa, re-download licensed raw occurrence records or fit effects.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from pathlib import Path

D=Path("results/ebms_denmark_2023_2024_WWTP_source_only")
SPEC=D/"DK2025_NEW399_ORIGINAL_65SPECIES_3VS5_FAMILY_DECOMPOSITION.json"
FAMILY=D/"DENMARK_2025_NEW399_SITE_LYCAENIDAE_FAMILY_POST_SOURCE_GATE_COUNTS.json"
GROUPS=("near5","remote10","intermediate")
FAMS=("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae")

def sha(b):return hashlib.sha256(b).hexdigest()

def oratio(a,b,fams):
    def odds(x):
        num=x["Lycaenidae"]
        den=sum(x[f] for f in fams if f!="Lycaenidae")
        return num/den if num and den else None
    x=odds(a);y=odds(b)
    return x/y if x and y else None

def check_close(a,b,label):
    if a is None or b is None:
        if a is not b:raise ValueError(label+" original undefined mismatch")
        return
    if not math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12):
        raise ValueError(f"{label} original mismatch {a} vs {b}")

def run(root):
    f1=root/SPEC;f2=root/FAMILY
    src=json.loads(f1.read_text())
    old=json.loads(f2.read_text())
    if src["source_JSON_SHA256"]["family_result_original_399_site"]!=sha(f2.read_bytes()):
        raise ValueError("Original species results use different family panel")
    if not src["reconciled_exact_399_source_2025_site_family_individualCounts"]:
        raise ValueError("Original 399 sites not reconciled")
    rows=src["all_accepted_five_family_species_2025_individual_count_by_source_band"]
    if len(rows)!=63 or len({x["gbif_accepted_usageKey"] for x in rows})!=63:
        raise ValueError("Expected 63 accepted canonical species usageIDs")
    if sum(x["original_scientificName_aliases_n"] for x in rows)!=65:
        raise ValueError("Original frozen 65 accepted name records not reconciled to 63 usage keys")
    if Counter(x["accepted_family"] for x in rows)!=Counter({
        "Lycaenidae":20,"Pieridae":7,"Papilionidae":1,"Nymphalidae":27,"Hesperiidae":8
    }):
        # If original GBIF taxonomy had different accepted synonym duplication,
        # audit only correct stable accepted source and do not silently drop.
        raise ValueError("Original canonical five-family source diversity changed")
    counts={g:Counter() for g in GROUPS}
    site_counts={g:0 for g in GROUPS}
    for r in rows:
        for g in GROUPS:
            v=r["per_group"][g]
            n=v["original_2025_eBMS_individualCount"]
            if type(n) is not int or n<0:raise ValueError("Original non-integer individualCount")
            counts[g][r["accepted_family"]]+=n
    for g in GROUPS:
        previously=old["2025_source_frozen_286_documented_site_counts"][g]
        if any(counts[g][f]!=previously["families_individualCount"][f] for f in FAMS):
            raise ValueError("Species aggregation fails independently published family totals")
        for f in FAMS:
            if counts[g][f]!=src["family_individualCount_by_source_group"][g][f]:
                raise ValueError("Species and original family source group disagree")
    original_or3=oratio(counts["near5"],counts["remote10"],FAMS[:3])
    original_or5=oratio(counts["near5"],counts["remote10"],FAMS)
    check_close(original_or3,src["original_near_vs_remote_three_family_pooled_observational_OR"],"Original 3family OR")
    check_close(original_or5,src["original_near_vs_remote_five_family_pooled_observational_OR"],"Original 5family OR")
    original_rate=counts["near5"]["Lycaenidae"]/539/(counts["remote10"]["Lycaenidae"]/757)
    check_close(original_rate,src["original_near_vs_remote_absolute_Lycaenidae_recorded_per_linked_day_ratio"],"Original absolute L rate")
    removals={}
    for row in rows:
        usage=row["gbif_accepted_usageKey"]
        family=row["accepted_family"]
        near=Counter(counts["near5"])
        remote=Counter(counts["remote10"])
        near[family]-=row["per_group"]["near5"]["original_2025_eBMS_individualCount"]
        remote[family]-=row["per_group"]["remote10"]["original_2025_eBMS_individualCount"]
        if any(v<0 for v in near.values()) or any(v<0 for v in remote.values()):
            raise ValueError("Removing recorded species exceeds original family population")
        removals[usage]={
            "or3":oratio(near,remote,FAMS[:3]),
            "or5":oratio(near,remote,FAMS),
            "species_name":row["accepted_GBIF_scientificName"],
            "original_near_n":row["per_group"]["near5"]["original_2025_eBMS_individualCount"],
            "original_remote_n":row["per_group"]["remote10"]["original_2025_eBMS_individualCount"]
        }
    orig_rem={x["accepted_GBIF_usageKey"]:x for x in src["species_sensitivity_original_pooled_three_vs_five_family_OR_after_one_species_removed"]}
    if set(removals)!=set(orig_rem):raise ValueError("Independent canonical accepted species list differs")
    for taxon,r in removals.items():
        check_close(r["or3"],orig_rem[taxon]["three_family_OR_after_deleting_species_in_both_groups"],taxon+" OR3")
        check_close(r["or5"],orig_rem[taxon]["five_family_OR_after_deleting_species_in_both_groups"],taxon+" OR5")
    selected={}
    for name,usage in (("Pieris brassicae","1920506"),
                       ("Maniola jurtina","5135513"),
                       ("Polyommatus icarus","5140214"),
                       ("Plebejus argus","5138561")):
        r=removals[usage]
        if name.lower() not in r["species_name"].lower():raise ValueError("GBIF usage key moved species")
        selected[name]=r
    source=src["largest_original_species_contributors_per_family_and_source_band"]
    for f in ("Pieridae","Nymphalidae","Lycaenidae"):
        for g in ("near5","remote10"):
            first=source[f][g][0]
            found=next(x for x in rows if x["gbif_accepted_usageKey"]==first["accepted_usageKey"])
            if first["original_individualCount"]!=found["per_group"][g]["original_2025_eBMS_individualCount"]:
                raise ValueError("Original contributor rank and species ID differ")
    return {
        "schema":"gesui_DK2025_new399_original_canonical_species_3_vs5_denominator_offline_independent_QA_v01",
        "source_species_results_SHA256":sha(f1.read_bytes()),
        "previous_source_family_results_SHA256":sha(f2.read_bytes()),
        "original_399_site_family_reconciliation_confirmed":True,
        "original_65_source_accepted_taxon_names_to_63_unique_GBIF_usage_ids":True,
        "canonical_species_accepted_family_diversity":dict(Counter(r["accepted_family"] for r in rows)),
        "source_class_five_family_totals_rebuilt":{g:dict(counts[g]) for g in GROUPS},
        "source_groups_absolute_Lycaenidae_recording_rate_ratio_rebuilt":original_rate,
        "source_groups_pooled_three_family_OR_rebuilt":original_or3,
        "source_groups_pooled_five_family_OR_rebuilt":original_or5,
        "all_63_species_leave_one_out_pooled_OR_both_denominators_rebuilt_identically":True,
        "selected_species_removal_sensitivity_observational_accounting_only":selected,
        "source_2025_event_dates_are_NOT_measurements_of_wastewater_dose":True,
        "original_2022_reported_EEA_WWTP_source_proximity_not_physical_larval_or_adult_exposure":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    x=run(root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(x,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"original_species":x["canonical_species_accepted_family_diversity"],
                       "or3":x["source_groups_pooled_three_family_OR_rebuilt"],
                       "or5":x["source_groups_pooled_five_family_OR_rebuilt"],
                       "removals":x["selected_species_removal_sensitivity_observational_accounting_only"]},indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
