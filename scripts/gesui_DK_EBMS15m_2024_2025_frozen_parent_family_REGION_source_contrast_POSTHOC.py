#!/usr/bin/env python3
"""POSTHOC region/source-unit analysis of frozen licensed EBMS15m parent-group counts.

Parent-group source QC, original Event/Occurrence association and original GBIF
accepted/synonym mapping were frozen independently before this calculation.
This computes observational composition contrasts, NOT downstream causal effects.
"""
from __future__ import annotations

import argparse,hashlib,json,math
from collections import defaultdict,Counter
from pathlib import Path

ROOT="results/ebms_denmark_15minute_2023_2025_WWTP_SOURCE_ONLY/"
PARENT=ROOT+"POSTHOC_ORIGINAL_PARENTGROUP_ANON_FAMILY_COUNT_REGION_AND_WWTP_QC.json"
BEFORE=ROOT+"ORIGINAL_DK_EBMS15M_2024_25_ACCEPTED_FAMILY_PARENT_SURVEY_COUNT_SUPPORT_BEFORE_EFFECT.json"
FAMS=("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae")
YEARS=(2023,2024,2025)
SOURCE="near5cumulative";CONTROL="remote10"

def sha(b):return hashlib.sha256(b).hexdigest()
def logit(l,others):
    if l+others<=0:return None
    return math.log((l+.5)/(others+.5))
def pct(a,b):return a/b if b else None

def compute(root):
    richp=root/PARENT
    src=json.loads(richp.read_text())
    old=json.loads((root/BEFORE).read_text())
    if not old["primary_2024_and_2025_after_original_parent_source_family_QC_gate"]:
        raise ValueError("Source/taxon/count original pre-effect sample quality gate failed")
    if not src["this_parent_group_export_added_after_original_pooled_family_QC"]:
        raise ValueError("Anonymized original survey group source missing")
    assert src["original_full_DwCA_v13_SHA256"]==old["original_full_DwCA_v13_SHA256"]
    records=src["original_anonymized_parent_group_family_counts_for_posthoc_region_and_outfall_audit"]
    if len({x["original_parent_survey_unit_SHA256"] for x in records})!=len(records):
        raise ValueError("Parent group repeated, not independent")
    make=lambda: {"groups":0,"linked_groups":0,"positive_threefamily_groups":0,"events":0,"occurrences":0,"counts":Counter()}
    yearband=defaultdict(make);regionband=defaultdict(make);sourceunit=defaultdict(make)
    pargroup_vals=defaultdict(list)
    for x in records:
        year=x["source_only_year"];band=x["original_qualifying_EEA2022_WWTP_distance_band"]
        reg=x["original_geographic_1degree_cell"]
        if year not in YEARS or band not in (SOURCE,CONTROL,"middle"):
            raise ValueError("Unexpected source/year classification")
        counts=Counter(x["accepted_exact_GBIF_species_family_individualCount_including_synonyms"])
        if any(n<0 for n in counts.values()):raise ValueError("Negative original count")
        target=sum(counts[f] for f in ("Lycaenidae","Pieridae","Papilionidae"))
        for acc in (yearband[(year,band)],regionband[(year,reg,band)]):
            acc["groups"]+=1
            acc["linked_groups"]+=bool(x["has_any_published_Original_Occurrence_record"])
            acc["positive_threefamily_groups"]+=bool(target>0)
            acc["events"]+=x["original_child_Event_core_rows"]
            acc["occurrences"]+=x["published_Original_Occurrence_rows"]
            acc["counts"].update(counts)
        if band==SOURCE:
            key=x["nearest_registered_EEA_WWTP_outfall_SHA256"]
            if not key:raise ValueError("Missing original near-source WWTP outfall ID")
            s=sourceunit[(year,key)]
            s["groups"]+=1
            s["linked_groups"]+=bool(x["has_any_published_Original_Occurrence_record"])
            s["positive_threefamily_groups"]+=bool(target>0)
            s["events"]+=x["original_child_Event_core_rows"]
            s["occurrences"]+=x["published_Original_Occurrence_rows"]
            s["counts"].update(counts)
        if target>0 and band in (SOURCE,CONTROL):
            pargroup_vals[(year,band)].append(counts["Lycaenidae"]/target)
    def obj(x):
        d=x["counts"];L=d["Lycaenidae"];P=d["Pieridae"]+d["Papilionidae"]
        return {"parent_groups":x["groups"],"documented_any_occurrence_parent_groups":x["linked_groups"],
                "parent_groups_positive_target_family":x["positive_threefamily_groups"],
                "child_Event_rows":x["events"],"linked_Occurrence_rows":x["occurrences"],
                "accepted_or_synonym_family_individuals":{f:d[f] for f in FAMS},
                "Lycaenidae_pooled_threefamily_share":pct(L,L+P),
                "Lycaenidae_site_equal_mean_targetpositive_parent_share":None,
                "Lycaenidae_pooled_logodds_if_positive_target":logit(L,P)}
    out={
        "schema":"gesui_2026_DK_EBMS15minute_original2024_25_parent_group_region_and_WWTP_source_vs_remote_POSTHOC_v01",
        "original_parent_qc_SHA256":sha(richp.read_bytes()),
        "original_before_effect_gate_SHA256":sha((root/BEFORE).read_bytes()),
        "original_source_taxonomy_and_parent_survey_group_QC_matches":True,
        "study_2024_2025_only_primary_as_predeclared":True,
        "study_2023_exploratory_not_original_parent_source_support":True,
        "local_near_source_is_proximity_NOT_actual_downstream_exposure":True,
        "analysis_compositional_individualCount_from_15min_Events_NOT_population_or_effect":True,
        "archived_only_frozen_anonymized_parent_group_counts":True,
        "year_by_proximity_parent_group_and_family_counts":{},
        "original_country_1degree_region_year_source_vs_control_within_region":{},
        "source_outfall_independent_units_by_year":{},
        "equal_region_results_by_year":{},
    }
    for year in YEARS:
        out["year_by_proximity_parent_group_and_family_counts"][str(year)]={}
        for band in (SOURCE,CONTROL,"middle"):
            z=obj(yearband[(year,band)])
            v=pargroup_vals[(year,band)]
            z["Lycaenidae_site_equal_mean_targetpositive_parent_share"]=sum(v)/len(v) if v else None
            out["year_by_proximity_parent_group_and_family_counts"][str(year)][band]=z
            original=old["source_band_by_year_parent_group_observation_and_family_count_support"][str(year)][band]
            if z["documented_any_occurrence_parent_groups"]!=original["original_parent_groups_having_any_linked_Occurrence"]:
                raise ValueError("Parent group observation count discrepancy with original sealed outcome QA")
            if z["accepted_or_synonym_family_individuals"]["Lycaenidae"]!=original["original_accepted_synonym_family_individualCount"].get("Lycaenidae",0):
                raise ValueError("Original accepted family Lycaenidae 2024/25 count mismatch")
        units={}
        sourcevalid=0
        for (y,s),data in sourceunit.items():
            if y!=year:continue
            rec=obj(data)
            rec["original_EEA_outfall_id_SHA256"]=s
            rec["source_outfall_has_ge2_original_occurrence_linked_parent_groups"]=rec["documented_any_occurrence_parent_groups"]>=2
            if rec["source_outfall_has_ge2_original_occurrence_linked_parent_groups"]:sourcevalid+=1
            units[s]=rec
        out["source_outfall_independent_units_by_year"][str(year)]={
            "original_distinct_source_outfall_nearby_count":len(units),
            "outfalls_with_two_or_more_original_record_published_parent_groups":sourcevalid,
            "original_qualified_outfall_family_counts":dict(sorted(units.items()))
        }
        assert sourcevalid==old["source_band_by_year_parent_group_observation_and_family_count_support"][str(year)][SOURCE]["source_outfalls_having_at_least_two_parent_groups"]
        cmp=[]
        regs={r for y,r,b in regionband if y==year}
        for reg in sorted(regs):
            if (year,reg,SOURCE) not in regionband or (year,reg,CONTROL) not in regionband:continue
            a=obj(regionband[(year,reg,SOURCE)]);b=obj(regionband[(year,reg,CONTROL)])
            x=a["Lycaenidae_pooled_logodds_if_positive_target"]
            z=b["Lycaenidae_pooled_logodds_if_positive_target"]
            ratio=a["Lycaenidae_pooled_threefamily_share"]
            vs=b["Lycaenidae_pooled_threefamily_share"]
            cmp.append({
                "original_country_1degree_region":reg,
                "near_original_source":a,"remote_registered_source_10km_or_more":b,
                "both_year_region_focal_3family_counts_positive":ratio is not None and vs is not None,
                "near_minus_remote_pooled_family_share_if_supported":ratio-vs if ratio is not None and vs is not None else None,
                "near_minus_remote_pooled_logodds_if_supported":x-z if x is not None and z is not None else None,
                "original_WWTP_source_systems_nearby_are_not_randomized_treatments":True
            })
        out["original_country_1degree_region_year_source_vs_control_within_region"][str(year)]=cmp
        differences=[x["near_minus_remote_pooled_logodds_if_supported"] for x in cmp if x["both_year_region_focal_3family_counts_positive"]]
        shar=[x["near_minus_remote_pooled_family_share_if_supported"] for x in cmp if x["both_year_region_focal_3family_counts_positive"]]
        rsum={
            "regions_containing_both_source_bands":len(cmp),
            "regions_both_source_bands_have_positive_threefamily_count":len(differences),
            "posthoc_equal_geographic_region_mean_Lycaenidae_share_near_minus_remote":sum(shar)/len(shar) if shar else None,
            "posthoc_equal_region_mean_logodds_near_minus_remote":sum(differences)/len(differences) if differences else None,
            "region_effect_signs_negative":sum(x<0 for x in differences),
            "region_effect_signs_positive":sum(x>0 for x in differences),
            "leave_one_region_out_region_equal_logodds_range":[min((sum(differences)-v)/(len(differences)-1) for v in differences),max((sum(differences)-v)/(len(differences)-1) for v in differences)] if len(differences)>1 else None,
            "not_randomized_treatment_effect_or_statistical_significance":True,
        }
        out["equal_region_results_by_year"][str(year)]=rsum
    for yr in (2024,2025):
        if not out["equal_region_results_by_year"][str(yr)]["regions_containing_both_source_bands"]>=3:
            raise ValueError("Not enough same 1-degree regions for original source study")
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    o=compute(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "primary_2024_2025":{
            str(y):{"pooled_family":o["year_by_proximity_parent_group_and_family_counts"][str(y)],
                    "same_geographical_region":o["equal_region_results_by_year"][str(y)],
                    "independent_outfalls":o["source_outfall_independent_units_by_year"][str(y)]["outfalls_with_two_or_more_original_record_published_parent_groups"]}
            for y in (2024,2025)},
        "2023_descriptive_only":o["equal_region_results_by_year"]["2023"]
    },indent=2))
if __name__=="__main__":main()
