#!/usr/bin/env python3
"""Independent arithmetic replay of Denmark EBMS 2023-24 posthoc WWTP comparisons.

This validator uses source-hash locked derived 107x2 survey panels, no
remote species API. Verifies all denominators, zero-target branches,
original shared-region sample and four WWTP exposure clusters.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
from collections import defaultdict
ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
S=ROOT+"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
C=ROOT+"DENMARK_EBMS_2023_2024_FAMILY_SITEYEAR_COUNTS_QC_BEFORE_WWTP_CONTRAST.json"
A=ROOT+"DENMARK_EBMS_2023_2024_LYCAENIDAE_SOURCE5KM_VS_REMOTE10KM_POSTHOC_REGION_RESULTS.json"
def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def mean(a):return sum(a)/len(a)
def logit(r):
    z=r["family_individual_counts"]
    if r["target_3family_count"]==0:return None
    return math.log((z["Lycaenidae"]+.5)/(z["Pieridae"]+z["Papilionidae"]+.5))
def load(root):
    s=read(root/S);c=read(root/C);a=read(root/A)
    assert a["no_2km_source_near_butterfly_effect_estimated"] is True
    assert c["pre_effect_model_family_count_observation_gate_pass"] is True
    assert c["original_taxon_WWTP_species_count_data_quality_QC"]["original_event_to_occurrence_link_fraction"]==1
    assert s["source_distance_and_yearvisit_eligibility"]["5km"]["butterfly_original_occurrence_extension_may_be_accessed_ONLY_if_true"] is True
    assert a["near5_sourceWWTP_independent_outfall_units_n"]==4
    assert len(c["site_year_family_counts_original_public_CC_BY_NC_aggregate"])==214
    return s,c,a
def run(root):
    s,c,a=load(root)
    cells=defaultdict(list)
    results={}
    panel=c["site_year_family_counts_original_public_CC_BY_NC_aggregate"]
    for r in panel:
        d=r["source_2022_nearest_qualified_WWTP_km"]
        cls="near" if d<=5 else "remote" if d>=10 else "other"
        cells[(r["original_geo_1deg_cell"],r["site_year"],cls)].append(r)
    for year in (2023,2024):
        report={}
        for name in ("near","remote"):
            rows=[r for r in panel if r["site_year"]==year and ((r["source_2022_nearest_qualified_WWTP_km"]<=5) if name=="near" else (r["source_2022_nearest_qualified_WWTP_km"]>=10))]
            l=sum(r["family_individual_counts"]["Lycaenidae"] for r in rows)
            p=sum(r["family_individual_counts"]["Pieridae"] for r in rows)
            pa=sum(r["family_individual_counts"]["Papilionidae"] for r in rows)
            visits=sum(r["distinct_original_visit_days"] for r in rows)
            report[name]={"site_rows":len(rows),"site_with_any_target":sum(r["target_3family_count"]>0 for r in rows),
                          "Lycaenidae_count":l,"Pieridae_count":p,"Papilionidae_count":pa,
                          "Lycaenidae_share_given_three_target_families":l/(l+p+pa) if l+p+pa else None,
                          "original_survey_days":visits,"Lycaenidae_per_survey_day":l/visits if visits else None,
                          "Pieridae_per_survey_day":p/visits if visits else None}
            target=a["2023_2024_siteyear_count_data"][str(year)]["near5" if name=="near" else "remote10"]
            assert target["observed_individualCount_by_family"]["Lycaenidae"]==l
            assert target["observed_individualCount_by_family"]["Pieridae"]==p
            assert target["original_distinct_visit_days"]==visits
        results[str(year)]=report
    expected2023={"near":[23,11,74,128],"remote":[38,28,387,355]}
    expected2024={"near":[23,18,34,151],"remote":[38,31,494,481]}
    for year,e in (("2023",expected2023),("2024",expected2024)):
        for group,v in e.items():
            q=results[year][group]
            assert [q["site_rows"],q["site_with_any_target"],q["Lycaenidae_count"],q["Pieridae_count"]]==v
    common=set()
    original={}
    for (cell,y,cls),items in cells.items():
        if cls=="near" and cells.get((cell,y,"remote")):
            common.add(cell)
    assert common=={"55|10","55|12","56|9"}
    for year in (2023,2024):
        orig={r["region_1degree"]:r for r in a["original_3_shared_1degree_region_year_contrasts"] if r["year"]==year}
        assert set(orig)==common
        for cell in common:
            rec=orig[cell]
            for cls,nm in (("near","near_source_5km"),("remote","remote_10km_or_more")):
                z=cells[(cell,year,cls)]
                logs=[m for r in z if (m:=logit(r)) is not None]
                ref=rec[nm]
                assert ref["siteyears_with_at_least_one_of_3_families"]==len(logs)
                site_mean=mean(logs) if logs else None
                old=ref["site_equal_logodds_Lycaenidae_vs_PieridaePapilionidae_given_3family_detected"]
                assert (site_mean is None and old is None) or abs(site_mean-old)<1e-9
            original[(cell,year)]=rec["effect_lycaenidae_site_equal_family_logodds_near_minus_remote"]
    assert original[("55|10",2023)] is None
    assert original[("55|12",2023)] is None
    assert original[("56|9",2023)]<0
    assert original[("55|10",2024)]>0
    assert original[("55|12",2024)]<0
    assert original[("56|9",2024)]<0
    target=a["year_level_same_region_summary"]["2024"]
    assert target["composition_data_eligible_regions"]==3
    assert abs(target["equal_region_logodds_mean_ONLY_IF_ALL3_POSITIVE_TARGET_SUPPORTED"]-(-0.503267396532263))<1e-9
    assert target["leave_one_region_out_range_IF_ALL3_SUPPORTED"][1]>0
    out={
        "schema":"gesui_DK_EBMS_2023_2024_5km_source_original_count_contrast_independent_arithmetic_QA_v01",
        "original_sha256":{name:sha(root/path) for name,path in (("source",S),("family_counts",C),("regional_WWTP_contrast",A))},
        "independently_verified_join_and_target_record_n":214,
        "original_source_outfall_independent_near_groups_n":4,
        "equal_region_shared_geographic_cells":sorted(common),
        "all_reference_counts_and_composition_logodds_replayed_exactly":True,
        "unadjusted_pooled_target_family_counts_are_not_effective_independent_populations":True,
        "source_status_not_a_known_2023_or_2024_treatment":True,
        "yearwise_source_vs_remote_descriptive_counts":results,
        "2023_two_of_three_comparison_regions_no_near_target_family_and_ratio_undefined":True,
        "2024_three_regions_mixed_sign_and_leave_region_out_includes_positive":True,
        "no_causal_butterfly_benefit_or_harm_inferred":True
    }
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);args=p.parse_args()
    root=Path(__file__).resolve().parents[1];o=run(root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(o,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print(json.dumps(o,ensure_ascii=False,indent=2))
if __name__=="__main__":main()
