#!/usr/bin/env python3
"""Independent offline arithmetic audit for source-locked Denmark same-month family counts."""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path

ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
SOURCE=ROOT+"DENMARK_2023_2024_SAME_SITE_MONTH_LYCAENIDAE_POSTHOC_FAMILY_SOURCE_CLUSTER_QC.json"
EARLIER=ROOT+"DENMARK_EBMS_2023_2024_LYCAENIDAE_SOURCE5KM_VS_REMOTE10KM_POSTHOC_REGION_RESULTS.json"
FAMS=("Lycaenidae","Pieridae","Nymphalidae","Hesperiidae")

def calc(data,fam):
    a=data["2023"];b=data["2024"]
    n=a["family_counts"][fam];m=b["family_counts"][fam]
    d=a["original_site_month_survey_days"];e=b["original_site_month_survey_days"]
    return (m/e)/(n/d) if n and d and e else None

def audit(root):
    p=root/SOURCE
    original=json.loads(p.read_text())
    old=json.loads((root/EARLIER).read_text())
    assert original["all_214_siteyear_original_family_sums_and_visit_days_exactly_replayed"] is True
    assert original["source_rows_qc"]["all_original_Event_rows"]==15606
    assert original["source_rows_qc"]["all_original_Occurrence_rows"]==21617
    assert original["source_rows_qc"]["original_accepted_family_counts_rows"]==4226
    assert original["licensed_original_Denmark_eBMS_v14_SHA256"]=="2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b"
    assert old["near5_sourceWWTP_independent_outfall_units_n"]==4
    results={}
    for key,study in original["season_results"].items():
        assert study["within_site_calendar_month_matched_sites_per_source_class"]=={
            "intermediate":46,"near5":23,"remote10":38
        }
        assert study["independent_registered_source_clusters_in_matched_sites"]==4
        near=study["matched_group_year"]["near5"]
        remote=study["matched_group_year"]["remote10"]
        for year in ("2023","2024"):
            source_summary=[x["match"][year] for x in study["four_source_outfall_near_group_2023_2024"]]
            assert sum(x["original_matched_source_site_n"] for x in source_summary)==23
            assert sum(x["original_site_month_survey_days"] for x in source_summary)==near[year]["original_site_month_survey_days"]
            for fam in (*FAMS,"Papilionidae"):
                assert sum(x["family_counts"][fam] for x in source_summary)==near[year]["family_counts"][fam]
            assert near[year]["original_matched_source_site_n"]==23
            assert remote[year]["original_matched_source_site_n"]==38
        metrics={}
        for f in FAMS:
            a=calc(near,f);b=calc(remote,f)
            expected=study["descriptive_near_vs_remote_ratio_of_2024_to_2023_visit_normalized_rates_by_family"][f]
            assert a is not None and b is not None
            assert abs(a-expected["near_2024_over_2023"])<1e-10
            assert abs(b-expected["remote_2024_over_2023"])<1e-10
            assert abs(a/b-expected["ratio_of_rate_ratios"])<1e-10
            metrics[f]={"relative_annual_rate_ratio_near_over_remote":a/b,
                        "near_2024_over_2023":a,"remote_2024_over_2023":b}
        l=metrics["Lycaenidae"]["relative_annual_rate_ratio_near_over_remote"]
        asserts=study["source_cluster_leave_one_out_5km_Lycaenidae_2023_2024_ratios"]
        assert len(asserts)==4
        for x in asserts:
            rec=[v for v in study["four_source_outfall_near_group_2023_2024"]
                 if v["outfall_id_sha256"]!=x["left_out_original_source_id_sha256"]]
            aa={"2023":{},"2024":{}}
            for yr in ("2023","2024"):
                aa[yr]={"family_counts":{"Lycaenidae":sum(u["match"][yr]["family_counts"]["Lycaenidae"] for u in rec)},
                        "original_site_month_survey_days":sum(u["match"][yr]["original_site_month_survey_days"] for u in rec)}
            nr=calc(aa,"Lycaenidae")
            assert abs(nr-metrics["Lycaenidae"]["remote_2024_over_2023"]*x["near_vs_remote_ratio_of_2024_2023_rate_ratios"])<1e-10
        results[key]={
            "source_status_never_known_randomized_or_changed_2023_2024":True,
            "four_independent_outfalls_behind_23_near_site_rows":4,
            "site_month_year_2023_2024":{yr:{
                g:{"matched_site_n":study["matched_group_year"][g][yr]["original_matched_source_site_n"],
                   "original_visit_days":study["matched_group_year"][g][yr]["original_site_month_survey_days"],
                   "L_count":study["matched_group_year"][g][yr]["family_counts"]["Lycaenidae"],
                   "P_count":study["matched_group_year"][g][yr]["family_counts"]["Pieridae"]}
                for g in ("near5","remote10")} for yr in ("2023","2024")},
            "family_year_rate_ratio_comparator":metrics,
            "source_outfall_leave_one_out_Lycaenidae_near_over_remote_rate_ratio_range":[
                min(x["near_vs_remote_ratio_of_2024_2023_rate_ratios"] for x in asserts),
                max(x["near_vs_remote_ratio_of_2024_2023_rate_ratios"] for x in asserts)],
            "4outfall_source_loo_not_statistically_independent_conclusive_test":True,
        }
    assert results["same_calendar_months_all"]["family_year_rate_ratio_comparator"]["Lycaenidae"]["relative_annual_rate_ratio_near_over_remote"]<.25
    assert results["same_calendar_months_May_Aug"]["family_year_rate_ratio_comparator"]["Lycaenidae"]["relative_annual_rate_ratio_near_over_remote"]<.25
    out={
        "schema":"gesui_DK_2023_2024_source_4outfall_monthmatched_independent_offline_replay_v01",
        "full_source_aggregate_sha256":hashlib.sha256(p.read_bytes()).hexdigest(),
        "all_mapped_family_214_siteyears_preserved":True,
        "same_site_same_month_source_replay_4near_outfalls_2years":results,
        "source_siteyear_data_are_licensed_CC_BY_NC_4_0":True,
        "observational_POSTHOC_not_registered_treatment_or_biochemical_dose":True,
        "no_original_individual_data_republished":True,
    }
    return out

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",type=Path,required=True)
    x=a.parse_args();root=Path(__file__).resolve().parents[1]
    v=audit(root)
    x.out.parent.mkdir(parents=True,exist_ok=True)
    x.out.write_text(json.dumps(v,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:{"ratio":v["same_site_same_month_source_replay_4near_outfalls_2years"][k]["family_year_rate_ratio_comparator"]["Lycaenidae"]["relative_annual_rate_ratio_near_over_remote"],
                         "source_loo":v["same_site_same_month_source_replay_4near_outfalls_2years"][k]["source_outfall_leave_one_out_Lycaenidae_near_over_remote_rate_ratio_range"]}
                       for k in v["same_site_same_month_source_replay_4near_outfalls_2years"]},indent=2))

if __name__=="__main__":
    main()
