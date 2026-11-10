#!/usr/bin/env python3
"""Denmark EBMS 2025 OUT-OF-TIME Lycaenidae count result on PRE-FROZEN 107-site cohort.

No new sites, source-distance radius, taxonomy or target family chosen from 2025
observations. Read 2025 original licensed individualCount only after source-only
original 107-site/2025-visit gate exists and passed. Reconcile earlier 214
site-year individualCount and visit-day totals EXACTLY before new result.
"""
from __future__ import annotations

import argparse, hashlib, json, math
from collections import Counter, defaultdict
from datetime import datetime,timezone
from pathlib import Path

import gesui_DK_2023_2024_ebms_same_site_calendar_month_FAMILY_seasoncheck as prior

ROOT=Path("results/ebms_denmark_2023_2024_WWTP_source_only")
GATE=ROOT/"DENMARK_EBMS_2025_107SITE_PREBUTTERFLY_HOLDOUT_SOURCE_GATE.json"
FAMILIES=("Lycaenidae","Pieridae","Nymphalidae","Hesperiidae","Papilionidae")
YRS=(2023,2024,2025)
NEAR="near5"; FAR="remote10";MID="intermediate"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_category(site):
    v=site["nearest_qualified_EEA2022_DK_FW_outfall_km"]
    return NEAR if v<=5 else FAR if v>=10 else MID


def proportion(counts):
    L=counts["Lycaenidae"]
    denom=L+counts["Pieridae"]+counts["Papilionidae"]
    return L/denom if denom>0 else None


def summarize_group(recs):
    result={}
    for year in YRS:
        rows=[x for x in recs if x["year"]==year]
        counts={f:sum(x["family_counts"][f] for x in rows) for f in FAMILIES}
        days=sum(x["visit_days"] for x in rows)
        out={
            "site_locations_n":len(rows),
            "original_distinct_survey_visit_days":days,
            "species_counts_by_accepted_family":counts,
            "Lycaenidae_share_of_three_families_given_positive_count":proportion(counts),
            "siteyears_with_three_family_count_positive":sum(x["family_counts"]["Lycaenidae"]+x["family_counts"]["Pieridae"]+x["family_counts"]["Papilionidae"]>0 for x in rows),
            "siteyears_with_three_family_count_zero":sum(x["family_counts"]["Lycaenidae"]+x["family_counts"]["Pieridae"]+x["family_counts"]["Papilionidae"]==0 for x in rows),
            "family_individualCount_per_original_survey_day":{f:counts[f]/days if days else None for f in FAMILIES},
        }
        result[str(year)]=out
    return result


def compare_rate(summary,year_old,year_new,fam):
    a=summary[str(year_old)]["family_individualCount_per_original_survey_day"][fam]
    b=summary[str(year_new)]["family_individualCount_per_original_survey_day"][fam]
    return b/a if a and b is not None else None


def rate_ratio_of_ratios(near,remote,year_old,year_new):
    o={}
    for f in FAMILIES:
        n=compare_rate(near,year_old,year_new,f)
        r=compare_rate(remote,year_old,year_new,f)
        o[f]={
            "near_time_rate_ratio":n,"remote_time_rate_ratio":r,
            "near_vs_remote_ratio_of_time_rate_ratios":n/r if n is not None and r is not None and r>0 else None,
        }
    return o


def survey_year_month_aggregation(locations, months, abundance, keep_years, predicate):
    recs=[]
    eligible_counts=Counter()
    for loc,site in locations.items():
        accepted=[m for m in range(1,13) if predicate(m) and all(months[(loc,y,m)] for y in keep_years)]
        if not accepted: continue
        accepted_counts=Counter()
        for y in keep_years:
            values={f:sum(abundance[(loc,y,m)][f] for m in accepted) for f in FAMILIES}
            days=sum(len(months[(loc,y,m)]) for m in accepted)
            assert days>0
            recs.append({
                "site_id_SHA256":site["site_id_SHA256"],
                "source_band":source_category(site),
                "source_geo_region":site["country_1deg_lat_lon_cell"],
                "year":y,"visit_days":days,"family_counts":values,
            })
        eligible_counts[source_category(site)]+=1
    return recs,dict(eligible_counts)


def read_exact_frozen_original(root):
    p=root/GATE
    gate=json.loads(p.read_text())
    if gate["schema"]!="gesui_denmark_2025_same107_original_source_ONLY_time_holdout_before_individualCount_v01":
        raise ValueError("Original independent 2025 source-only site cohort was not frozen")
    if not gate["source_sample_gate_prespecified"]["2025_family_outcome_access_permitted"]:
        raise ValueError("STOP. 2025 original source site visit gate failed")
    if gate["2025_butterfly_occurrence_data_rows_opened"]!=0:
        raise ValueError("Original source-only stage included 2025 outcome data")
    sites,families,prior_panel,near_ids,prior_sha=prior.yearsources(root)
    sitehash={s["site_id_SHA256"]:s for s in sites.values()}
    rows=gate["original_2025_repeated_site_year_metadata_fixed"]
    if len(rows)!=107 or len(sitehash)!=107 or {s["site_id_SHA256"] for s in rows}!=set(sitehash):
        raise ValueError("Source gate and original 107 site cohort are not identical")
    for site in rows:
        st=sitehash[site["site_id_SHA256"]]
        if st["visit_days_2025_descriptive"]!=site["2025_original_visit_days"]:
            raise ValueError("2025 source-visit cohort was not fixed")
        if abs(st["nearest_qualified_EEA2022_DK_FW_outfall_km"]-site["2022_source_distance_km"])>1e-10:
            raise ValueError("WWTP source changed")
    if len(set(near_ids.values()))!=4:
        raise ValueError("Exactly four source outfall IDs must be retained")
    for name,path in (
        ("original_site_and_source",prior.SOURCE),
        ("frozen_taxonomy",prior.first.TAX),
        ("original_2023_24_siteyear_counts",str(ROOT/"DENMARK_EBMS_2023_2024_FAMILY_SITEYEAR_COUNTS_QC_BEFORE_WWTP_CONTRAST.json")),
        ("earlier_season_analysis",str(ROOT/"DENMARK_2023_2024_SAME_SITE_MONTH_LYCAENIDAE_POSTHOC_FAMILY_SOURCE_CLUSTER_QC.json"))
    ):
        if sha(root/path)!=gate["original_source_file_SHA256"][name]:
            raise ValueError("Previously frozen study file changed: "+name)
    return gate,sites,families,prior_panel,near_ids


def independent_compare_prior(months,abundance,locations,baseline):
    # Original 2023/24 accepted family/visits were already released, but the new
    # 2025 analysis MUST byte-for-byte reproduce them from source.
    for loc,site in locations.items():
        for year in (2023,2024):
            prior_result=baseline[(site["site_id_SHA256"],year)]
            actual_visit_dates=set()
            for mon in range(1,13):
                actual_visit_dates.update(months[(loc,year,mon)])
            if len(actual_visit_dates)!=prior_result["distinct_original_visit_days"]:
                raise ValueError("Previously audited annual survey-day counts not reproduced")
            for fam in FAMILIES:
                total=sum(abundance[(loc,year,mon)][fam] for mon in range(1,13))
                if total!=prior_result["family_individual_counts"][fam]:
                    raise ValueError("Previously audited family counts not reproduced: "+fam)
    return True


def all_full_year_rows(months,abundance,locations,prior_year_panel):
    rows=[]
    for loc,site in locations.items():
        for yr in YRS:
            all_days=set()
            for mon in range(1,13):
                all_days.update(months[(loc,yr,mon)])
            vals={f:sum(abundance[(loc,yr,m)][f] for m in range(1,13)) for f in FAMILIES}
            if yr<2025 and prior_year_panel:
                old=prior_year_panel[(site["site_id_SHA256"],yr)]
                assert old["distinct_original_visit_days"]==len(all_days)
                assert all(old["family_individual_counts"][f]==vals[f] for f in FAMILIES)
            if yr==2025 and len(all_days)!=site["visit_days_2025_descriptive"]:
                raise ValueError("New 2025 actual visit days do not match prior source-only count")
            rows.append({
                "source_site_id_SHA256":site["site_id_SHA256"],
                "source_band":source_category(site),
                "original_1deg_geographic_cell":site["country_1deg_lat_lon_cell"],
                "year":yr,
                "visit_days":len(all_days),
                "family_counts":vals,
                "original_WWTP_outfall_2022_distance_km":site["nearest_qualified_EEA2022_DK_FW_outfall_km"],
            })
    if len(rows)!=321:
        raise ValueError("Expected 107 sites times 3 observed source years")
    return rows


def run(root):
    gate,sites,families,original_panel,near_ids=read_exact_frozen_original(root)
    # Make the immutable older monthly code read the 2025 observation year as
    # an explicitly predeclared extension; original family/source code unchanged.
    prior.YEARS=YRS
    archive=prior.first.download()
    if hashlib.sha256(archive).hexdigest()!=prior.first.ORIGINAL_SHA:
        raise ValueError("Published eBMS licensed archive SHA differs")
    month,abund,qc=prior.counts_by_month(archive,sites,families)
    assert qc["all_original_Event_rows"]==15606
    assert qc["all_original_Occurrence_rows"]==21617
    if not independent_compare_prior(month,abund,sites,original_panel):
        raise ValueError("2023/24 independent original audit failed")
    allrows=all_full_year_rows(month,abund,sites,original_panel)
    bands={name:summarize_group([r for r in allrows if r["source_band"]==name]) for name in (NEAR,FAR,MID)}
    # Check the site-year repetition and group-year totals against original
    # cohort already frozen before 2025; 2025 outcomes are newly opened now.
    for name,expected in ((NEAR,23),(FAR,38),(MID,46)):
        assert all(bands[name][str(y)]["site_locations_n"]==expected for y in YRS)
        assert bands[name]["2025"]["original_distinct_survey_visit_days"]==gate["repeated_2025_source_only_group_sites_and_days"][name]["2025_original_distinct_survey_days"]

    current=bands[NEAR]["2025"];far=bands[FAR]["2025"]
    direction={
        "predeclared_near_2025_Lycaenidae_per_visit_below_remote":
            current["family_individualCount_per_original_survey_day"]["Lycaenidae"]<far["family_individualCount_per_original_survey_day"]["Lycaenidae"],
        "predeclared_near_2025_Lycaenidae_share_below_remote":
            (current["Lycaenidae_share_of_three_families_given_positive_count"] is not None
             and far["Lycaenidae_share_of_three_families_given_positive_count"] is not None
             and current["Lycaenidae_share_of_three_families_given_positive_count"]<far["Lycaenidae_share_of_three_families_given_positive_count"]),
    }
    direction["both_predeclared_2025_negative_proximity_directions_replicate"]=all(direction.values())

    # Four true near registered-source systems, NOT 23 independent effluent units.
    by_outfall=[]
    for src in sorted(set(near_ids.values())):
        sset={sites[loc]["site_id_SHA256"] for loc,k in near_ids.items() if k==src}
        rs=[x for x in allrows if x["source_site_id_SHA256"] in sset]
        assert len(rs)==len(sset)*3
        by_outfall.append({
            "original_outfall_id_SHA256":src,
            "2022_reported_outfall_status_static_all_three_years":True,
            "distinct_near_WWTP_public_monitor_site_locations":len(sset),
            "source_near_group_2023_2024_2025":summarize_group(rs)})
    if len(by_outfall)!=4 or sum(x["distinct_near_WWTP_public_monitor_site_locations"] for x in by_outfall)!=23:
        raise ValueError("Original four-outfall-source cohort changed")
    loo=[]
    for excluded in by_outfall:
        ids={sites[loc]["site_id_SHA256"] for loc,k in near_ids.items() if k!=excluded["original_outfall_id_SHA256"]}
        sub=summarize_group([x for x in allrows if x["source_site_id_SHA256"] in ids])
        nf=sub["2025"]["family_individualCount_per_original_survey_day"]["Lycaenidae"]
        rf=far["family_individualCount_per_original_survey_day"]["Lycaenidae"]
        ns=sub["2025"]["Lycaenidae_share_of_three_families_given_positive_count"]
        fs=far["Lycaenidae_share_of_three_families_given_positive_count"]
        loo.append({
            "left_out_original_source_outfall_SHA256":excluded["original_outfall_id_SHA256"],
            "remaining_outfall_systems":3,
            "remaining_near_2025_Lycaenidae_rate_vs_remote_ratio":nf/rf if rf>0 else None,
            "remaining_near_2025_Lycaenidae_threefamily_share_minus_remote":ns-fs if ns is not None and fs is not None else None,
            "outfall_clustering_warning_not_confirmatory_significance":True,
        })

    regional={}
    for cell in gate["shared_1degree_near_remote_cells"]:
        regional[cell]={}
        for y in YRS:
            items=[x for x in allrows if x["original_1deg_geographic_cell"]==cell and x["year"]==y]
            n=summarize_group([x for x in items if x["source_band"]==NEAR])["%d"%y]
            f=summarize_group([x for x in items if x["source_band"]==FAR])["%d"%y]
            ns=n["Lycaenidae_share_of_three_families_given_positive_count"]
            fs=f["Lycaenidae_share_of_three_families_given_positive_count"]
            regional[cell][str(y)]={
                "near5":n,"remote10":f,
                "source_response_conditioned_composition_supported_on_both_sides":ns is not None and fs is not None,
                "source_near_minus_remote_pooled_logodds_only_if_both_have_target_families":math.log(ns/(1-ns))-math.log(fs/(1-fs)) if ns is not None and fs is not None and 0<ns<1 and 0<fs<1 else None,
                "source_region_effect_not_independent_WWTP_causality":True
            }

    seasons={}
    configs={
        "same_2024_and_2025_calendar_months":((2024,2025),tuple(range(1,13))),
        "same_2024_and_2025_May_Aug":((2024,2025),(5,6,7,8)),
        "same_2023_2024_2025_calendar_months":(YRS,tuple(range(1,13))),
        "same_2023_2024_2025_May_Aug":(YRS,(5,6,7,8)),
    }
    for name,(years,specified_months) in configs.items():
        rows,elig=survey_year_month_aggregation(sites,month,abund,years,lambda m:m in specified_months)
        subset={group:summarize_group([x for x in rows if x["source_band"]==group]) for group in (NEAR,FAR,MID)}
        info={"original_subset_years":list(years),
              "calendar_months_allowed_prior_to_2025_species_outcomes":list(specified_months),
              "sampled_same_site_same_month_site_count_by_source_band":elig,
              "group_original_survey_day_and_butterfly_family_counts":{g:{str(y):v[str(y)] for y in years} for g,v in subset.items()}}
        if 2024 in years and 2025 in years:
            info["2024_to_2025_family_rate_ratio_of_ratios"]=rate_ratio_of_ratios(subset[NEAR],subset[FAR],2024,2025)
        seasons[name]=info

    out={
        "schema":"gesui_Denmark_2025_107_site_eBMS_outoftime_source_locked_Lycaenidae_individualCount_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "preoutcome_original_2025_gate_SHA256":sha(root/GATE),
        "licensed_official_eBMS_v14_original_ZIP_SHA256":prior.first.ORIGINAL_SHA,
        "original_Event_and_Occurrence_rows_in_full_archive":dict(qc),
        "all_prior_214_2023_2024_family_individual_sums_and_days_reproduced":True,
        "original_2025_new_site_year_accepted_aggregate_individualCounts":[x for x in allrows if x["year"]==2025],
        "source_static_2022_EEA_never_measured_2025_actual_downstream_exposure":True,
        "original_107site_full_year_three_year_source_band_counts":bands,
        "2025_predeclared_proximity_directional_outoftime_test":direction,
        "cross_year_family_2024_vs_2025_near_remote_rate_ratio_of_ratios":rate_ratio_of_ratios(bands[NEAR],bands[FAR],2024,2025),
        "cross_year_family_2023_vs_2025_near_remote_rate_ratio_of_ratios":rate_ratio_of_ratios(bands[NEAR],bands[FAR],2023,2025),
        "true_registered_WWTP_near_outfall_cluster_count":len(by_outfall),
        "original_4_source_near_annual_group_counts":by_outfall,
        "original_2025_leave_one_WWTP_outfall_out":loo,
        "original_shared_1degree_region_annual_family_outcomes":regional,
        "same_site_same_calendar_months_precommitted_sensitivities":seasons,
        "no_2025_regulatory_treatment_year_or_downstream_nutrient_dose_measured":True,
        "no_individual_butterfly_personal_ID_or_raw_DwCA_extension_repotted":True,
        "outcome_is_public_licensed_eBMS_individualCount_not_complete_census":True,
        "report_all_signs_not_only_positive_WWTP_butterfly_results":True,
    }
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    value=run(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(value,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "source_2025_pass_before_species":True,
        "original_214_previous_family_counts_match":value["all_prior_214_2023_2024_family_individual_sums_and_days_reproduced"],
        "2025_pooled_near_and_remote":{g:value["original_107site_full_year_three_year_source_band_counts"][g]["2025"] for g in (NEAR,FAR)},
        "2025_outoftime_result":value["2025_predeclared_proximity_directional_outoftime_test"],
        "2024_to_2025_family_change":value["cross_year_family_2024_vs_2025_near_remote_rate_ratio_of_ratios"],
        "2024_25_season_matched":value["same_site_same_calendar_months_precommitted_sensitivities"]["same_2024_and_2025_calendar_months"]["2024_to_2025_family_rate_ratio_of_ratios"],
        "near_WWTP_source_systems":value["true_registered_WWTP_near_outfall_cluster_count"]
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
