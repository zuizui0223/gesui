#!/usr/bin/env python3
"""Audit Simpson-like aggregation reversal in a NEW399 Denmark 2025 EBMS sample.

Original fully source-locked accepted family individualCount vs a static 2022
registered WWTP outfall proximity proxy, before this script is fit.
Posthoc effect interpretation: compositional ratios are NOT butterfly abundance
or evidence of a wastewater causal effect. Six shared 1-degree regions, 17
5-km connected original WWTP point source units; no treatment randomization.
"""
from __future__ import annotations

import argparse,hashlib,json,math,random,statistics
from collections import Counter,defaultdict
from pathlib import Path

ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
SOURCE=ROOT+"DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json"
FAMILY=ROOT+"DENMARK_2025_NEW399_SOURCE_LOCKED_ACCEPTED_BUTTERFLY_FAMILY_COUNTS.json"
EVENT=ROOT+"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
TARGET=("Lycaenidae","Pieridae","Papilionidae")
FAMILIES=("Lycaenidae","Pieridae","Nymphalidae","Hesperiidae","Papilionidae")
SEED=20261110


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def met(items):
    n=len(items)
    nday=sum(x["original_any_Occurrence_linked_distinct_date_count"] for x in items)
    counts={f:sum(x["family_individual_count_original"][f] for x in items) for f in FAMILIES}
    L=counts["Lycaenidae"];other=counts["Pieridae"]+counts["Papilionidae"]
    valid=[x for x in items if x["target_3family_positive_individuals"]>0]
    site_logs=[math.log((x["family_individual_count_original"]["Lycaenidae"]+.5)/
        (x["family_individual_count_original"]["Pieridae"]+x["family_individual_count_original"]["Papilionidae"]+.5))
        for x in valid]
    return {
        "site_rows":n,"original_any_published_Occurrence_site_dates":nday,
        "recorded_counts":counts,
        "site_rows_with_any_of_three_families":len(valid),
        "Lycaenidae_share_L_over_3_families":L/(L+other) if L+other else None,
        "Lycaenidae_per_any_published_Occurrence_date":L/nday if nday else None,
        "site_equal_mean_3family_positive_log_Lyca_over_other":statistics.mean(site_logs) if site_logs else None,
        "pooled_log_Lyca_over_other_continuity_0_5":math.log((L+.5)/(other+.5)) if L+other else None,
        "zero_target_site_not_assigned_arbitrary_compositional_logit":True,
    }


def h(a,b):
    r=math.pi/180
    dx=(a["lat"]-b["lat"])*r;dy=(a["lon"]-b["lon"])*r
    z=math.sin(dx/2)**2+math.cos(a["lat"]*r)*math.cos(b["lat"]*r)*math.sin(dy/2)**2
    return 12742*math.asin(min(1.,math.sqrt(z)))


def components(stations,threshold):
    """Spatially connected source/remote locations, not independent populations."""
    k=list(stations)
    pa={key:key for key in k}
    def find(x):
        while x!=pa[x]:
            pa[x]=pa[pa[x]]
            x=pa[x]
        return x
    for i,a in enumerate(k):
        for b in k[i+1:]:
            if h(stations[a],stations[b])<=threshold:
                pa[find(a)]=find(b)
    size=Counter(find(x) for x in k)
    return {"original_location_points":len(k),
            "source_independent_km_is_heuristic_not_true_river_independence":threshold,
            "geodesic_connected_components_n":len(size),
            "component_site_count_descending":sorted(size.values(),reverse=True)}


def contrasts(near,remote):
    n=met(near);f=met(remote)
    ll=n["pooled_log_Lyca_over_other_continuity_0_5"]
    ff=f["pooled_log_Lyca_over_other_continuity_0_5"]
    return {
        "original_near5":n,"original_remote10":f,
        "near_minus_remote_pooled_share_percentage_points":
            100*(n["Lycaenidae_share_L_over_3_families"]-f["Lycaenidae_share_L_over_3_families"])
            if n["Lycaenidae_share_L_over_3_families"] is not None and f["Lycaenidae_share_L_over_3_families"] is not None else None,
        "pooled_family_logodds_near_minus_remote":ll-ff if ll is not None and ff is not None else None,
        "pooled_family_oddsratio_near_divided_remote":math.exp(ll-ff) if ll is not None and ff is not None else None,
        "Lycaenidae_per_linked_date_rate_ratio":n["Lycaenidae_per_any_published_Occurrence_date"]/
            f["Lycaenidae_per_any_published_Occurrence_date"] if f["Lycaenidae_per_any_published_Occurrence_date"] else None,
        "site_equal_target_family_logodds_near_minus_remote":
            n["site_equal_mean_3family_positive_log_Lyca_over_other"]-
            f["site_equal_mean_3family_positive_log_Lyca_over_other"]
            if n["site_equal_mean_3family_positive_log_Lyca_over_other"] is not None and f["site_equal_mean_3family_positive_log_Lyca_over_other"] is not None else None,
    }


def percent(a,p):
    v=sorted(a)
    return v[int((len(v)-1)*p)]


def analyze(root):
    paths={"original_independent_NEW399_preoutcome_source_manifest":root/SOURCE,
           "original_licensed_family_count_for_NEW399":root/FAMILY,
           "original_preoutcome_Denmark_eBMS_event_source":root/EVENT}
    src=json.loads(paths["original_independent_NEW399_preoutcome_source_manifest"].read_text())
    dat=json.loads(paths["original_licensed_family_count_for_NEW399"].read_text())
    sourceSite=json.loads(paths["original_preoutcome_Denmark_eBMS_event_source"].read_text())
    assert src["frozen_before_2025_NEW399_butterfly_counts_gate"]["2025_NEW399_butterfly_family_count_outcome_permitted"]
    assert dat["focal_family_outcome_support_GATE"]["original_new399_2025_family_effect_descriptive_contrasts_permitted"]
    assert dat["original_source_399_2025_event_any_occurrence_gate_SHA256"]==sha(paths["original_independent_NEW399_preoutcome_source_manifest"])
    assert len(dat["all399_source_site_original_family_counts_private_hash_aggregate_only"])==399
    eligible=[x for x in dat["all399_source_site_original_family_counts_private_hash_aggregate_only"]
              if x["pre_2025_family_any_occurrence_linked_site_source_gate_pass"]]
    assert len(eligible)==286
    near=[r for r in eligible if r["source_band"]=="near5"]
    remote=[r for r in eligible if r["source_band"]=="remote10"]
    assert len(near)==70 and len(remote)==92
    original=contrasts(near,remote)
    byregion=defaultdict(lambda:defaultdict(list))
    for x in eligible:byregion[x["original_1degree_region"]][x["source_band"]].append(x)
    shared=sorted(r for r,groups in byregion.items() if groups["near5"] and groups["remote10"])
    assert len(shared)==6 and shared==sorted(src["shared_source_regions"])
    cellresults=[]
    for cell in shared:
        a=contrasts(byregion[cell]["near5"],byregion[cell]["remote10"])
        cellresults.append({
            "region_1degree":cell,
            "near_registered_WWTP_outfall_point_IDs":len(set(x["2022_nearest_registered_EEA_FW_outfall_SHA256"] for x in byregion[cell]["near5"])),
            "source_near_and_remote_comparison":a,
            "local_family_share_direction_positive":
                a["near_minus_remote_pooled_share_percentage_points"]>0,
            "local_absolute_Lycaenidae_recorded_rate_direction_positive":
                a["Lycaenidae_per_linked_date_rate_ratio"]>1,
        })
    region_log=[x["source_near_and_remote_comparison"]["pooled_family_logodds_near_minus_remote"] for x in cellresults]
    region_site_log=[x["source_near_and_remote_comparison"]["site_equal_target_family_logodds_near_minus_remote"] for x in cellresults]
    region_rate=[math.log(x["source_near_and_remote_comparison"]["Lycaenidae_per_linked_date_rate_ratio"])
                 for x in cellresults if x["source_near_and_remote_comparison"]["Lycaenidae_per_linked_date_rate_ratio"]>0]
    rmean=statistics.mean(region_log)
    site_mean=statistics.mean(region_site_log)
    ratemean=statistics.mean(region_rate)
    rand=random.Random(SEED)
    draws=[statistics.mean(rand.choice(region_log) for j in range(len(region_log))) for i in range(10000)]
    site_draws=[statistics.mean(rand.choice(region_site_log) for j in range(len(region_site_log))) for i in range(10000)]
    nearcomponents=defaultdict(list)
    for x in near:nearcomponents[x["original_WWTP_spatial_component_5km_SHA256"]].append(x)
    leavecomp=[]
    for uid,rows in sorted(nearcomponents.items()):
        n=[x for x in near if x["original_WWTP_spatial_component_5km_SHA256"]!=uid]
        result=contrasts(n,remote)
        leavecomp.append({"removed_component_SHA256":uid,
                          "removed_near_source_site_rows":len(rows),
                          "remove_original_source_comp_remaining_n":len(nearcomponents)-1,
                          "near_vs_remote_3family_pooled_OR":result["pooled_family_oddsratio_near_divided_remote"],
                          "near_vs_remote_recorded_Lycaenidae_per_linked_day_rate_ratio":result["Lycaenidae_per_linked_date_rate_ratio"]})
    lool=[]
    for i,c in enumerate(shared):
        keep=[x for x in cellresults if x["region_1degree"]!=c]
        m=statistics.mean(x["source_near_and_remote_comparison"]["pooled_family_logodds_near_minus_remote"] for x in keep)
        lool.append({"excluded_region":c,"mean_remaining_region_equal_logOR":m,"mean_remaining_region_equal_OR":math.exp(m)})
    nsum=met(near);rsum=met(remote)
    rawnear={x["site_id_SHA256"] for x in near}
    rawfar={x["site_id_SHA256"] for x in remote}
    allsites=sourceSite["original_2023_2024_repeated_event_site_source_metadata"]
    newsourcecoords={s["site_id_SHA256"]:s for s in allsites if s["site_id_SHA256"] in rawnear or s["site_id_SHA256"] in rawfar}
    if len(newsourcecoords)!=len(near)+len(remote):
        raise ValueError("Original WWTP source-only 2025 site metadata missing")
    remote_sites={sid:{"lat":p["lat"],"lon":p["lon"]}
                  for sid,p in newsourcecoords.items() if sid in rawfar}
    nearby_sites={sid:{"lat":p["lat"],"lon":p["lon"]}
                  for sid,p in newsourcecoords.items() if sid in rawnear}
    source_and_geography={
        "original_near5_registered_outfall_site_count":len(nearcomponents),
        "remotely_monitored_geodesic_components":{
            str(d)+"km":components(remote_sites,d) for d in (0.5,2,5,10)},
        "near_site_geodesic_components":{
            str(d)+"km":components(nearby_sites,d) for d in (0.5,2,5,10)},
        "source_near_WWTP_cluster_5km_sourceID_grouping":{
            uid:met(items) for uid,items in sorted(nearcomponents.items())},
    }
    nL=nsum["recorded_counts"]["Lycaenidae"]
    remotePierTop=sorted(remote,key=lambda x:x["family_individual_count_original"]["Pieridae"],reverse=True)
    remotePierCount=rsum["recorded_counts"]["Pieridae"]
    deleted={}
    for n in (1,3,5,6,10):
        subset=remotePierTop[n:]
        a=contrasts(near,subset)
        deleted[str(n)]={
            "excluded_largest_Pieridae_recorded_sites":n,
            "removed_Pieridae_individualCount":remotePierCount-met(subset)["recorded_counts"]["Pieridae"],
            "removed_Pieridae_fraction_of_original_remote":(remotePierCount-met(subset)["recorded_counts"]["Pieridae"])/remotePierCount,
            "pooled_family_oddsratio_near_divided_remote":a["pooled_family_oddsratio_near_divided_remote"],
            "pooled_share_difference_pp":a["near_minus_remote_pooled_share_percentage_points"],
            "NOT_predeclared_evidence_selection_on_response":True,
        }

    output={
        "schema":"gesui_DK_2025_NEW399_near_source_3family_SIMPSOn_scale_sensitivity_posthoc_v01",
        "source_original_sha256":{k:sha(v) for k,v in paths.items()},
        "all_effect_checks_posthoc_on_known_new399_2025_direction":True,
        "accepted_source_first_2025_new_site_sample_not_2023_24_old107":True,
        "new399_total_2025_Event_source_candidates":399,
        "2025_family_eligibility_2plus_days_linked_any_original_occurrence_n":len(eligible),
        "registered_WWTP_source_point_IDs_near_2025_are_not_hydrological_exposure_proof":True,
        "pooled_new399_2025_butterfly_family_count_vs_daily_recorded_rate":original,
        "preoriginal_original_1degree_region_stratified_analyses":cellresults,
        "same_1degree_region_n":len(shared),
        "same_1degree_regions_positive_family_share_relation_n":sum(x>0 for x in region_log),
        "same_1degree_regions_positive_abs_rate_n":sum(x>0 for x in region_rate),
        "equal_1degree_region_unweighted_mean_pooled_family_log_odds_near_minus_remote":rmean,
        "equal_1degree_region_geometric_mean_pooled_family_OR_near_divided_remote":math.exp(rmean),
        "equal_region_site_equal_mean_local_family_log_odds_near_minus_remote":site_mean,
        "equal_region_site_equal_OR_near_divided_remote":math.exp(site_mean),
        "equal_region_geometric_Lycaenidae_rate_near_vs_remote_ratio":math.exp(ratemean),
        "10000_1degree_region_descriptive_bootstrap_cluster95_logOR":[percent(draws,.025),percent(draws,.975)],
        "10000_1degree_region_descriptive_bootstrap_cluster95_site_equal_logOR":[percent(site_draws,.025),percent(site_draws,.975)],
        "region_bootstrap_seed":SEED,
        "six_region_leave_one_out_unweighted_OR":lool,
        "5km_near_registered_WWTP_spatial_component_leave_one_out":leavecomp,
        "source_near_and_remote_geodesic_pseudoreplication":source_and_geography,
        "posthoc_top_remote_Pieridae_site_leverage":deleted,
        "2025_new399_absolute_Lycaenidae_per_original_recorded_survey_date_NOT_total_true_abundance":True,
        "Simpson_LIKE_not_strict_every_subgroup_direction":True,
        "no_2025_actual_organic_N_sodium_phosphate_receiving_water_butterfly_exposure_measured":True,
        "no_sewage_population_effect_inferred":True,
    }
    assert abs(original["pooled_family_oddsratio_near_divided_remote"]-2.26)<0.2
    assert len(cellresults)==6
    assert sum(x>0 for x in region_log)==2
    assert abs(rmean+0.09652783660906716)<1e-9
    assert abs(site_mean+0.6031584942459333)<1e-9
    return output


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    v=analyze(root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(v,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"pooled":v["pooled_new399_2025_butterfly_family_count_vs_daily_recorded_rate"],
        "n_regions":v["same_1degree_region_n"],
        "share_regions_positive":v["same_1degree_regions_positive_family_share_relation_n"],
        "equal_region_OR":v["equal_1degree_region_geometric_mean_pooled_family_OR_near_divided_remote"],
        "site_equal_region_OR":v["equal_region_site_equal_OR_near_divided_remote"],
        "region_bootstrap_CI":v["10000_1degree_region_descriptive_bootstrap_cluster95_logOR"],
        "top_far_Pieridae_leverage":v["posthoc_top_remote_Pieridae_site_leverage"]},indent=2))


if __name__=="__main__":
    main()
