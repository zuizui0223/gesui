#!/usr/bin/env python3
"""Independent math/provenance QC for source-frozen 2025 Danish new399.

Recomputes Lycaenidae vs Pieridae+Papilionidae source proportions and
geographic signs from published aggregate site-year counts. Never
downloads or discloses licensed individual DwCA Event/Occurrence rows.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random,statistics
from collections import Counter,defaultdict
from pathlib import Path

DIR="results/ebms_denmark_2023_2024_WWTP_source_only/"
FILE_GATE=DIR+"DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json"
FILE_FAMILY=DIR+"DENMARK_2025_NEW399_SOURCE_LOCKED_ACCEPTED_BUTTERFLY_FAMILY_COUNTS.json"
FILE_SCALE=DIR+"DENMARK_2025_NEW399_LYCAENIDAE_2022_WWTP_AGGREGATE_REGION_COMPONENT_SENSITIVITY.json"
FILE_OLD=DIR+"DENMARK_2025_ORIGINAL107_LYCAENIDAE_OUTOFTIME_FAMILY_AND_MONTH_RESULTS.json"
FILE_ORIGINAL_EVENT=DIR+"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def met(rows):
    L=sum(r["family_individual_count_original"]["Lycaenidae"] for r in rows)
    P=sum(r["family_individual_count_original"]["Pieridae"]+r["family_individual_count_original"]["Papilionidae"] for r in rows)
    days=sum(r["original_any_Occurrence_linked_distinct_date_count"] for r in rows)
    return {"sites":len(rows),"L":L,"others":P,"linkedDays":days,
            "share":L/(L+P) if L+P else None,
            "rate":L/days if days else None,
            "logit":math.log((L+.5)/(P+.5)) if L+P else None}


def check(root):
    path={key:root/name for key,name in
          (("gate",FILE_GATE),("family",FILE_FAMILY),("spatial",FILE_SCALE),
           ("old107",FILE_OLD),("event",FILE_ORIGINAL_EVENT))}
    d={key:json.loads(f.read_text()) for key,f in path.items()}
    gate,family,sc,old,event=(d[x] for x in ("gate","family","spatial","old107","event"))
    assert gate["NEW2025_Lycaenidae_or_any_butterfly_family_individualCount_NOT_READ"]
    assert gate["frozen_before_2025_NEW399_butterfly_counts_gate"]["2025_NEW399_butterfly_family_count_outcome_permitted"]
    assert len(gate["all_original_399_source_only_new_site_2025_documented_visit_metadata"])==399
    assert family["original_source_399_2025_event_any_occurrence_gate_SHA256"]==sha(path["gate"])
    assert sc["source_original_sha256"]["original_licensed_family_count_for_NEW399"]==sha(path["family"])
    records=family["all399_source_site_original_family_counts_private_hash_aggregate_only"]
    assert len(records)==399
    ids={s["site_id_SHA256"] for s in records}
    oldids={s["site_id_SHA256"] for s in event["original_2023_2024_repeated_event_site_source_metadata"]
            if s["repeat_2023_2024_dates_gate"] and s["precise_500m_explicit_uncertainty_gate"]
               and s["rounded5_site_coord_scatter_le2m"]}
    assert len(oldids)==107 and not(ids & oldids)
    sample=[r for r in records if r["pre_2025_family_any_occurrence_linked_site_source_gate_pass"]]
    assert len(sample)==286
    near=[r for r in sample if r["source_band"]=="near5"]
    far=[r for r in sample if r["source_band"]=="remote10"]
    mid=[r for r in sample if r["source_band"]=="intermediate"]
    assert (len(near),len(far),len(mid))==(70,92,124)
    totalN,totalF=met(near),met(far)
    assert (totalN["L"],totalN["others"],totalN["linkedDays"])==(382,680,539)
    assert (totalF["L"],totalF["others"],totalF["linkedDays"])==(644,2602,757)
    oldn=old["original_107site_full_year_three_year_source_band_counts"]["near5"]["2025"]
    oldf=old["original_107site_full_year_three_year_source_band_counts"]["remote10"]["2025"]
    assert oldn["species_counts_by_accepted_family"]["Lycaenidae"]==95
    assert oldn["species_counts_by_accepted_family"]["Pieridae"]==148
    assert oldf["species_counts_by_accepted_family"]["Lycaenidae"]==585
    assert oldf["species_counts_by_accepted_family"]["Pieridae"]==324
    ratio=math.exp(totalN["logit"]-totalF["logit"])
    assert abs(ratio-sc["pooled_new399_2025_butterfly_family_count_vs_daily_recorded_rate"]["pooled_family_oddsratio_near_divided_remote"])<1e-12
    assert abs(totalN["rate"]/totalF["rate"]-sc["pooled_new399_2025_butterfly_family_count_vs_daily_recorded_rate"]["Lycaenidae_per_linked_date_rate_ratio"])<1e-12
    cell=defaultdict(lambda:defaultdict(list))
    for r in sample:cell[r["original_1degree_region"]][r["source_band"]].append(r)
    shared=sorted(x for x,v in cell.items() if v["near5"] and v["remote10"])
    assert len(shared)==6 and shared==sorted(gate["shared_source_regions"])
    region=[]
    equal_site=[]
    for name in shared:
        n=cell[name]["near5"];f=cell[name]["remote10"]
        mn,mf=met(n),met(f)
        l=mn["logit"]-mf["logit"]
        region.append(l)
        def site_mean(a):
            z=[]
            for x in a:
                ly=x["family_individual_count_original"]["Lycaenidae"]
                other=x["family_individual_count_original"]["Pieridae"]+x["family_individual_count_original"]["Papilionidae"]
                if ly+other>0:z.append(math.log((ly+.5)/(other+.5)))
            return statistics.mean(z)
        equal_site.append(site_mean(n)-site_mean(f))
    rmean=statistics.mean(region)
    site_equal=statistics.mean(equal_site)
    assert abs(rmean-sc["equal_1degree_region_unweighted_mean_pooled_family_log_odds_near_minus_remote"])<1e-12
    assert abs(site_equal-sc["equal_region_site_equal_mean_local_family_log_odds_near_minus_remote"])<1e-12
    assert sum(v>0 for v in region)==2
    rng=random.Random(20261110)
    draws=[statistics.mean(rng.choice(region) for _ in region) for _ in range(10000)]
    sd=[statistics.mean(rng.choice(equal_site) for _ in equal_site) for _ in range(10000)]
    ci=lambda a:[sorted(a)[int((len(a)-1)*.025)],sorted(a)[int((len(a)-1)*.975)]]
    assert all(abs(a-b)<1e-12 for a,b in zip(ci(draws),sc["10000_1degree_region_descriptive_bootstrap_cluster95_logOR"]))
    assert all(abs(a-b)<1e-12 for a,b in zip(ci(sd),sc["10000_1degree_region_descriptive_bootstrap_cluster95_site_equal_logOR"]))
    original_outfalls={x["2022_nearest_registered_EEA_FW_outfall_SHA256"] for x in near}
    comp5={x["original_WWTP_spatial_component_5km_SHA256"] for x in near}
    assert len(comp5)==17 and len(original_outfalls)==18
    remote_top=sorted(far,key=lambda x:x["family_individual_count_original"]["Pieridae"],reverse=True)
    top5=sum(x["family_individual_count_original"]["Pieridae"] for x in remote_top[:5])
    assert top5==1321
    return {
        "schema":"gesui_DK_2025_NEW399_original_source_integrity_and_six_area_families_INDEPENDENT_QA_v01",
        "original_SHA256":{k:sha(v) for k,v in path.items()},
        "original_107_sites_excluded_before_new2025_target_family_counts":True,
        "new2025_399_sites_source_locked_before_species":True,
        "2025_2plus_ANY_published_occurrence_dates_n":len(sample),
        "new2025_source_near5":totalN,"new2025_source_remote10":totalF,
        "new2025_pooled_Lycaenidae_family_OR_near_divided_remote":ratio,
        "new2025_Lycaenidae_absolute_rate_near_divided_remote":totalN["rate"]/totalF["rate"],
        "new2025_6_shared_original_1degree_area_pooled_family_logOR":dict(zip(shared,region)),
        "new2025_equal_region_pooled_family_logOR":rmean,
        "new2025_equal_region_family_OR":math.exp(rmean),
        "new2025_equal_region_site_equal_logOR":site_equal,
        "new2025_equal_region_site_equal_OR":math.exp(site_equal),
        "new2025_region_positive_local_family_logOR_n":sum(x>0 for x in region),
        "new2025_region_based_descriptive_bootstrap95_pooled_logOR":ci(draws),
        "new2025_region_based_descriptive_bootstrap95_site_equal_logOR":ci(sd),
        "new2025_near_source_registered_outfall_ID_n":len(original_outfalls),
        "new2025_near_source_5km_geodesic_component_count":len(comp5),
        "new2025_remote_Pieridae_top5_site_individualCount":top5,
        "new2025_remote_Pieridae_top5_fraction":top5/totalF["others"],
        "old_107site_cohort": {
            "original_near_2025_Lycaenidae_share":95/(95+148),
            "original_remote_2025_Lycaenidae_share":585/(585+324),
            "different_sample_selection_and_published_occurrence_coverage_not_causal":True,
        },
        "new399_2025_family_share_direction_differs_from_original107_even_with_same_eBMS_year":True,
        "no_2025_WWTP_downstream_receiving_water_exposure_or_random_treatment_measured":True,
        "Simpson_TYPE_aggregation_not_strict_all_cells_sign_reversal":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    v=check(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(v,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "new_new399_positive_aggregate_family_odds":v["new2025_pooled_Lycaenidae_family_OR_near_divided_remote"],
        "recorded_rate_ratio":v["new2025_Lycaenidae_absolute_rate_near_divided_remote"],
        "equal_1degree_region_odds":v["new2025_equal_region_family_OR"],
        "equal_region_site_odds":v["new2025_equal_region_site_equal_OR"],
        "source_geodesic_WWTP_components":v["new2025_near_source_5km_geodesic_component_count"],
        "remote_Pieridae_top5_contribution":v["new2025_remote_Pieridae_top5_fraction"]},indent=2))

if __name__=="__main__":
    main()
