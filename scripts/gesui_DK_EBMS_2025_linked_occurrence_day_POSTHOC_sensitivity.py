#!/usr/bin/env python3
"""Independent posthoc quality sensitivity: 2025 original EBMS Event without Occurrence.

Do not equate an original eBMS Event lacking an Occurrence extension row with
zero butterflies. Use ANY published occurrence ID linkage (independent of
target family) as an imperfect indicator of documented observation effort.
This sensitivity is posthoc AFTER seeing original 2025 butterfly counts.
"""
from __future__ import annotations

import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path

ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
SOURCE=ROOT+"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
PANEL=ROOT+"DENMARK_EBMS_2023_2024_FAMILY_SITEYEAR_COUNTS_QC_BEFORE_WWTP_CONTRAST.json"
NEW2025=ROOT+"DENMARK_2025_ORIGINAL107_LYCAENIDAE_OUTOFTIME_FAMILY_AND_MONTH_RESULTS.json"
COMPLETE=ROOT+"DK_EBMS_2025_ORIGINAL_EVENT_OCCURRENCE_PUBLICATION_MISSINGNESS_QA.json"
FAMILIES=("Lycaenidae","Pieridae","Nymphalidae","Hesperiidae","Papilionidae")
YEARS=(2023,2024,2025)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def hav(a,b):
    r=math.pi/180
    dl=(a["lat"]-b["lat"])*r; dg=(a["lon"]-b["lon"])*r
    v=math.sin(dl/2)**2+math.cos(a["lat"]*r)*math.cos(b["lat"]*r)*math.sin(dg/2)**2
    return 12742*math.asin(min(1,math.sqrt(v)))

def cls(s):
    d=s["nearest_qualified_EEA2022_DK_FW_outfall_km"]
    return "near5" if d<=5 else "remote10" if d>=10 else "intermediate"

def summarize(rows):
    counts={fam:sum(x["counts"][fam] for x in rows) for fam in FAMILIES}
    days=sum(x["link_days"] for x in rows)
    eventdays=sum(x["event_days"] for x in rows)
    L=counts["Lycaenidae"];P=counts["Pieridae"];Pa=counts["Papilionidae"]
    return {
        "n_source_linked_siteyears":len(rows),
        "published_Occurrence_linked_distinct_survey_days":days,
        "original_unfiltered_Event_distinct_days":eventdays,
        "recorded_individualCount_family":counts,
        "Lycaenidae_per_Occurrence_linked_date":L/days if days else None,
        "Lycaenidae_share_of_target_three_families":L/(L+P+Pa) if L+P+Pa else None,
        "family_individualCount_per_Occurrence_linked_date":{fam:counts[fam]/days if days else None for fam in FAMILIES},
        "zero_target_family_count_siteyears":sum(x["counts"]["Lycaenidae"]+x["counts"]["Pieridae"]+x["counts"]["Papilionidae"]==0 for x in rows)
    }

def rate_ratio(a,b,f):
    y=a["family_individualCount_per_Occurrence_linked_date"][f]
    z=b["family_individualCount_per_Occurrence_linked_date"][f]
    return z/y if y and z is not None else None

def perform(root):
    files={"original_sites":root/SOURCE,"2023_2024_family":root/PANEL,
           "2025_original_new_family":root/NEW2025,"2023_2025_Event_Occ_linkage_QA":root/COMPLETE}
    obj={k:json.loads(v.read_text(encoding="utf-8")) for k,v in files.items()}
    s=obj["original_sites"];old=obj["2023_2024_family"];new=obj["2025_original_new_family"];qa=obj["2023_2025_Event_Occ_linkage_QA"]
    assert qa["only_event_and_occurrence_ID_join_not_family_counts_read_in_this_stage"]
    assert len(old["site_year_family_counts_original_public_CC_BY_NC_aggregate"])==214
    assert len(new["original_2025_new_site_year_accepted_aggregate_individualCounts"])==107
    assert len(qa["all_source_site_year_published_record_link_support"])==321
    assert new["all_prior_214_2023_2024_family_individual_sums_and_days_reproduced"]
    assert new["preoutcome_original_2025_gate_SHA256"]
    sites={x["site_id_SHA256"]:x for x in s["original_2023_2024_repeated_event_site_source_metadata"]
           if x["repeat_2023_2024_dates_gate"] and x["precise_500m_explicit_uncertainty_gate"] and x["rounded5_site_coord_scatter_le2m"]}
    assert len(sites)==107
    linked={(x["source_rounded_site_hash"],x["original_year"]):x for x in qa["all_source_site_year_published_record_link_support"]}
    raw={}
    for x in old["site_year_family_counts_original_public_CC_BY_NC_aggregate"]:
        k=(x["source_rounded5_site_hash"],x["site_year"])
        raw[k]={"counts":x["family_individual_counts"],"original_days":x["distinct_original_visit_days"]}
    for x in new["original_2025_new_site_year_accepted_aggregate_individualCounts"]:
        k=(x["source_site_id_SHA256"],x["year"])
        if k in raw or x["year"]!=2025:raise ValueError("Unexpected overlapping original 2025 counted site/year")
        raw[k]={"counts":x["family_counts"],"original_days":x["visit_days"]}
    assert len(raw)==321 and set(raw)==set(linked)
    combined={}
    for (site,year),v in raw.items():
        q=linked[(site,year)]
        assert q["original_Event_distinct_dates_n"]==v["original_days"]
        assert 0<=q["original_distinct_dates_with_at_least_one_Occurrence_n"]<=v["original_days"]
        if q["original_distinct_dates_with_at_least_one_Occurrence_n"]==0:
            assert sum(v["counts"][f] for f in FAMILIES)==0
        combined[(site,year)]={"id":site,"year":year,"source_band":cls(sites[site]),
                               "counts":v["counts"],"event_days":v["original_days"],
                               "link_days":q["original_distinct_dates_with_at_least_one_Occurrence_n"],
                               "original_Occurrence_rows":q["original_Occurrence_extension_rows_linked"]}
    out={
        "schema":"gesui_DK_EBMS_2025_source_band_quality_corrected_LINKED_occurrence_survey_days_posthoc_v01",
        "original_source_SHA256":{k:sha(v) for k,v in files.items()},
        "analysis_POSTHOC_after_2025_source_proximity_family_result_seen":True,
        "linked_Occurrence_rows_are_only_a_proxy_of_reported_survey_completeness":True,
        "Event_without_Occurrence_is_NA_in_sensitivity_not_confirmed_biological_zero":True,
        "source_static_2022_WWTP_near_distance_not_treatment_or_biochemical_dose":True,
        "original_full_siteyears":321,
        "2025_original_source_outfall_count":4,
        "yearbands_orig_All107_siteyears":{},
        "quality_matched_2024_2025_linked_2plus_same_sites":{},
        "quality_matched_2023_2024_2025_linked_2plus_same_sites":{},
        "original_2025_sourcecluster_reporting_completeness":[],
        "report_limited_subsets_not_population_effect":True
    }
    for yr in YEARS:
        for g in ("near5","remote10","intermediate"):
            vals=[r for r in combined.values() if r["year"]==yr and r["source_band"]==g]
            eligible=[x for x in vals if x["link_days"]>=2]
            out["yearbands_orig_All107_siteyears"][f"{yr}_{g}"]={
                "original_total_site_years":len(vals),
                "siteyears_2plus_any_occurrence_linked_days":len(eligible),
                "unlinked_2plus_event_rows_not_false_zero":len(vals)-len(eligible),
                "all_linked_day_denominator_descriptive":summarize(vals),
                "posthoc_2plus_linked_day_denominator_subset":summarize(eligible)
            }
    for yrs,key in [((2024,2025),"quality_matched_2024_2025_linked_2plus_same_sites"),
                    ((2023,2024,2025),"quality_matched_2023_2024_2025_linked_2plus_same_sites")]:
        matches={g:[site for site,s in sites.items() if cls(s)==g and
                    all(combined[(site,y)]["link_days"]>=2 for y in yrs)]
                 for g in ("near5","remote10","intermediate")}
        out[key]={"same_original_geographic_site_required_all_source_years":True,
                  "quality_gate_at_least_2_separate_dates_with_any_original_Occurrence_record_per_year":True,
                  "eligible_same_source_site_n":{g:len(v) for g,v in matches.items()},
                  "group_year":{},"relative_Lycaenidae_change_2024_to_2025":None}
        for g in matches:
            out[key]["group_year"][g]={}
            for yr in yrs:
                selected=[combined[(site,yr)] for site in matches[g]]
                out[key]["group_year"][g][str(yr)]=summarize(selected)
        if 2024 in yrs:
            a=out[key]["group_year"]["near5"];b=out[key]["group_year"]["remote10"]
            effects={}
            for f in FAMILIES:
                rn=rate_ratio(a["2024"],a["2025"],f)
                rf=rate_ratio(b["2024"],b["2025"],f)
                effects[f]={"near_2025_over_2024":rn,"remote_2025_over_2024":rf,
                            "rate_ratio_of_ratios":rn/rf if rn is not None and rf is not None and rf else None}
            out[key]["relative_Lycaenidae_change_2024_to_2025"]=effects
    sources=s["original_2022_EEA_DK_qualifying_Fw_outfalls"]
    cluster=defaultdict(list)
    for site,x in sites.items():
        if cls(x)!="near5":continue
        ranked=sorted((hav(x,z),z["outfall_id_SHA256"]) for z in sources)
        if abs(ranked[0][0]-x["nearest_qualified_EEA2022_DK_FW_outfall_km"])>1e-8:
            raise ValueError("Original source location and outlet membership changed")
        cluster[ranked[0][1]].append(site)
    assert len(cluster)==4
    for outfall,ids in sorted(cluster.items()):
        d={}
        for yr in YEARS:
            all_recs=[combined[(s,yr)] for s in ids]
            active=[x for x in all_recs if x["link_days"]>=2]
            d[str(yr)]={"original_event_site_count":len(all_recs),
                        "original_sites_with_2plus_published_Occurrence_dates":len(active),
                        "total_original_Events":sum(linked[(s,yr)]["original_Events_n"] for s in ids),
                        "total_original_Events_with_any_Occurrence":sum(linked[(s,yr)]["original_Events_with_at_least_one_Occurrence_n"] for s in ids),
                        "Lycaenidae_individual_count_original_recorded":sum(x["counts"]["Lycaenidae"] for x in all_recs)}
        out["original_2025_sourcecluster_reporting_completeness"].append({"outfall_id_SHA256":outfall,"original_near_site_rows":len(ids),"years":d})
    assert sum(x["original_near_site_rows"] for x in out["original_2025_sourcecluster_reporting_completeness"])==23
    out["warning_entire_WWTP_cluster_without_published_2025_occurrences"]=[
        x["outfall_id_SHA256"] for x in out["original_2025_sourcecluster_reporting_completeness"]
        if x["years"]["2025"]["total_original_Events_with_any_Occurrence"]==0]
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=perform(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"sameSite_2024_25":out["quality_matched_2024_2025_linked_2plus_same_sites"],
                      "sameSite_2023_25_N":out["quality_matched_2023_2024_2025_linked_2plus_same_sites"]["eligible_same_source_site_n"],
                      "WWTP_2025_cluster_completeness":out["original_2025_sourcecluster_reporting_completeness"],
                      "zero_occurrence_source_groups":len(out["warning_entire_WWTP_cluster_without_published_2025_occurrences"])},indent=2))

if __name__=="__main__":
    main()
