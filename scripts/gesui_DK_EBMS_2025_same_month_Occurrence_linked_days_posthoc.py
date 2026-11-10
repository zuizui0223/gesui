#!/usr/bin/env python3
"""Posthoc original Danish eBMS 2023-25 SAME SITE × SAME CALENDAR MONTH ×
ANY OCCURRENCE-LINKED survey-day audit for Lycaenidae and butterfly comparators.

Original 2025 source-year design already completed. This reinterprets Event
rows without attached published Occurrence as unknown, not true-zero.
Reuses original licensed DWCA SHA + 107 sites + blinded GBIF family mapping;
source-group data are not randomized wastewater exposures.
"""
from __future__ import annotations

import argparse,hashlib,io,json,math,zipfile
import xml.etree.ElementTree as ET
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import gesui_DK_2023_2024_ebms_same_site_calendar_month_FAMILY_seasoncheck as monthly

SOURCE_ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
GATE=SOURCE_ROOT+"DENMARK_EBMS_2025_107SITE_PREBUTTERFLY_HOLDOUT_SOURCE_GATE.json"
2025FAMILY=SOURCE_ROOT+"DENMARK_2025_ORIGINAL107_LYCAENIDAE_OUTOFTIME_FAMILY_AND_MONTH_RESULTS.json"
LINK_QC=SOURCE_ROOT+"DK_EBMS_2025_ORIGINAL_EVENT_OCCURRENCE_PUBLICATION_MISSINGNESS_QA.json"
YEARS=(2023,2024,2025)
FAMILIES=("Lycaenidae","Pieridae","Nymphalidae","Hesperiidae","Papilionidae")
CONFIG={
    "2024_2025_all_common_months":((2024,2025),range(1,13)),
    "2024_2025_May_Aug":((2024,2025),(5,6,7,8)),
    "2023_2024_2025_all_common_months":(YEARS,range(1,13)),
    "2023_2024_2025_May_Aug":(YEARS,(5,6,7,8)),
}

def sha(b):return hashlib.sha256(b).hexdigest()

def source_band(site):
    d=site["nearest_qualified_EEA2022_DK_FW_outfall_km"]
    return "near5" if d<=5 else "remote10" if d>=10 else "intermediate"

def ratio(x,y):
    return (x/y) if y else None

def aggregate(rows):
    counts={f:sum(x["family_counts"][f] for x in rows) for f in FAMILIES}
    original_days=sum(x["original_all_event_days_for_kept_months"] for x in rows)
    linked_days=sum(x["occurrence_linked_days_for_kept_months"] for x in rows)
    L=counts["Lycaenidae"];P=counts["Pieridae"];Pa=counts["Papilionidae"]
    return {
        "qualified_site_rows":len(rows),
        "original_all_event_days_in_matched_months":original_days,
        "any_Occurrence_linked_distinct_survey_days_in_matched_months":linked_days,
        "family_individualCount_original":counts,
        "Lycaenidae_count_per_any_original_Occurrence_linked_day":ratio(L,linked_days),
        "Lycaenidae_target_three_family_share":ratio(L,L+P+Pa),
        "family_per_any_original_Occurrence_linked_day":{f:ratio(counts[f],linked_days) for f in FAMILIES},
        "sites_with_no_Lycaenidae_count_in_matched_months":sum(x["family_counts"]["Lycaenidae"]==0 for x in rows),
        "sites_with_no_target_three_family_counts":sum(x["family_counts"]["Lycaenidae"]+x["family_counts"]["Pieridae"]+x["family_counts"]["Papilionidae"]==0 for x in rows),
    }

def month_linked_days(raw,sites):
    """Event core and Occurrence extension JOIN only, no taxonomy/abundance read here."""
    event_lookup={}
    events=defaultdict(set)
    linked=defaultdict(set)
    qc=Counter()
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        x=ET.fromstring(z.read("meta.xml"))
        ec=x.find(monthly.NS+"core")
        ext=[u for u in x.findall(monthly.NS+"extension") if u.attrib.get("rowType")=="http://rs.tdwg.org/dwc/terms/Occurrence"]
        if ec is None or len(ext)!=1:raise ValueError("Unexpected eBMS archive shape")
        ef,efile,es,efield=monthly.first.core_terms(ec)
        of,ofile,os,ofield=monthly.first.core_terms(ext[0])
        idindex=efield("id");occidindex=ofield("coreid")
        for line in monthly.first.stream_csv(z,efile,es):
            qc["all_original_Event_rows"]+=1
            if len(line)<=max(idindex,ef["eventDate"],ef["decimalLatitude"],ef["decimalLongitude"]):
                raise ValueError("Original Event archive has truncated row")
            event_id=line[idindex].strip()
            if not event_id or event_id in event_lookup:raise ValueError("Duplicate original Event key")
            date=line[ef["eventDate"]].strip()[:10]
            try:
                yr=int(date[:4]);mon=int(date[5:7]);day=int(date[8:10])
                lat=float(line[ef["decimalLatitude"]]);lon=float(line[ef["decimalLongitude"]])
            except (ValueError,TypeError):continue
            if yr not in YEARS or not(1<=mon<=12 and 1<=day<=31) or not(math.isfinite(lat) and math.isfinite(lon)):
                continue
            loc=f"{lat:.5f}|{lon:.5f}"
            if loc not in sites:continue
            event_lookup[event_id]=(loc,yr,mon,date)
            events[(loc,yr,mon)].add(date)
            qc["event_rows_107original_sites_2023_25"]+=1
        for line in monthly.first.stream_csv(z,ofile,os):
            qc["all_original_Occurrence_rows"]+=1
            if len(line)<=occidindex:raise ValueError("Original Occurrence extension row truncated")
            identifier=line[occidindex].strip()
            if identifier in event_lookup:
                loc,year,mon,date=event_lookup[identifier]
                linked[(loc,year,mon)].add(date)
                qc["original_occurrence_join_rows_at_107_sites"]+=1
    if qc["all_original_Event_rows"]!=15606 or qc["all_original_Occurrence_rows"]!=21617:
        raise ValueError("Changed original Danish source rows")
    return events,linked,qc

def validate_join_from_prior(events,linked,sites,root):
    qa=json.loads((root/LINK_QC).read_text(encoding="utf-8"))
    lookup={(r["source_rounded_site_hash"],r["original_year"]):r for r in qa["all_source_site_year_published_record_link_support"]}
    if len(lookup)!=321:raise ValueError("Original 2023/24/25 completeness QA missing")
    for loc,s in sites.items():
        for yr in YEARS:
            valid_all=set()
            valid_link=set()
            for m in range(1,13):
                valid_all.update(events[(loc,yr,m)])
                valid_link.update(linked[(loc,yr,m)])
            x=lookup[(s["site_id_SHA256"],yr)]
            if len(valid_all)!=x["original_Event_distinct_dates_n"] or len(valid_link)!=x["original_distinct_dates_with_at_least_one_Occurrence_n"]:
                raise ValueError("Independent Event/Occurrence linked survey days do not match previous source QA")
    return qa

def run(root):
    gate=json.loads((root/GATE).read_text(encoding="utf-8"))
    old=json.loads((root/2025FAMILY).read_text(encoding="utf-8"))
    if gate["2025_butterfly_occurrence_data_rows_opened"]!=0 or not gate["source_sample_gate_prespecified"]["2025_family_outcome_access_permitted"]:
        raise ValueError("Pre-outcome 2025 source-only gate not available")
    sites,families,previous,clusters,source_sha=monthly.yearsources(root)
    if len(sites)!=107 or len(set(clusters.values()))!=4:raise ValueError("Original source units/site frame changed")
    monthly.YEARS=YEARS
    raw=monthly.first.download()
    if sha(raw)!=monthly.first.ORIGINAL_SHA:raise ValueError("Licensed v1.4 ZIP changed")
    all_days,link_days,join_qc=month_linked_days(raw,sites)
    completeness=validate_join_from_prior(all_days,link_days,sites,root)
    _unfiltered,abund,_qc=monthly.counts_by_month(raw,sites,families)
    assert len(old["original_2025_new_site_year_accepted_aggregate_individualCounts"])==107
    prev2025={x["source_site_id_SHA256"]:x for x in old["original_2025_new_site_year_accepted_aggregate_individualCounts"]}
    for loc,site in sites.items():
        report=prev2025[site["site_id_SHA256"]]
        for fam in FAMILIES:
            n=sum(abund[(loc,2025,m)][fam] for m in range(1,13))
            if n!=report["family_counts"][fam]:
                raise ValueError("2025 source-locked family count sum changed")
    out={
        "schema":"gesui_DK_eBMS_2025_posthoc_same_month_any_occurrence_linked_day_Lycaenidae_v01",
        "date_utc":datetime.now(timezone.utc).isoformat(),
        "posthoc_analysis_after_original_2025_outoftime_pattern_known":True,
        "exact_original_DWCA_zip_SHA256":monthly.first.ORIGINAL_SHA,
        "pre2025_species_original_source_gate_SHA256":sha((root/GATE).read_bytes()),
        "previous_2025_individualCount_original_results_SHA256":sha((root/2025FAMILY).read_bytes()),
        "prior_independent_Event_Occ_link_QA_SHA256":sha((root/LINK_QC).read_bytes()),
        "original_107_sites_and_four_registered_source_units_verified":True,
        "original_2025_family_individualCount_reproduced_exactly":True,
        "Event_without_Occurrence_rows_not_claimed_true_zero":True,
        "linked_days_selection_may_condition_on_total_biota_reporting":True,
        "three_years_do_not_imply_WWTP_changed_treatment":True,
        "no_causal_or_confirmatory_WWTP_butterfly_effect":True,
        "original_source_join_QC":dict(join_qc),
        "configurations":{},
    }
    for label,(years,months_allowed) in CONFIG.items():
        rows=[]
        nsites=Counter()
        for loc,s in sites.items():
            keep=[m for m in months_allowed if all(link_days[(loc,y,m)] for y in years)]
            if not keep:continue
            band=source_band(s)
            nsites[band]+=1
            for y in years:
                vals={f:sum(abund[(loc,y,m)][f] for m in keep) for f in FAMILIES}
                nday=sum(len(link_days[(loc,y,m)]) for m in keep)
                e=sum(len(all_days[(loc,y,m)]) for m in keep)
                if nday<=0 or e<nday:raise ValueError("Matched month without actual linked observation")
                rows.append({
                    "site_id_SHA256":s["site_id_SHA256"],
                    "source_band":band,
                    "1deg_geo_region":s["country_1deg_lat_lon_cell"],
                    "nearest_registered_source_outfall_id_SHA256":clusters.get(loc),
                    "year":y,
                    "original_all_event_days_for_kept_months":e,
                    "occurrence_linked_days_for_kept_months":nday,
                    "family_counts":vals,
                })
        result={"matched_years":list(years),
                "matching_rule":"Same original public site, same calendar month in every selected year, at least one published ANY Occurrence on at least one survey date in each year; numerator original accepted families, denominator distinct dates with any Occurrence only.",
                "matched_site_counts_by_source_band":dict(nsites),
                "all_years":{}, "relative_2024_2025_change":{},
                "WWTP_near_cluster_count_with_reported_individuals":None,
                "outfall_specific_near_source_2025":[],
                "four_outfall_leave_one_out_2024_2025_ratio_sensitivity":[]}
        for g in ("near5","remote10","intermediate"):
            result["all_years"][g]={}
            for y in years:
                result["all_years"][g][str(y)]=aggregate([z for z in rows if z["source_band"]==g and z["year"]==y])
        if 2024 in years and 2025 in years:
            near=result["all_years"]["near5"]
            remote=result["all_years"]["remote10"]
            for fam in FAMILIES:
                nr=ratio(near["2025"]["family_per_any_original_Occurrence_linked_day"][fam],
                         near["2024"]["family_per_any_original_Occurrence_linked_day"][fam])
                fr=ratio(remote["2025"]["family_per_any_original_Occurrence_linked_day"][fam],
                         remote["2024"]["family_per_any_original_Occurrence_linked_day"][fam])
                result["relative_2024_2025_change"][fam]={
                    "near_rate_ratio":nr,"remote_rate_ratio":fr,
                    "near_over_remote_rate_ratio_ratio":ratio(nr,fr) if nr is not None and fr is not None else None}
            nearrows=[z for z in rows if z["source_band"]=="near5"]
            clusters_all=sorted(set(z["nearest_registered_source_outfall_id_SHA256"] for z in nearrows))
            for group in clusters_all:
                c=aggregate([z for z in nearrows if z["year"]==2025 and z["nearest_registered_source_outfall_id_SHA256"]==group])
                result["outfall_specific_near_source_2025"].append({
                    "outfall_id_SHA256":group,"2025_group":c,
                    "cannot_count_0_published_individuals_as_demonstrated_no_butterflies":True})
                rest=[z for z in nearrows if z["nearest_registered_source_outfall_id_SHA256"]!=group]
                if not rest:continue
                a=aggregate([z for z in rest if z["year"]==2024])
                b=aggregate([z for z in rest if z["year"]==2025])
                rr=ratio(b["Lycaenidae_count_per_any_original_Occurrence_linked_day"],
                         a["Lycaenidae_count_per_any_original_Occurrence_linked_day"])
                remote_rr=result["relative_2024_2025_change"]["Lycaenidae"]["remote_rate_ratio"]
                result["four_outfall_leave_one_out_2024_2025_ratio_sensitivity"].append({
                    "left_out_source_SHA256":group,
                    "remaining_near_to_remote_2025_2024_Lycaenidae_rate_ratio_ratio":
                        ratio(rr,remote_rr) if rr is not None else None,
                    "not_an_independently_replicated_significance_test":True})
            result["WWTP_near_cluster_count_with_reported_individuals"]=sum(x["2025_group"]["family_counts_individualCount_original"]["Lycaenidae"]>0 for x in result["outfall_specific_near_source_2025"]) if False else sum(x["2025_group"]["family_individualCount_original"]["Lycaenidae"]>0 for x in result["outfall_specific_near_source_2025"])
        out["configurations"][label]=result
    return out


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    x=run(root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(x,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:{"sites":v["matched_site_counts_by_source_band"],
                         "2024_25_Lycaenidae":v["relative_2024_2025_change"].get("Lycaenidae"),
                         "source_reporting_groups":v["WWTP_near_cluster_count_with_reported_individuals"]}
                      for k,v in x["configurations"].items()},indent=2))

if __name__=="__main__":
    main()
