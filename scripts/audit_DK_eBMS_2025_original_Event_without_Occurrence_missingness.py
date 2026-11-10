#!/usr/bin/env python3
"""POSTHOC publication-completeness audit of Denmark original 2023–2025 eBMS.

Crucial distinction: original Event rows are not guaranteed to have a linked
Occurrence row. Event with no Occurrence cannot automatically be called
true-zero butterfly abundance. Uses ONLY original Event ID/Occurrence core ID
joins (not taxon names, family effect or individualCount) for the support gate.
"""
from __future__ import annotations

import argparse,hashlib,io,json,math,zipfile
import xml.etree.ElementTree as ET
from collections import Counter,defaultdict
from pathlib import Path

import gesui_DK_EBMS_V14_2023_2024_licensed_family_count_panel_after_taxonomy as first

ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
SOURCE=ROOT+"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
PANEL=ROOT+"DENMARK_2025_ORIGINAL107_LYCAENIDAE_OUTOFTIME_FAMILY_AND_MONTH_RESULTS.json"
YEARS=(2023,2024,2025)

def sha(b):return hashlib.sha256(b).hexdigest()

def analyze(root):
    src=json.loads((root/SOURCE).read_text(encoding="utf-8"))
    family=json.loads((root/PANEL).read_text(encoding="utf-8"))
    assert family["licensed_official_eBMS_v14_original_ZIP_SHA256"]==first.ORIGINAL_SHA
    sites={f"{s['lat']:.5f}|{s['lon']:.5f}":s for s in src["original_2023_2024_repeated_event_site_source_metadata"]
           if s["repeat_2023_2024_dates_gate"] and s["precise_500m_explicit_uncertainty_gate"] and s["rounded5_site_coord_scatter_le2m"]}
    assert len(sites)==107
    orig=first.download()
    if sha(orig)!=first.ORIGINAL_SHA:raise ValueError("Original licensed eBMS v1.4 source changed")
    events={}
    eventid_seen=set()
    siteyear=defaultdict(lambda:{"events":set(),"days":set(),"linked_events":set(),"linked_days":set(),"occurrence_rows":0})
    counts=Counter()
    with zipfile.ZipFile(io.BytesIO(orig)) as zipf:
        root_xml=ET.fromstring(zipf.read("meta.xml"))
        ec=root_xml.find(first.NS+"core")
        extensions=[z for z in root_xml.findall(first.NS+"extension") if z.attrib.get("rowType")=="http://rs.tdwg.org/dwc/terms/Occurrence"]
        if ec is None or len(extensions)!=1:raise ValueError("DwCA Event + Occurrence extension expected")
        ef,efile,es,efield=first.core_terms(ec)
        of,ofile,os,ofield=first.core_terms(extensions[0])
        ek=efield("id");ok=ofield("coreid")
        for row in first.stream_csv(zipf,efile,es):
            counts["full_original_Event_core_rows"]+=1
            if len(row)<=max(ek,ef["eventDate"],ef["decimalLatitude"],ef["decimalLongitude"]):
                raise ValueError("Truncated original Event row")
            key=row[ek].strip()
            if not key or key in eventid_seen:raise ValueError("Duplicate original Event IDs")
            eventid_seen.add(key)
            dt=row[ef["eventDate"]].strip()[:10]
            if len(dt)!=10 or not dt[:4].isdigit():continue
            year=int(dt[:4])
            if year not in YEARS:continue
            try:lat=float(row[ef["decimalLatitude"]]);lon=float(row[ef["decimalLongitude"]])
            except (ValueError,TypeError):continue
            if not math.isfinite(lat) or not math.isfinite(lon):continue
            loc=f"{lat:.5f}|{lon:.5f}"
            if loc not in sites:continue
            events[key]=(loc,year,dt)
            group=(loc,year)
            siteyear[group]["events"].add(key)
            siteyear[group]["days"].add(dt)
            counts["eligible_original_Event_core_rows"]+=1
        for row in first.stream_csv(zipf,ofile,os):
            counts["full_original_Occurrence_extension_rows"]+=1
            if len(row)<=ok:raise ValueError("Truncated original Occurrence row")
            key=row[ok].strip()
            if key not in events:continue
            loc,yr,date=events[key]
            siteyear[(loc,yr)]["linked_events"].add(key)
            siteyear[(loc,yr)]["linked_days"].add(date)
            siteyear[(loc,yr)]["occurrence_rows"]+=1
            counts["eligible_original_occurrence_join_rows"]+=1
    assert counts["full_original_Event_core_rows"]==15606
    assert counts["full_original_Occurrence_extension_rows"]==21617
    assert len(siteyear)==107*3, len(siteyear)
    details=[]
    for loc,site in sites.items():
        d=site["nearest_qualified_EEA2022_DK_FW_outfall_km"]
        c="near5" if d<=5 else "remote10" if d>=10 else "intermediate"
        for year in YEARS:
            x=siteyear[(loc,year)]
            expected=site[f"visit_days_{year}" if year<2025 else "visit_days_2025_descriptive"]
            if len(x["days"])!=expected:
                raise ValueError("Original visit days do not match previously frozen source metadata")
            details.append({
                "source_rounded_site_hash":site["site_id_SHA256"],
                "source_band":c,"region_1degree":site["country_1deg_lat_lon_cell"],
                "original_year":year,
                "original_Events_n":len(x["events"]),
                "original_Event_distinct_dates_n":len(x["days"]),
                "original_Events_with_at_least_one_Occurrence_n":len(x["linked_events"]),
                "original_distinct_dates_with_at_least_one_Occurrence_n":len(x["linked_days"]),
                "original_Events_without_any_published_Occurrence_n":len(x["events"]-x["linked_events"]),
                "original_Occurrence_extension_rows_linked":x["occurrence_rows"],
                "contains_original_Event_without_occurrence_not_treated_as_confirmed_zero":
                    len(x["events"]-x["linked_events"])>0,
            })
    result={
        "schema":"gesui_DK_2023_2025_107site_original_Event_Occurrence_publication_completeness_posthoc_QA_v01",
        "posthoc_2025_quality_inspection_AFTER_2025_count_pattern_seen":True,
        "original_full_licensed_DK_eBMS_v14_ZIP_SHA256":first.ORIGINAL_SHA,
        "original_2025_family_count_panel_SHA256":sha((root/PANEL).read_bytes()),
        "original_DK_source_and_sourceyear_frame_SHA256":sha((root/SOURCE).read_bytes()),
        "only_event_and_occurrence_ID_join_not_family_counts_read_in_this_stage":True,
        "missing_Occurrence_extension_for_Event_cannot_identify_true_0_butterflies":True,
        "n_sites":len(sites),"year_scope":list(YEARS),
        "source_original_row_qc":dict(counts),
        "all_source_site_year_published_record_link_support":details,
        "year_and_exposure_group_publication_quality":{},
    }
    for year in YEARS:
        for group in ("near5","remote10","intermediate"):
            a=[x for x in details if x["original_year"]==year and x["source_band"]==group]
            result["year_and_exposure_group_publication_quality"][f"{year}_{group}"]={
                "original_eligible_sites":len(a),
                "siteyears_WITH_ZERO_published_occurrence_rows":sum(x["original_Occurrence_extension_rows_linked"]==0 for x in a),
                "original_Events_n":sum(x["original_Events_n"] for x in a),
                "original_Events_with_any_published_Occurrence_n":sum(x["original_Events_with_at_least_one_Occurrence_n"] for x in a),
                "original_occurrence_rows_published":sum(x["original_Occurrence_extension_rows_linked"] for x in a),
                "original_all_site_event_distinct_dates":sum(x["original_Event_distinct_dates_n"] for x in a),
                "original_linked_occurrence_site_dates_n":sum(x["original_distinct_dates_with_at_least_one_Occurrence_n"] for x in a),
                "siteyears_with_2plus_appearance_days":sum(x["original_distinct_dates_with_at_least_one_Occurrence_n"]>=2 for x in a),
                "fraction_of_Event_rows_having_Any_Occurrence":sum(x["original_Events_with_at_least_one_Occurrence_n"] for x in a)/sum(x["original_Events_n"] for x in a),
            }
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    v=analyze(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(v,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"sourceRows":v["source_original_row_qc"],"yearGroup":v["year_and_exposure_group_publication_quality"]},indent=2))
if __name__=="__main__":
    main()
