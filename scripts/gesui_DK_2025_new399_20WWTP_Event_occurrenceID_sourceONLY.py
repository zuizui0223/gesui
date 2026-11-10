#!/usr/bin/env python3
"""Freeze 399 NEW Denmark 2025 eBMS sites + original 2022 EEA WWTP sources.

SOURCE/effort ONLY. Reuses original independently captured DwC-A v1.4 and
pre-existing 735 source-only site metadata. From original Occurrence extension
read ONLY the event core-ID linkage, never scientificName or individualCount.
No 2025 Lycaenidae family counts are opened by this script.
"""
from __future__ import annotations

import argparse,hashlib,io,json,math,zipfile
import xml.etree.ElementTree as ET
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import gesui_DK_EBMS_V14_2023_2024_licensed_family_count_panel_after_taxonomy as prior

ORIGINAL_SOURCE="results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
EXPECTED_ZIP="2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b"
SOURCE_NEAR="near5"
SOURCE_REMOTE="remote10"
PREVIOUS="original107"
YEAR=2025


def sha(content):
    return hashlib.sha256(content).hexdigest()


def h(a,b):
    r=math.pi/180
    dy=(a["lat"]-b["lat"])*r;dx=(a["lon"]-b["lon"])*r
    q=math.sin(dy/2)**2+math.cos(a["lat"]*r)*math.cos(b["lat"]*r)*math.sin(dx/2)**2
    return 12742*math.asin(min(1.0,math.sqrt(q)))


def classification(v):
    return "near5" if v<=5 else "remote10" if v>=10 else "intermediate"


def component_map(source_points,dist):
    keys=sorted(source_points)
    parent={key:key for key in keys}
    def find(key):
        while key!=parent[key]:
            parent[key]=parent[parent[key]]
            key=parent[key]
        return key
    for i,ki in enumerate(keys):
        for kj in keys[i+1:]:
            if h(source_points[ki],source_points[kj])<=dist:
                parent[find(ki)]=find(kj)
    groups=defaultdict(list)
    for key in keys:groups[find(key)].append(key)
    result={}
    for ids in groups.values():
        component=sha(("WWTP_spatial_component_"+str(dist)+"|"+"|".join(sorted(ids))).encode())
        for key in ids:result[key]=component
    return result,len(groups),sorted([len(g) for g in groups.values()],reverse=True)


def run(root):
    path=root/ORIGINAL_SOURCE
    data=json.loads(path.read_text())
    original=data["original_2023_2024_repeated_event_site_source_metadata"]
    outlets=data["original_2022_EEA_DK_qualifying_Fw_outfalls"]
    if data["butterfly_species_2023_2024_occurrence_extension_READ"]:
        raise ValueError("Only Event-only original site source qualifies for new spatial frame")
    if len(outlets)!=203 or len(original)!=735:
        raise ValueError("Official source count changed. No butterfly outcomes permitted.")
    oldids={x["site_id_SHA256"] for x in original
        if x["repeat_2023_2024_dates_gate"]
        and x["precise_500m_explicit_uncertainty_gate"]
        and x["rounded5_site_coord_scatter_le2m"]}
    if len(oldids)!=107:raise ValueError("Original 107 benchmark sites not frozen")
    base=[x for x in original if x["visit_days_2025_descriptive"]>=2
        and x["precise_500m_explicit_uncertainty_gate"]
        and x["rounded5_site_coord_scatter_le2m"]]
    if len(base)!=506:raise ValueError("2025 original public Event source cohort changed")
    selected=[x for x in base if x["site_id_SHA256"] not in oldids]
    if len(selected)!=399:raise ValueError("Expected 399 new public Event sites before outcomes")
    sources={x["outfall_id_SHA256"]:x for x in outlets}
    stations={}
    for s in selected:
        key=s["site_id_SHA256"]
        if key in stations:raise ValueError("Repeated site original hash")
        best=min((h(s,out),out["outfall_id_SHA256"]) for out in outlets)
        if abs(best[0]-s["nearest_qualified_EEA2022_DK_FW_outfall_km"])>1e-7:
            raise ValueError("Original source/outfall geometry changed")
        stations[key]={
            "site_id_SHA256":key,"original_1degree_region":s["country_1deg_lat_lon_cell"],
            "original_2025_Event_dates_expected":int(s["visit_days_2025_descriptive"]),
            "source_band":classification(best[0]),
            "nearest_original_2022_EEA_outfall_id_SHA256":best[1],
            "source_distance_km":best[0],
            "original_event_date_count":0,"original_event_rows":0,
            "original_linked_occurrence_dates":0,
            "original_linked_event_rows":0,
            "2025_event_date_counts_by_month":{},
            "2025_occurrence_linked_date_counts_by_month":{},
        }
    assert Counter(v["source_band"] for v in stations.values())=={
        "near5":100,"remote10":127,"intermediate":172}
    nearkeys={x["nearest_original_2022_EEA_outfall_id_SHA256"]
              for x in stations.values() if x["source_band"]=="near5"}
    if len(nearkeys)!=20:raise ValueError("Original near outfall-source support changed")
    sub={key:sources[key] for key in nearkeys}
    cmp5,n5,nsize5=component_map(sub,5)
    cmp10,n10,nsize10=component_map(sub,10)
    if n5!=19 or n10!=16:
        raise ValueError("Unexpected original source-to-source geography 5/10km")
    for s in stations.values():
        if s["source_band"]=="near5":
            s["5km_original_WWTP_spatial_component_SHA256"]=cmp5[s["nearest_original_2022_EEA_outfall_id_SHA256"]]
            s["10km_original_WWTP_spatial_component_SHA256"]=cmp10[s["nearest_original_2022_EEA_outfall_id_SHA256"]]

    # Download original licensed archive, no occurrence taxa nor individualCount.
    raw=prior.download()
    if sha(raw)!=EXPECTED_ZIP:raise ValueError("DwCA original source byte SHA changed")
    meta={}
    event_index={}
    site_dates=defaultdict(set)
    eventrows=Counter()
    linked_dates=defaultdict(set)
    linked_eventids=set()
    stat=Counter()
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        meta["zip_original_file_members"]=z.namelist()
        rootxml=ET.fromstring(z.read("meta.xml"))
        core=rootxml.find(prior.NS+"core")
        ext=[e for e in rootxml.findall(prior.NS+"extension")
             if e.attrib.get("rowType")=="http://rs.tdwg.org/dwc/terms/Occurrence"]
        if core is None or len(ext)!=1:raise ValueError("Expected Event/Occurrence original DwCA")
        ef,efilename,es,efield=prior.core_terms(core)
        of,ofilename,os,ofield=prior.core_terms(ext[0])
        id_index=efield("id")
        occ_id=ofield("coreid")
        if "eventDate" not in ef or "decimalLatitude" not in ef or "decimalLongitude" not in ef:
            raise ValueError("Original Event terms unavailable")
        seen_ids=set()
        for row in prior.stream_csv(z,efilename,es):
            stat["original_full_event_rows"]+=1
            if len(row)<=max(id_index,ef["eventDate"],ef["decimalLatitude"],ef["decimalLongitude"]):
                raise ValueError("Original Event source short row")
            eid=row[id_index].strip()
            if not eid or eid in seen_ids:raise ValueError("Original Event IDs not uniquely usable")
            seen_ids.add(eid)
            date=row[ef["eventDate"]].strip()[:10]
            try:
                yr=int(date[:4]);month=int(date[5:7]);day=int(date[8:10])
                lat=float(row[ef["decimalLatitude"]]);lon=float(row[ef["decimalLongitude"]])
            except (TypeError,ValueError):continue
            if yr!=YEAR or month<1 or month>12 or day<1 or day>31 or not(math.isfinite(lat) and math.isfinite(lon)):continue
            loc=f"{lat:.5f}|{lon:.5f}"
            sid=sha(loc.encode())
            if sid not in stations:continue
            site_dates[sid].add(date)
            eventrows[sid]+=1
            event_index[eid]=(sid,date,month)
            stat["original_2025_399site_event_rows"]+=1
        for row in prior.stream_csv(z,ofilename,os):
            stat["original_full_occurrence_rows"]+=1
            if len(row)<=occ_id:raise ValueError("Original Occurrence source short row")
            eid=row[occ_id].strip()
            if eid not in event_index:continue
            sid,date,month=event_index[eid]
            linked_dates[sid].add(date)
            linked_eventids.add(eid)
            stat["original_2025_399site_linked_any_occurrence_rows"]+=1
    assert stat["original_full_event_rows"]==15606
    assert stat["original_full_occurrence_rows"]==21617
    for sid,s in stations.items():
        s["original_event_date_count"]=len(site_dates[sid])
        s["original_event_rows"]=eventrows[sid]
        s["original_linked_occurrence_dates"]=len(linked_dates[sid])
        s["original_linked_event_rows"]=sum(k[0]==sid for key,k in event_index.items() if key in linked_eventids)
        if s["original_event_date_count"]!=s["original_2025_Event_dates_expected"]:
            raise ValueError("Original 2025 Event-only site reported days mismatch source")
        for month in range(1,13):
            a=sum(date.startswith(f"2025-{month:02d}-") for date in site_dates[sid])
            b=sum(date.startswith(f"2025-{month:02d}-") for date in linked_dates[sid])
            if a:s["2025_event_date_counts_by_month"][str(month)]=a
            if b:s["2025_occurrence_linked_date_counts_by_month"][str(month)]=b
        s["pre_2025_family_any_occurrence_linked_site_source_gate_pass"]=s["original_linked_occurrence_dates"]>=2
    eligible=[s for s in stations.values() if s["pre_2025_family_any_occurrence_linked_site_source_gate_pass"]]
    sizes=Counter(s["source_band"] for s in eligible)
    near=[s for s in eligible if s["source_band"]=="near5"]
    nearbyid=Counter(x["nearest_original_2022_EEA_outfall_id_SHA256"] for x in near)
    nearbycomponent=Counter(x["5km_original_WWTP_spatial_component_SHA256"] for x in near)
    geo=defaultdict(lambda:Counter())
    for s in eligible:
        geo[s["original_1degree_region"]][s["source_band"]]+=1
    shared=[k for k,v in geo.items() if v["near5"]>0 and v["remote10"]>0]
    near_outfall_with2=sum(n>=2 for n in nearbyid.values())
    passed=(sizes["near5"]>=40 and sizes["remote10"]>=60
            and near_outfall_with2>=8 and len(nearbycomponent)>=6 and len(shared)>=4)
    output={
        "schema":"gesui_DK_eBMS_2025_NEW399_original_event_and_any_occurrence_link_source_only_v01",
        "source_freeze_utc":datetime.now(timezone.utc).isoformat(),
        "NEW2025_Lycaenidae_or_any_butterfly_family_individualCount_NOT_READ":True,
        "original_2025_family_taxon_or_scientificName_NOT_READ_in_this_stage":True,
        "original_denmark_eBMS_v14_ZIP_SHA256":EXPECTED_ZIP,
        "original_2023_24_Event_source_only_735site_SHA256":sha(path.read_bytes()),
        "original_2025_399_new_site_source_all_original_days_both_prior_family_panels_excluded":True,
        "original_2022_EEA_DK_fullFW_reporter_source_count":len(outlets),
        "original_public_site_Event_source_year_2025_total_all_qualifying":506,
        "original_previous_2023_2024_butterfly_outcome_examined_site_count_EXCLUDED":len(oldids),
        "new_2025_year_sourced_station_n":len(stations),
        "new_source_band_counts_prior_to_individual_butterfly_outcomes":dict(Counter(s["source_band"] for s in stations.values())),
        "2025_EBMS_event_to_any_occurrence_ID_source_QC":dict(stat),
        "source_near_registered_FW_outfall_ID_count_before_butterfly":len(nearkeys),
        "source_registered_FW_outfall_components_5km":n5,
        "source_registered_FW_outfall_components_10km":n10,
        "source_component_member_sizes":{"5km":nsize5,"10km":nsize10},
        "source_any_occurrence_date_gate_minimum":2,
        "source_eligible_new_sites_with_2plus_any_published_occurrence_dates_by_band":dict(sizes),
        "source_near_outfall_IDs_with_2plus_eligible_site_locations":near_outfall_with2,
        "source_near_5km_spatial_components_with_eligible_sites":len(nearbycomponent),
        "country_1degree_regions_with_both_source_bands":len(shared),
        "shared_source_regions":sorted(shared),
        "country_1deg_source_band_counts":{key:dict(value) for key,value in sorted(geo.items())},
        "frozen_before_2025_NEW399_butterfly_counts_gate":{
            "min_source_near_sites_with_2plus_ANY_published_occurrence_dates":40,
            "min_remote_sites_with_2plus_ANY_published_occurrence_dates":60,
            "min_distinct_near_registered_outfall_IDs_with_2plus_sites":8,
            "min_source_near_outfall_connected_5km_components":6,
            "min_1degree_regions_with_both_near_and_remote_eligible_sites":4,
            "2025_NEW399_butterfly_family_count_outcome_permitted":bool(passed)
        },
        "all_original_399_source_only_new_site_2025_documented_visit_metadata":sorted(stations.values(),key=lambda s:s["site_id_SHA256"]),
        "source_outfall_ID_to_5km_spatial_component":cmp5,
        "source_outfall_ID_to_10km_spatial_component":cmp10,
        "only_original_Event_and_Occurrence_core_id_fields_read":True,
        "no_2025_butterfly_family_response_accessed_for_new_sites":True,
        "any_occurrence_link_date_may_select_on_other_biota_reporting":True,
        "WWTP_source_proximity_not_hydrological_downstream_butterfly_exposure":True,
        "country_src_Fw_reported_2022_not_dose_year_2025":True,
    }
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    o=run(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"new_sites":o["new_2025_year_sourced_station_n"],
         "original_near_outfalls":o["source_near_registered_FW_outfall_ID_count_before_butterfly"],
         "old_site_ids_excluded":o["original_previous_2023_2024_butterfly_outcome_examined_site_count_EXCLUDED"],
         "documented_new_sites_by_band":o["source_eligible_new_sites_with_2plus_any_published_occurrence_dates_by_band"],
         "outfall_systems_with_2plus_sites":o["source_near_outfall_IDs_with_2plus_eligible_site_locations"],
         "shared_geo_regions":o["country_1degree_regions_with_both_source_bands"],
         "can_open_new_2025_butterfly_outcomes":o["frozen_before_2025_NEW399_butterfly_counts_gate"]["2025_NEW399_butterfly_family_count_outcome_permitted"]},indent=2))


if __name__=="__main__":
    main()
