#!/usr/bin/env python3
"""Original Danish eBMS 2025 NEW399 accepted-family counts BEFORE WWTP contrasts.

Immutable source-only 399-site/2025 published-Occurrence opportunity gate was
committed on this branch before opening any original family individual counts
at these NEW locations. All 107 earlier-outcome sites are excluded. Original
CC-BY-NC eBMS v1.4 Event/Occurrence join, taxonomy frozen by GBIF beforehand.
Exports only accepted family counts by hashed site and source-only cohorts.
"""
from __future__ import annotations

import argparse, hashlib, io, json, math, zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import gesui_DK_EBMS_V14_2023_2024_licensed_family_count_panel_after_taxonomy as previous

ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
GATE=ROOT+"DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json"
TAXONOMY=ROOT+"DENMARK_EBMS_V14_ORIGINAL_SCIENTIFICNAME_GBIF_FAMILY_BLIND_TAXONOMY.json"
ORIGINAL_SOURCE=ROOT+"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
YEAR=2025
FAMILIES=("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae","Riodinidae")

def sha(data):
    return hashlib.sha256(data).hexdigest()

def build(root):
    files={name:root/path for name,path in (
       ("frozen_preoutcome_NEW399_site_source",GATE),
       ("frozen_original_taxonomy_BLIND_to_WWTP",TAXONOMY),
       ("original_2022_EEA_and_eBMS_Event_source",ORIGINAL_SOURCE))}
    frozen={name:json.loads(path.read_text(encoding="utf-8")) for name,path in files.items()}
    pre=frozen["frozen_preoutcome_NEW399_site_source"]
    tax=frozen["frozen_original_taxonomy_BLIND_to_WWTP"]
    if not pre["frozen_before_2025_NEW399_butterfly_counts_gate"]["2025_NEW399_butterfly_family_count_outcome_permitted"]:
        raise ValueError("Original source-only NEW399 pre-family gate failed. STOP.")
    if not pre["NEW2025_Lycaenidae_or_any_butterfly_family_individualCount_NOT_READ"]:
        raise ValueError("Pre-outcome source manifest contaminated")
    if not pre["original_2025_family_taxon_or_scientificName_NOT_READ_in_this_stage"]:
        raise ValueError("Pre-outcome source stage had species data")
    if pre["original_denmark_eBMS_v14_ZIP_SHA256"]!=previous.ORIGINAL_SHA:
        raise ValueError("Underlying eBMS archive changed")
    if tax["original_eBMS_DwCA_SHA256"]!=previous.ORIGINAL_SHA or not tax["original_individualCount_values_NOT_READ"]:
        raise ValueError("Taxonomic family mapping not blindly frozen")
    if not tax["original_EEA_WWTP_near_or_remote_exposure_NOT_ACCESSED"]:
        raise ValueError("Taxonomy previously viewed source exposure")
    if pre["new_2025_year_sourced_station_n"]!=399 or pre["original_previous_2023_2024_butterfly_outcome_examined_site_count_EXCLUDED"]!=107:
        raise ValueError("Original new 2025 cohort not independent of old sites")
    if pre["original_2022_EEA_DK_fullFW_reporter_source_count"]!=203:
        raise ValueError("Original 2022 Danish source cohort changed")
    sourceSites=pre["all_original_399_source_only_new_site_2025_documented_visit_metadata"]
    if len(sourceSites)!=399 or len({x["site_id_SHA256"] for x in sourceSites})!=399:
        raise ValueError("Frozen 399 site IDs not unique")
    sites={x["site_id_SHA256"]:x for x in sourceSites}
    originalTaxa=tax["original_scientificName_GBIF_taxonomic_source_match_audit"]
    if len(originalTaxa)!=112:
        raise ValueError("Original accepted species registry length changed")
    names={x["original_scientificName"]:x for x in originalTaxa}
    if len(names)!=len(originalTaxa):
        raise ValueError("Original taxon string duplicated")
    mapping={k:v["gbif_family"] for k,v in names.items() if v["authoritative_family_accepted"]}
    if not set(mapping.values()).issubset(set(FAMILIES)):
        raise ValueError("Unexpected blindly accepted non-target taxonomy")
    raw=previous.download()
    if sha(raw)!=previous.ORIGINAL_SHA:
        raise ValueError("Original eBMS licensed source ZIP SHA mismatch")
    event_map={}
    seen=set()
    visitdays=defaultdict(set)
    eventrows=Counter()
    original_2025_399_occurrence_linked_event_ids=set()
    linked_days=defaultdict(set)
    accepted=defaultdict(Counter)
    diagnostics=Counter()
    statuses=Counter()
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        parsed=ET.fromstring(z.read("meta.xml"))
        eventCore=parsed.find(previous.NS+"core")
        extensions=[x for x in parsed.findall(previous.NS+"extension")
                   if x.attrib.get("rowType")=="http://rs.tdwg.org/dwc/terms/Occurrence"]
        if eventCore is None or len(extensions)!=1 or eventCore.attrib.get("rowType")!="http://rs.tdwg.org/dwc/terms/Event":
            raise ValueError("Original licensed DwCA schema not recognized")
        ef,efilename,es,eidField=previous.core_terms(eventCore)
        of,ofilename,os,oidField=previous.core_terms(extensions[0])
        eventID=eidField("id")
        occEventID=oidField("coreid")
        if not {"scientificName","individualCount","occurrenceStatus"}.issubset(of):
            raise ValueError("Original occurrence fields changed")
        if not {"eventDate","decimalLatitude","decimalLongitude"}.issubset(ef):
            raise ValueError("Original Event source fields changed")

        for row in previous.stream_csv(z,efilename,es):
            diagnostics["all_original_Event_rows"]+=1
            if len(row)<=max(eventID,ef["eventDate"],ef["decimalLatitude"],ef["decimalLongitude"]):
                raise ValueError("Truncated original Event row")
            eid=row[eventID].strip()
            if not eid or eid in seen:
                raise ValueError("Unusable original Event core ID")
            seen.add(eid)
            date=row[ef["eventDate"]].strip()[:10]
            if not date.startswith("2025-"):continue
            try:
                lat=float(row[ef["decimalLatitude"]])
                lon=float(row[ef["decimalLongitude"]])
                year=int(date[:4]);month=int(date[5:7]);day=int(date[8:10])
            except (ValueError,TypeError):continue
            if year!=YEAR or not(1<=month<=12 and 1<=day<=31 and math.isfinite(lat) and math.isfinite(lon)):continue
            loc=f"{lat:.5f}|{lon:.5f}"
            key=sha(loc.encode())
            if key not in sites:continue
            event_map[eid]=(key,date)
            visitdays[key].add(date)
            eventrows[key]+=1
            diagnostics["original_2025_NEW399_Event_rows"]+=1

        for row in previous.stream_csv(z,ofilename,os):
            diagnostics["all_original_Occurrence_rows"]+=1
            if len(row)<=max(occEventID,of["scientificName"],of["individualCount"],of["occurrenceStatus"]):
                raise ValueError("Truncated original Occurrence row")
            eid=row[occEventID].strip()
            if eid not in event_map:continue
            key,date=event_map[eid]
            diagnostics["original_2025_NEW399_Event_linked_Occurrence_rows"]+=1
            original_2025_399_occurrence_linked_event_ids.add(eid)
            linked_days[key].add(date)
            status=row[of["occurrenceStatus"]].strip().lower()
            statuses[status or "(blank)"]+=1
            if status not in ("", "present","presence"):
                diagnostics["explicit_nonpresent_or_ambiguous_status"]+=1
                continue
            count=previous.floatcount(row[of["individualCount"]])
            if count is None:
                diagnostics["missing_or_nonpositive_noninteger_individualCount"]+=1
                continue
            scientific=row[of["scientificName"]].strip()
            if scientific not in names:
                raise ValueError("2025 scientificName absent from original 112-name frozen taxonomy")
            fam=mapping.get(scientific)
            if fam in FAMILIES:
                accepted[key][fam]+=count
                diagnostics["accepted_family_positive_occurrence_rows"]+=1
            else:
                diagnostics["other_or_unresolved_species_occurrence_rows"]+=1

    if diagnostics["all_original_Event_rows"]!=15606 or diagnostics["all_original_Occurrence_rows"]!=21617:
        raise ValueError("Original full licensed dataset counts changed")
    if diagnostics["original_2025_NEW399_Event_rows"]!=pre["2025_EBMS_event_to_any_occurrence_ID_source_QC"]["original_2025_399site_event_rows"]:
        raise ValueError("Original source-only original 2025 event rows not reproduced")
    if diagnostics["original_2025_NEW399_Event_linked_Occurrence_rows"]!=pre["2025_EBMS_event_to_any_occurrence_ID_source_QC"]["original_2025_399site_linked_any_occurrence_rows"]:
        raise ValueError("Original source-only linked occurrence row count mismatch")
    for key,s in sites.items():
        if len(visitdays[key])!=s["original_event_date_count"] or eventrows[key]!=s["original_event_rows"]:
            raise ValueError("Original NEW399 site Event/Date QA failed")
        if len(linked_days[key])!=s["original_linked_occurrence_dates"]:
            raise ValueError("Pre-outcome source any-Occurrence linked days mismatch")
        if s["pre_2025_family_any_occurrence_linked_site_source_gate_pass"]!=(len(linked_days[key])>=2):
            raise ValueError("Original pre-outcome any-Occurrence site eligibility was not frozen correctly")
    # Full independent original 399 site-year panel: includes source-only sites with
    # zero published Occurrence (unknown completeness, NOT verified zero fauna).
    entries=[]
    for key,s in sorted(sites.items()):
        n={f:int(accepted[key][f]) for f in FAMILIES}
        total=n["Lycaenidae"]+n["Pieridae"]+n["Papilionidae"]
        entries.append({
            "site_SHA256":key,
            "source_original_2022_registered_outfall_band":s["source_band"],
            "original_country_1deg_region":s["original_1degree_region"],
            "nearest_original_2022_WWTP_outfall_SHA256":s["nearest_original_2022_EEA_outfall_id_SHA256"],
            "near_outfall_5km_spatial_cluster_SHA256":s.get("5km_original_WWTP_spatial_component_SHA256"),
            "near_outfall_10km_spatial_cluster_SHA256":s.get("10km_original_WWTP_spatial_component_SHA256"),
            "original_2025_source_distance_km":s["source_distance_km"],
            "original_2025_distinct_Event_dates":len(visitdays[key]),
            "original_2025_Event_rows":int(eventrows[key]),
            "original_2025_distinct_dates_with_any_published_Occurrence":len(linked_days[key]),
            "original_2025_2plus_any_Occurrence_date_source_gate":s["pre_2025_family_any_occurrence_linked_site_source_gate_pass"],
            "original_2025_family_individualCount":n,
            "original_2025_target_three_family_count":total,
            "original_2025_Lycaenidae_fraction_given_3family_count":
                (n["Lycaenidae"]/total) if total>0 else None,
            "original_2025_no_published_Occurrence_CANNOT_be_confirmed_no_butterflies":
                len(linked_days[key])==0,
        })
    if len(entries)!=399 or len({r["site_SHA256"] for r in entries})!=399:
        raise ValueError("New independent 399 locations changed")
    by_class={}
    for cls in ("near5","remote10","intermediate"):
        vals=[r for r in entries if r["source_original_2022_registered_outfall_band"]==cls]
        supported=[r for r in vals if r["original_2025_2plus_any_Occurrence_date_source_gate"]]
        def sample(z):
            return {
                "sites":len(z),
                "total_original_Event_dates":sum(r["original_2025_distinct_Event_dates"] for r in z),
                "total_any_Occurrence_linked_dates":sum(r["original_2025_distinct_dates_with_any_published_Occurrence"] for r in z),
                "sites_with_target_threefamily_positive":sum(r["original_2025_target_three_family_count"]>0 for r in z),
                "sites_with_zero_published_Occurrence":sum(r["original_2025_no_published_Occurrence_CANNOT_be_confirmed_no_butterflies"] for r in z),
                "family_recorded_individualCount":{f:sum(r["original_2025_family_individualCount"][f] for r in z) for f in FAMILIES},
            }
        by_class[cls]={"all_original_site_rows":sample(vals),"eligible_prior_source_reporting_sites":sample(supported)}
        if cls=="near5" and len(supported)!=70:raise ValueError("Pre-frozen source-near sample not 70")
        if cls=="remote10" and len(supported)!=92:raise ValueError("Pre-frozen source-far sample not 92")
        if cls=="intermediate" and len(supported)!=124:raise ValueError("Pre-frozen source-intermediate sample not 124")
    model_quality={
        "source_prospectively_frozen_documented_new_sites_qualify":True,
        "near5_prior_any_occurrence_supported_sites":70,
        "remote10_prior_any_occurrence_supported_sites":92,
        "min_near_source_outfall_5km_spatial_clusters_with_target_butterflies":10,
        "min_near_2plus_documented_sites":40,
        "min_remote_2plus_documented_sites":60,
        "min_shared_1deg_near_remote_cells_with_any_threefamily_butterflies":3,
        "independent_WWTP_source_units_are_unique_reported_outfall_IDs_not_causal_treatment_systems":True,
    }
    supported=[r for r in entries if r["original_2025_2plus_any_Occurrence_date_source_gate"]]
    near=[r for r in supported if r["source_original_2022_registered_outfall_band"]=="near5"]
    far=[r for r in supported if r["source_original_2022_registered_outfall_band"]=="remote10"]
    source_clusters={r["near_outfall_5km_spatial_cluster_SHA256"]
                     for r in near if r["original_2025_target_three_family_count"]>0}
    near_outfall_ids={r["nearest_original_2022_WWTP_outfall_SHA256"]
                     for r in near if r["original_2025_target_three_family_count"]>0}
    cells={}
    for r in supported:
        cell=r["original_country_1deg_region"]
        g=r["source_original_2022_registered_outfall_band"]
        cells.setdefault(cell,{}).setdefault(g,[]).append(r)
    shared_positive=[cell for cell,by in cells.items()
         if "near5" in by and "remote10" in by
         and sum(x["original_2025_target_three_family_count"] for x in by["near5"])>0
         and sum(x["original_2025_target_three_family_count"] for x in by["remote10"])>0]
    model_quality["observed_precontrast_family_recording_supported_near5_outfalls"]=len(near_outfall_ids)
    model_quality["observed_near_source_5km_spatial_clusters_with_target_butterflies"]=len(source_clusters)
    model_quality["observed_shared_1deg_cells_with_target_threefamily_near_and_far"]=len(shared_positive)
    model_quality["post_family_data_quality_gate_supports_descriptive_comparison"]=bool(
        len(near)>=40 and len(far)>=60 and len(source_clusters)>=10 and len(shared_positive)>=3)
    output={
        "schema":"gesui_2025_DK_eBMS_NEW399_source_locked_taxonomy_original_family_count_before_distance_contrasts_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_source_SHA256":{k:sha(p.read_bytes()) for k,p in files.items()},
        "original_eBMS_DwCA_v1_4_SHA256":previous.ORIGINAL_SHA,
        "original_eBMS_DwCA_license":"CC-BY-NC 4.0",
        "original_GBIF_dataset":"caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c",
        "all_prior_107_2023_24_2025_outcome_sites_EXCLUDED":True,
        "2025_original_Event_to_Occurrence_and_taxon_count_QC":dict(diagnostics),
        "original_2025_occurrence_status_counts_at_NEW399":dict(statuses),
        "source_stage_original_Fw_outfalls_2022_reported_not_actual_exposure":True,
        "original_2025_NEW399_individualCount_source_site_records":entries,
        "pre_effect_original_source_year_region_family_data_quality_gate":model_quality,
        "original_2025_source_band_family_count_coverage_before_effect":by_class,
        "no_near_vs_remote_butterfly_family_effect_or_sign_computed_in_this_stage":True,
        "no_butterfly_causal_wastewater_effect_claim":True,
        "raw_original_personal_observerIDs_and_occurrence_rows_NOT_EXPORTED":True,
    }
    return output

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=build(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "source_frame_NEW399_exact":len(out["original_2025_NEW399_individualCount_source_site_records"]),
        "source_original_taxonomy_observations_QC":out["2025_original_Event_to_Occurrence_and_taxon_count_QC"],
        "source_family_quality_gate":out["pre_effect_original_source_year_region_family_data_quality_gate"],
        "family_photo_count_source_by_band":out["original_2025_source_band_family_count_coverage_before_effect"],
    },indent=2))

if __name__=="__main__":
    main()
