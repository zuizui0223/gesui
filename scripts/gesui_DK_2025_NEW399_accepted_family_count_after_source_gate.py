#!/usr/bin/env python3
"""POSTHOC 2025 Denmark independent NEW399 eBMS sites x 2022 registered WWTP.

Read butterfly families ONLY AFTER new399 2025 Event -> ANY Occurrence source
opportunity was source-hash-frozen. New399 excludes the 107 previously
examined 2023/24 repeat sites. Family taxonomy frozen before all WWTP contrasts.
No real recipient WWTP dose, no randomized treatment, no true population census.
"""
from __future__ import annotations

import argparse,hashlib,io,json,math,zipfile
import xml.etree.ElementTree as ET
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import gesui_DK_EBMS_V14_2023_2024_licensed_family_count_panel_after_taxonomy as first

ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
MANIFEST=ROOT+"DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json"
TAX=ROOT+"DENMARK_EBMS_V14_ORIGINAL_SCIENTIFICNAME_GBIF_FAMILY_BLIND_TAXONOMY.json"
FAMILIES=("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae","Riodinidae")
TARGET=("Lycaenidae","Pieridae","Papilionidae")
YEAR=2025


def sha(v):return hashlib.sha256(v).hexdigest()


def stats(items):
    n=len(items)
    days=sum(x["original_any_Occurrence_linked_distinct_date_count"] for x in items)
    originals=sum(x["original_all_Event_distinct_date_count"] for x in items)
    counts={f:sum(x["family_individual_count_original"][f] for x in items) for f in FAMILIES}
    t=sum(counts[f] for f in TARGET)
    l=counts["Lycaenidae"]
    return {
        "qualified_site_coordinate_rows":n,
        "original_any_Occurrence_linked_distinct_survey_days":days,
        "original_all_Event_distinct_survey_days":originals,
        "accepted_family_original_individual_counts":counts,
        "original_3family_positive_recorded_sites":sum(sum(x["family_individual_count_original"][f] for f in TARGET)>0 for x in items),
        "original_Lycaenidae_positive_recorded_sites":sum(x["family_individual_count_original"]["Lycaenidae"]>0 for x in items),
        "Lycaenidae_share_of_three_target_families_if_positive_total":l/t if t else None,
        "Lycaenidae_per_any_published_Occurrence_linked_date":l/days if days else None,
        "accepted_families_per_any_published_Occurrence_linked_date":{
            f:counts[f]/days if days else None for f in FAMILIES},
        "threefamily_denominator_is_not_filled_if_all_target_absent":True,
    }


def family_mapping(raw):
    tax=json.loads(raw)
    if tax["original_eBMS_DwCA_SHA256"]!=first.ORIGINAL_SHA or not tax["original_individualCount_values_NOT_READ"]:
        raise ValueError("Previous accepted GBIF taxonomy not frozen before 2025 family counts")
    records=tax["original_scientificName_GBIF_taxonomic_source_match_audit"]
    if len(records)!=112 or len({r["original_scientificName"] for r in records})!=112:
        raise ValueError("Original blind taxonomy source altered")
    return {x["original_scientificName"]:x["gbif_family"] for x in records
            if x["authoritative_family_accepted"] and x.get("gbif_family") in FAMILIES}


def process_archive(raw,selected,tax):
    eligible={x["site_id_SHA256"]:x for x in selected if x["source_gate_documented_2plus_any_occurrences"]}
    all_candidates={x["site_id_SHA256"]:x for x in selected}
    bysite=defaultdict(Counter)
    seen_events=set()
    eventmap={}
    everdays=defaultdict(set)
    raw_eventrows=Counter()
    occ_eventids=set()
    linked=defaultdict(set)
    qc=Counter()
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        source=ET.fromstring(z.read("meta.xml"))
        e=source.find(first.NS+"core")
        extensions=[x for x in source.findall(first.NS+"extension")
             if x.attrib.get("rowType")=="http://rs.tdwg.org/dwc/terms/Occurrence"]
        if e is None or len(extensions)!=1:raise ValueError("Not exactly one original Event core / Occurrence extension")
        ef,efile,es,efield=first.core_terms(e)
        of,ofile,os,ofield=first.core_terms(extensions[0])
        ei=efield("id");oi=ofield("coreid")
        for f in ("scientificName","individualCount","occurrenceStatus"):
            if f not in of:raise ValueError("Original licensed accepted butterfly source field missing")
        for row in first.stream_csv(z,efile,es):
            qc["full_original_Event_rows"]+=1
            if len(row)<=max(ei,ef["eventDate"],ef["decimalLatitude"],ef["decimalLongitude"]):
                raise ValueError("Original Event table truncated")
            eid=row[ei].strip()
            if not eid or eid in seen_events:raise ValueError("Original EventID collision")
            seen_events.add(eid)
            date=row[ef["eventDate"]].strip()[:10]
            try:
                yr=int(date[:4]);month=int(date[5:7]);day=int(date[8:10])
                lat=float(row[ef["decimalLatitude"]]);lon=float(row[ef["decimalLongitude"]])
            except (TypeError,ValueError):continue
            if yr!=YEAR or not(1<=month<=12 and 1<=day<=31) or not(math.isfinite(lat) and math.isfinite(lon)):continue
            key=sha(f"{lat:.5f}|{lon:.5f}".encode())
            if key not in all_candidates:continue
            eventmap[eid]=(key,date,month)
            everdays[key].add(date)
            raw_eventrows[key]+=1
            qc["original_new399_2025_Event_rows"]+=1
        for row in first.stream_csv(z,ofile,os):
            qc["full_original_Occurrence_rows"]+=1
            if len(row)<=max(oi,of["scientificName"],of["individualCount"],of["occurrenceStatus"]):
                raise ValueError("Original Occurrence table truncated")
            eid=row[oi].strip()
            if eid not in eventmap:continue
            key,date,month=eventmap[eid]
            occ_eventids.add(eid)
            linked[key].add(date)
            qc["original_new399_any_taxa_published_occurrence_rows"]+=1
            status=row[of["occurrenceStatus"]].strip().lower()
            if status and status not in ("present","presence"):
                qc["occurrence_not_present_or_unexpected_status"]+=1
                continue
            n=first.floatcount(row[of["individualCount"]])
            if n is None:
                qc["invalid_individualCount_not_substituted_zero"]+=1
                continue
            name=row[of["scientificName"]].strip()
            if name not in tax:
                qc["non_accepted_or_non_butterfly_family_species_rows"]+=1
                continue
            bysite[key][tax[name]]+=n
            qc["source_accepted_frozen_butterfly_family_occurrence_rows"]+=1
    if qc["full_original_Event_rows"]!=15606 or qc["full_original_Occurrence_rows"]!=21617:
        raise ValueError("eBMS licensed ZIP original source no longer byte/row coherent")
    original=[]
    for key,s in sorted(all_candidates.items()):
        if len(everdays[key])!=s["original_event_date_count"]:
            raise ValueError("Source-locked original 2025 Event visit days changed")
        if len(linked[key])!=s["original_linked_occurrence_dates"]:
            raise ValueError("Pre-outcome original Event-to-ANY-occurrence linked date support changed")
        if raw_eventrows[key]!=s["original_event_rows"]:
            raise ValueError("Pre-outcome original 2025 Event fragments changed")
        values={f:int(bysite[key].get(f,0)) for f in FAMILIES}
        original.append({
            "site_id_SHA256":key,
            "source_band":s["source_band"],
            "original_1degree_region":s["original_1degree_region"],
            "2022_nearest_registered_EEA_FW_outfall_SHA256":s["nearest_original_2022_EEA_outfall_id_SHA256"],
            "2022_nearest_qualifying_WWTP_source_distance_km":s["source_distance_km"],
            "original_WWTP_spatial_component_5km_SHA256":s.get("5km_original_WWTP_spatial_component_SHA256"),
            "original_WWTP_spatial_component_10km_SHA256":s.get("10km_original_WWTP_spatial_component_SHA256"),
            "original_all_Event_distinct_date_count":len(everdays[key]),
            "original_any_Occurrence_linked_distinct_date_count":len(linked[key]),
            "source_gate_documented_2plus_any_occurrences":s["source_gate_documented_2plus_any_occurrences"],
            "family_individual_count_original":values,
            "target_3family_positive_individuals":sum(values[f] for f in TARGET),
            "original_2025_monthly_ANY_Occurrence_date_counts":s["2025_occurrence_linked_date_counts_by_month"],
            "not_claimed_unbiased_population_census":True,
        })
    if len(original)!=399:
        raise ValueError("Some selected new original sites absent from source audit")
    if sum(x["source_gate_documented_2plus_any_occurrences"] for x in original)!=len(eligible):
        raise ValueError("Source-gate sample count changed")
    return original,dict(qc)


def run(root):
    mp=root/MANIFEST;tp=root/TAX
    m=json.loads(mp.read_text(encoding="utf-8"))
    if m["schema"]!="gesui_DK_eBMS_2025_NEW399_original_event_and_any_occurrence_link_source_only_v01":
        raise ValueError("Pre-outcome new399 source manifest not frozen")
    if not m["frozen_before_2025_NEW399_butterfly_counts_gate"]["2025_NEW399_butterfly_family_count_outcome_permitted"]:
        raise ValueError("Original source site/effort gate FAILED: cannot open new399 family counts")
    if not m["NEW2025_Lycaenidae_or_any_butterfly_family_individualCount_NOT_READ"]:
        raise ValueError("Source gate was not independently sampled before family result")
    sites=m["all_original_399_source_only_new_site_2025_documented_visit_metadata"]
    if len(sites)!=399:raise ValueError("Source frame lost one or more original 2025 sites")
    accepted=family_mapping(tp.read_text(encoding="utf-8"))
    if len(accepted)<30:raise ValueError("Original fixed family taxonomy unexpectedly small")
    raw=first.download()
    if sha(raw)!=first.ORIGINAL_SHA:raise ValueError("Licensed source ZIP mismatch")
    allsites,qc=process_archive(raw,sites,accepted)
    supported=[x for x in allsites if x["source_gate_documented_2plus_any_occurrences"]]
    cats=Counter(x["source_band"] for x in supported)
    if cats!=Counter(m["source_eligible_new_sites_with_2plus_any_published_occurrence_dates_by_band"]):
        raise ValueError("Pre-outcome support classes changed after opening family counts")
    region=defaultdict(lambda:defaultdict(list))
    for r in supported:region[r["original_1degree_region"]][r["source_band"]].append(r)
    both=sorted([r for r,x in region.items() if x["near5"] and x["remote10"]])
    near=[x for x in supported if x["source_band"]=="near5"]
    bysource=defaultdict(list)
    bycomponent5=defaultdict(list)
    for x in near:
        bysource[x["2022_nearest_registered_EEA_FW_outfall_SHA256"]].append(x)
        bycomponent5[x["original_WWTP_spatial_component_5km_SHA256"]].append(x)
    n_withtarget=sum(
        sum(z["target_3family_positive_individuals"] for z in arr)>0 for arr in bycomponent5.values())
    both_and_target=[
        cell for cell in both if
        sum(x["target_3family_positive_individuals"] for x in region[cell]["near5"])>0 and
        sum(x["target_3family_positive_individuals"] for x in region[cell]["remote10"])>0
    ]
    outcome_supported=bool(cats["near5"]>=40 and cats["remote10"]>=60 and
        n_withtarget>=10 and len(both_and_target)>=3)
    outcome={
        "schema":"gesui_DK_eBMS_2025_NEW399_source_locked_licensed_individualCount_family_panel_v01",
        "original_source_399_2025_event_any_occurrence_gate_SHA256":sha(mp.read_bytes()),
        "original_2023_2025_blinded_accepted_species_taxonomy_SHA256":sha(tp.read_bytes()),
        "original_2025_eBMS_v14_zip_SHA256":first.ORIGINAL_SHA,
        "new_2025_butterfly_family_OUTCOME_READ_only_AFTER_source_gate":True,
        "original_source_year_2022_EEA_registered_WWTP_is_NOT_measured_2025_biochemical_dose":True,
        "original_new399_not_previously_examined_107_sites":True,
        "2025_original_source_and_outcome_QA_counts":qc,
        "original_new399_source_with_published_any_occurrence_date_support_site_n":len(supported),
        "source_band_2025_documented_site_counts":dict(cats),
        "all399_source_site_original_family_counts_private_hash_aggregate_only":allsites,
        "source_band_2025_family_original_counts_and_documented_visit_days":{
            typ:stats([x for x in supported if x["source_band"]==typ])
            for typ in ("near5","remote10","intermediate")
        },
        "source_component5_and_outfall_original_source_count":len(bycomponent5),
        "original_registered_outfall_IDs_with_eligible_near_sites":len(bysource),
        "source_5km_components_with_any_2025_original_target_family_individualCount":n_withtarget,
        "country_1degree_regions_with_near_remote_eligible_sites":both,
        "regions_with_both_source_bands_and_2025_3family_recorded":both_and_target,
        "focal_family_outcome_support_GATE":{
            "min_near_source_components_with_at_least_one_three_family_butterfly_record":10,
            "min_source_near_sites":40,
            "min_registered_remote_sites":60,
            "min_1degree_regions_with_both_source_status_and_target_3family_positive":3,
            "original_new399_2025_family_effect_descriptive_contrasts_permitted":outcome_supported,
        },
        "source_region_original_2025_families":{
            cell:{typ:stats(x) for typ,x in types.items()}
            for cell,types in sorted(region.items())
        },
        "registered_source_outfall_original_2025_family_counts":{
            sid:stats(items) for sid,items in sorted(bysource.items())
        },
        "registered_source_5km_spatial_component_original_2025_family_counts":{
            sid:stats(items) for sid,items in sorted(bycomponent5.items())
        },
        "original_event_without_any_occurrence_not_treated_as_verified_zero_butterflies":True,
        "source_site_rows_not_independent_WWTP_treatment_replicates":True,
        "2025_all_observation_site_family_counts_are_public_licensed_CC_BY_NC_4_0_aggregates":True,
        "no_raw_DWCA_occurrence_or_original_personal_observer_record_republished":True,
    }
    if outcome_supported:
        nearstat=outcome["source_band_2025_family_original_counts_and_documented_visit_days"]["near5"]
        remstat=outcome["source_band_2025_family_original_counts_and_documented_visit_days"]["remote10"]
        outcome["posthoc_registered_WWTP_near_vs_remote_source_2025_descriptive_RECORDING_only"]={
            "original_near5_Lycaenidae_per_any_occurrence_linked_day":nearstat["Lycaenidae_per_any_published_Occurrence_linked_date"],
            "original_remote10_Lycaenidae_per_any_occurrence_linked_day":remstat["Lycaenidae_per_any_published_Occurrence_linked_date"],
            "near_divided_remote_Lycaenidae_recording_per_linked_day_ratio":
                nearstat["Lycaenidae_per_any_published_Occurrence_linked_date"]/
                remstat["Lycaenidae_per_any_published_Occurrence_linked_date"]
                if remstat["Lycaenidae_per_any_published_Occurrence_linked_date"] else None,
            "original_near5_family_photo_or_survey_count_threefamily_share":nearstat["Lycaenidae_share_of_three_target_families_if_positive_total"],
            "original_remote10_family_photo_or_survey_count_threefamily_share":remstat["Lycaenidae_share_of_three_target_families_if_positive_total"],
            "no_source_causal_claim_or_confirmatory_p_value":True,
        }
    return outcome


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    d=run(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "source_only_gate_checked_before_family":d["new_2025_butterfly_family_OUTCOME_READ_only_AFTER_source_gate"],
        "source_groups":d["source_band_2025_family_original_counts_and_documented_visit_days"],
        "published_family_group_support":d["focal_family_outcome_support_GATE"],
        "new399_near_remote_comparison":d.get("posthoc_registered_WWTP_near_vs_remote_source_2025_descriptive_RECORDING_only"),
        "effective_outfall_component_count":d["source_component5_and_outfall_original_source_count"]
    },indent=2))


if __name__=="__main__":
    main()
