#!/usr/bin/env python3
"""POSTHOC source-blind-to-response Denmark 2025 NEW399 eBMS family audit.

Reads species/individualCount ONLY after the original NEW399 Event-to-Occurrence
ID-only support gate has been SHA-frozen and passed, with 107 previously
analyzed 2023–25 coordinate sites excluded. Reports all signs. A 2022
registered WWTP source 5km from a site is NOT confirmed hydrological exposure.
"""
from __future__ import annotations
import argparse,hashlib,io,json,math,zipfile
import xml.etree.ElementTree as ET
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import gesui_DK_EBMS_V14_2023_2024_licensed_family_count_panel_after_taxonomy as prior

DIR=Path("results/ebms_denmark_2023_2024_WWTP_source_only")
GATE=DIR/"DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json"
ORIGINAL=DIR/"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
FAMILIES=("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae")
TARGET=("Lycaenidae","Pieridae","Papilionidae")

def checksum(data):
    return hashlib.sha256(data).hexdigest()

def group_report(rows):
    counts={f:sum(r["family_individualCount"][f] for r in rows) for f in FAMILIES}
    L=counts["Lycaenidae"];target=sum(counts[f] for f in TARGET)
    days=sum(r["linked_visit_days"] for r in rows)
    original_days=sum(r["original_event_visit_days"] for r in rows)
    return {
        "site_rows":len(rows),
        "original_Event_visit_dates":original_days,
        "dates_with_any_published_Occurrence":days,
        "families_individualCount":counts,
        "sites_with_any_target_family_individuals":sum(r["total_target_3families"]>0 for r in rows),
        "sites_with_Lycaenidae_recorded":sum(r["family_individualCount"]["Lycaenidae"]>0 for r in rows),
        "Lycaenidae_per_any_published_Occurrence_date":L/days if days else None,
        "Lycaenidae_share_among_three_family_counts":L/target if target else None,
        "family_recorded_per_any_published_Occurrence_date":{f:counts[f]/days if days else None for f in FAMILIES},
        "not_population_survey_or_WWTP_response":True,
    }

def nearest_band_summary(rows,group):
    return group_report([r for r in rows if r["source_band"]==group])

def run(root):
    gate_file=root/GATE
    gate=json.loads(gate_file.read_text(encoding="utf-8"))
    if gate["schema"]!="gesui_DK_eBMS_2025_NEW399_original_event_and_any_occurrence_link_source_only_v01":
        raise ValueError("Unexpected pre-outcome 2025 NEW399 source cohort")
    if gate["no_2025_butterfly_family_response_accessed_for_new_sites"] is not True:
        raise ValueError("Original NEW399 gate included butterfly outcomes")
    if gate["frozen_before_2025_NEW399_butterfly_counts_gate"]["2025_NEW399_butterfly_family_count_outcome_permitted"] is not True:
        raise ValueError("ORIGINAL new-site source/effort gate failed; STOP outcome read")
    if gate["original_denmark_eBMS_v14_ZIP_SHA256"]!=prior.ORIGINAL_SHA:
        raise ValueError("Frozen source archive differs from taxonomy source")
    if gate["original_2023_24_Event_source_only_735site_SHA256"]!=checksum((root/ORIGINAL).read_bytes()):
        raise ValueError("Underlying original independently fixed Event sites source changed")
    if gate["new_2025_year_sourced_station_n"]!=399:
        raise ValueError("NEW399 original independent coordinate frame changed")
    m=gate["all_original_399_source_only_new_site_2025_documented_visit_metadata"]
    if len(m)!=399:raise ValueError("Bad source-only site metadata cohort")
    st={x["site_id_SHA256"]:x for x in m}
    if len(st)!=399:raise ValueError("Duplicate frozen station SHA")
    if Counter(x["source_band"] for x in st.values())!=Counter({"near5":100,"remote10":127,"intermediate":172}):
        raise ValueError("Original source exposure categories changed")
    if gate["source_eligible_new_sites_with_2plus_any_published_occurrence_dates_by_band"]!={"intermediate":124,"near5":70,"remote10":92}:
        raise ValueError("Previous recording-support class counts unexpectedly changed")

    original=json.loads((root/ORIGINAL).read_text(encoding="utf-8"))
    sitebyid={x["site_id_SHA256"]:x for x in original["original_2023_2024_repeated_event_site_source_metadata"]}
    site_coords={f"{sitebyid[k]['lat']:.5f}|{sitebyid[k]['lon']:.5f}":k for k in st}
    if len(site_coords)!=399:raise ValueError("Published exact 5-decimal coordinates not unique")
    assert all(checksum(loc.encode())==sid for loc,sid in site_coords.items())

    __,taxon_manifest,__,old_sites,taxa,families,source_hashes=prior.load_sources(root)
    if len(old_sites)!=107 or len(taxa)!=112:
        raise ValueError("Original blinded taxonomy mapping is incomplete")
    if not taxon_manifest["original_individualCount_values_NOT_READ"]:
        raise ValueError("Taxonomy family mapping selected from butterfly counts")
    raw=prior.download()
    if checksum(raw)!=prior.ORIGINAL_SHA:raise ValueError("Licensed original DwCA ZIP checksum changed")

    event_lookup={}
    src_day=defaultdict(set)
    linked_day=defaultdict(set)
    counts_by_site=defaultdict(Counter)
    name_missing=Counter()
    status=Counter()
    audit=Counter()
    with zipfile.ZipFile(io.BytesIO(raw)) as zipf:
        meta=ET.fromstring(zipf.read("meta.xml"))
        core=meta.find(prior.NS+"core")
        ext=[x for x in meta.findall(prior.NS+"extension")
             if x.attrib.get("rowType")=="http://rs.tdwg.org/dwc/terms/Occurrence"]
        if core is None or len(ext)!=1:raise ValueError("Expected single original Event core and Occurrence extension")
        ef,ename,eskip,eindex=prior.core_terms(core)
        of,oname,oskip,oindex=prior.core_terms(ext[0])
        idf=eindex("id");ido=oindex("coreid")
        for field in ("eventDate","decimalLatitude","decimalLongitude"):
            if field not in ef:raise ValueError("Missing original Event column")
        for field in ("scientificName","individualCount","occurrenceStatus"):
            if field not in of:raise ValueError("Missing original Occurrence field")
        unique_event=set()
        for e in prior.stream_csv(zipf,ename,eskip):
            audit["original_Event_rows"]+=1
            if len(e)<=max(idf,ef["eventDate"],ef["decimalLatitude"],ef["decimalLongitude"]):
                raise ValueError("Truncated original Event row")
            eid=e[idf].strip()
            if not eid or eid in unique_event:raise ValueError("Duplicate or blank Event ID")
            unique_event.add(eid)
            date=e[ef["eventDate"]].strip()[:10]
            if len(date)!=10 or not date[:4].isdigit() or int(date[:4])!=2025:continue
            try:
                la=float(e[ef["decimalLatitude"]]);lo=float(e[ef["decimalLongitude"]])
            except (ValueError,TypeError):continue
            if not math.isfinite(la) or not math.isfinite(lo):continue
            sid=site_coords.get(f"{la:.5f}|{lo:.5f}")
            if sid is None:continue
            src_day[sid].add(date)
            event_lookup[eid]=(sid,date)
            audit["original_NEW399_2025_Event_rows"]+=1

        for v in prior.stream_csv(zipf,oname,oskip):
            audit["original_Occurrence_rows"]+=1
            if len(v)<=max(ido,of["scientificName"],of["individualCount"],of["occurrenceStatus"]):
                raise ValueError("Truncated original Occurrence row")
            evt=event_lookup.get(v[ido].strip())
            if not evt:continue
            sid,date=evt
            linked_day[sid].add(date)
            audit["NEW399_2025_linked_Occurrence_rows"]+=1
            s=v[of["occurrenceStatus"]].strip()
            status[s or "(blank)"]+=1
            if s and s.lower() not in ("present","presence"):continue
            n=prior.floatcount(v[of["individualCount"]])
            if n is None:
                audit["noninteger_or_nonpositive_individualCount"]+=1
                continue
            scientific=v[of["scientificName"]].strip()
            if scientific not in taxa:
                name_missing[scientific]+=1
                continue
            family=families.get(scientific)
            if family not in FAMILIES:
                audit["other_or_unresolved_taxon_rows"]+=1
                continue
            counts_by_site[sid][family]+=n
            audit["accepted_original_family_occurrence_rows"]+=1

    if audit["original_Event_rows"]!=15606 or audit["original_Occurrence_rows"]!=21617:
        raise ValueError("Original official eBMS archive row totals differ")
    if name_missing:raise ValueError("New taxon scientificName absent from independently frozen GBIF mapping")
    outcome_rows=[]
    for sid in sorted(st):
        r=st[sid]
        if r["original_event_date_count"]!=len(src_day[sid]):
            raise ValueError("Original site Event days no longer match source freeze")
        if r["original_linked_occurrence_dates"]!=len(linked_day[sid]):
            raise ValueError("Original site Event-to-Occurrence date support changed")
        family={f:int(counts_by_site[sid][f]) for f in FAMILIES}
        row={
            "original_public_site_SHA256":sid,
            "source_band":r["source_band"],
            "original_geo_1deg_cell":r["original_1degree_region"],
            "original_registered_WWTP_id_SHA256":r["nearest_original_2022_EEA_outfall_id_SHA256"] if r["source_band"]=="near5" else None,
            "nearest_WWTP_5km_spatial_component_SHA256":r.get("5km_original_WWTP_spatial_component_SHA256"),
            "nearest_WWTP_10km_spatial_component_SHA256":r.get("10km_original_WWTP_spatial_component_SHA256"),
            "source_distance_km":r["source_distance_km"],
            "original_event_visit_days":len(src_day[sid]),
            "linked_visit_days":len(linked_day[sid]),
            "source_gate_2plus_any_Occurrence_dates":r["pre_2025_family_any_occurrence_linked_site_source_gate_pass"],
            "family_individualCount":family,
            "total_target_3families":sum(family[f] for f in TARGET),
            "no_original_Event_occurrence_personal_identifiers_exported":True,
        }
        outcome_rows.append(row)
    assert len(outcome_rows)==399
    sampled=[x for x in outcome_rows if x["source_gate_2plus_any_Occurrence_dates"]]
    assert Counter(x["source_band"] for x in sampled)==Counter({"near5":70,"remote10":92,"intermediate":124})

    near=[x for x in sampled if x["source_band"]=="near5"]
    remote=[x for x in sampled if x["source_band"]=="remote10"]
    near_by_source=defaultdict(list)
    near_by_component=defaultdict(list)
    for r in near:
        near_by_source[r["original_registered_WWTP_id_SHA256"]].append(r)
        near_by_component[r["nearest_WWTP_5km_spatial_component_SHA256"]].append(r)
    counted_components=set(
        r["nearest_WWTP_5km_spatial_component_SHA256"]
        for r in near if r["total_target_3families"]>0)
    cells=defaultdict(lambda:defaultdict(list))
    for row in sampled:
        cells[row["original_geo_1deg_cell"]][row["source_band"]].append(row)
    supported_cells=[
        k for k,v in cells.items()
        if any(r["total_target_3families"]>0 for r in v["near5"])
        and any(r["total_target_3families"]>0 for r in v["remote10"])
    ]
    biological_support={
        "prespecified_source_only_sample_supported":True,
        "min_registered_near_WWTP_5km_components_with_any_3family_record":10,
        "actual_near_WWTP_5km_components_with_any_3family_record":len(counted_components),
        "min_source_near_site_rows_with_2plus_linked_occ_dates":40,
        "min_remote_site_rows_with_2plus_linked_occ_dates":60,
        "actual_near_site_rows":len(near),
        "actual_remote_site_rows":len(remote),
        "min_1degree_regions_with_both_statuses_and_target_family":3,
        "actual_1degree_regions_with_both_and_any_target":len(supported_cells),
        "original_2025_family_contrast_reporting_gate_PASS":(
            len(counted_components)>=10 and len(near)>=40 and len(remote)>=60
            and len(supported_cells)>=3
        ),
        "source_family_observation_gate_selected_after_count_availability_not_sign":True,
        "geographical_WWTP_sources_are_not_independent_causal_treatments":True,
    }
    group={g:nearest_band_summary(sampled,g) for g in ("near5","remote10","intermediate")}
    out={
        "schema":"gesui_DK_eBMS_2025_new399_sites_286_documented_EEA22_WWTP_Lycaenidae_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "source_blinded_to_2025_family_frozen_PREOUTCOME_JSON_SHA256":checksum(gate_file.read_bytes()),
        "original_eBMS_v14_CC_BY_NC_licensed_DWCA_SHA256":prior.ORIGINAL_SHA,
        "original_Event_to_Occurrence_source_data_rows":dict(audit),
        "2025_new399_taxon_family_accepted_count":len(families),
        "2025_source_original_399_site_year_counts_all_before_quality_subset":{g:nearest_band_summary(outcome_rows,g) for g in ("near5","remote10","intermediate")},
        "2025_source_frozen_286_documented_site_counts":group,
        "source_family_model_sample_gate":biological_support,
        "2025_original_399_site_accepted_family_counts":outcome_rows,
        "qualifying_near_2022_registered_WWTP_outfall_group_aggregates":[
            {"original_WWTP_outfall_id_SHA256":k,
             "source_group":group_report(v),
             "2022_original_source_static_not_a_measured_2025_nutrient_dose":True}
            for k,v in sorted(near_by_source.items())],
        "qualifying_near_2022_registered_WWTP_5km_connected_source_group_aggregates":[
            {"WWTP_5km_connected_component_SHA256":k,"WWTP_outfall_IDs_n":len(set(r["original_registered_WWTP_id_SHA256"] for r in v)),"source_group":group_report(v)}
            for k,v in sorted(near_by_component.items())],
        "1degree_source_near_vs_remote_annual_family_raw_group_counts":[
            {"original_1degree_cell":k,
             "near5":group_report(v["near5"]),
             "remote10":group_report(v["remote10"]),
             "intermediate":group_report(v["intermediate"]),
             "near_remote_families_both_observed":k in supported_cells}
            for k,v in sorted(cells.items())],
        "new2025_site_independent_from_prior_107_censored_coordinates":True,
        "only_original_2025_standardized_individualCount_not_PN_water_or_butterfly_biomass":True,
        "2022_register_source_distance_NOT_actual_downstream_exposure":True,
        "recorded_butterflies_on_published_date_are_not_complete_faunal_census":True,
        "current_posthoc_positive_or_negative_response_cannot_establish_WWTP_causal_effect":True,
    }
    # Numeric source-distance comparison permitted only if original source/family
    # observation gates are both already passed. All signs are retained.
    if biological_support["original_2025_family_contrast_reporting_gate_PASS"]:
        n=group["near5"];r=group["remote10"]
        lnr=n["Lycaenidae_per_any_published_Occurrence_date"]
        lrr=r["Lycaenidae_per_any_published_Occurrence_date"]
        sn=n["Lycaenidae_share_among_three_family_counts"]
        sr=r["Lycaenidae_share_among_three_family_counts"]
        out["2025_new_sites_source_near_vs_remote_descriptive"]={
            "Lycaenidae_recording_rate_ratio_near_vs_remote":lnr/lrr if lrr else None,
            "Lycaenidae_per_any_linked_day_near":lnr,
            "Lycaenidae_per_any_linked_day_remote":lrr,
            "Lycaenidae_threefamily_share_near":sn,
            "Lycaenidae_threefamily_share_remote":sr,
            "Lycaenidae_share_percentage_point_difference_near_minus_remote":100*(sn-sr) if sn is not None and sr is not None else None,
            "ratio_of_recording_rates_by_family_near_remote":{
                f:(n["family_recorded_per_any_published_Occurrence_date"][f]/r["family_recorded_per_any_published_Occurrence_date"][f]
                   if r["family_recorded_per_any_published_Occurrence_date"][f] else None)
                for f in FAMILIES
            },
            "only_observational_recorded_individuals_no_confirmed_WWTP_receptor_connection":True,
        }
    else:
        out["2025_new_sites_source_near_vs_remote_descriptive"]=None
        out["stop_reason"]="PREDECLARED_BUTTERFLY_FAMILY_SUPPORT_FAIL_no_spatial_comparison_interpretation"
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    res=run(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(res,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "gate":res["source_family_model_sample_gate"],
        "near":res["2025_source_frozen_286_documented_site_counts"]["near5"],
        "remote":res["2025_source_frozen_286_documented_site_counts"]["remote10"],
        "comparison":res["2025_new_sites_source_near_vs_remote_descriptive"],
        "sample_source_rows":res["original_Event_to_Occurrence_source_data_rows"]
    },indent=2))

if __name__=="__main__":
    main()
