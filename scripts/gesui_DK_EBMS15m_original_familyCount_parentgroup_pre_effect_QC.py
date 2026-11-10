#!/usr/bin/env python3
"""Licensed DK EBMS15m v1.3 original 2023-25 Event/Occurrence count integrity, grouped by parent surveys.

This is a PRE-SPATIAL-EFFECT gate. It does not fit near-vs-remote coefficients,
p-values or selected species effects. Frozen original species taxonomy and
source-only parent/outfall gates are checked BEFORE opening individualCount.
"""
from __future__ import annotations

import argparse,csv,hashlib,io,json,math,urllib.request,zipfile
import xml.etree.ElementTree as ET
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

NS="{http://rs.tdwg.org/dwc/text/}"
FOLDER="results/ebms_denmark_15minute_2023_2025_WWTP_SOURCE_ONLY/"
EEA="results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
EVENT=FOLDER+"ORIGINAL_DK_EBMS15M_2023_25_EEA2022_203_OUTFALL_EVENT_ONLY_QA.json"
PARENT=FOLDER+"ORIGINAL_DK_EBMS15MIN_EVENT_PARENT_SOURCE_SYSTEM_NOCOUNTS_QA.json"
TAX=FOLDER+"ORIGINAL_DK_EBMS15MIN_BLIND_GBIF_ACCEPTED_FAMILIES_NO_COUNTS.json"
RAW="https://ipt.brc.ac.uk/archive.do?r=ebms_denmark_15-minute_counts&v=1.3"
PINNED="c32be29200134fd2ee06b7dd34d76e1420bd6a7810574fb8ffb345c1451169cb"
FAMILIES=("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae","Riodinidae")
YEARS=(2023,2024,2025)

def sha(b):return hashlib.sha256(b).hexdigest()
def nfloat(v):
    try:
        z=float(v)
        return z if math.isfinite(z) else None
    except (TypeError,ValueError):return None
def hav(a,b,c,d):
    r=math.pi/180
    dl=(a-c)*r;dg=(b-d)*r
    v=math.sin(dl/2)**2+math.cos(a*r)*math.cos(c*r)*math.sin(dg/2)**2
    return 12742*math.asin(min(1,math.sqrt(v)))
def openarchive():
    with urllib.request.urlopen(urllib.request.Request(RAW,headers={"User-Agent":"gesui-EBMS15m-parent-group-count-locked-v13/0.1"}),timeout=140) as r:b=r.read()
    if sha(b)!=PINNED:raise ValueError("Original licensed v1.3 bytes changed")
    return b
def source(root):
    files={"event":EVENT,"taxonomy":TAX,"hierarchy":PARENT,"EEA2022":EEA}
    data={k:json.loads((root/v).read_text(encoding="utf-8")) for k,v in files.items()}
    if not data["event"]["strict_500m_primary_5km_butterfly_count_access_gate"]:
        raise ValueError("Original source Event count gate STOP")
    if not data["hierarchy"]["2024_2025_both_parent_cluster_support_gate_pass"]:
        raise ValueError("Parent independent source gate STOP")
    if data["taxonomy"]["original_source_archive_SHA256"]!=PINNED or not data["taxonomy"]["no_original_individualCount_values_read"]:
        raise ValueError("No original source-blind taxonomy")
    if not data["taxonomy"]["primary_species_status_accepted_or_synonym_fixed_before_any_individualCount"]:
        raise ValueError("No original GBIF accepted/synonym family classification policy")
    if data["hierarchy"]["original_archive_sha256"]!=PINNED:raise ValueError("Different parent Event source")
    if data["taxonomy"]["original_Event_only_source_gate_file_SHA256"]!=sha((root/EVENT).read_bytes()):
        raise ValueError("Original taxonomy frozen against different Event source")
    outfalls=data["EEA2022"]["original_2022_EEA_DK_qualifying_Fw_outfalls"]
    if len(outfalls)!=203:raise ValueError("Original 203 outfalls changed")
    byname={x["original_scientificName"]:x for x in data["taxonomy"]["original_scientificName_GBIF_family_mapping_frozen"]}
    if len(byname)!=81:raise ValueError("Original taxonomy names changed")
    return data,outfalls,byname,{k:sha((root/v).read_bytes()) for k,v in files.items()}

def meta_tables(zipfile):
    m=ET.fromstring(zipfile.read("meta.xml"))
    e=m.find(NS+"core")
    o=[x for x in m.findall(NS+"extension") if x.attrib.get("rowType")=="http://rs.tdwg.org/dwc/terms/Occurrence"]
    if e is None or len(o)!=1:raise ValueError("No canonical Event + Occurrence tables")
    return e,o[0]

def mapping(node):
    terms={x.attrib["term"].rsplit("/",1)[-1]:int(x.attrib["index"]) for x in node.findall(NS+"field")}
    f=node.find(".//"+NS+"files/"+NS+"location")
    if f is None or not f.text:raise ValueError("No source table path")
    sep=node.attrib.get("fieldsTerminatedBy",r"\t")
    if sep==r"\t":sep="\t"
    if sep!="\t":raise ValueError("Unexpected original separator")
    skip=int(node.attrib.get("ignoreHeaderLines","1"))
    return terms,f.text,skip

def rows(zipfile,path,skip):
    with io.TextIOWrapper(zipfile.open(path),encoding="utf-8",newline="") as f:
        r=csv.reader(f,delimiter="\t")
        for _ in range(skip):next(r,None)
        for x in r:yield x

def positive_integer(x):
    try:
        v=float(x)
        if math.isfinite(v) and v>0 and v.is_integer():return int(v)
    except (ValueError,TypeError,OverflowError):pass
    return None

def run(root):
    data,outfalls,taxa,source_sha=source(root)
    eventmaps={}
    parent={}
    totals=Counter()
    originalids=set()
    archive=openarchive()
    outpoints=[(float(s["lat"]),float(s["lon"]),s["outfall_id_SHA256"]) for s in outfalls]
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        e,o=meta_tables(z)
        ef,ep,es=mapping(e);of,op,os=mapping(o)
        eidx=e.find(NS+"id");oidx=o.find(NS+"coreid")
        if eidx is None or oidx is None:raise ValueError("Event id or Occurrence coreid unavailable")
        ei=int(eidx.attrib["index"]);oi=int(oidx.attrib["index"])
        for name in ("eventDate","decimalLatitude","decimalLongitude"):
            if name not in ef:raise ValueError("original Event critical source term missing")
        for name in ("scientificName","individualCount","occurrenceStatus"):
            if name not in of:raise ValueError("Original licensed Occurrence source term missing")
        for a in rows(z,ep,es):
            totals["full_original_Event_rows"]+=1
            if len(a)<=max(ei,*ef.values()):raise ValueError("Truncated Event")
            key=a[ei].strip()
            if not key or key in originalids:raise ValueError("Duplicated Event core source IDs")
            originalids.add(key)
            date=a[ef["eventDate"]].strip()[:10]
            year=date[:4]
            if year not in ("2023","2024","2025"):continue
            lat=nfloat(a[ef["decimalLatitude"]]);lon=nfloat(a[ef["decimalLongitude"]])
            unc=nfloat(a[ef["coordinateUncertaintyInMeters"]]) if "coordinateUncertaintyInMeters" in ef else None
            if lat is None or lon is None or unc is None or unc>500 or not(54<=lat<=58 and 7<=lon<=16):continue
            yr=int(year)
            par=a[ef["parentEventID"]].strip() if "parentEventID" in ef else ""
            groupID=par or key
            nearest=min(((hav(lat,lon,aa,bb),sid) for aa,bb,sid in outpoints),key=lambda t:(t[0],t[1]))
            distance,sid=nearest
            category="near2" if distance<=2 else "near5" if distance<=5 else "remote10" if distance>=10 else "middle"
            band="near5cumulative" if category in ("near2","near5") else category
            region=str(math.floor(lat))+"|"+str(math.floor(lon))
            groupkey=(yr,groupID)
            eventmaps[key]=groupkey
            p=parent.setdefault(groupkey,{
                "year":yr,"Event_child_n":0,"band":set(),"region":set(),
                "near_outfall_sha256":set(),"has_any_original_Occurrence_link":False,
                "original_Occurrence_rows":0,
                "individual_counts_accepted_or_synonym":Counter(),
                "individual_counts_accepted_only":Counter(),
            })
            p["Event_child_n"]+=1;p["band"].add(band);p["region"].add(region)
            if band=="near5cumulative":p["near_outfall_sha256"].add(sid)
            totals["source_DK_strict_coordinate_Event_rows"]+=1
        if totals["full_original_Event_rows"]!=44884:raise ValueError("Original original Event count 44884 changed")
        for a in rows(z,op,os):
            totals["full_original_Occurrence_rows"]+=1
            if len(a)<=max(oi,of["scientificName"],of["individualCount"],of["occurrenceStatus"]):
                raise ValueError("Truncated Occurrence")
            event=a[oi].strip()
            if event not in eventmaps:
                totals["original_occurrence_on_other_eventyears_or_geographies"]+=1
                continue
            g=parent[eventmaps[event]]
            g["has_any_original_Occurrence_link"]=True
            g["original_Occurrence_rows"]+=1
            totals["original_occurrence_rows_linked_eligible_source_year_quality"]+=1
            name=a[of["scientificName"]].strip()
            status=a[of["occurrenceStatus"]].strip()
            if status and status.lower() not in ("present","presence"):
                totals["non_present_explicit_original_occurrence_status_rows"]+=1
                continue
            qty=positive_integer(a[of["individualCount"]].strip())
            if qty is None:
                totals["invalid_or_zero_unverified_original_individualCount_rows"]+=1
                continue
            tax=taxa.get(name)
            if tax is None:
                totals["source_scientificName_absent_from_frozen_81_name_taxonomy"]+=1
                raise ValueError("Original recorded scientificName not in frozen taxonomy")
            fam=tax.get("original_GBIF_family")
            if tax["family_supported_exact_species_accepted_or_synonym"]:
                g["individual_counts_accepted_or_synonym"][fam]+=qty
                totals["family_accepted_or_synonym_original_Occurrence_rows"]+=1
            else:totals["strict_taxonomy_unclassified_positive_original_count_rows"]+=1
            if tax["accepted_strict_target_butterfly_species"]:
                g["individual_counts_accepted_only"][fam]+=qty
        if totals["full_original_Occurrence_rows"]!=39905:raise ValueError("Original 39905 licensed Occurrence count changed")
    group_summary=defaultdict(lambda:{"groups":0,"groups_with_occurrence":0,"event_children":0,"occurrence_records":0,
                                      "family_counts":Counter(),"family_counts_strict_accepted":Counter(),
                                      "near_outfall_group_ids":defaultdict(set),
                                      "country1deg_regions":set(),"groups_with_target_positive":0,"ambiguous":0})
    original_unassigned=Counter()
    unit_rows=[]
    for (yr,groupid),g in parent.items():
        if len(g["band"])!=1 or len(g["region"])!=1 or len(g["near_outfall_sha256"])>1:
            original_unassigned[f"ambiguous_parent_group_{yr}"]+=1
            continue
        band=next(iter(g["band"]))
        region=next(iter(g["region"]))
        k=(yr,band)
        t=group_summary[k]
        t["groups"]+=1;t["event_children"]+=g["Event_child_n"]
        t["groups_with_occurrence"]+=bool(g["has_any_original_Occurrence_link"])
        t["occurrence_records"]+=g["original_Occurrence_rows"]
        t["country1deg_regions"].add(region)
        t["family_counts"].update(g["individual_counts_accepted_or_synonym"])
        t["family_counts_strict_accepted"].update(g["individual_counts_accepted_only"])
        t["groups_with_target_positive"]+=sum(g["individual_counts_accepted_or_synonym"][f] for f in ("Lycaenidae","Pieridae","Papilionidae"))>0
        sourceID=next(iter(g["near_outfall_sha256"])) if band=="near5cumulative" else None
        unit_rows.append({
            "original_parent_survey_unit_SHA256":sha((str(yr)+"|"+groupid).encode()),
            "source_only_year":yr,
            "original_geographic_1degree_cell":region,
            "original_qualifying_EEA2022_WWTP_distance_band":band,
            "nearest_registered_EEA_WWTP_outfall_SHA256":sourceID,
            "original_child_Event_core_rows":g["Event_child_n"],
            "published_Original_Occurrence_rows":g["original_Occurrence_rows"],
            "has_any_published_Original_Occurrence_record":g["has_any_original_Occurrence_link"],
            "accepted_exact_GBIF_species_family_individualCount_including_synonyms":dict(g["individual_counts_accepted_or_synonym"]),
            "strict_accepted_only_species_family_counts":dict(g["individual_counts_accepted_only"]),
            "no_site_lat_lon_or_original_personal_event_id_republished":True
        })
        if band=="near5cumulative":
            sid=sourceID
            if g["has_any_original_Occurrence_link"]:
                t["near_outfall_group_ids"][sid].add(groupid)
    stats={}
    for yr in YEARS:
        y={}
        for band in ("near5cumulative","remote10","middle"):
            t=group_summary[(yr,band)]
            ratios=t["family_counts"]
            target=sum(ratios[f] for f in ("Lycaenidae","Pieridae","Papilionidae"))
            sourceN=sum(len(x)>=2 for x in t["near_outfall_group_ids"].values()) if band=="near5cumulative" else None
            y[band]={
                "original_eligible_parent_group_count":t["groups"],
                "original_child_Event_row_count":t["event_children"],
                "original_parent_groups_having_any_linked_Occurrence":t["groups_with_occurrence"],
                "original_parent_groups_with_any_focal_threefamily_count":t["groups_with_target_positive"],
                "original_occurrence_records_linked_to_group":t["occurrence_records"],
                "original_accepted_synonym_family_individualCount":dict(ratios),
                "original_strict_accepted_only_family_individualCount":dict(t["family_counts_strict_accepted"]),
                "target_threefamily_individualCount_total_for_quality":target,
                "source_outfalls_having_at_least_two_parent_groups":sourceN,
                "original_source_1degree_regions":len(t["country1deg_regions"]),
                "spatially_ambiguous_parent_groups_excluded":original_unassigned[f"ambiguous_parent_group_{yr}"],
            }
        near=y["near5cumulative"];far=y["remote10"]
        regionsN=group_summary[(yr,"near5cumulative")]["country1deg_regions"]
        regionsF=group_summary[(yr,"remote10")]["country1deg_regions"]
        shared=len(regionsN & regionsF)
        passed=bool(
            near["original_parent_groups_having_any_linked_Occurrence"]>=50
            and far["original_parent_groups_having_any_linked_Occurrence"]>=100
            and near["source_outfalls_having_at_least_two_parent_groups"]>=10
            and shared>=3
            and near["target_threefamily_individualCount_total_for_quality"]>=50
            and far["target_threefamily_individualCount_total_for_quality"]>=50
        )
        y["matched_1degree_regions_with_both_source_bands"]=shared
        y["precommitted_parent_outfall_family_count_support_gate"]=passed
        stats[str(yr)]=y
    # 2023 was source-only parent unit fail and shall not be a primary inferential year.
    overall=stats["2024"]["precommitted_parent_outfall_family_count_support_gate"] and stats["2025"]["precommitted_parent_outfall_family_count_support_gate"]
    return {
        "schema":"gesui_DK_EBMS15m_2023_25_2022WWTP_licenced_original_event_occurrence_parentgroup_countQC_BEFORE_spatial_effect_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_full_DwCA_v13_SHA256":PINNED,
        "original_source_input_file_SHA256":source_sha,
        "source_taxonomy_family_accepted_synonym_exact_gbif_frozen_before_individualCount":True,
        "original_row_and_taxonomic_quality":dict(totals),
        "original_parent_group_unstable_source_or_geography_excluded":dict(original_unassigned),
        "source_band_by_year_parent_group_observation_and_family_count_support":stats,
        "original_anonymized_parent_group_family_counts_for_posthoc_region_and_outfall_audit":sorted(unit_rows,key=lambda x:(x["source_only_year"],x["original_qualifying_EEA2022_WWTP_distance_band"],x["original_parent_survey_unit_SHA256"])),
        "this_parent_group_export_added_after_original_pooled_family_QC":True,
        "primary_2024_and_2025_after_original_parent_source_family_QC_gate":overall,
        "original_2023_descriptive_only_not_independent_parent_source_support":True,
        "full_original_DwCA_nonpublic_rows_not_republished":True,
        "group_count_with_no_published_Occurrence_UNCERTAIN_not_confirmed_zero_butterflies":True,
        "no_spatial_WWTP_vs_butterfly_effect_computed_at_this_stage":True,
        "license":"CC-BY-NC 4.0"
    }

def main():
    a=argparse.ArgumentParser();a.add_argument("--out",required=True,type=Path);z=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    result=run(root)
    z.out.parent.mkdir(parents=True,exist_ok=True)
    z.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "quality":result["original_row_and_taxonomic_quality"],
        "year":result["source_band_by_year_parent_group_observation_and_family_count_support"],
        "primary_go":result["primary_2024_and_2025_after_original_parent_source_family_QC_gate"]
    },indent=2))
if __name__=="__main__":main()
