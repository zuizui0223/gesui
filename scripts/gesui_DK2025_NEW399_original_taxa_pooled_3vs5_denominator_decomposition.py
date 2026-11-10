#!/usr/bin/env python3
"""Original Danish eBMS 2025 species-level denominator decomposition (POSTHOC).

Re-reads CC-BY-NC v1.4 original Event/Occurrence exact SHA and pre-frozen GBIF
strict species→family mapping. Excludes prior analyzed 107 sites and keeps
original NEW399, 286 publishing-supported sites, 70 near/92 remote. Strictly
reconciles every 399-site×five-family earlier frozen aggregate before releasing
any species-level result. No larval host or causal effluent exposure measured.
"""
from __future__ import annotations
import argparse,hashlib,io,json,math,zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict,Counter
from datetime import datetime,timezone
from pathlib import Path

import gesui_DK_EBMS_V14_2023_2024_licensed_family_count_panel_after_taxonomy as prior

D=Path("results/ebms_denmark_2023_2024_WWTP_source_only")
SOURCE=D/"DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json"
FAMILY=D/"DENMARK_2025_NEW399_SITE_LYCAENIDAE_FAMILY_POST_SOURCE_GATE_COUNTS.json"
ORIGINAL=D/"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
TAXON=D/"DENMARK_EBMS_V14_ORIGINAL_SCIENTIFICNAME_GBIF_FAMILY_BLIND_TAXONOMY.json"
FAMS=("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae")
GROUPS=("near5","remote10","intermediate")
TARGET3=("Lycaenidae","Pieridae","Papilionidae")
TARGET5=FAMS
MIN_SPECIES_RECORDS=20
MIN_SPECIES_SITES=5

def sha(raw):return hashlib.sha256(raw).hexdigest()

def family_totals(counts):
    out={g:Counter() for g in GROUPS}
    for (sid,family),n in counts.items():
        if family not in FAMS:raise ValueError("Unexpected original family")
    return out

def odds_from_counts(group,include_families):
    L=group.get("Lycaenidae",0)
    rest=sum(group.get(f,0) for f in include_families if f!="Lycaenidae")
    return L/rest if L>0 and rest>0 else None

def effect(group_by_species,species):
    # Effect of REMOVING one original species from both near and remote,
    # NOT treatment, or evidence of direct interspecific competition.
    result={}
    for g in ("near5","remote10"):
        x=Counter()
        for taxon,z in group_by_species[g].items():
            if taxon==species:continue
            x[z["family"]]+=z["individuals"]
        result[g]=x
    a3=odds_from_counts(result["near5"],TARGET3)
    b3=odds_from_counts(result["remote10"],TARGET3)
    a5=odds_from_counts(result["near5"],TARGET5)
    b5=odds_from_counts(result["remote10"],TARGET5)
    return {
        "leave_taxon_out_pooled_three_family_OR":a3/b3 if a3 is not None and b3 is not None else None,
        "leave_taxon_out_pooled_five_family_OR":a5/b5 if a5 is not None and b5 is not None else None,
        "hypothetical_species_removal_not_observed_or_treatment":True,
    }

def calculate(root):
    srcpath=root/SOURCE;fampath=root/FAMILY;eventpath=root/ORIGINAL;taxpath=root/TAXON
    src=json.loads(srcpath.read_text(encoding="utf-8"))
    original=json.loads(eventpath.read_text(encoding="utf-8"))
    panel=json.loads(fampath.read_text(encoding="utf-8"))
    tax=json.loads(taxpath.read_text(encoding="utf-8"))
    if not src["frozen_before_2025_NEW399_butterfly_counts_gate"]["2025_NEW399_butterfly_family_count_outcome_permitted"]:
        raise ValueError("Original pre-new399 source-only gate failed")
    if not src["NEW2025_Lycaenidae_or_any_butterfly_family_individualCount_NOT_READ"]:
        raise ValueError("Original site frame source is not response-blind")
    if not panel["source_family_model_sample_gate"]["original_2025_family_contrast_reporting_gate_PASS"]:
        raise ValueError("Previously established family sample gate failed")
    if panel["source_blinded_to_2025_family_frozen_PREOUTCOME_JSON_SHA256"]!=sha(srcpath.read_bytes()):
        raise ValueError("Original NEW399 family source hash changed")
    if src["original_2023_24_Event_source_only_735site_SHA256"]!=sha(eventpath.read_bytes()):
        raise ValueError("Original 735 Event-only source hash changed")
    if src["original_denmark_eBMS_v14_ZIP_SHA256"]!=prior.ORIGINAL_SHA or tax["original_eBMS_DwCA_SHA256"]!=prior.ORIGINAL_SHA:
        raise ValueError("Original licensed Danish v1.4 archive SHA changed")
    if not tax["original_individualCount_values_NOT_READ"]:
        raise ValueError("Taxonomic mapping was not made blind to butterfly responses")
    original_sites={x["site_id_SHA256"]:x for x in original["original_2023_2024_repeated_event_site_source_metadata"]}
    src_sites={x["site_id_SHA256"]:x for x in src["all_original_399_source_only_new_site_2025_documented_visit_metadata"]}
    old_panel={x["original_public_site_SHA256"]:x for x in panel["2025_original_399_site_accepted_family_counts"]}
    if not(len(src_sites)==len(old_panel)==399 and set(src_sites)==set(old_panel)):
        raise ValueError("Original 399 site identities changed")
    qualified=[key for key in sorted(src_sites) if src_sites[key]["pre_2025_family_any_occurrence_linked_site_source_gate_pass"]]
    if len(qualified)!=286:
        raise ValueError("Original published-date support must be 286")
    if Counter(src_sites[k]["source_band"] for k in qualified)!=Counter({"near5":70,"remote10":92,"intermediate":124}):
        raise ValueError("Original 2025 group support changed")

    taxons=tax["original_scientificName_GBIF_taxonomic_source_match_audit"]
    mapping={x["original_scientificName"]:x for x in taxons}
    if len(mapping)!=112:raise ValueError("Original GBIF scientificName cohort changed")
    accepted={}
    for name,t in mapping.items():
        if not t["authoritative_family_accepted"] or t["gbif_family"] not in FAMS:continue
        accepted[name]=(str(t["gbif_usageKey"]),t["gbif_family"],t["gbif_scientificName"])
    expected_map={f:sum(v[1]==f for v in accepted.values()) for f in FAMS}
    if expected_map!={"Lycaenidae":21,"Pieridae":7,"Papilionidae":1,"Nymphalidae":28,"Hesperiidae":8}:
        raise ValueError("Independent accepted GBIF classification was mutated")
    taxid_index={}
    for source_name,(usage,family,canonical) in accepted.items():
        if usage in taxid_index:
            if taxid_index[usage]["family"]!=family:
                raise ValueError("Original GBIF useKey has multiple incompatible families")
            taxid_index[usage]["source_names"].append(source_name)
        else:
            taxid_index[usage]={"family":family,"canonical_species":canonical,"source_names":[source_name]}

    site_coords={}
    for siteid in src_sites:
        r=original_sites[siteid]
        pt=f"{r['lat']:.5f}|{r['lon']:.5f}"
        if sha(pt.encode())!=siteid:raise ValueError("Public site coordinate key mismatch")
        if pt in site_coords:raise ValueError("Rounding joined two original sites")
        site_coords[pt]=siteid
    fam_site=defaultdict(Counter)
    species_site=defaultdict(Counter)
    eventlookup={}
    date_by_site=defaultdict(set)
    linked_by_site=defaultdict(set)
    qc=Counter()
    raw=prior.download()
    if sha(raw)!=prior.ORIGINAL_SHA:raise ValueError("Official licensed Danish DwC-A changed")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        rootxml=ET.fromstring(z.read("meta.xml"))
        core=rootxml.find(prior.NS+"core")
        ext=[x for x in rootxml.findall(prior.NS+"extension") if x.attrib.get("rowType")=="http://rs.tdwg.org/dwc/terms/Occurrence"]
        if core is None or len(ext)!=1:raise ValueError("Original DwCA Event/Occurrence morphology differs")
        ef,ename,eskip,eidx=prior.core_terms(core)
        of,oname,oskip,oidx=prior.core_terms(ext[0])
        eid=eidx("id");oid=oidx("coreid")
        for req in ("eventDate","decimalLatitude","decimalLongitude"):
            if req not in ef:raise ValueError("Original Event source lacks "+req)
        for req in ("scientificName","individualCount","occurrenceStatus"):
            if req not in of:raise ValueError("Original Occurrence source lacks "+req)
        seen=set()
        for row in prior.stream_csv(z,ename,eskip):
            qc["full_original_Event_records"]+=1
            if len(row)<=max(eid,ef["eventDate"],ef["decimalLatitude"],ef["decimalLongitude"]):
                raise ValueError("Truncated original Event record")
            id=row[eid].strip()
            if not id or id in seen:raise ValueError("Duplicate original Event core key")
            seen.add(id)
            date=row[ef["eventDate"]].strip()[:10]
            if date[:4]!="2025":continue
            try:
                lat=float(row[ef["decimalLatitude"]]);lon=float(row[ef["decimalLongitude"]])
            except (ValueError,TypeError):continue
            if not math.isfinite(lat) or not math.isfinite(lon):continue
            siteid=site_coords.get(f"{lat:.5f}|{lon:.5f}")
            if siteid is None:continue
            eventlookup[id]=(siteid,date)
            date_by_site[siteid].add(date)
            qc["original_2025_new399_Event_rows"]+=1
        for row in prior.stream_csv(z,oname,oskip):
            qc["full_original_Occurrence_records"]+=1
            if len(row)<=max(oid,of["scientificName"],of["individualCount"],of["occurrenceStatus"]):
                raise ValueError("Truncated original Occurrence record")
            evt=eventlookup.get(row[oid].strip())
            if not evt:continue
            siteid,date=evt
            qc["original_2025_new399_linked_Occurrence_rows"]+=1
            linked_by_site[siteid].add(date)
            status=row[of["occurrenceStatus"]].strip()
            if status and status.lower() not in ("present","presence"):continue
            count=prior.floatcount(row[of["individualCount"]])
            if count is None:continue
            name=row[of["scientificName"]].strip()
            if name not in mapping:
                raise ValueError("New unrecognized scientificName; cannot reconcile source-blind taxonomy")
            if name not in accepted:continue
            key,family,_=accepted[name]
            fam_site[siteid][family]+=count
            species_site[siteid][key]+=count
            qc["accepted_five_family_original_Occurrence_rows"]+=1
    if qc["full_original_Event_records"]!=15606 or qc["full_original_Occurrence_records"]!=21617:
        raise ValueError("Official licensed archive counts don't match original source")
    for siteid,source in src_sites.items():
        original=old_panel[siteid]
        if len(linked_by_site[siteid])!=original["linked_visit_days"]:
            raise ValueError("Old source-linked reporting dates mismatch")
        if len(date_by_site[siteid])!=original["original_event_visit_days"]:
            raise ValueError("Original Event visit dates mismatch")
        if len(linked_by_site[siteid])!=source["original_linked_occurrence_dates"]:
            raise ValueError("Original source-only before-outcome event linkage changed")
        for family in FAMS:
            if fam_site[siteid][family]!=original["family_individualCount"][family]:
                raise ValueError(f"Original 399 source site familyCount mismatch {siteid} {family}")
        rebuilt=Counter()
        for usage,n in species_site[siteid].items():
            rebuilt[taxid_index[usage]["family"]]+=n
        if any(rebuilt[f]!=fam_site[siteid][f] for f in FAMS):
            raise ValueError("GBIF accepted species totals do not sum to original five families")

    grp={g:defaultdict(lambda:{"individuals":0,"sites":set(),"regions":set()}) for g in GROUPS}
    famgrp={g:Counter() for g in GROUPS}
    region_fam=defaultdict(lambda:defaultdict(Counter))
    for sid in qualified:
        s=src_sites[sid]
        g=s["source_band"]
        for usage,n in species_site[sid].items():
            grp[g][usage]["individuals"]+=n
            grp[g][usage]["sites"].add(sid)
            grp[g][usage]["regions"].add(s["original_1degree_region"])
        for f in FAMS:
            famgrp[g][f]+=fam_site[sid][f]
            region_fam[s["original_1degree_region"]][g][f]+=fam_site[sid][f]
    expected_totals={
        "near5":{"Lycaenidae":382,"Pieridae":680,"Papilionidae":0,"Nymphalidae":3822,"Hesperiidae":372},
        "remote10":{"Lycaenidae":644,"Pieridae":2602,"Papilionidae":0,"Nymphalidae":3544,"Hesperiidae":192},
        "intermediate":{"Lycaenidae":2268,"Pieridae":3595,"Papilionidae":1,"Nymphalidae":9465,"Hesperiidae":866}}
    if any(dict(famgrp[g])!=expected_totals[g] for g in GROUPS):
        raise ValueError("Independent source-locked GROUP family count is not reproduced")
    docs={g:sum(src_sites[sid]["original_linked_occurrence_dates"] for sid in qualified if src_sites[sid]["source_band"]==g) for g in GROUPS}
    if docs!={"near5":539,"remote10":757,"intermediate":1286}:
        raise ValueError("Original source opportunity date count mismatch")

    species_by_group={}
    species_lookup={}
    for usage,info in sorted(taxid_index.items(),key=lambda z:z[1]["canonical_species"]):
        f=info["family"]
        row={
            "gbif_accepted_usageKey":usage,
            "accepted_GBIF_scientificName":info["canonical_species"],
            "accepted_family":f,
            "original_scientificName_aliases_n":len(info["source_names"]),
            "per_group":{}
        }
        for g in GROUPS:
            r=grp[g][usage]
            n=r["individuals"]
            row["per_group"][g]={
                "original_2025_eBMS_individualCount":n,
                "original_distinct_monitored_site_centroids_with_record":len(r["sites"]),
                "original_distinct_1degree_regions_with_record":len(r["regions"]),
                "proportion_of_2025_family_individualCount":n/famgrp[g][f] if famgrp[g][f] else None,
                "individualCount_per_any_published_Occurrence_linked_day":n/docs[g]
            }
        v=row["per_group"]
        nn=v["near5"]["original_2025_eBMS_individualCount"]
        rr=v["remote10"]["original_2025_eBMS_individualCount"]
        count=nn+rr
        sites=v["near5"]["original_distinct_monitored_site_centroids_with_record"]+v["remote10"]["original_distinct_monitored_site_centroids_with_record"]
        row["cross_group_support_min20_individuals_min5_sites"]=count>=MIN_SPECIES_RECORDS and sites>=MIN_SPECIES_SITES
        row["near_remote_relative_individualCount_per_day"]=(nn/docs["near5"])/(rr/docs["remote10"]) if rr>0 else None
        species_lookup[usage]=row
        species_by_group[usage]=row

    def pooled_or(spec_group):
        a=odds_from_counts(spec_group["near5"],TARGET3)
        b=odds_from_counts(spec_group["remote10"],TARGET3)
        aa=odds_from_counts(spec_group["near5"],TARGET5)
        bb=odds_from_counts(spec_group["remote10"],TARGET5)
        return {"OR3":a/b if a and b else None,"OR5":aa/bb if aa and bb else None}
    pooled=pooled_or(famgrp)
    if abs(pooled["OR3"]-(382/680)/(644/2602))>1e-12:
        raise ValueError("Original pooled three family source OR differs")
    if abs(pooled["OR5"]-(382/(680+3822+372))/(644/(2602+3544+192)))>1e-12:
        raise ValueError("Original pooled five family source OR differs")
    # Each species removal is a bookkeeping counterfactual on the SAME data.
    group_species={g:{} for g in ("near5","remote10")}
    for usage,info in taxid_index.items():
        for g in group_species:
            group_species[g][usage]={"family":info["family"],"individuals":grp[g][usage]["individuals"]}
    removal=[]
    for usage,row in species_lookup.items():
        v=effect(group_species,usage)
        r3=v["leave_taxon_out_pooled_three_family_OR"]
        r5=v["leave_taxon_out_pooled_five_family_OR"]
        removal.append({
            "accepted_GBIF_usageKey":usage,
            "scientificName":row["accepted_GBIF_scientificName"],
            "family":row["accepted_family"],
            "near_recorded_individuals":row["per_group"]["near5"]["original_2025_eBMS_individualCount"],
            "remote_recorded_individuals":row["per_group"]["remote10"]["original_2025_eBMS_individualCount"],
            "three_family_OR_after_deleting_species_in_both_groups":r3,
            "five_family_OR_after_deleting_species_in_both_groups":r5,
            "delta_log_pooled_three_OR":math.log(r3/pooled["OR3"]) if r3 else None,
            "delta_log_pooled_five_OR":math.log(r5/pooled["OR5"]) if r5 else None,
            "counterfactual_only_not_population_or_treatment":True
        })
    removal.sort(key=lambda x:(-max(abs(x["delta_log_pooled_three_OR"] or 0),abs(x["delta_log_pooled_five_OR"] or 0)),x["scientificName"]))
    ranked={}
    for family in FAMS:
        ranked[family]={}
        for g in ("near5","remote10"):
            rows=[row for row in species_lookup.values() if row["accepted_family"]==family and row["per_group"][g]["original_2025_eBMS_individualCount"]>0]
            rows.sort(key=lambda r:(-r["per_group"][g]["original_2025_eBMS_individualCount"],r["accepted_GBIF_scientificName"]))
            ranked[family][g]=[{
                "scientificName":r["accepted_GBIF_scientificName"],
                "accepted_usageKey":r["gbif_accepted_usageKey"],
                "original_individualCount":r["per_group"][g]["original_2025_eBMS_individualCount"],
                "source_reporting_site_centroids":r["per_group"][g]["original_distinct_monitored_site_centroids_with_record"],
                "percentage_of_recorded_family_individuals":100*r["per_group"][g]["proportion_of_2025_family_individualCount"]
            } for r in rows[:8]]
    metrics={}
    for family in FAMS:
        metrics[family]={}
        for g in ("near5","remote10"):
            xs=[r["per_group"][g]["original_2025_eBMS_individualCount"] for r in species_lookup.values() if r["accepted_family"]==family]
            tot=sum(xs)
            prop=[x/tot for x in xs if x>0] if tot else []
            metrics[family][g]={
                "species_with_original_recorded_individuals_n":len(prop),
                "original_total_family_individualCount":tot,
                "species_Hill_N2_effective_count":1/sum(z*z for z in prop) if prop else None,
                "largest_species_fraction_of_family":max(prop) if prop else None
            }
    out={
        "schema":"gesui_DK2025_new399_exact_species_denominator_decomposition_source_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "posthoc_species_result_inspection_after_2025_family_3vs5_sign_known":True,
        "licensed_original_Denmark_eBMS_v14_source_SHA256":prior.ORIGINAL_SHA,
        "source_JSON_SHA256":{"source_only_NEW399_before_families":sha(srcpath.read_bytes()),
                             "family_result_original_399_site":sha(fampath.read_bytes()),
                             "source_blind_112taxonomy":sha(taxpath.read_bytes())},
        "source_license":"CC-BY-NC 4.0",
        "original_event_and_occurrence_QC":dict(qc),
        "reconciled_exact_399_source_2025_site_family_individualCounts":True,
        "source_original_five_family_taxonomy_n":len(taxid_index),
        "original_2025_qualifying_source_site_year_counts":dict(Counter(src_sites[k]["source_band"] for k in qualified)),
        "original_2025_qualifying_source_group_published_dates":docs,
        "family_individualCount_by_source_group":{g:dict(famgrp[g]) for g in GROUPS},
        "original_near_vs_remote_three_family_pooled_observational_OR":pooled["OR3"],
        "original_near_vs_remote_five_family_pooled_observational_OR":pooled["OR5"],
        "original_near_vs_remote_absolute_Lycaenidae_recorded_per_linked_day_ratio":
            (famgrp["near5"]["Lycaenidae"]/docs["near5"])/(famgrp["remote10"]["Lycaenidae"]/docs["remote10"]),
        "all_accepted_five_family_species_2025_individual_count_by_source_band":list(species_lookup.values()),
        "species_sensitivity_original_pooled_three_vs_five_family_OR_after_one_species_removed":removal,
        "largest_original_species_contributors_per_family_and_source_band":ranked,
        "original_site_species_recorded_concentration_by_family":metrics,
        "original_shared_1degree_region_five_family_counts":{
            region:{g:{f:int(q[g][f]) for f in FAMS} for g in GROUPS}
            for region,q in region_fam.items()},
        "species_display_min_total_20_individuals_min_5_reporting_sites_only_not_family_estimate_filter":True,
        "outcome_all_species_included_even_sparse_and_zero_original_recorded":True,
        "raw_original_licensed_species_site_individual_rows_not_redistributed":True,
        "recording_composition_not_population_survey_or_causal_WWTP_exposure":True,
        "non_focal_EEA_source_year_2022_not_actual_butterfly_2025_wastewater_dose":True,
        "original_HydroRIVERS_2025_river_directed_source_support_gate_fails_only_two_outfall_reaches":True,
    }
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=calculate(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    major={}
    for f in ("Lycaenidae","Pieridae","Nymphalidae"):
        major[f]={g:out["largest_original_species_contributors_per_family_and_source_band"][f][g][:4] for g in ("near5","remote10")}
    print(json.dumps({
        "source":out["source_JSON_SHA256"],
        "families":out["family_individualCount_by_source_group"],
        "source_sample":out["original_2025_qualifying_source_site_year_counts"],
        "OR3":out["original_near_vs_remote_three_family_pooled_observational_OR"],
        "OR5":out["original_near_vs_remote_five_family_pooled_observational_OR"],
        "top":major,
        "largest_leave_species_out_effect":out["species_sensitivity_original_pooled_three_vs_five_family_OR_after_one_species_removed"][:8]
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
