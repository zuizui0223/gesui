#!/usr/bin/env python3
"""Denmark EBMS original licensed DwC-A 2023/2024 butterfly family count panel.

Measured butterfly outcome data are opened ONLY after prior frozen EEA
source/site support and blinded taxonomic mapping pass. This script
audits source Event<->Occurrence linkage and constructs a site-year panel.
It does NOT compute any WWTP-near vs remote family effect or p-value.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, math, zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path
from datetime import datetime,timezone
from urllib.request import urlopen,Request

URL="https://ipt.brc.ac.uk/archive.do?r=ebms_denmark_transects&v=1.4"
ORIGINAL_SHA="2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b"
NS="{http://rs.tdwg.org/dwc/text/}"
FAMILIES=("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae","Riodinidae")
FIRST=("Lycaenidae","Pieridae","Papilionidae")
YRS=(2023,2024)
ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
SITE=ROOT+"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
TAX=ROOT+"DENMARK_EBMS_V14_ORIGINAL_SCIENTIFICNAME_GBIF_FAMILY_BLIND_TAXONOMY.json"
CLUST=ROOT+"DENMARK_2023_24_EVENT_SITE_WWTP_CLUSTER_PSEUDOREPLICATION_PREOUTCOME_QA.json"

def sha(b):return hashlib.sha256(b).hexdigest()

def core_terms(node):
    fields={x.attrib["term"].rsplit("/",1)[-1]:int(x.attrib["index"]) for x in node.findall(NS+"field")}
    def index(element):
        n=node.find(NS+element)
        if n is None:raise ValueError("Source lacks "+element)
        return int(n.attrib["index"])
    f=node.find(".//"+NS+"files/"+NS+"location")
    if f is None or not f.text:raise ValueError("Missing table location")
    skip=int(node.attrib.get("ignoreHeaderLines","1"))
    delim=node.attrib.get("fieldsTerminatedBy",r"\t")
    if delim==r"\t":delim="\t"
    if delim!="\t" or skip not in (0,1):raise ValueError("Unexpected DwCA delimiter/header")
    return fields,f.text,skip,index

def stream_csv(z,name,skip):
    with io.TextIOWrapper(z.open(name),encoding="utf-8",newline="") as f:
        r=csv.reader(f,delimiter="\t")
        for _ in range(skip):next(r,None)
        for row in r:yield row

def floatcount(raw):
    try:
        n=float(raw)
        if not math.isfinite(n) or n<=0 or not n.is_integer():return None
        return int(n)
    except (ValueError,TypeError,OverflowError):return None

def load_sources(root):
    sb=(root/SITE).read_bytes();tb=(root/TAX).read_bytes();cb=(root/CLUST).read_bytes()
    s=json.loads(sb);t=json.loads(tb);c=json.loads(cb)
    if not s["source_distance_and_yearvisit_eligibility"]["5km"]["butterfly_original_occurrence_extension_may_be_accessed_ONLY_if_true"]:
        raise ValueError("No pre-outcome 5km source gate; cannot open original species counts")
    if s["source_distance_and_yearvisit_eligibility"]["2km"]["butterfly_original_occurrence_extension_may_be_accessed_ONLY_if_true"]:
        raise ValueError("Unexpected newly permitted 2km source study")
    if t["original_eBMS_DwCA_SHA256"]!=ORIGINAL_SHA or not t["original_individualCount_values_NOT_READ"]:
        raise ValueError("Taxonomic family mapping was not source/blinded frozen")
    if c["distinct_nearest_WWTP_outfalls_with_ge5km_near_sites"]!=4 or not c["no_butterfly_outcomes_read"]:
        raise ValueError("Physical independent WWTP source-unit QC not frozen")
    valid=s["original_2023_2024_repeated_event_site_source_metadata"]
    valid=[x for x in valid if x["repeat_2023_2024_dates_gate"]
           and x["precise_500m_explicit_uncertainty_gate"]
           and x["rounded5_site_coord_scatter_le2m"]]
    if len(valid)!=107:raise ValueError("Changed pre-outcome site sample")
    sites={f"{x['lat']:.5f}|{x['lon']:.5f}":x for x in valid}
    if len(sites)!=len(valid):raise ValueError("Duplicate rounded locations in source frame")
    taxa={x["original_scientificName"]:x for x in t["original_scientificName_GBIF_taxonomic_source_match_audit"]}
    if len(taxa)!=112:raise ValueError("Frozen unique source taxonomy not matching original")
    families={k:v.get("gbif_family") for k,v in taxa.items() if v["authoritative_family_accepted"]}
    return s,t,c,sites,taxa,families,{"event_source_SHA256":sha(sb),"taxonomy_SHA256":sha(tb),"source_independence_SHA256":sha(cb)}

def download():
    with urlopen(Request(URL,headers={"User-Agent":"gesui-DK-EBMS-v14-standardized-count-locked-source/0.1"}),timeout=120) as r:
        b=r.read()
    if sha(b)!=ORIGINAL_SHA:raise ValueError("Original Denmark v1.4 zip SHA changed; STOP")
    return b

def make_panel(raw,sites,taxa,families):
    count=Counter()
    stats=defaultdict(Counter)
    occ_status=Counter()
    family_rows=Counter()
    family_counts=Counter()
    event_key_map={}
    event_keys_seen=set()
    all_eventIDs=set()
    event_qualified_days=defaultdict(set)
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        xml=ET.fromstring(z.read("meta.xml"))
        core=xml.find(NS+"core")
        ext=[x for x in xml.findall(NS+"extension") if x.attrib.get("rowType")=="http://rs.tdwg.org/dwc/terms/Occurrence"]
        if core is None or len(ext)!=1 or core.attrib.get("rowType")!="http://rs.tdwg.org/dwc/terms/Event":
            raise ValueError("Wrong original DwCA Event/Occurrence structure")
        ef,efile,es,efield=core_terms(core)
        of,ofile,os,ofield=core_terms(ext[0])
        event_id_index=efield("id");occ_coreid_index=ofield("coreid")
        for x in ["eventDate","decimalLatitude","decimalLongitude"]:
            if x not in ef:raise ValueError("Missing Event source term "+x)
        for x in ["scientificName","individualCount","occurrenceStatus"]:
            if x not in of:raise ValueError("Missing Occurrence source term "+x)
        for row in stream_csv(z,efile,es):
            count["all_event_rows"]+=1
            if len(row)<=max(event_id_index,ef["eventDate"],ef["decimalLatitude"],ef["decimalLongitude"]):
                count["invalid_short_event_rows"]+=1;continue
            corekey=row[event_id_index].strip()
            if not corekey or corekey in event_keys_seen:
                count["event_missing_or_duplicate_corekey"]+=1;continue
            event_keys_seen.add(corekey)
            dt=row[ef["eventDate"]].strip()[:10]
            if len(dt)<10 or not dt[:4].isdigit():
                count["unusable_event_date"]+=1;continue
            yr=int(dt[:4])
            try:
                lat=float(row[ef["decimalLatitude"]]);lon=float(row[ef["decimalLongitude"]])
            except ValueError:
                count["invalid_event_coordinates"]+=1;continue
            if not math.isfinite(lat) or not math.isfinite(lon):continue
            loc=f"{lat:.5f}|{lon:.5f}"
            event_key_map[corekey]=(loc,yr,dt)
            if yr in YRS and loc in sites:
                event_qualified_days[(loc,yr)].add(dt)
                count["qualified_siteyear_event_core_rows"]+=1
        for row in stream_csv(z,ofile,os):
            count["all_occurrence_rows"]+=1
            if len(row)<=max(occ_coreid_index,of["scientificName"],of["individualCount"],of["occurrenceStatus"]):
                count["invalid_short_occurrence_rows"]+=1;continue
            coreid=row[occ_coreid_index].strip()
            if coreid not in event_key_map:
                count["unlinked_occurrence_core_rows"]+=1
                continue
            count["successfully_event_joined_occurrence_rows"]+=1
            loc,yr,dt=event_key_map[coreid]
            if yr not in YRS or loc not in sites:continue
            count["eligible_siteyear_occurrence_rows"]+=1
            name=row[of["scientificName"]].strip()
            status=row[of["occurrenceStatus"]].strip()
            occ_status[status or "(blank)"]+=1
            if status and status.lower() not in ("present","presence"):
                count["explicit_absent_or_other_status"]+=1
                continue
            n=floatcount(row[of["individualCount"]])
            if n is None:
                count["missing_zero_noninteger_or_invalid_individualCount"]+=1
                continue
            count["positive_individualCount_occurrence_rows"]+=1
            if name not in taxa:
                count["species_name_not_in_frozen_taxonomy"]+=1;continue
            fam=families.get(name)
            if fam not in FAMILIES:
                count["non_target_or_unresolved_family_positive_rows"]+=1
                if taxa[name].get("gbif_rank") in ("GENUS","FAMILY"):
                    count["family_unresolved_due_higher_rank"]+=1
                continue
            count["accepted_confident_butterfly_family_rows"]+=1
            stats[(loc,yr)][fam]+=n
            family_rows[fam]+=1
            family_counts[fam]+=n
    if count["all_event_rows"]!=15606 or count["all_occurrence_rows"]!=21617:
        raise ValueError("Original event/occurrence record totals changed")
    valid=count["successfully_event_joined_occurrence_rows"]/(count["all_occurrence_rows"]-count["invalid_short_occurrence_rows"])
    # Wrong/ambiguous source joins invalidate the outcome model, not treated as zeros.
    if valid<.95:count["sample_source_join_quality_GATE_FAILED"]=1
    panel=[]
    for loc,site in sites.items():
        for yr in YRS:
            visited=len(event_qualified_days[(loc,yr)])
            if visited<2 or visited!=site[f"visit_days_{yr}"]:
                count["siteyear_visit_days_event_source_mismatch"]+=1
            vals={f:int(stats[(loc,yr)].get(f,0)) for f in FAMILIES}
            tot=sum(vals[f] for f in FIRST)
            panel.append({
                "source_rounded5_site_hash":site["site_id_SHA256"],
                "original_geo_1deg_cell":site["country_1deg_lat_lon_cell"],
                "source_2022_nearest_qualified_WWTP_km":site["nearest_qualified_EEA2022_DK_FW_outfall_km"],
                "site_year":yr,"distinct_original_visit_days":visited,
                "source_eventrows_site_year":site[f"Event_rows_{yr}"],
                "family_individual_counts":vals,
                "target_3family_count":tot,
                "target_family_recorded":tot>0,
                "target_3family_count_per_visit":tot/visited if visited else None,
                "Lycaenidae_count_per_visit":vals["Lycaenidae"]/visited if visited else None,
                "no_individual_observer_or_occurrence_ids_exported":True
            })
    if len(panel)!=214:raise ValueError("Not exactly 107 repeated sites times two years")
    diagnostics={
        "source_rows":dict(count),
        "original_occurrence_status_distribution_on_eligible_sites":dict(occ_status),
        "accepted_family_original_rows":dict(family_rows),
        "accepted_family_original_individual_sums":dict(family_counts),
        "original_event_to_occurrence_link_fraction":valid,
        "distinct_valid_siteyear_2023_2024":len(panel),
        "qualified_site_year_target_family_recorded":sum(x["target_family_recorded"] for x in panel),
        "year_target_family_site_count":{
            str(y):sum(x["target_family_recorded"] for x in panel if x["site_year"]==y) for y in YRS
        },
        "year_target_family_individual_sums":{
            str(y):{
                f:sum(x["family_individual_counts"][f] for x in panel if x["site_year"]==y)
                for f in FIRST} for y in YRS
        },
        "visit_days_mismatch":count["siteyear_visit_days_event_source_mismatch"]
    }
    return panel,diagnostics

def stopgates(panel,diagnostics):
    # These are data-quality/observation gates only, no association effect is computed.
    if diagnostics["original_event_to_occurrence_link_fraction"]<0.95:
        return False,"EVENT_OCCURRENCE_JOIN_LT95PCT"
    if diagnostics["visit_days_mismatch"]>0:
        return False,"FROZEN_SITE_VISIT_DAYS_MISMATCH"
    if diagnostics["source_rows"].get("species_name_not_in_frozen_taxonomy",0):
        return False,"UNFROZEN_TAXON_NAME_PRESENT"
    per={}
    for yr in YRS:
        group={}
        for name,rule in (("near5",lambda d:d<=5),("remote10",lambda d:d>=10)):
            x=[r for r in panel if r["site_year"]==yr and rule(r["source_2022_nearest_qualified_WWTP_km"])]
            n=sum(r["target_family_recorded"] for r in x)
            group[name]={"all_source_qualified_sites":len(x),"target_3family_recorded_sites":n}
        per[str(yr)]=group
    enough=all(per[str(y)][name]["target_3family_recorded_sites"]>=5 for y in YRS for name in ("near5","remote10"))
    return enough,("FAMILY_RESPONSE_SUPPORTED" if enough else "INSUFFICIENT_TARGET_FAMILY_SITES_PER_GROUP_YEAR"),per

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    site,tax,clust,sites,taxa,fams,sha_files=load_sources(root)
    raw=download()
    panel,qc=make_panel(raw,sites,taxa,fams)
    gate=stopgates(panel,qc)
    if len(gate)==2:
        allowed,reason=gate
        byyr={}
    else:allowed,reason,byyr=gate
    qc["precommitted_family_outcome_basic_support_gate_pass"]=allowed
    qc["family_outcome_model_gate_reason"]=reason
    qc["family_outcome_group_year_support"]=byyr
    out={
        "schema":"gesui_DK_ebms_v14_2023_2024_individualCount_family_panel_SOURCELOCKED_quality_before_effect_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_licensed_archive":"eBMS Denmark v1.4 CC-BY-NC 4.0 GBIF caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c",
        "original_archive_SHA256":ORIGINAL_SHA,
        "exact_blinded_taxonomy_and_source_site_file_sha256":sha_files,
        "source_n_original_repeated_survey_sites":107,
        "source_near5_sites":23,"registered_source_remote10_sites":38,
        "original_independent_registered_WWTP_outfall_exposure_groups":4,
        "site_rows_not_independent_WWTP_exposure_replicates":True,
        "year_2023_and_2024_are_repeated_sites_not_independent_WWTP_treatments":True,
        "only_frozen_taxonomy_accepted_species_butterfly_families_counted":True,
        "unmatched_odonata_or_other_non_butterfly_taxa_not_mislabeled_as_failed_butterfly_identification":True,
        "original_taxon_WWTP_species_count_data_quality_QC":qc,
        "site_year_family_counts_original_public_CC_BY_NC_aggregate":panel,
        "pre_effect_model_family_count_observation_gate_pass":allowed,
        "pre_effect_model_gate_reason":reason,
        "no_near_vs_remote_family_effect_coefficient_computed_in_this_stage":True,
        "not_causal_WWTP_biochemical_or_puddling_result":True
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"event_occurrence_QC":qc,"family_model_may_run":allowed,"gate_reason":reason},indent=2))
if __name__=="__main__":main()
