#!/usr/bin/env python3
"""BLIND species-to-family taxonomy for Danish 15min EBMS v1.3.

Only Occurrence extension scientificName is read. No original individualCount,
Occurrence coreID, event location/date, source-near counts or WWTP effect read.
Freeze mapping BEFORE the new counts model. Original GBIF family authority.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,time,urllib.parse,urllib.request,zipfile
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

FOLDER="results/ebms_denmark_15minute_2023_2025_WWTP_SOURCE_ONLY/"
SOURCE=FOLDER+"ORIGINAL_DK_EBMS15M_2023_25_EEA2022_203_OUTFALL_EVENT_ONLY_QA.json"
OLD_TAX="results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_EBMS_V14_ORIGINAL_SCIENTIFICNAME_GBIF_FAMILY_BLIND_TAXONOMY.json"
ZIP="https://ipt.brc.ac.uk/archive.do?r=ebms_denmark_15-minute_counts&v=1.3"
PINNED="c32be29200134fd2ee06b7dd34d76e1420bd6a7810574fb8ffb345c1451169cb"
FAMILIES=("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae","Riodinidae")
NS="{http://rs.tdwg.org/dwc/text/}"

def sha(b):return hashlib.sha256(b).hexdigest()

def download():
    with urllib.request.urlopen(urllib.request.Request(ZIP,headers={"User-Agent":"gesui-DK-ebms15m-source-blind-taxonomy-20261110/0.1"}),timeout=140) as f:
        raw=f.read()
    if sha(raw)!=PINNED:raise ValueError("Original published EBMS15m archive changed")
    return raw

def unique_names(raw):
    counter=Counter()
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        meta=ET.fromstring(z.read("meta.xml"))
        ext=[x for x in meta.findall(NS+"extension") if x.attrib.get("rowType")=="http://rs.tdwg.org/dwc/terms/Occurrence"]
        if len(ext)!=1:raise ValueError("Expected exactly one original Occurrence extension")
        e=ext[0]
        terms={x.attrib["term"].rsplit("/",1)[-1]:int(x.attrib["index"]) for x in e.findall(NS+"field")}
        if "scientificName" not in terms:raise ValueError("No scientificName source")
        filename=e.find(".//"+NS+"files/"+NS+"location")
        if filename is None or not filename.text:raise ValueError("No Occurrence table")
        sep=e.attrib.get("fieldsTerminatedBy",r"\t")
        if sep==r"\t":sep="\t"
        if sep!="\t":raise ValueError("Unexpected delimiter")
        skip=int(e.attrib.get("ignoreHeaderLines","1"))
        with io.TextIOWrapper(z.open(filename.text),encoding="utf-8",newline="") as f:
            reader=csv.reader(f,delimiter=sep)
            for _ in range(skip):next(reader,None)
            for row in reader:
                counter["total_rows_IN_OCCURRENCE_EXT"]=counter["total_rows_IN_OCCURRENCE_EXT"]+1
                name=row[terms["scientificName"]].strip() if len(row)>terms["scientificName"] else ""
                if not name:counter["blank_scientificName_rows"]+=1
                else:counter["unique_name_source_only_"+name]=1
    if counter["total_rows_IN_OCCURRENCE_EXT"]!=39905:
        raise ValueError("Original published source 39905 records mismatch")
    names=sorted(x[len("unique_name_source_only_"):] for x in counter if x.startswith("unique_name_source_only_"))
    return names,{"original_occurrence_extension_rows":counter["total_rows_IN_OCCURRENCE_EXT"],
                  "missing_original_scientificName":counter["blank_scientificName_rows"],
                  "unique_original_species_or_higher_names":len(names)}

def gbif(name):
    url="https://api.gbif.org/v1/species/match?"+urllib.parse.urlencode({"name":name,"strict":"true"})
    e=None
    for attempt in range(4):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"gesui-DK-ebms15m-taxonomy-only/0.1"})
            with urllib.request.urlopen(req,timeout=30) as f: b=f.read()
            x=json.loads(b)
            return {
                "gbif_match_URL":url,
                "GBIF_original_response_SHA256":sha(b),
                "original_GBIF_status":x.get("status"),
                "original_GBIF_matchType":x.get("matchType"),
                "original_GBIF_rank":x.get("rank"),
                "original_GBIF_family":x.get("family"),
                "original_GBIF_usageKey":x.get("usageKey"),
                "original_GBIF_confidence":x.get("confidence"),
                "match_status":"GBIF_FROZEN",
            }
        except Exception as z:
            e=str(z)
            if attempt<3:time.sleep((attempt+1)*4)
    return {"match_status":"GBIF_FAILED","error":e,"gbif_match_URL":url}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    sourcepath=root/SOURCE
    audit=json.loads(sourcepath.read_text())
    if audit["source_archive_receipt"]["original_ZIP_sha256"]!=PINNED:
        raise ValueError("Original Event-only source hash not frozen")
    if not audit["strict_500m_primary_5km_butterfly_count_access_gate"]:
        raise ValueError("STOP source gate did not pass; do not open original Occurrence scientificName")
    assert audit["original_Occurrence_extension_READ_in_this_stage"] is False
    oldpath=root/OLD_TAX
    original=json.loads(oldpath.read_text())
    if not original["original_individualCount_values_NOT_READ"]:raise ValueError("Prior taxonomy was count contaminated")
    names,meta=unique_names(download())
    priornames={x["original_scientificName"]:x for x in original["original_scientificName_GBIF_taxonomic_source_match_audit"]}
    result=[]
    for k,name in enumerate(names,1):
        prev=priornames.get(name)
        if prev is not None:
            entry={
                "original_scientificName":name,
                "GBIF_original_response_SHA256":prev.get("source_GBIF_JSON_SHA256"),
                "gbif_match_URL":prev.get("source_GBIF_name_match_URL"),
                "original_GBIF_status":prev.get("gbif_status"),
                "original_GBIF_matchType":prev.get("gbif_matchType"),
                "original_GBIF_rank":prev.get("gbif_rank"),
                "original_GBIF_family":prev.get("gbif_family"),
                "original_GBIF_usageKey":prev.get("gbif_usageKey"),
                "original_GBIF_confidence":prev.get("gbif_confidence"),
                "match_status":"REUSED_PRIOR_FROZEN_BLIND_GBIF",
            }
        else:
            entry={"original_scientificName":name,**gbif(name)}
            time.sleep(0.3)
        entry["accepted_strict_target_butterfly_species"] = bool(
            entry["original_GBIF_family"] in FAMILIES
            and entry["original_GBIF_rank"]=="SPECIES"
            and entry["original_GBIF_status"]=="ACCEPTED"
            and entry["original_GBIF_matchType"]=="EXACT"
            and isinstance(entry["original_GBIF_confidence"],(int,float))
            and entry["original_GBIF_confidence"]>=95
        )
        result.append(entry)
    if len(result)!=meta["unique_original_species_or_higher_names"]:
        raise ValueError("Species name mapping not complete")
    byfamily=Counter(x["original_GBIF_family"] for x in result if x["accepted_strict_target_butterfly_species"])
    out={
        "schema":"gesui_2026_DK_EBMS15m_original_species_blind_GBIF_family_taxonomy_before_individualCount_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_source_archive_SHA256":PINNED,
        "original_Event_only_source_gate_file_SHA256":sha(sourcepath.read_bytes()),
        "original_prior_taxonomy_SHA256":sha(oldpath.read_bytes()),
        "original_Occurrence_only_scientificName_unique_names_QC":meta,
        "no_original_individualCount_values_read":True,
        "no_original_event_coordinates_or_WWTP_proximity_or_species_count_effect_READ_in_this_stage":True,
        "original_scientificName_GBIF_family_mapping_frozen":result,
        "target_taxonomic_families_source_only":list(FAMILIES),
        "source_accepted_species_by_family":dict(byfamily),
        "unresolved_names":sum(not x["accepted_strict_target_butterfly_species"] for x in result),
        "GBIF_source_lookup_failures":sum(x["match_status"]=="GBIF_FAILED" for x in result),
        "original_15min_count_data_not_read_before_source_taxonomy_freeze":True
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"source":meta,"accepted_species_by_family":dict(byfamily),
                      "failed_GBIF":out["GBIF_source_lookup_failures"],
                      "original_Occurrence_individualCount_values_read":False},indent=2))

if __name__=="__main__":
    main()
