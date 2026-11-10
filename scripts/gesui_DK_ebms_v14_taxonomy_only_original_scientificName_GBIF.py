#!/usr/bin/env python3
"""Denmark eBMS v1.4: freeze original scientificName -> GBIF family taxonomy ONLY.

Original individualCount, occurrenceStatus, site, year and exposure values are
not read in this stage; no source-near versus remote butterfly response computed.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,time,urllib.parse,urllib.request,zipfile
import xml.etree.ElementTree as ET
from datetime import datetime,timezone
from pathlib import Path
SOURCE="https://ipt.brc.ac.uk/archive.do?r=ebms_denmark_transects&v=1.4"
SHA="2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b"
NS="{http://rs.tdwg.org/dwc/text/}"
EXPECTED_DATASET="caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c"
TARGET=("Lycaenidae","Pieridae","Papilionidae")
SECONDARY=("Nymphalidae","Hesperiidae","Riodinidae")
VALID_FAMILY=set(TARGET+SECONDARY)
def sha(data):return hashlib.sha256(data).hexdigest()
def get(url):
    last=None
    for attempt in range(4):
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"gesui-DK-eBMS-v14-scientificName-GBIF-family-QA/0.1",
                "Accept":"application/json,application/zip"})
            with urllib.request.urlopen(req,timeout=90) as r:
                return r.read()
        except Exception as e:
            last=str(e)
            if attempt<3:time.sleep(2*(attempt+1))
    raise RuntimeError("Original source unavailable: "+str(last))
def terms(node):
    return {f.attrib["term"].rsplit("/",1)[-1]:int(f.attrib["index"]) for f in node.findall(NS+"field")}
def filenames(node):
    loc=node.find(".//"+NS+"files/"+NS+"location")
    if loc is None or not loc.text:raise ValueError("Missing extension file")
    return loc.text
def values(z,node):
    delim=node.attrib.get("fieldsTerminatedBy",r"\t")
    if delim==r"\t":delim="\t"
    if delim!="\t":raise ValueError("Unexpected source field delimiter")
    skip=int(node.attrib.get("ignoreHeaderLines","1"))
    with io.TextIOWrapper(z.open(filenames(node)),encoding="utf-8",newline="") as stream:
        reader=csv.reader(stream,delimiter=delim)
        for _ in range(skip):next(reader,None)
        for row in reader:yield row
def match_species(name):
    query=urllib.parse.urlencode({"name":name,"strict":"true"})
    url="https://api.gbif.org/v1/species/match?"+query
    err=None
    for attempt in range(4):
        try:
            raw=get(url)
            d=json.loads(raw)
            confidence=d.get("confidence",0)
            family=d.get("family")
            matchtype=d.get("matchType")
            rank=d.get("rank")
            accepted=bool(family in VALID_FAMILY and matchtype in ("EXACT","FUZZY") and
                          confidence>=90 and rank in ("SPECIES","SUBSPECIES","VARIETY"))
            return {"original_scientificName":name,"source_GBIF_name_match_URL":url,
                    "source_GBIF_JSON_SHA256":sha(raw),"gbif_usageKey":d.get("usageKey"),
                    "gbif_rank":rank,"gbif_family":family,"gbif_matchType":matchtype,
                    "gbif_confidence":confidence,"authoritative_family_accepted":accepted,
                    "gbif_status":d.get("status"),"gbif_scientificName":d.get("scientificName")}
        except Exception as e:
            err=str(e)
            if attempt<3:time.sleep(3*(attempt+1))
    return {"original_scientificName":name,"family_match_lookup_failed":True,"error":err,
            "authoritative_family_accepted":False}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    src=get(SOURCE)
    if sha(src)!=SHA:raise ValueError("Public original eBMS archive SHA changed; STOP")
    b=io.BytesIO(src)
    with zipfile.ZipFile(b) as z:
        root=ET.fromstring(z.read("meta.xml"))
        ex=[x for x in root.findall(NS+"extension") if x.attrib.get("rowType")=="http://rs.tdwg.org/dwc/terms/Occurrence"]
        if len(ex)!=1:raise ValueError("Original Occurrence extension absent or nonunique")
        node=ex[0];cols=terms(node)
        if not {"scientificName","taxonID"}.issubset(cols):raise ValueError("Missing taxonomy source columns")
        keys=set();row_n=0;empty_n=0
        for row in values(z,node):
            row_n+=1
            if len(row)<=max(cols["scientificName"],cols["taxonID"]):
                raise ValueError("Original species source row shorter than taxonomy fields")
            n=row[cols["scientificName"]].strip()
            if not n:empty_n+=1;continue
            keys.add(n)
        # Never access individualCount / eventID / coordinates / occurrenceStatus in this stage.
    matches=[]
    for k,name in enumerate(sorted(keys),1):
        result=match_species(name)
        matches.append(result)
        if k%20==0:print(f"Taxonomy-only {k}/{len(keys)} original scientific names",flush=True)
        time.sleep(0.15)
    out={"schema":"gesui_DK_eBMS_v14_original_2023_2024_taxon_name_GBIF_family_before_count_contrasts_v01",
         "utc":datetime.now(timezone.utc).isoformat(),
         "original_GBIF_dataset":EXPECTED_DATASET,
         "original_eBMS_DwCA_SHA256":SHA,"original_metaXML_SHA256":sha(zipfile.ZipFile(io.BytesIO(src)).read("meta.xml")),
         "original_license":"CC-BY-NC 4.0",
         "original_Occurrence_extension_species_values_read_only_for_taxonomy":True,
         "original_individualCount_values_NOT_READ":True,
         "original_occurrenceStatus_or_locality_or_year_values_NOT_READ":True,
         "original_EEA_WWTP_near_or_remote_exposure_NOT_ACCESSED":True,
         "original_occurrence_row_n":row_n,"empty_scientificName_n":empty_n,
         "unique_original_scientificName_n":len(keys),
         "strict_matching_rules":{"families_primary":list(TARGET),"secondary":list(SECONDARY),
             "GBIF_species_match_confidence_min":90,"matchType_allowed":["EXACT","FUZZY"],
             "ranks_allowed":["SPECIES","SUBSPECIES","VARIETY"],
             "unresolved_genus_or_species_group_not_assigned_to_family":True},
         "original_scientificName_GBIF_taxonomic_source_match_audit":matches,
         "accepted_family_mapping_counts_by_family":dict(collections.Counter(x.get("gbif_family") for x in matches if x["authoritative_family_accepted"])),
         "source_taxon_classification_failures_n":sum(not x["authoritative_family_accepted"] for x in matches),
         "no_butterfly_near_vs_remote_count_effect_calculated":True}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print(json.dumps({"taxa":len(keys),"accepted":out["accepted_family_mapping_counts_by_family"],
                      "unresolved":out["source_taxon_classification_failures_n"],"no_outcome_values":True}))
if __name__=="__main__":main()
