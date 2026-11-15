#!/usr/bin/env python3
"""Source-only GBIF SeBMS v1.13 2022 event-core feasibility.

NEVER touch occurrence.txt, species names/keys, HNS trait scores, source-treatment
geography, or butterfly response. Only examine Event-core schema, 2022 dates,
sampling-event protocol, coordinates and uncertainty. A centroid of a 0.5-3km
transect does not locate individual butterflies relative to a river.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import io
import json
import re
import time
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

DATASET="be77e203-486c-4651-91b9-8347968b728c"
VERSION="1.13"
URL="https://www.gbif.se/ipt/archive.do?r=lu_sebms&v=1.13"
FROZEN_PROTOCOL_GIT_BLOB="b9ce415c36d398eee77d8318a20aef4873e0bf44"
EXPECTED_SOURCE_EVENT_ROWS=39466
MAX_ARCHIVE_BYTES=50*1024*1024
MAX_SITE_TRUNCATION_SAMPLE=0


def git_blob(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()


def get_source_url(url):
    if url!=URL:
        raise ValueError("Only frozen canonical source archive URL allowed")
    err=None
    for attempt in range(3):
        try:
            req=urllib.request.Request(url,headers={
              "User-Agent":"gesui-sourceonly-20261115/1.0 (academic reproducibility; no occurrence rows)",
              "Accept":"application/zip,application/octet-stream,*/*",
            })
            with urllib.request.urlopen(req,timeout=65) as resp:
                mime=resp.headers.get("Content-Type","")
                buf=resp.read(MAX_ARCHIVE_BYTES+1)
                if len(buf)>MAX_ARCHIVE_BYTES:raise ValueError("Archive >50MB, STOP")
                if not zipfile.is_zipfile(io.BytesIO(buf)):
                    raise RuntimeError(f"Official source response not DWCA ZIP (HTTP {resp.status}; content-type {mime}; first bytes={buf[:20]!r})")
                return buf,mime
        except Exception as exc:
            err=str(exc)
            if attempt<2:time.sleep(2*(attempt+1))
    raise RuntimeError("Official latest SeBMS v1.13 source archive inaccessible: "+str(err))


def attr(x,name,default=""):
    return x.attrib.get(name,default)


def delimit(raw):
    # Darwin Core meta.xml uses literal strings '	tab','\t', ',' etc
    if raw in [r"\t","\\t","tab","TAB"]:return "\t"
    if raw in [r"\n",r"\r\n"]:return "\n"
    if raw=='': return ','
    return raw


def parse_event_core(zipdata):
    with zipfile.ZipFile(io.BytesIO(zipdata)) as z:
        names=z.namelist()
        meta=next((n for n in names if n.rsplit("/",1)[-1]=="meta.xml"),None)
        if meta is None:raise ValueError("No Darwin Core meta.xml found")
        root=ET.fromstring(z.read(meta))
        ns={"a":"http://rs.tdwg.org/dwc/text/"}
        core=root.find("a:core",ns)
        if core is None:raise ValueError("No core table in Darwin Core Archive")
        if not attr(core,"rowType").lower().endswith("/event"):
            raise ValueError("GBIF source core not Event")
        loc=core.find("a:files/a:location",ns)
        if loc is None or not (loc.text or "").strip():raise ValueError("No core location")
        locname=loc.text.strip()
        base=meta.rsplit("/",1)[0] if "/" in meta else ""
        eventname=locname if locname in names else (base+"/"+locname).lstrip("/")
        if eventname not in names:raise ValueError("Missing event core member "+locname)
        enc=attr(core,"encoding","utf-8")
        sep=delimit(attr(core,"fieldsTerminatedBy",r"\t"))
        q=attr(core,"fieldsEnclosedBy",'"')
        header=int(attr(core,"ignoreHeaderLines","0"))
        idx={}
        for field in core.findall("a:field",ns):
            name=attr(field,"term").rsplit("/",1)[-1].rsplit("#",1)[-1]
            idx[name]=int(attr(field,"index"))
        id_elem=core.find("a:id",ns)
        if id_elem is not None:idx.setdefault("_coreID",int(attr(id_elem,"index")))
        if "eventDate" not in idx:raise ValueError("No real eventDate in event core")
        if "decimalLatitude" not in idx or "decimalLongitude" not in idx:
            raise ValueError("Event core lacks decimal coordinate columns")
        if "eventID" not in idx and "_coreID" not in idx:
            raise ValueError("Event core missing event ID")
        # Only open event core, NEVER read any occurrence extension or species trait table.
        raw=z.read(eventname)
        tf=io.TextIOWrapper(io.BytesIO(raw),encoding=enc,errors="strict",newline="")
        rdr=csv.reader(tf,delimiter=sep,quotechar=q if len(q)==1 else None)
        for _ in range(header):next(rdr)
        records=[]
        for line in rdr:
            records.append({k:line[i] if i<len(line) else "" for k,i in idx.items()})
        return records,{"member":eventname,"meta_xml_path":meta,"encoding":enc,
                        "field_names":sorted(idx),"n_original_zip_members":len(names),
                        "file_sha256":hashlib.sha256(raw).hexdigest(),
                        "event_core_only":True}


def num(v):
    try:return float(v)
    except (ValueError,TypeError):return None


def datekey(s):
    s=(s or "").strip()
    m=re.match(r"^(20\d\d)-(\d\d)-(\d\d)",s)
    if not m:return None
    y,mo,d=map(int,m.groups())
    try:
        import datetime
        return datetime.date(y,mo,d)
    except ValueError:return None


def audit_events(events,meta,archive_sha,mime):
    if len(events)!=EXPECTED_SOURCE_EVENT_ROWS:
        raise ValueError(f"Source event core count {len(events)} differs from published 2026 v1.13 {EXPECTED_SOURCE_EVENT_ROWS}; cannot silently mix versions")
    selected=[]
    all_2022=0
    all_years=collections.Counter()
    for v in events:
        d=datekey(v.get("eventDate",""))
        if d is None:continue
        all_years[str(d.year)]+=1
        if d.year==2022:
            all_2022+=1
            if 5<=d.month<=9:selected.append((v,d))
    protocols=collections.Counter()
    quality=collections.Counter()
    uncertain_values=[]
    sites=collections.defaultdict(lambda:{"date_set":set(),"events":0,"protocols":set()})
    for v,d in selected:
        prot=(v.get("samplingProtocol") or "").strip() or "__MISSING__"
        protocols[prot]+=1
        eid=(v.get("eventID") or v.get("_coreID") or "").strip()
        if not eid:quality["missing_event_id"]+=1
        lat=num(v.get("decimalLatitude"))
        lon=num(v.get("decimalLongitude"))
        u=num(v.get("coordinateUncertaintyInMeters"))
        if lat is None or lon is None or not (-90<=lat<=90 and -180<=lon<=180):
            quality["coords_missing_invalid"]+=1
            continue
        quality["coords_valid"]+=1
        if u is None:quality["uncertainty_missing"]+=1
        elif u<0:quality["uncertainty_invalid_negative"]+=1
        else:
            quality["uncertainty_known_nonnegative"]+=1
            uncertain_values.append(u)
            if u<=100:quality["uncertainty_le100m"]+=1
            else:quality["uncertainty_gt100m"]+=1
        # Sites identified from original locationID if populated; fallback to exact coordinates
        site=(v.get("locationID") or "").strip() or (f"coords:{lat:.8f}|{lon:.8f}")
        sites[site]["date_set"].add(d.isoformat())
        sites[site]["events"]+=1
        sites[site]["protocols"].add(prot)
        if (v.get("footprintWKT") or "").strip():
            quality["footprintWKT_present"]+=1
        else:quality["footprintWKT_missing"]+=1
        if (v.get("samplingEffort") or "").strip(): quality["samplingEffort_present"]+=1
        if (v.get("samplingSizeValue") or "").strip(): quality["samplingSizeValue_present"]+=1
        if (v.get("eventTime") or "").strip(): quality["eventTime_present"]+=1
    eligible_sites_5=sum(len(v["date_set"])>=5 for v in sites.values())
    # When event-site geometry is a whole transect centroid, the site is not a butterfly point
    return {
      "schema":"GESUI_SEBMS_2026_113_2022_EVENT_CORE_SOURCE_ONLY_V1",
      "status":"2022_EVENT_OPPORTUNITY_ONLY_NOT_RIPARIAN_8CELL_QUALIFICATION",
      "source_dataset":DATASET,
      "source_version":VERSION,
      "publisher_archive_url":URL,
      "official_source_archive_sha256":archive_sha,
      "downloaded_mime":mime,
      "source_meta":meta,
      "n_all_event_rows":len(events),
      "n_events_with_parsed_date_by_year":dict(sorted(all_years.items())),
      "n_2022_full_year_event_rows":all_2022,
      "n_2022_May_to_September_event_rows":len(selected),
      "original_summer_event_samplingProtocol_counts":dict(sorted(protocols.items())),
      "original_summer_event_source_field_status":dict(sorted(quality.items())),
      "summer_distinct_event_site_ids_with_valid_coords":len(sites),
      "summer_sites_at_least_5_distinct_event_dates":eligible_sites_5,
      "coord_uncertainty_numeric_min_max":(
          [min(uncertain_values),max(uncertain_values)] if uncertain_values else None),
      "date_2022_May_to_September_in_frozen_window":True,
      "site_id_fallback_coordinates_not_proof_location_unit":True,
      "transect_centroids_do_not_locate_species":True,
      "point_count_25m_sampling_radius_not_identical_to_uncertainty":True,
      "source_event_core_only_occurrence_extension_not_read":True,
      "no_taxon_or_species_scores_read":True,
      "no_lycaenid_HNS_calculated":True,
      "no_river_directions_matched":True,
      "no_outfall_exposure_assigned":True,
      "no_8cell_source_reference_pairs_verified":True,
      "frozen_primary_outcome_test_ran":False,
      "frozen_primary_stage":"STOP_NOT_YET_SOURCE_QUALIFIED",
      "source_metadata_2026_version_not_assumed_same_as_2025_index":True
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--freeze",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    if git_blob(a.freeze.read_bytes())!=FROZEN_PROTOCOL_GIT_BLOB:
        raise ValueError("Frozen 2026-11-15 protocol Git blob changed")
    blob,mime=get_source_url(URL)
    events,meta=parse_event_core(blob)
    result=audit_events(events,meta,hashlib.sha256(blob).hexdigest(),mime)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print("GESUI_SEBMS_SOURCE_ONLY="+json.dumps(result,sort_keys=True,ensure_ascii=False,separators=(",",":")),flush=True)


if __name__=="__main__": main()
