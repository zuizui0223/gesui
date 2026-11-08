#!/usr/bin/env python3
"""Audit freely published national eBMS DwC-A event cores WITHOUT butterfly counts.

National archives: Sweden CC0, Denmark CC-BY-NC. Central eBMS data remain
controlled by coordinator approval; never access/rehost them.
No occurrence extension, species identities or count values are opened.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,re,urllib.request,urllib.error,zipfile
import xml.etree.ElementTree as ET
from datetime import datetime,timezone
from pathlib import Path

SOURCES={
 "Sweden_SeBMS":{
    "url":"https://www.gbif.se/ipt/archive.do?r=lu_sebms",
    "datasetKey":"be77e203-486c-4651-91b9-8347968b728c",
    "license":"CC0 1.0",
    "official_publisher":"Lund University, Sweden",
 },
 "Denmark_eBMS":{
    "url":"https://ipt.brc.ac.uk/archive.do?r=ebms_denmark_transects",
    "datasetKey":"caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c",
    "license":"CC-BY-NC 4.0",
    "official_publisher":"Biological Records Centre (UK CEH), Danish EBMS",
 }}
EVENT="http://rs.tdwg.org/dwc/terms/Event"
NAMESPACE="{http://rs.tdwg.org/dwc/text/}"
FIELDS=["eventID","parentEventID","locationID","eventDate","year","country",
    "countryCode","decimalLatitude","decimalLongitude",
    "coordinateUncertaintyInMeters","samplingProtocol","samplingEffort",
    "sampleSizeValue","sampleSizeUnit","eventTime"]
MAX_DOWNLOAD=80_000_000

def obtain(source,path):
    req=urllib.request.Request(source["url"],headers={
        "User-Agent":"gesui-eBMS-official-event-only-inventory/0.1",
        "Accept":"application/zip,*/*"})
    h=hashlib.sha256();size=0
    with urllib.request.urlopen(req,timeout=140) as r,path.open("wb") as f:
        if r.status!=200:raise ValueError("Official source archive HTTP failed")
        for b in iter(lambda:r.read(1048576),b""):
            size+=len(b)
            if size>MAX_DOWNLOAD:raise ValueError("Archive unexpectedly >80MB")
            h.update(b);f.write(b)
    if not zipfile.is_zipfile(path):raise ValueError("Source is not publicly available DwC-A ZIP")
    return {"download_byte_size":size,"download_sha256":h.hexdigest()}

def nonnull(v):
    return bool(v and str(v).strip() and str(v).lower() not in ("nan","null","na"))
def summarize(core, archive):
    fpath=core.find(".//"+NAMESPACE+"files/"+NAMESPACE+"location")
    rowtype=core.attrib.get("rowType","")
    terms={}
    for x in core.findall(NAMESPACE+"field"):
        i=x.attrib.get("index");term=x.attrib.get("term","").rsplit("/",1)[-1]
        if i is not None:terms[term]=int(i)
    id_field=core.find(NAMESPACE+"id")
    if id_field is not None and "eventID" not in terms:
        terms["eventID"]=int(id_field.attrib.get("index","0"))
    if rowtype!=EVENT:
        return {"rowType":rowtype,"event_core_compatible":False,
           "NO_SPECIES_OR_COUNTS_READ":True,
           "field_names":list(terms),"reason":"Core is not Event; do not read occurrence/species counts"},None
    if fpath is None or not fpath.text or fpath.text not in archive.namelist():
        raise ValueError("DwC-A Event core file location missing")
    enc=core.attrib.get("encoding","UTF-8")
    delim=core.attrib.get("fieldsTerminatedBy","\\t")
    if delim=="\\t":delim="\t"
    elif delim=="\\,":delim=","
    elif delim=="\\;":delim=";"
    enclosed=core.attrib.get("fieldsEnclosedBy",'"')
    skip=int(core.attrib.get("ignoreHeaderLines","1"))
    if skip<0 or skip>2:raise ValueError("Unexpected Event core header count")
    n=0;year=collections.Counter();country=collections.Counter()
    coords=collections.Counter();identified_events=set();ids=set()
    n2022=0;coord2022=0;good2022=0;visit2022=0
    methods=collections.Counter();valid_years=collections.Counter()
    maxcols=max(terms.values())+1 if terms else 0
    f=archive.open(fpath.text)
    with io.TextIOWrapper(f,encoding=enc,errors="replace",newline="") as stream:
        reader=csv.reader(stream,delimiter=delim,quotechar=enclosed or '"')
        for _ in range(skip):next(reader,None)
        for row in reader:
            if len(row)<maxcols:continue
            def field(key):return row[terms[key]].strip() if key in terms else ""
            n+=1
            date=field("eventDate");yr=field("year")
            if not yr and len(date)>=4 and date[:4].isdigit():yr=date[:4]
            year[yr or "MISSING"]+=1
            lat=field("decimalLatitude");lon=field("decimalLongitude")
            country[field("countryCode") or field("country") or "MISSING"]+=1
            method=field("samplingProtocol")
            if method:methods[method[:80]]+=1
            if yr=="2022":
                n2022+=1
                evt=field("eventID")
                if evt:identified_events.add(evt)
                if lat and lon:
                    coord2022+=1
                    try:
                        f_lat=float(lat);f_lon=float(lon)
                        if -90<=f_lat<=90 and -180<=f_lon<=180:
                            coord=(round(f_lat,5),round(f_lon,5))
                            coords[coord]+=1
                            if 44<=f_lat<=64 and -12<=f_lon<=26:good2022+=1
                    except ValueError:pass
                e=field("samplingEffort")
                if e:visit2022+=1
            if n>3_000_000:raise ValueError("Event core exceeds unexpectedly safe row cap")
    return {
       "rowType":rowtype,"event_core_compatible":True,
       "event_core_relative_path":fpath.text,
       "source_event_field_names":sorted(terms),
       "source_event_rows":n,
       "source_year_event_counts":dict(year),
       "source_country_code_counts":dict(country),
       "original_sampling_protocol_counts":dict(methods.most_common(15)),
       "2022_event_rows":n2022,
       "2022_unique_eventIDs":len(identified_events),
       "2022_event_rows_with_lat_lon":coord2022,
       "2022_independent_unique_coordinate_points":len(coords),
       "2022_events_in_Europe_study_bbox":good2022,
       "2022_events_with_explicit_samplingEffort_field":visit2022,
       "event_row_locations_not_redistributed":True,
       "NO_OCCURRENCE_SPECIES_DATA_OPENED":True,
       "NO_BUTTERFLY_COUNTS_OPENED":True,
    },None

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--tmpdir",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    a=parser.parse_args()
    a.tmpdir.mkdir(parents=True,exist_ok=True)
    out={"schema":"gesui_eBMS_official_open_national_sampling_event_source_only_QC_v01",
        "date_utc":datetime.now(timezone.utc).isoformat(),
        "source_method":"original IPT archive.do official public national DwC-A only",
        "no_occurrence_taxa_or_count_data_opened":True,"source_reports":{}}
    for name,info in SOURCES.items():
        file=a.tmpdir/(name+".zip")
        data={"gbif_datasetKey":info["datasetKey"],
              "publication_license":info["license"],
              "publisher":info["official_publisher"],
              "original_publisher_URL":info["url"]}
        try:
            data["original_archive"]=obtain(info,file)
            with zipfile.ZipFile(file) as z:
                meta=next((x for x in z.namelist() if x.endswith("meta.xml")),None)
                if meta is None:raise ValueError("No DwC-A meta.xml available")
                mb=z.read(meta);data["meta_xml_SHA256"]=hashlib.sha256(mb).hexdigest()
                root=ET.fromstring(mb)
                core=root.find(NAMESPACE+"core")
                data["ZIP_original_entry_names"]=[x for x in z.namelist() if not ("occurrence" in x.lower() or "taxon" in x.lower())][:80]
                if core is None:raise ValueError("No DwC-A core")
                qc,_=summarize(core,z)
                data["event_source_only_quality"]=qc
        except Exception as ex:
            data["unavailable_source_reason"]=type(ex).__name__+": "+str(ex)
        out["source_reports"][name]=data
        summary={"dataset":name,"downloaded":bool(data.get("original_archive")),
                 "eventQC":data.get("event_source_only_quality",{}).get("event_core_compatible"),
                 "events2022":data.get("event_source_only_quality",{}).get("2022_event_rows"),
                 "error":data.get("unavailable_source_reason")}
        print("GESUI_EBMS_SOURCE_ONLY_PUBLIC_EVENT_SOURCE_GATES",json.dumps(summary,sort_keys=True),flush=True)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_EBMS_NATIONAL_SOURCE_ONLY_PUBLIC_ARCHIVE_QC_COMPLETE",flush=True)
if __name__=="__main__":main()
