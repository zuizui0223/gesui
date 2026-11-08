#!/usr/bin/env python3
"""Source-only DwC Event core 2022 coverage for open Norway and Flanders butterfly schemes.

Reads ONLY publisher-original public archive meta.xml + Event core, NEVER
Occurrence/measurement/species/individualCount, no EEA or butterfly outcomes.
Reports independent-site confidence separately from GPS/source support.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,re,urllib.request,zipfile
from collections import Counter,defaultdict
from pathlib import Path
from datetime import datetime, timezone
from xml.etree import ElementTree as ET

SOURCES={
 "Norway_NINA":{
   "url":"https://ipt.nina.no/archive.do?r=butterflies_bumblebees2020",
   "uuid":"aea17af8-5578-4b04-b5d3-7adf0c5a1e60",
   "license":"CC BY 4.0",
   "country":"NO"},
 "Flanders_Meetnetten":{
   "url":"https://ipt.inbo.be/archive.do?r=meetnetten-butterflies-occurrences",
   "uuid":"63938753-1fec-4c08-ae39-e9f8a6576521",
   "license":"CC0 1.0",
   "country":"BE"},
 "Flanders_historical_1991_2020":{
   "url":"https://ipt.inbo.be/archive.do?r=bms-occurrences",
   "uuid":"81b909c7-5612-4148-8e70-99fea2cdb572",
   "license":"CC BY-NC 4.0",
   "country":"BE"},
}
NAMESPACE="{http://rs.tdwg.org/dwc/text/}"
YEAR="2022"
MAX_DOWNLOAD=140_000_000
MAX_EVENTS=1_500_000
MAX_ACCURACY_M=500.0
NEED_VISITS=2

def hash_file(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()

def read_zip(info,path):
    request=urllib.request.Request(info["url"],headers={
      "User-Agent":"gesui-independent-2022-Norway-Flanders-public-Event-only-QC/1.0",
      "Accept":"application/zip,*/*"})
    n=0
    with urllib.request.urlopen(request,timeout=180) as r,path.open("wb") as f:
        if r.status!=200:raise ValueError("Official source returned HTTP "+str(r.status))
        for data in iter(lambda:r.read(1024*1024),b""):
            n+=len(data)
            if n>MAX_DOWNLOAD:raise ValueError("Archive exceeds source safety bound")
            f.write(data)
    if not zipfile.is_zipfile(path):raise ValueError("Publisher archive not ZIP")
    return {"source_url":info["url"],"download_bytes":n,"archive_sha256":hash_file(path)}

def valid_coord(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except (TypeError,ValueError):return None

def delim(v):
    if v in (r"\t",r"\\t","&#9;","&#x09;"):return "\t"
    if v in (r"\,",r"\\,"):return ","
    if v in (r"\;",r"\\;"):return ";"
    return v or "\t"

def event_archive(path,info):
    result={"original_zip_sha256":hash_file(path),"dataset_key":info["uuid"],
             "licence":info["license"],"source_country":info["country"],
             "no_2022_butterfly_taxa_or_counts_read":True,
             "strict_2022_site_coordinate_uncertainty_max_m":MAX_ACCURACY_M,
             "strict_independent_site_min_distinct_visits":NEED_VISITS}
    with zipfile.ZipFile(path) as z:
        metadata=next((n for n in z.namelist() if n.endswith("meta.xml")),None)
        if metadata is None:raise ValueError("DwC-A no meta.xml")
        mb=z.read(metadata);result["meta_sha256"]=hashlib.sha256(mb).hexdigest()
        root=ET.fromstring(mb)
        core=root.find(NAMESPACE+"core")
        if core is None or core.get("rowType")!="http://rs.tdwg.org/dwc/terms/Event":
            result["event_core_available"]=False
            result["reason"]="Original DwC core is not Event; no occurrence table read"
            return result,[]
        result["event_core_available"]=True
        location=core.find(".//"+NAMESPACE+"files/"+NAMESPACE+"location")
        if location is None or not location.text:raise ValueError("Event core no file")
        source_name=location.text
        if source_name not in z.namelist():
            maybe=[n for n in z.namelist() if n.endswith("/"+source_name)]
            if len(maybe)!=1:raise ValueError("Cannot find Event core file")
            source_name=maybe[0]
        if any(v in source_name.lower() for v in ("occurrence","taxon","measurementorfact")):
            raise ValueError("Refuse to open butterfly outcome source")
        result["event_only_core_path"]=source_name
        fields={}
        for fld in core.findall(NAMESPACE+"field"):
            ix=fld.get("index");term=fld.get("term","").rsplit("/",1)[-1]
            if ix is not None: fields[term]=int(ix)
        keytag=core.find(NAMESPACE+"id")
        if keytag is not None and "eventID" not in fields:
            fields["eventID"]=int(keytag.get("index","0"))
        result["original_event_field_names"]=sorted(fields)
        result["original_dwc_extensions_present_BUT_NOT_OPENED"]=[
            (e.get("rowType","").rsplit("/",1)[-1]) for e in root.findall(NAMESPACE+"extension")]
        sep=delim(core.get("fieldsTerminatedBy",r"\t"))
        enc=core.get("encoding","utf-8")
        quote=core.get("fieldsEnclosedBy",'"')
        header=int(core.get("ignoreHeaderLines","1"))
        if sep not in ("\t",",",";","|"):raise ValueError("Unsupported Event delimiter "+repr(sep))
        if header<0 or header>2:raise ValueError("Unexpected Event header count")
        maxfield=max(fields.values()) if fields else 0
        counts=Counter();grouped=defaultdict(list); site_formats=Counter()
        sample_events=[]
        with io.TextIOWrapper(z.open(source_name),encoding=enc,errors="replace",newline="") as f:
            rd=csv.reader(f,delimiter=sep,quotechar=quote if len(quote)==1 else '"')
            for _ in range(header): next(rd,None)
            for row in rd:
                counts["all_events"]+=1
                if counts["all_events"]>MAX_EVENTS:raise ValueError("Unexpected enormous Event core")
                if len(row)<=maxfield:
                    counts["short_row"]+=1;continue
                def fld(k):
                    i=fields.get(k)
                    return row[i].strip() if i is not None and i<len(row) else ""
                date=fld("eventDate")
                yr=fld("year")
                if not yr and len(date)>=4 and date[:4].isdigit():yr=date[:4]
                counts["year_"+str(yr or "unknown")]+=1
                if yr!=YEAR:continue
                counts["events_targetyear"]+=1
                valid_date=bool(re.match(r"^2022-\d{2}-\d{2}$",date))
                if not valid_date:
                    counts["non_dayresolved_targetyear_event"]+=1
                    continue
                counts["dated_events_targetyear"]+=1
                lat=valid_coord(fld("decimalLatitude"));lon=valid_coord(fld("decimalLongitude"))
                if lat is None or lon is None or not(-90<=lat<=90 and -180<=lon<=180):
                    counts["no_finite_global_coord"]+=1;continue
                counts["2022_events_geo"]+=1
                country=fld("countryCode") or info["country"]
                if country!=info["country"]:
                    counts["mismatched_country_event"]+=1;continue
                acc=valid_coord(fld("coordinateUncertaintyInMeters"))
                if acc is None:counts["missing_coordinate_uncertainty_2022"]+=1
                elif acc>MAX_ACCURACY_M:counts["over_500m_coordinate_uncertainty_2022"]+=1
                else:counts["coordinate_uncertainty_le_500m_2022"]+=1
                protocol=fld("samplingProtocol")
                group=(
                  fld("locationID") or fld("samplingLocationID") or
                  fld("locality") or fld("parentEventID"))
                if group:
                    kind="original_site_metadata_key"
                    sourcekey=group+"|"+protocol
                else:
                    # Provisional sites cannot imply independent Lucas grids.
                    kind="provisional_public_coord_plus_protocol"
                    sourcekey=f"PROVISIONAL:{lat:.4f}:{lon:.4f}|{protocol}"
                site_formats[kind]+=1
                ev={
                  "date":date,
                  "lat":lat,"lon":lon,"coord_uncertainty_m":acc,
                  "protocol":protocol,"parentEventID":fld("parentEventID"),
                  "original_locationID_present":bool(fld("locationID")),
                  "original_locality_present":bool(fld("locality")),
                  "original_parentID_present":bool(fld("parentEventID")),
                  "site_key":sourcekey,"site_key_kind":kind,
                  "event_id_sha256":hashlib.sha256(fld("eventID").encode()).hexdigest() if fld("eventID") else None,
                  "country":country
                }
                grouped[sourcekey].append(ev)
                if len(sample_events)<4:
                    sample_events.append({k:v for k,v in ev.items() if k in (
                      "protocol","coord_uncertainty_m","date","site_key_kind",
                      "original_locationID_present","original_locality_present","original_parentID_present")})
        result["event_counts"]=dict(counts)
        result["site_key_kind_counts"]=dict(site_formats)
        result["first_four_2022_event_metadata_examples_NO_BUTTERFLY"]=sample_events
        result["2022_distinct_source_derived_group_n"]=len(grouped)
        site_rows=[];fail=Counter()
        for k,events in sorted(grouped.items()):
            dates=set(e["date"] for e in events)
            unique_positions={(round(e["lat"],5),round(e["lon"],5)) for e in events}
            first=events[0]
            # For metadata-only fuzzy group (e.g. 1.5x1.5km Lucas grid),
            # the supplied coordinates should match within 500m to avoid
            # attributing separately located sections to the same point.
            ref=first
            if len(dates)<NEED_VISITS:fail["fewer_than_two_distinct_dates"]+=1
            if any(e["coord_uncertainty_m"] is None or e["coord_uncertainty_m"]<0 or
                   e["coord_uncertainty_m"]>MAX_ACCURACY_M for e in events):
                fail["missing_or_excessive_geolocation_uncertainty"]+=1
            if len(unique_positions)>1:fail["multiple_public_coordinates_for_site_group"]+=1
            if len(dates)<NEED_VISITS or any(
                   e["coord_uncertainty_m"] is None or e["coord_uncertainty_m"]<0 or
                   e["coord_uncertainty_m"]>MAX_ACCURACY_M for e in events) or len(unique_positions)>1:
                continue
            site_rows.append({
               "site_key_sha256":hashlib.sha256(k.encode()).hexdigest(),
               "site_key_kind":first["site_key_kind"],
               "countryCode":first["country"],
               "lat":first["lat"],"lon":first["lon"],
               "samplingProtocol":first["protocol"],
               "2022_unique_visit_dates":len(dates),
               "2022_event_rows":len(events),
               "coordinateUncertaintyInMeters":max(e["coord_uncertainty_m"] for e in events),
               "site_independence_not_confirmed_if_provisional":first["site_key_kind"]!="original_site_metadata_key"
            })
        result["2022_site_group_rejections"]=dict(fail)
        result["2022_repeat_visit_public_geo_500m_sites_n"]=len(site_rows)
        result["2022_repeat_visit_public_geo_500m_provisional_n"]=sum(
            x["site_independence_not_confirmed_if_provisional"] for x in site_rows)
        result["original_independent_2022_monitoring_grid_units_verified"]=False
        result["notes"]=[
            "Only Event core read; never infer butterfly identities/counts from event metadata.",
            "Original locationID, locality and parentEventID may identify visits/routes rather than independent Lucas grid units.",
            "Public 1/5/10km generalized Belgium locations must not be reverse-geocoded.",
            "Missing coordinate uncertainty is unknown, not zero.",
            "Strictly count independent monitoring grid units before any butterfly abundance modeling."
        ]
        return result,site_rows

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--temporary-dir",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    args.temporary_dir.mkdir(parents=True,exist_ok=True)
    reports={};allrows={}
    for name,info in SOURCES.items():
        path=args.temporary_dir/(name+".zip")
        item={
          "source_dataset":info["uuid"],"licence":info["license"],
          "original_url":info["url"],"expected_country":info["country"],
          "no_occurrence_extension_opened":True,
        }
        try:
            item["archive_source"]=read_zip(info,path)
            diag,rows=event_archive(path,info)
            item["original_event_diagnostics"]=diag
            allrows[name]=rows
        except Exception as e:
            item["source_access_failure"]=type(e).__name__+": "+str(e)
            allrows[name]=[]
        reports[name]=item
        d=item.get("original_event_diagnostics",{})
        print("GESUI_NORWAY_FLANDERS_PUBLIC_2022_EVENT_SOURCE_GATE",
          json.dumps({"dataset":name,"archive_size":item.get("archive_source",{}).get("download_bytes"),
                      "2022_events":d.get("event_counts",{}).get("events_targetyear"),
                      "2022_dated":d.get("event_counts",{}).get("dated_events_targetyear"),
                      "strict_site_groups":d.get("2022_repeat_visit_public_geo_500m_sites_n"),
                      "event_fields":d.get("original_event_field_names"),
                      "error":item.get("source_access_failure")},sort_keys=True),flush=True)
    output={
       "schema":"gesui_2022_Norway_Flanders_original_public_BMS_event_site_source_only_v01",
       "utc":datetime.now(timezone.utc).isoformat(),
       "source_only_country_and_site_population_before_EEA_WWTP_comparison":True,
       "NO_BUTTERFLY_OR_EQR_RESPONSES_OPENED":True,
       "archive_reports":reports,
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(output,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    for name,rows in allrows.items():
        p=args.out.parent/(name+"_2022_STRICT_PUBLIC_EVENT_SITE_ONLY.csv")
        with p.open("w",newline="") as f:
            cols=["site_key_sha256","site_key_kind","countryCode","lat","lon",
                  "samplingProtocol","2022_unique_visit_dates","2022_event_rows",
                  "coordinateUncertaintyInMeters","site_independence_not_confirmed_if_provisional"]
            wr=csv.DictWriter(f,fieldnames=cols);wr.writeheader();wr.writerows(rows)
        output["archive_reports"][name]["strict_public_site_CSV"]=p.name
        output["archive_reports"][name]["strict_public_site_CSV_SHA256"]=hash_file(p)
    args.out.write_text(json.dumps(output,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_NORWAY_FLANDERS_OFFICIAL_EVENT_ONLY_ARCHIVE_NO_SPECIES_QC_COMPLETE",
          json.dumps({"strict_sites":{k:len(v) for k,v in allrows.items()},"no_butterfly_results":True},
           sort_keys=True),flush=True)

if __name__=="__main__":main()
