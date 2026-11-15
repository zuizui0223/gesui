#!/usr/bin/env python3
"""Original 2022 Event-core source-only gate: Denmark and Flanders.

No read of Occurrence, MeasurementOrFact, species, counts or HNS. Avoid using
event-core public generalized positions as fine riverbank butterfly positions.
No ecological outcome or primary 8-cell source/reference treatment is fitted.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,re,time,urllib.request,zipfile
from collections import Counter,defaultdict
from datetime import date
from pathlib import Path
import xml.etree.ElementTree as ET

FROZEN_PROTOCOL_GIT_BLOB="b9ce415c36d398eee77d8318a20aef4873e0bf44"
SOURCES=[
  {"country":"DK","label":"Danish eBMS transects v1.4", "gbif_dataset":"caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c",
   "url":"https://ipt.brc.ac.uk/archive.do?r=ebms_denmark_transects&v=1.4",
   "expected_SHA256":"2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b",
   "licence":"CC-BY-NC-4.0"},
  {"country":"BE","label":"Flanders priority butterfly monitoring event archive",
   "gbif_dataset":"63938753-1fec-4c08-ae39-e9f8a6576521",
   "url":"https://ipt.inbo.be/archive.do?r=meetnetten-butterflies-occurrences",
   "expected_SHA256":None,"licence":"GBIF registry CC0-1.0; original archive coordinate generalization caution"},
]
MAX_ZIP=45*1024*1024
MAX_EVENT_CORE=125*1024*1024
NS={"d":"http://rs.tdwg.org/dwc/text/"}
def git_sha(data:bytes)->str:
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
def fetch_source(url):
    if url not in [x["url"] for x in SOURCES]:
        raise ValueError("Unexpected dataset; only preselected event source allowed")
    error=None
    for retry in range(3):
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"gesui-20261115-BE-DK-event-only-source-gate",
                "Accept":"application/zip,application/octet-stream,*/*"})
            with urllib.request.urlopen(req,timeout=110) as conn:
                raw=conn.read(MAX_ZIP+1)
                if len(raw)>MAX_ZIP: raise RuntimeError("Source archive > 45MB")
            if not zipfile.is_zipfile(io.BytesIO(raw)):
                raise RuntimeError(f"Publisher did not supply Darwin Core ZIP: {raw[:35]!r}")
            return raw,None
        except Exception as e:
            error=str(e)
            if retry<2:time.sleep(retry+2)
    return None,error

def delim(x):
    if x in ("\\t","\\\\t"):return "\t"
    if x in ("\\n","\\\\n"):return "\n"
    if x in ("\\r","\\\\r"):return "\r"
    return x or ","
def fnum(x):
    try:
        y=float(x)
        return y if math.isfinite(y) else None
    except (ValueError,TypeError):return None
def dateval(s):
    m=re.match(r"^(20\d\d)-(\d\d)-(\d\d)",(s or "").strip())
    if not m:return None
    try:return date(*map(int,m.groups()))
    except ValueError:return None

def open_only_event_core(blob):
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        names=z.namelist()
        meta=next((n for n in names if n.rsplit("/",1)[-1]=="meta.xml"),None)
        if meta is None:raise RuntimeError("Darwin Core meta.xml missing")
        root=ET.fromstring(z.read(meta))
        core=root.find("d:core",NS)
        if core is None or not core.attrib.get("rowType","").lower().endswith("/event"):
            raise RuntimeError("Refuse to read non-Event core")
        loc=core.find("d:files/d:location",NS)
        if loc is None or not loc.text:raise RuntimeError("No official Event core member")
        corefile=loc.text.strip()
        basename=meta.rsplit("/",1)[0] if "/" in meta else ""
        path=corefile if corefile in names else (basename+"/"+corefile).lstrip("/")
        if path not in names:raise RuntimeError("Original Event core file not in archive "+corefile)
        if "occurrence" in path.lower():raise RuntimeError("Refuse to read taxon outcomes")
        members=z.getinfo(path)
        if members.file_size>MAX_EVENT_CORE:raise RuntimeError("Event core too large")
        fm={}
        for t in core.findall("d:field",NS):
            term=t.attrib.get("term","").rsplit("/",1)[-1].rsplit("#",1)[-1]
            fm[term]=int(t.attrib.get("index","0"))
        elem=core.find("d:id",NS)
        if elem is not None:fm.setdefault("_coreID",int(elem.attrib["index"]))
        if "eventDate" not in fm:raise RuntimeError("No sampling event date")
        q=core.attrib.get("fieldsEnclosedBy",'"')
        fmt=delim(core.attrib.get("fieldsTerminatedBy",r"\t"))
        if len(fmt)!=1:raise RuntimeError("Unsupported delimiter: "+repr(fmt))
        raw=z.read(path)
        wrapped=io.TextIOWrapper(io.BytesIO(raw),encoding=core.attrib.get("encoding","utf-8"),
                                 errors="strict",newline="")
        reader=csv.reader(wrapped,delimiter=fmt,quotechar=q if len(q)==1 else '"')
        ignored=int(core.attrib.get("ignoreHeaderLines","0"))
        for i in range(ignored):next(reader)
        needed=["eventDate","decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
                "samplingProtocol","eventID","locationID","parentEventID","samplingEffort",
                "samplingSizeValue","samplingSizeUnit","footprintWKT","dataGeneralizations",
                "informationWithheld","countryCode","eventTime","locality","verbatimLocality"]
        records=[]; totals=Counter()
        years=Counter()
        for line in reader:
            totals["all_source_event_core_rows"]+=1
            if len(line)<=max(fm.values()):
                totals["short_source_event_rows"]+=1
                continue
            get=lambda k: line[fm[k]].strip() if k in fm else ""
            when=dateval(get("eventDate"))
            if when:years[str(when.year)]+=1
            if when is None or when.year!=2022:continue
            totals["2022_full_year_rows"]+=1
            if not 5<=when.month<=9:continue
            totals["2022_May_Sep_rows"]+=1
            row={f:get(f) for f in needed}
            row["_valid_date"]=when.isoformat()
            records.append(row)
        return records,{
            "event_meta":meta,"event_member":path,
            "event_core_source_SHA256":hashlib.sha256(raw).hexdigest(),
            "event_core_available_columns":sorted(fm),
            "all_event_rows":totals["all_source_event_core_rows"],
            "source_event_short_rows":totals["short_source_event_rows"],
            "events_2022_full_year":totals["2022_full_year_rows"],
            "events_2022_May_Sep":totals["2022_May_Sep_rows"],
            "parsed_year_row_counts":dict(sorted(years.items())),
            "archive_member_count":len(names),
            "occurrence_extension_was_NOT_opened":True,
            "frozen_source_year_2022":True,
        }

def counts(events):
    qual=Counter();protocol=Counter();uncertainty=Counter();generalization=Counter()
    effort=Counter();days=defaultdict(set);precise_days=defaultdict(set)
    for x in events:
        prot=x.get("samplingProtocol") or "__MISSING__";protocol[prot]+=1
        lat=fnum(x.get("decimalLatitude"));lon=fnum(x.get("decimalLongitude"))
        unc=fnum(x.get("coordinateUncertaintyInMeters"))
        unctxt=x.get("coordinateUncertaintyInMeters") or "__MISSING__"
        uncertainty[unctxt]+=1
        if x.get("dataGeneralizations"):generalization["dataGeneralizations_present"]+=1
        if x.get("informationWithheld"):generalization["informationWithheld_present"]+=1
        if x.get("footprintWKT"):generalization["footprintWKT_present"]+=1
        if x.get("samplingEffort"):effort["samplingEffort_present"]+=1
        if x.get("samplingSizeValue"):effort["samplingSizeValue_present"]+=1
        if x.get("eventTime"):effort["eventTime_present"]+=1
        if x.get("eventID"):effort["eventID_present"]+=1
        if lat is None or lon is None or abs(lat)>90 or abs(lon)>180:
            qual["coordinates_missing_or_invalid"]+=1
            continue
        qual["coordinates_valid"]+=1
        if unc is None:qual["uncertainty_missing"]+=1
        elif unc>100 or unc<0:qual["uncertainty_over100_or_invalid"]+=1
        else:qual["uncertainty_le100"]+=1
        # This key is an approximate *metadata location* only and NEVER an
        # eligible original physical site if survey sections are generalized.
        sk=(x.get("locationID") or "")+"|"+f"{lat:.8f}|{lon:.8f}|{prot}"
        days[sk].add(x["_valid_date"])
        if unc is not None and 0<=unc<=100:
            precise_days[sk].add(x["_valid_date"])
    return {
      "protocol_counts":dict(protocol.most_common(20)),
      "uncertainty_m_values":dict(uncertainty.most_common(15)),
      "location_quality_counts":dict(qual),
      "generalization_flags":dict(generalization),
      "sampling_event_fields_populated":dict(effort),
      "approx_coordinate_location_groups_all_methods":len(days),
      "approx_coordinate_location_groups_with_5days":sum(len(d)>=5 for d in days.values()),
      "approx_coordinate_location_groups_with_5days_and_explicit_uncertainty_le100":sum(len(d)>=5 for d in precise_days.values()),
      "site_groups_are_not_confirmed_original_transect_section_positions":True,
      "no_species_or_five_family_counts_opened":True,
    }
def run(args):
    if git_sha(args.freeze.read_bytes())!=FROZEN_PROTOCOL_GIT_BLOB:
        raise RuntimeError("Frozen one-shot protocol byte identity changed")
    result={"schema":"GESUI_SOURCEONLY_DK_BE_2022_ARCHIVE_EVENT_GATE_V1",
        "frozen_protocol_git_blob":FROZEN_PROTOCOL_GIT_BLOB,
        "candidate_dataset_count":len(SOURCES),"results":[],
        "no_occurrence_extension_loaded":True,
        "no_taxon_or_HNS_scored":True,
        "no_registered_WWTP_outfall_exposure_tested":True,
        "no_2022_eightcell_pairs_qualified":True,
        "primary_HNS_test_ran":False,
        "frozen_primary_status":"STOP_SOURCE_OPPORTUNITY_NOT_VERIFIED"}
    for s in SOURCES:
        data,err=fetch_source(s["url"])
        if data is None:
            result["results"].append({**s,"status":"SOURCE_ARCHIVE_UNAVAILABLE","error":err})
            continue
        digest=hashlib.sha256(data).hexdigest()
        if s["expected_SHA256"] and digest!=s["expected_SHA256"]:
            raise RuntimeError("Immutable publisher source SHA mismatch for "+s["label"])
        rows,metadata=open_only_event_core(data)
        record={**{k:v for k,v in s.items() if k!="expected_SHA256"},
            "source_original_zip_SHA256":digest,"source_archive_bytes":len(data),
            "event_core":metadata,"source_only_event_quality":counts(rows),
            "status":"ORIGINAL_ARCHIVE_EVENT_ONLY_VERIFIED",
            "eligible_2022_exact_100m_river_8cell_pairs":0,
            "zero_pairs_means_NOT_YET_GEOMETRICALLY_CONFIRMED":True}
        result["results"].append(record)
        print("ORIGINAL_2022_SOURCE_EVENT_ONLY="+json.dumps({
            "country":record["country"],"archive_sha256":digest,
            "n_2022_full":metadata["events_2022_full_year"],
            "n_2022_summer":metadata["events_2022_May_Sep"],
            "source_quality":record["source_only_event_quality"],
            "no_biological_test":True},sort_keys=True),flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_TWOCOUNTRY_SOURCE_EVENT_AUDIT="+json.dumps({
      "statuses":[r["status"] for r in result["results"]],
      "countries":len(result["results"]),"outcome_test_ran":False,
      "frozen_stage":result["frozen_primary_status"]},sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--freeze",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
