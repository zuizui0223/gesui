#!/usr/bin/env python3
"""Convert public 2022 Norway NINA Event sections to independent Lucas grid units.

No butterfly occurrence/species/count extension accessed. The source Event
core records 50m segments and repeated dates, NOT independent butterfly sites.
Grid grouping is accepted ONLY if an original public locality+county identifier
groups compatible segments in a compact <=3km footprint. No spatial guessing
or recombination of withheld sensitive locations.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,re,urllib.request,zipfile
from collections import Counter,defaultdict
from pathlib import Path
from datetime import datetime,timezone
from xml.etree import ElementTree as ET

URL="https://ipt.nina.no/archive.do?r=butterflies_bumblebees2020"
EXPECTED_SHA="44928a2ad02ca85e58325a3f2025e2ea5edea0fb8afefd15bdca52a0a1cf45af"
NS="{http://rs.tdwg.org/dwc/text/}"
YEAR="2022"
MAX_POSITION_ERROR_M=500
MAX_RADIUS_KM=2.5
MIN_SEGMENT_COUNT=5
MAX_SEGMENT_COUNT=25
MIN_VISIT_DATES=2

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def getzip(path):
    req=urllib.request.Request(URL,headers={"User-Agent":"gesui-public-Norway-2022-BMS-Lucas-grid-Event-only-QC/1.0"})
    with urllib.request.urlopen(req,timeout=155) as r,path.open("wb") as f:
        for b in iter(lambda:r.read(1024*1024),b""):f.write(b)
    if sha(path)!=EXPECTED_SHA:raise ValueError("Original NINA Event source ZIP drift")
def dkm(la,lo,lb,lon):
    z=math.pi/180
    u=math.sin((la-lb)*z/2)**2+math.cos(la*z)*math.cos(lb*z)*math.sin((lo-lon)*z/2)**2
    return 6371*2*math.atan2(math.sqrt(u),math.sqrt(max(0,1-u)))

def read_event(p):
    with zipfile.ZipFile(p) as z:
        meta=ET.fromstring(z.read("meta.xml"))
        core=meta.find(NS+"core")
        if core is None or core.get("rowType")!="http://rs.tdwg.org/dwc/terms/Event":
            raise ValueError("Cannot safely read original Event-only core")
        fields={f.get("term","").rsplit("/",1)[-1]:int(f.get("index")) for f in core.findall(NS+"field") if f.get("index") is not None}
        loc=core.find(".//"+NS+"files/"+NS+"location").text
        if "occurrence" in loc.lower():raise ValueError("Refuse to read taxon/occurrence source")
        sep=core.get("fieldsTerminatedBy",r"\t")
        sep="\t" if sep==r"\t" else sep
        if sep not in ("\t",","):raise ValueError("Unknown core delimiter")
        skip=int(core.get("ignoreHeaderLines","1"))
        records=[];count=Counter()
        with io.TextIOWrapper(z.open(loc),encoding=core.get("encoding","UTF-8"),newline="",errors="replace") as f:
            reader=csv.reader(f,delimiter=sep)
            for _ in range(skip):next(reader,None)
            for row in reader:
                if len(row)<=max(fields.values()):continue
                def value(k):return row[fields[k]].strip() if k in fields else ""
                date=value("eventDate")
                if not re.match(r"^2022-\d{2}-\d{2}$",date):continue
                try:
                    lat=float(value("decimalLatitude"));lon=float(value("decimalLongitude"))
                    err=float(value("coordinateUncertaintyInMeters"))
                except (ValueError,TypeError):
                    count["2022_missing_coordinates_or_accuracy"]+=1;continue
                if not (math.isfinite(lat) and math.isfinite(lon) and math.isfinite(err)):continue
                site=value("locationID");locality=value("locality")
                county=value("stateProvince");municipality=value("municipality")
                parent=value("parentEventID")
                protocol=value("samplingProtocol")
                if not locality:
                    count["2022_missing_original_grid_locality"]+=1
                    continue
                gridid=county+"|"+municipality+"|"+locality
                records.append({
                    "grid_key":gridid,"locationID":site,"date":date,
                    "protocol":protocol,"parentEventID":parent,"county":county,
                    "municipality":municipality,"lat":lat,"lon":lon,
                    "accuracy":err,
                })
                count["2022_dated_precise_section_events"]+=1
        return records,dict(count),fields

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--temporary-dir",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    a.temporary_dir.mkdir(parents=True,exist_ok=True)
    source=a.temporary_dir/"Norway_BMS_original_2022_events.zip"
    getzip(source)
    records,loadqc,fieldmap=read_event(source)
    groups=defaultdict(list)
    for v in records:groups[v["grid_key"]].append(v)
    summary=Counter();gridrows=[];examples=[]
    for key,part in sorted(groups.items()):
        sections=set(p["locationID"] for p in part if p["locationID"])
        dates=set(p["date"] for p in part)
        protocols=set(p["protocol"] for p in part)
        locs={(round(p["lat"],7),round(p["lon"],7)) for p in part}
        lat=sum(x[0] for x in locs)/len(locs);lon=sum(x[1] for x in locs)/len(locs)
        radius=max(dkm(lat,lon,la,lo) for la,lo in locs)
        err=max(p["accuracy"] for p in part)
        summary["original_locality_grid_groups"]+=1
        summary["grid_groups_with_20_sections"]+=(len(sections)==20)
        summary["grid_groups_with_3_distinct_dates"]+=(len(dates)==3)
        if len(examples)<6:
            examples.append({"original_grid_key_SHA256":hashlib.sha256(key.encode()).hexdigest(),
             "sections":len(sections),"dates":len(dates),"protocols":sorted(protocols),
             "max_segment_radius_km":radius,"max_coord_uncertainty_m":err})
        if len(sections)<MIN_SEGMENT_COUNT or len(sections)>MAX_SEGMENT_COUNT:
            summary["excluded_section_count_not_plausible"]+=1;continue
        if len(dates)<MIN_VISIT_DATES:
            summary["excluded_too_few_2022_visit_days"]+=1;continue
        if len(protocols)!=1:
            summary["excluded_multiple_sampling_protocols_per_grid"]+=1;continue
        if radius>MAX_RADIUS_KM:
            summary["excluded_grid_footprint_over2_5km"]+=1;continue
        if err>MAX_POSITION_ERROR_M or err<0:
            summary["excluded_missing_or_large_coordinate_uncertainty"]+=1;continue
        code=hashlib.sha256(key.encode()).hexdigest()
        gridrows.append({
         "original_locality_grid_key_SHA256":code,
         "countryCode":"NO",
         "county":part[0]["county"],"municipality":part[0]["municipality"],
         "grid_centroid_lat":round(lat,8),"grid_centroid_lon":round(lon,8),
         "grid_max_section_distance_from_center_km":round(radius,6),
         "coordinateUncertaintyInMeters_max":round(err,3),
         "unique_source_section_locationIDs":len(sections),
         "original_event_rows_2022":len(part),
         "original_2022_unique_visit_dates":len(dates),
         "original_2022_protocol":part[0]["protocol"],
         "source_event_data_only":True,
         "grid_independence_not_guaranteed_by_locality_only":True,
        })
    # Verify original source label is stable and groups plausibly match ~20 sections per Lucas square.
    plausible_ratio=summary["grid_groups_with_20_sections"]/max(1,summary["original_locality_grid_groups"])
    note="Source grid grouping requires >=80% of locality groups to have exactly 20 locationIDs"
    verified=bool(plausible_ratio>=0.8 and len(gridrows)>0)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    csvp=a.out.parent/"NORWAY_2022_LUCAS_GRID_SOURCE_ONLY_SITES.csv"
    cols=["original_locality_grid_key_SHA256","countryCode","county","municipality",
          "grid_centroid_lat","grid_centroid_lon",
          "grid_max_section_distance_from_center_km",
          "coordinateUncertaintyInMeters_max",
          "unique_source_section_locationIDs","original_event_rows_2022",
          "original_2022_unique_visit_dates","original_2022_protocol",
          "source_event_data_only","grid_independence_not_guaranteed_by_locality_only"]
    with csvp.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(gridrows if verified else [])
    report={
     "schema":"gesui_2022_Norway_NINA_source_only_locality_Lucas_grid_independence_QC_v01",
     "timestamp":datetime.now(timezone.utc).isoformat(),
     "original_national_source_archive_SHA256":EXPECTED_SHA,
     "event_fields_used":sorted(fieldmap),
     "2022_event_metadata_count":loadqc,
     "source_independence_summary":dict(summary),
     "original_2022_locality_grid_group_count":len(groups),
     "original_2022_source_only_site_count_before_grid_independence_verification":len(gridrows),
     "20_section_grid_group_share":plausible_ratio,
     "grid_from_original_locality_20_segment_criterion_pass":verified,
     "grid_independence_criterion":note,
     "source_only_2022_verified_grid_csv_path":csvp.name,
     "source_only_2022_verified_grid_csv_sha256":sha(csvp),
     "first_six_anonymized_site_structure_examples":examples,
     "max_grid_radius_km":MAX_RADIUS_KM,
     "grid_footprint_not_a_single_butterfly_point":True,
     "no_butterfly_or_EQR_species_or_count_data_read":True,
     "no_2022_qualifying_WWTP_data_used_to_select_or_redefine_grid_units":True,
     "grid_sites_should_be_treated_as_independent_sampling_units":verified,
     "no_association_results":True
    }
    a.out.write_text(json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_NORWAY_2022_LUCAS_GRID_UNIT_INDEPENDENCE",
     json.dumps({"original_groups":len(groups),"20_segment_share":plausible_ratio,
                "sites_after_source_only_repeated_visits_and_precision_gate":len(gridrows),
                "verified":verified,"count":dict(summary)},sort_keys=True),flush=True)
if __name__=="__main__":main()
