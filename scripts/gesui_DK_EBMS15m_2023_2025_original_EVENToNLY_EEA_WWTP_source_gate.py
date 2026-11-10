#!/usr/bin/env python3
"""Denmark EBMS 15min v1.3 EVENT-ONLY WWTP proximity feasibility.

No original Occurrence table, taxon names, family counts, individualCount
or licensed event-level observer IDs read/re-published in this stage.
"""
from __future__ import annotations

import argparse,csv,hashlib,io,json,math,urllib.request,zipfile
import xml.etree.ElementTree as ET
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

SOURCE_URL="https://ipt.brc.ac.uk/archive.do?r=ebms_denmark_15-minute_counts&v=1.3"
DATASET="99f9976a-dfdd-4359-bf4a-5538f29458cf"
SOURCE_FILE="results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
NS="{http://rs.tdwg.org/dwc/text/}"
YEARSET=(2023,2024,2025)
BANDS=("near2","near5","intermediate","remote10")

def sha(data):return hashlib.sha256(data).hexdigest()

def source_download(path):
    request=urllib.request.Request(
        SOURCE_URL,
        headers={"User-Agent":"gesui-DK-ebms15m-original-Event-only-source-QA-20261110/1.0","Accept":"application/zip"})
    with urllib.request.urlopen(request,timeout=160) as u,path.open("wb") as f:
        while True:
            b=u.read(1048576)
            if not b:break
            f.write(b)
    return {"archive_version":"1.3","archive_url":SOURCE_URL,"original_ZIP_sha256":sha(path.read_bytes()),"ZIP_bytes":path.stat().st_size}

def strict_number(t):
    try:
        x=float(t)
        if math.isfinite(x):return x
    except (TypeError,ValueError):pass
    return None

def haversine(a,b,c,d):
    rad=math.pi/180
    dl=(a-c)*rad
    dg=(b-d)*rad
    z=math.sin(dl/2)**2+math.cos(a*rad)*math.cos(c*rad)*math.sin(dg/2)**2
    return 12742*math.asin(min(1,math.sqrt(z)))

def parse_events(path,source_outfalls):
    source_positions=[(float(s["lat"]),float(s["lon"]),s["outfall_id_SHA256"]) for s in source_outfalls]
    if len(source_positions)!=203 or len({x[2] for x in source_positions})!=203:
        raise ValueError("Official source frame must remain exactly 203 unique DK WWTP outfalls")
    count=Counter()
    sitecells=defaultdict(set)
    nearest_sources=defaultdict(Counter)
    region_bands=defaultdict(lambda:defaultdict(set))
    region_event_counts=defaultdict(Counter)
    cluster_by_year=defaultdict(Counter)
    raw_parent_event_counts=defaultdict(Counter)
    with zipfile.ZipFile(path) as archive:
        meta=archive.read("meta.xml")
        root=ET.fromstring(meta)
        event=root.find(NS+"core")
        if event is None or event.attrib.get("rowType")!="http://rs.tdwg.org/dwc/terms/Event":
            raise ValueError("Not a DarwinCore Event-core archive; never open species counts")
        where=event.find(".//"+NS+"files/"+NS+"location")
        if where is None or not where.text or "occurrence" in where.text.lower() or where.text not in archive.namelist():
            raise ValueError("Source Event table unavailable")
        term={f.attrib["term"].rsplit("/",1)[-1]:int(f.attrib["index"]) for f in event.findall(NS+"field")}
        if not {"eventDate","decimalLatitude","decimalLongitude"}.issubset(term):
            raise ValueError("Original 15min Event core lacks date/coordinates")
        delim=event.attrib.get("fieldsTerminatedBy",r"\t")
        if delim==r"\t":delim="\t"
        if delim!="\t":raise ValueError("Unexpected Event delimiter; do not infer parsing")
        skip=int(event.attrib.get("ignoreHeaderLines","1"))
        if skip not in (0,1):raise ValueError("Unexpected original row header")
        names=archive.namelist()
        extra={
            "original_event_table":where.text,
            "original_event_field_names":sorted(term),
            "source_original_meta_xml_SHA256":sha(meta),
            "original_archive_member_names_not_opened_except_event_and_meta":names,
            "original_occurrence_extension_never_opened":True
        }
        with io.TextIOWrapper(archive.open(where.text),encoding="utf-8",newline="") as file:
            rows=csv.reader(file,delimiter=delim)
            for _ in range(skip):next(rows,None)
            for row in rows:
                count["full_original_Event_rows"]+=1
                if len(row)<=max(term.values()):
                    count["event_rows_short_or_unparseable"]+=1
                    continue
                get=lambda k:row[term[k]].strip() if k in term else ""
                yr=get("eventDate")[:4]
                count["all_event_year_"+yr]+=1
                if yr not in ("2023","2024","2025"):continue
                year=int(yr)
                lat=strict_number(get("decimalLatitude"))
                lon=strict_number(get("decimalLongitude"))
                if lat is None or lon is None or not (54<=lat<=58 and 7<=lon<=16):
                    count["targetyear_outside_DK_or_invalid_coordinates"]+=1
                    continue
                count["targetyear_inside_DK_with_coordinate_event_rows"]+=1
                uncertainty=strict_number(get("coordinateUncertaintyInMeters"))
                if uncertainty is None:count["original_coord_uncertainty_missing"]+=1
                elif uncertainty>500:count["original_coord_uncertainty_gt500"]+=1
                parent=get("parentEventID")
                if parent:count["original_parentEventID_nonempty_for_targetyear"]+=1
                else:count["original_parentEventID_empty_for_targetyear"]+=1
                for quality in ("strict500","missing_allowed_sensitivity"):
                    if quality=="strict500" and (uncertainty is None or uncertainty>500):continue
                    if quality=="missing_allowed_sensitivity" and uncertainty is not None and uncertainty>500:continue
                    count[f"{quality}_eligible_Event_rows_{year}"]+=1
                    nearest=min(
                        ((haversine(lat,lon,la,lo),outid) for la,lo,outid in source_positions),
                        key=lambda p:(p[0],p[1]))
                    dist,nearestID=nearest
                    band="near2" if dist<=2 else "near5" if dist<=5 else "remote10" if dist>=10 else "intermediate"
                    cell=f"{math.floor(lat*100)}|{math.floor(lon*100)}"
                    region=f"{math.floor(lat)}|{math.floor(lon)}"
                    count[f"{quality}_year{year}_{band}_events"]+=1
                    sitecells[(quality,year,band)].add(cell)
                    region_bands[(quality,year)][region].add(band)
                    region_event_counts[(quality,year,region)][band]+=1
                    if band in ("near2","near5"):
                        # Cumulative <=5km near independently from <=2km primary.
                        nearest_sources[(quality,year,band)][nearestID]+=1
                        nearest_sources[(quality,year,"near5CUMULATIVE")][nearestID]+=1
                        cluster_by_year[(quality,year)][nearestID]+=1
                    if band=="near2":
                        count[f"{quality}_year{year}_near5cumulative_events"]+=1
                    elif band=="near5":
                        count[f"{quality}_year{year}_near5cumulative_events"]+=1
    if count["full_original_Event_rows"]!=44884:
        raise ValueError("Original EBMS 15min v1.3 44884 Event rows changed; do not proceed to species outcome")
    stats={}
    for quality in ("strict500","missing_allowed_sensitivity"):
        stats[quality]={}
        for year in YEARSET:
            rawNear2=count[f"{quality}_year{year}_near2_events"]
            rawNear5CUM=count[f"{quality}_year{year}_near5cumulative_events"]
            remote=count[f"{quality}_year{year}_remote10_events"]
            near2cells=sitecells[(quality,year,"near2")]
            near5cells=sitecells[(quality,year,"near2")] | sitecells[(quality,year,"near5")]
            remotecells=sitecells[(quality,year,"remote10")]
            near2sources=nearest_sources[(quality,year,"near2")]
            near5sources=nearest_sources[(quality,year,"near5CUMULATIVE")]
            regions=region_bands[(quality,year)]
            joint2=sum("near2" in b and "remote10" in b for b in regions.values())
            joint5=sum(("near2" in b or "near5" in b) and "remote10" in b for b in regions.values())
            primary=rawNear5CUM>=100 and remote>=200 and len(near5sources)>=10 and len(near5cells)>=20 and len(remotecells)>=30 and joint5>=3
            narrow=rawNear2>=50 and remote>=200 and len(near2sources)>=8 and len(near2cells)>=10 and len(remotecells)>=30 and joint2>=3
            stats[quality][str(year)]={
                "original_15min_Event_rows_inside_valid_Denmark_and_accuracy":count[f"{quality}_eligible_Event_rows_{year}"],
                "outfall_near2_Events":rawNear2,
                "outfall_near5_cumulative_Events":rawNear5CUM,
                "registered_source_remote10_Events":remote,
                "intermediate_Events":count[f"{quality}_year{year}_intermediate_events"],
                "unique_original_WWTP_outfall_ids_nearest_near2":len(near2sources),
                "unique_original_WWTP_outfall_ids_nearest_near5":len(near5sources),
                "near2_spatial_approx_1km_grid_cells":len(near2cells),
                "near5_spatial_approx_1km_grid_cells":len(near5cells),
                "remote10_spatial_approx_1km_grid_cells":len(remotecells),
                "original_1degree_regions_containing_near2_and_remote":joint2,
                "original_1degree_regions_containing_near5_and_remote":joint5,
                "predeclared_year_2km_source_support":narrow,
                "predeclared_year_5km_source_support":primary,
                "near5_per_outfall_original_event_counts":dict(sorted(near5sources.items())),
                "near2_per_outfall_original_event_counts":dict(sorted(near2sources.items())),
                "Event_without_Occurrence_not_zero_diagnosed_in_this_stage":True
            }
        stats[quality]["two_or_more_source_qualified_years_for_2km"]=sum(stats[quality][str(y)]["predeclared_year_2km_source_support"] for y in YEARSET)>=2
        stats[quality]["two_or_more_source_qualified_years_for_5km"]=sum(stats[quality][str(y)]["predeclared_year_5km_source_support"] for y in YEARSET)>=2
    return count,extra,stats

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--tmpdir",type=Path,required=True)
    args=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    sourcepath=root/SOURCE_FILE
    src=json.loads(sourcepath.read_text(encoding="utf-8"))
    if src["schema"]!="gesui_2023_2024_Denmark_eBMS_2022_EEA_WWTP_source_only_no_butterfly_outcomes_v01":
        raise ValueError("Original 2022 source frame is not the already frozen Denmark EEA data")
    original=src["original_2022_EEA_DK_qualifying_Fw_outfalls"]
    if len(original)!=203:raise ValueError("Original source count is not 203")
    args.tmpdir.mkdir(parents=True,exist_ok=True)
    tmp=args.tmpdir/"dk-ebms15m-v1.3-source-only.zip"
    receipt=source_download(tmp)
    counts,extra,stats=parse_events(tmp,original)
    out={
        "schema":"gesui_DK_EBMS_15min_2023_2025_2022WWTP_presource_EVENToNLY_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_DK_EBMS_15min_GBIF_dataset":DATASET,
        "original_publisher":"Biological Records Centre",
        "source_license":"CC-BY-NC 4.0",
        "source_archive_receipt":receipt,
        "original_EEA2022_203_outfall_source_file_SHA256":sha(sourcepath.read_bytes()),
        "original_EEA_source_reporter_freshwater_bothNP":203,
        "source_Event_only_raw_QC":dict(counts),
        "source_original_DwCA_event_metadata":extra,
        "source_grouped_year_and_accuracy_support":stats,
        "strict_500m_primary_2km_butterfly_count_access_gate":stats["strict500"]["two_or_more_source_qualified_years_for_2km"],
        "strict_500m_primary_5km_butterfly_count_access_gate":stats["strict500"]["two_or_more_source_qualified_years_for_5km"],
        "original_Occurrence_extension_READ_in_this_stage":False,
        "no_butterfly_taxon_families_or_individualCount_READ":True,
        "original_EBMS_15min_is_not_fixed_site_transect_or_independent_event_replicates":True,
        "near_outfall_geographic_proximity_not_downstream_WWTP_exposure":True,
        "stage0_fail_does_NOT_establish_sewage_no_effect":True
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "original_archive_SHA256":receipt["original_ZIP_sha256"],
        "original_Event_rows":counts["full_original_Event_rows"],
        "year_counts":{str(y):counts["all_event_year_"+str(y)] for y in YEARSET},
        "quality_missing_events":counts["original_coord_uncertainty_missing"],
        "primary_500m_source_gate":{y:{k:v for k,v in stats["strict500"][str(y)].items() if k!="near5_per_outfall_original_event_counts" and k!="near2_per_outfall_original_event_counts"} for y in YEARSET},
        "primary_5km_go":out["strict_500m_primary_5km_butterfly_count_access_gate"],
        "primary_2km_go":out["strict_500m_primary_2km_butterfly_count_access_gate"]
    },indent=2))

if __name__=="__main__":
    main()
