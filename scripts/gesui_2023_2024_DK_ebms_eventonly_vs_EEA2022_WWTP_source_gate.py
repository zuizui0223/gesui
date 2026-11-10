#!/usr/bin/env python3
"""Official Denmark EBMS EVENT-only monitoring vs 2022 EEA WWTP source feasibility.

Do not read any butterfly occurrence/family or abundance records.
The fixed 2022 EEA reported source is NOT known 2023/2024 wastewater exposure.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,re,time,urllib.request,zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict,Counter
from datetime import datetime,timezone
from pathlib import Path

import pandas as pd
import gesui_2022_sebms_events_to_eea_wwtp_source_only as eea

URL="https://ipt.brc.ac.uk/archive.do?r=ebms_denmark_transects&v=1.4"
EXPECTED_ARCHIVE_SHA="2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b"
ORIGINAL_GBIF_DATASET="caa38f94-ddeb-43e7-9a27-1ac46d7ddc5c"
SOURCE_YEAR=2022
YEARS=(2023,2024)
DATE_BOX=(54.,58.,7.,16.)
SOURCE_REP="DK2022"
MIN_DISTINCT_DAYS=2
MAX_UNCERTAINTY_M=500.0
PRIMARY_KM=2.0
SECONDARY_KM=5.0
REMOTE_KM=10.0
NS="{http://rs.tdwg.org/dwc/text/}"

def sha(data):
    return hashlib.sha256(data).hexdigest()

def nfloat(val):
    try:
        z=float(val)
        if math.isfinite(z):return z
    except (TypeError,ValueError):pass
    return None

def coord_box(lat,lon):
    return 54<=lat<=58 and 7<=lon<=16

def haversine(a,b,c,d):
    r=math.pi/180
    x=(a-c)*r;y=(b-d)*r
    q=math.sin(x/2)**2+math.cos(a*r)*math.cos(c*r)*math.sin(y/2)**2
    return 12742*math.asin(min(1,math.sqrt(q)))

def download_archive(path):
    req=urllib.request.Request(URL,headers={"User-Agent":"gesui-EBMS-DK-event-only-2026-11-10/0.1","Accept":"application/zip"})
    with urllib.request.urlopen(req,timeout=140) as f, path.open("wb") as o:
        while True:
            chunk=f.read(1048576)
            if not chunk:break
            o.write(chunk)
    checksum=sha(path.read_bytes())
    if checksum!=EXPECTED_ARCHIVE_SHA:
        raise RuntimeError("Original exact Denmark EBMS v1.4 source ZIP SHA changed; do NOT proceed to species counts: "+checksum)
    return {"source_URL":URL,"source_archive_SHA256":checksum,"source_archive_bytes":path.stat().st_size}

def read_events_only(archive):
    data=defaultdict(lambda:{"date_by_year":defaultdict(set),"eventrows_by_year":Counter(),
                              "max_uncertainty":None,"missing_uncertainty":False,
                              "localities":set(),"parents":set(),"exact_latlon":set()})
    totals=Counter()
    with zipfile.ZipFile(archive) as z:
        members=z.namelist()
        root=ET.fromstring(z.read("meta.xml"))
        core=root.find(NS+"core")
        if core is None or core.attrib.get("rowType")!="http://rs.tdwg.org/dwc/terms/Event":
            raise ValueError("No original Darwin Core EVENT source core; never open butterfly outcomes")
        loc=core.find(".//"+NS+"files/"+NS+"location")
        if loc is None or not loc.text or "occurrence" in loc.text.lower() or loc.text not in members:
            raise ValueError("Invalid original Event core location")
        fields={}
        for f in core.findall(NS+"field"):
            fields[f.attrib["term"].rsplit("/",1)[-1]]=int(f.attrib["index"])
        missing=["eventDate","decimalLatitude","decimalLongitude"]
        if any(x not in fields for x in missing):
            raise ValueError("Original Event core lacks required source-only fields")
        skip=int(core.attrib.get("ignoreHeaderLines","1"))
        if skip not in (0,1):raise ValueError("Unrecognized original Event header handling")
        delim=core.attrib.get("fieldsTerminatedBy","\t")
        if delim==r"\t":delim="\t"
        if len(delim)!=1:raise ValueError("Unexpected source delimiter")
        with io.TextIOWrapper(z.open(loc.text),encoding="utf-8",newline="") as stream:
            reader=csv.reader(stream,delimiter=delim)
            for __ in range(skip):next(reader,None)
            for record in reader:
                totals["source_original_Event_rows"]+=1
                if len(record)<=max(fields.values()):
                    totals["short_or_corrupted_event_rows"]+=1
                    continue
                get=lambda k:record[fields[k]].strip() if k in fields else ""
                date=get("eventDate")
                year=date[:4]
                totals["Event_year_"+year]+=1
                if year not in ("2022","2023","2024","2025"):continue
                lat=nfloat(get("decimalLatitude"));lon=nfloat(get("decimalLongitude"))
                if lat is None or lon is None or not coord_box(lat,lon):
                    totals["non_Denmark_usable_public_coordinate_events"]+=1
                    continue
                uncertainty=nfloat(get("coordinateUncertaintyInMeters"))
                if uncertainty is None:totals["coordinate_uncertainty_missing_events"]+=1
                elif uncertainty>500:totals["coordinate_uncertainty_gt500_events"]+=1
                # Country centroid 5-decimal public coordinate -> site grouping.
                point=f"{lat:.5f}|{lon:.5f}"
                x=data[point]
                x["eventrows_by_year"][year]+=1
                visit=date[:10]
                if re.fullmatch(r"\d{4}-\d{2}-\d{2}",visit):
                    x["date_by_year"][year].add(visit)
                else:totals["unparseable_source_visit_dates"]+=1
                x["exact_latlon"].add((lat,lon))
                if get("locality"):x["localities"].add(get("locality"))
                if get("parentEventID"):x["parents"].add(get("parentEventID"))
                if uncertainty is None:x["missing_uncertainty"]=True
                else:x["max_uncertainty"]=max(uncertainty,x["max_uncertainty"] or 0.)
        qc={"archive_file_names":members,
            "meta_XML_SHA256":sha(z.read("meta.xml")),
            "core_Event_file":loc.text,
            "Event_term_names":sorted(fields),
            "has_original_DwCA_occurrence_extension":any("occurrence" in n.lower() for n in members),
            "no_original_species_occurrence_or_abundance_rows_opened":True,
            "source_event_row_counts":dict(totals),
            "event_site_distinct_public_points":len(data)}
    sites=[]
    for k,x in data.items():
        lat,lon=map(float,k.split("|"))
        days={y:len(x["date_by_year"].get(str(y),set())) for y in (*YEARS,2025)}
        visits={y:int(x["eventrows_by_year"].get(str(y),0)) for y in (*YEARS,2025)}
        repeated=all(days[y]>=MIN_DISTINCT_DAYS for y in YEARS)
        precision=not x["missing_uncertainty"] and x["max_uncertainty"] is not None and x["max_uncertainty"]<=MAX_UNCERTAINTY_M
        sites.append({"site_id_SHA256":sha(k.encode()),"lat":lat,"lon":lon,
                      "repeat_2023_2024_dates_gate":repeated,
                      "precise_500m_explicit_uncertainty_gate":precision,
                      "no_separate_event_row_independence_assumed":True,
                      "visit_days_2023":days[2023],"visit_days_2024":days[2024],
                      "visit_days_2025_descriptive":days[2025],
                      "Event_rows_2023":visits[2023],"Event_rows_2024":visits[2024],
                      "largest_original_coordinate_uncertainty_m":x["max_uncertainty"],
                      "any_missing_coordinate_uncertainty":x["missing_uncertainty"],
                      "distinct_published_locality_strings":len(x["localities"]),
                      "distinct_original_parentEventIDs":len(x["parents"]),
                      "source_original_coords_consistent_after_5decimal_round":len(x["exact_latlon"])==1})
    qc["source_event_site_with_2days_each_year"]=sum(x["repeat_2023_2024_dates_gate"] for x in sites)
    qc["source_event_site_with_explicit_uncertainty_le500"]=sum(x["precise_500m_explicit_uncertainty_gate"] for x in sites)
    qc["source_event_site_with_both_visit_and_accuracy_gates"]=sum(x["repeat_2023_2024_dates_gate"] and x["precise_500m_explicit_uncertainty_gate"] and x["source_original_coords_consistent_after_5decimal_round"] for x in sites)
    return sites,qc

def read_eea_sources():
    # Source-only uses existing independently audited EEA core methods.
    p,prec=eea.download_layer(eea.PLANT,eea.BOTH_WHERE,eea.PLANT_FIELDS,5492)
    d,drec=eea.download_layer(eea.OUTFALL,"1=1",eea.OUT_FIELDS,24416)
    p=eea.normalize_key(p);d=eea.normalize_key(d)
    p=p[p.compound_join_key.notna()].copy()
    d=d[d.compound_join_key.notna()].copy()
    pc=p.compound_join_key.value_counts();dc=d.compound_join_key.value_counts()
    p=p[p.compound_join_key.map(pc).eq(1)]
    p=p[p.compound_join_key.map(dc).fillna(0).eq(1)]
    m=p.merge(d.drop_duplicates("compound_join_key"),on="compound_join_key",suffixes=("_plant","_outfall"),validate="one_to_one")
    repcol="repCode_plant" if "repCode_plant" in m else "repCode"
    m=m[m[repcol].astype("string").str.strip().str.upper().eq(SOURCE_REP)].copy()
    m["dcpLatitude"]=pd.to_numeric(m.dcpLatitude,errors="coerce")
    m["dcpLongitude"]=pd.to_numeric(m.dcpLongitude,errors="coerce")
    m=m[m.dcpWaterBodyType.astype("string").str.strip().str.upper().eq("FW")]
    m=m[m.dcpLatitude.between(54,58)&m.dcpLongitude.between(7,16)]
    out=[]
    for item in m.itertuples(index=False):
        key=str(item.compound_join_key)
        out.append({"outfall_id_SHA256":sha(key.encode()),
                    "lat":float(item.dcpLatitude),"lon":float(item.dcpLongitude),
                    "original_2022_reported_N_t_year":float(item.uwwNDischargeMeasured),
                    "original_2022_reported_P_t_year":float(item.uwwPDischargeMeasured)})
    if len(out)!=69:
        raise ValueError(f"Unexpected original EEA DK2022 qualified 1-FW reporter outfalls ({len(out)} != 69)")
    qc={"original_2022_EEA_measured_bothNP_plants_5492":True,
        "original_2022_EEA_discharge_points_24416":True,
        "strict_reporter_code":SOURCE_REP,
        "one_full_reporter_plant_key_only":True,
        "one_direct_FreshWater_outfall_only":True,
        "unique_Danish_qualified_2022_FW_bothNP_outfalls_n":len(out),
        "all_original_EEA_source_pages_SHA_receipts":{
            "plants":prec,"discharge_points":drec},
        "annual_reported_source_mass_not_instream_concentration_or_butterfly_dose":True}
    return out,qc

def test_source_gate(sites,sources):
    for x in sites:
        distances=sorted((haversine(x["lat"],x["lon"],y["lat"],y["lon"]) for y in sources))
        x["nearest_qualified_EEA2022_DK_FW_outfall_km"]=distances[0] if distances else None
        x["country_1deg_lat_lon_cell"]=f"{math.floor(x['lat'])}|{math.floor(x['lon'])}"
    qualified=[x for x in sites if x["repeat_2023_2024_dates_gate"] and x["precise_500m_explicit_uncertainty_gate"] and x["source_original_coords_consistent_after_5decimal_round"]]
    report={}
    for rad,minimum in ((PRIMARY_KM,8),(SECONDARY_KM,12)):
        near=[x for x in qualified if x["nearest_qualified_EEA2022_DK_FW_outfall_km"] is not None and x["nearest_qualified_EEA2022_DK_FW_outfall_km"]<=rad]
        far=[x for x in qualified if x["nearest_qualified_EEA2022_DK_FW_outfall_km"] is None or x["nearest_qualified_EEA2022_DK_FW_outfall_km"]>=REMOTE_KM]
        gn={x["country_1deg_lat_lon_cell"] for x in near}
        gf={x["country_1deg_lat_lon_cell"] for x in far}
        shared=gn&gf
        passed=len(near)>=minimum and len(far)>=minimum and len(shared)>=3
        report[f"{int(rad)}km"]={"all_original_2023_2024_repeated_public_500m_eligible_sites":len(qualified),
              "sites_near_registered_2022_EEA_NP_FW_source":len(near),
              "sites_remote_ge10km_registered_qualifying_source":len(far),
              "country_1deg_cells_with_both_near_and_remote":len(shared),
              "sites_in_intermediate_gap_not_a_control":sum(rad<x["nearest_qualified_EEA2022_DK_FW_outfall_km"]<REMOTE_KM for x in qualified if x["nearest_qualified_EEA2022_DK_FW_outfall_km"] is not None),
              "precommitted_source_gate_min_near_and_far_each":minimum,
              "precommitted_min_1deg_cells_both_classes":3,
              "butterfly_original_occurrence_extension_may_be_accessed_ONLY_if_true":passed}
    return report

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--tmpdir",type=Path,required=True)
    args=ap.parse_args()
    args.tmpdir.mkdir(parents=True,exist_ok=True)
    path=args.tmpdir/"Denmark_eBMS_transect_v1.4_DwCA.zip"
    receipt=download_archive(path)
    sites,qc=read_events_only(path)
    src,sqc=read_eea_sources()
    gates=test_source_gate(sites,src)
    out={"schema":"gesui_2023_2024_Denmark_eBMS_2022_EEA_WWTP_source_only_no_butterfly_outcomes_v01",
         "utc":datetime.now(timezone.utc).isoformat(),
         "original_eBMS_2026_v1_4_source":receipt,
         "eBMS_dataset_key":ORIGINAL_GBIF_DATASET,
         "eBMS_original_license":"CC-BY-NC 4.0",
         "butterfly_species_2023_2024_occurrence_extension_READ":False,
         "original_Event_core_source_QC":qc,
         "original_EEA2022_outfall_source_QC":sqc,
         "original_2023_2024_repeated_event_site_source_metadata":sites,
         "original_2022_EEA_DK_qualifying_Fw_outfalls":src,
         "source_distance_and_yearvisit_eligibility":gates,
         "WWTP_source_year_not_dose_or_year_of_butterfly_treatment":True,
         "no_causal_butterfly_effect_claim":True,
         "source_localization_not_proof_of_physical_downstream_exposure":True,
         "exact_original_external_DWCA_source_licence_must_be_respected":True}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"source":"Denmark v1.4 standardized Event core, EEA DK2022 outfalls, NO butterfly species records",
                      "year_counts":{k:v for k,v in qc["source_event_row_counts"].items() if k.startswith("Event_year")},
                      "event_sites":qc["event_site_distinct_public_points"],
                      "repeated_2023_2024_sites":qc["source_event_site_with_2days_each_year"],
                      "accurate_repeat_sites":qc["source_event_site_with_both_visit_and_accuracy_gates"],
                      "qualifying_EEA2022_DK_outfalls":sqc["unique_Danish_qualified_2022_FW_bothNP_outfalls_n"],
                      "source_gate":gates},indent=2))
if __name__=="__main__":
    main()
