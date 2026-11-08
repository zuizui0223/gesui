#!/usr/bin/env python3
"""2022 open SeBMS visit SITE metadata near reported 2022 EEA WWTP outfalls.

NO original butterfly occurrence/species/count file opened.
NO original aquatic BQE numeric score opened.
This is a spatial data feasibility test; physical 2/5km distance to a
FW WWTP outlet is not confirmed butterfly effluent or nutrition exposure.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,urllib.request,zipfile,xml.etree.ElementTree as ET
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np,pandas as pd
import gesui_eea_nutrient_load_river_crosswalk as eea

SE_SOURCE="https://www.gbif.se/ipt/archive.do?r=lu_sebms"
SE_ARCHIVE_SHA="384a2e5145df29268fe2de79fde9e196fc113e68adfc922a81acd58069437348"
YEAR="2022"
BBOX=(-12.,44.,26.,64.)
NS="{http://rs.tdwg.org/dwc/text/}"
EXPECTED_PLANTS=5492
EXPECTED_DCP=24416
MAX_UNCERTAINTY_M=500.0
MIN_DISTINCT_VISITS=2
SOURCE_NEAR_KM=2.0
SOURCE_SECONDARY_KM=5.0
SOURCE_REMOTE_KM=10.0
PRIMARY_PROTOCOL="Fixed route Pollard walk transect"

def sha(file):
    h=hashlib.sha256()
    with open(file,"rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def numeric(x):
    try:
        y=float(x)
        return y if math.isfinite(y) else None
    except (TypeError,ValueError):return None

def download(path):
    req=urllib.request.Request(SE_SOURCE,headers={"User-Agent":"gesui-sourceonly-SEBMS-EEA-WWTP-nearby-site-QC/0.1"})
    with urllib.request.urlopen(req,timeout=125) as rsp,path.open("wb") as f:
        while True:
            b=rsp.read(1048576)
            if not b:break
            f.write(b)
    if sha(path)!=SE_ARCHIVE_SHA:raise ValueError("Original public SeBMS 2022 full source archive changed")
    return SE_ARCHIVE_SHA

def event_metadata(path):
    sites=defaultdict(list);summary=Counter();field_set={}
    with zipfile.ZipFile(path) as z:
        mb=z.read("meta.xml")
        root=ET.fromstring(mb)
        core=root.find(NS+"core")
        if core is None or core.attrib.get("rowType")!="http://rs.tdwg.org/dwc/terms/Event":
            raise ValueError("No verified Event source core")
        loc=core.find(".//"+NS+"files/"+NS+"location")
        if loc is None or not loc.text:raise ValueError("No source event file")
        terms={}
        for field in core.findall(NS+"field"):
            terms[field.attrib["term"].rsplit("/",1)[-1]]=int(field.attrib["index"])
        if loc.text.lower().find("occurrence")>=0:raise ValueError("Refusing to read butterfly occurrence/count data")
        text=io.TextIOWrapper(z.open(loc.text),encoding="utf-8",newline="")
        reader=csv.reader(text,delimiter="\t")
        for _ in range(int(core.attrib.get("ignoreHeaderLines","1"))):next(reader,None)
        for r in reader:
            if len(r)<=max(terms.values()):continue
            def val(k):return r[terms[k]].strip() if k in terms else ""
            date=val("eventDate")
            yr=val("year") or date[:4]
            if yr!=YEAR:continue
            summary["2022_event_rows"]+=1
            pr=val("samplingProtocol")
            summary["protocol_"+(pr or "MISSING")]+=1
            lat=numeric(val("decimalLatitude"));lon=numeric(val("decimalLongitude"))
            if lat is None or lon is None or not(BBOX[1]<=lat<=BBOX[3] and BBOX[0]<=lon<=BBOX[2]):
                summary["coordinates_unusable_or_outside_study_extent"]+=1;continue
            unc=numeric(val("coordinateUncertaintyInMeters"))
            if unc is None:summary["uncertainty_missing"]+=1
            elif unc>MAX_UNCERTAINTY_M:summary["uncertainty_over_500m"]+=1
            locid=val("locationID")
            sitekey=(locid if locid else f"COORD:{lat:.5f}:{lon:.5f}")+"|"+pr
            sites[sitekey].append({"date":date,"eventID":val("eventID"),
                    "protocol":pr,"country":val("countryCode"),
                    "county":val("stateProvince"),"lat":lat,"lon":lon,
                    "uncertainty_m":unc,"locationID":locid})
        field_set={x:terms.get(x) for x in ("eventID","eventDate","locationID","samplingProtocol",
             "countryCode","stateProvince","decimalLatitude","decimalLongitude",
             "coordinateUncertaintyInMeters")}
    out=[];reject=Counter()
    for key,evs in sites.items():
        dates=set(x["date"] for x in evs if x["date"])
        coords=set((round(x["lat"],5),round(x["lon"],5)) for x in evs)
        first=evs[0]
        if len(coords)!=1:reject["multiple_public_site_coordinates"]+=1;continue
        if len(dates)<MIN_DISTINCT_VISITS:reject["less_than_2_distinct_dates"]+=1;continue
        if any(e["uncertainty_m"] is None or e["uncertainty_m"]>MAX_UNCERTAINTY_M for e in evs):
            reject["missing_or_large_coordinate_uncertainty"]+=1;continue
        out.append({"site_key":key,"locationID":first["locationID"],
                    "protocol":first["protocol"],"country":first["country"],
                    "county":first["county"],"lat":first["lat"],"lon":first["lon"],
                    "2022_events_n":len(evs),"2022_distinct_visit_days":len(dates),
                    "coordinate_uncertainty_m":max(e["uncertainty_m"] for e in evs)})
    return out,{"event_metadata_fields":field_set,"event_summary":dict(summary),
                "site_unique_id_groups_before_qc":len(sites),
                "site_rejection_reasons":dict(reject),
                "qc_eligible_sites":len(out),
                "no_original_occurrence_species_or_counts_read":True}

def haversine(lat1,lon1,lat2,lon2):
    a=math.pi/180
    dlat=(lat2-lat1)*a;dlon=(lon2-lon1)*a
    q=math.sin(dlat/2)**2+math.cos(lat1*a)*math.cos(lat2*a)*math.sin(dlon/2)**2
    return 6371*2*math.asin(min(1,math.sqrt(q)))

def source_eea():
    fields=eea.PLANT_FIELDS
    w,wp=eea.download_layer(eea.PLANT,eea.BOTH_WHERE,fields,EXPECTED_PLANTS)
    d,dp=eea.download_layer(eea.OUTFALL,"1=1",eea.OUT_FIELDS,EXPECTED_DCP)
    w=eea.normalize_key(w);d=eea.normalize_key(d)
    # Source-only reporter+facility full ID, do not substitute 2018/2020 sources.
    w=w[w.compound_join_key.notna()].copy()
    d=d[d.compound_join_key.notna()].copy()
    ws=w.compound_join_key.value_counts()
    dc=d.compound_join_key.value_counts()
    w=w[w.compound_join_key.map(ws).eq(1)].copy()
    w["n_dcp"]=w.compound_join_key.map(dc).fillna(0).astype(int)
    eligible=w[w.n_dcp.eq(1)].merge(
        d.drop_duplicates("compound_join_key"),on="compound_join_key",
        suffixes=("_p","_d"),validate="one_to_one")
    for c in ("dcpLatitude","dcpLongitude","uwwNDischargeMeasured","uwwPDischargeMeasured"):
        eligible[c]=pd.to_numeric(eligible[c],errors="coerce")
    eligible=eligible[eligible.dcpWaterBodyType.astype("string").str.strip().str.upper().eq("FW")]
    eligible=eligible[eligible.dcpLatitude.between(BBOX[1],BBOX[3])&
                      eligible.dcpLongitude.between(BBOX[0],BBOX[2])]
    # Prefer Sweden reporter code, not an arbitrary geographically neighboring country's plant.
    repkey=eligible["repCode_p"] if "repCode_p" in eligible.columns else eligible["repCode"]
    source=eligible[repkey.astype("string").str.strip().str.upper().str.startswith(("SE_","SE-"))].copy()
    qc={
        "original_2022_both_measured_positive_reported_WWTP_count":EXPECTED_PLANTS,
        "original_2022_all_public_discharge_points_count":EXPECTED_DCP,
        "Sweden_qualifying_unique_FW_outfalls_in_44_64N_bbox":len(source),
        "EEA_source_query_page_receipts_SHA256":{
            "WWTP":hashlib.sha256(json.dumps(wp,sort_keys=True).encode()).hexdigest(),
            "DCP":hashlib.sha256(json.dumps(dp,sort_keys=True).encode()).hexdigest()},
        "outlets_not_any_nonFW_or_multi_outlet":True,
        "Sweden_source_report_year_2022":True}
    return source,qc

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--tmpdir",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    a.tmpdir.mkdir(parents=True,exist_ok=True)
    archive=a.tmpdir/"SeBMS_original_public_DwCA.zip"
    download(archive)
    sites,bqc=event_metadata(archive)
    sources,sqc=source_eea()
    src=[(float(q.dcpLatitude),float(q.dcpLongitude)) for q in sources.itertuples(index=False)]
    if not src:
        print("GESUI_EBMS_2022_SWEDEN_ZERO_ELIGIBLE_QUALIFYING_REGISTERED_WWTP_SOURCE",flush=True)
    out=[]
    for site in sites:
        d=min((haversine(site["lat"],site["lon"],la,lo) for la,lo in src),default=None)
        s={**site,"nearest_qualified_SWEDEN_2022_FW_WWTP_km":d,
           "source_within2km":bool(d is not None and d<=SOURCE_NEAR_KM),
           "source_within5km":bool(d is not None and d<=SOURCE_SECONDARY_KM),
           "registered_source_remote10km":bool(d is None or d>=SOURCE_REMOTE_KM)}
        out.append(s)
    summary={}
    for pr in (PRIMARY_PROTOCOL,"Point site count"):
        obs=[s for s in out if s["protocol"]==pr]
        near2=[s for s in obs if s["source_within2km"]]
        near5=[s for s in obs if s["source_within5km"]]
        far=[s for s in obs if s["registered_source_remote10km"]]
        countiesNear2=set(s["county"] for s in near2)
        countiesNear5=set(s["county"] for s in near5)
        countiesFar=set(s["county"] for s in far)
        summary[pr]={
          "sites_eligible_public_coordinates_2_2022_visits":len(obs),
          "2022_event_visits_n":sum(s["2022_events_n"] for s in obs),
          "source_near_within2km_n":len(near2),
          "source_near_within5km_n":len(near5),
          "registered_source_remote_10km_or_more_n":len(far),
          "counties_near2km":len(countiesNear2),
          "counties_near5km":len(countiesNear5),
          "counties_remote":len(countiesFar),
          "county_overlap_primary_near2_and_far_n":len(countiesNear2&countiesFar),
          "county_overlap_secondary_near5_and_far_n":len(countiesNear5&countiesFar),
          "predeclared_primary_feasible":len(near2)>=20 and len(far)>=20 and
                   len(countiesNear2)>=2 and len(countiesFar)>=2 and len(countiesNear2&countiesFar)>=2,
          "predeclared_secondary_feasible":len(near5)>=20 and len(far)>=20 and
                   len(countiesNear5)>=2 and len(countiesFar)>=2 and len(countiesNear5&countiesFar)>=2,
        }
        print("GESUI_SEBMS_2022_STANDARDIZED_TRANSECT_WWTP_SOURCE_ONLY",json.dumps({
            "method":pr,**summary[pr]},sort_keys=True),flush=True)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    csvpath=a.out.parent/"SWEDEN_2022_SEBMS_EVENT_SITE_METADATA_OUTFALL_DISTANCES_SOURCE_ONLY.csv"
    if out:
        with csvpath.open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(out[0]))
            writer.writeheader();writer.writerows(out)
    report={
       "schema":"gesui_2022_SeBMS_public_event_sites_vs_EEA_reported_FW_WWTP_source_only_v01",
       "original_SeBMS_archive_SHA256":SE_ARCHIVE_SHA,
       "original_GBIF_dwc_dataset":"be77e203-486c-4651-91b9-8347968b728c",
       "SEbms_license":"CC0 1.0",
       "no_2022_butterfly_occurrence_species_counts_loaded":True,
       "original_EEA2022_WWTP_source":sqc,
       "original_SeBMS_event_metadata":bqc,
       "study_protocol_primary":PRIMARY_PROTOCOL,
       "source_distance_bands_km":{"near_primary":SOURCE_NEAR_KM,
          "near_secondary":SOURCE_SECONDARY_KM,"remote":SOURCE_REMOTE_KM},
       "min_repeated_2022_visits":MIN_DISTINCT_VISITS,
       "max_original_position_uncertainty_m":MAX_UNCERTAINTY_M,
       "source_only_protocol_support":summary,
       "public_source_only_site_CSV":csvpath.name,
       "public_source_only_site_CSV_SHA256":sha(csvpath) if out else None,
       "source_proximity_is_not_measured_butterfly_exposure":True,
       "not_sewage_causality":True,
       "limits":[
           "2022 reported measured WWTP N/P mass is not recipient organic N/sodium/protein or butterfly puddling exposure",
           "BMS transects are terrestrial and do not physically follow directed river contaminant paths",
           "Within-source-distance near controls may differ in land use, effort and native host plants",
           "Only one publicly reported direct freshwater outfall included per plant",
           "Missing registered sources do not mean no pollution",
           "Exact original 2022 site coordinates might have location uncertainty and distribution sensitivity; >500m excluded"
       ]
    }
    a.out.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_2022_SEBMS_EEA_WWTP_NEARBY_OR_REMOTE_SOURCE_ONLY_GATE_COMPLETE",
        json.dumps({"bms_site_n":len(out),"swedish_WWTP_n":len(sources),
          "primary_butterfly_outcome_permission":summary[PRIMARY_PROTOCOL]["predeclared_primary_feasible"],
          "secondary_butterfly_outcome_permission":summary[PRIMARY_PROTOCOL]["predeclared_secondary_feasible"]},sort_keys=True),flush=True)
if __name__=="__main__":main()
