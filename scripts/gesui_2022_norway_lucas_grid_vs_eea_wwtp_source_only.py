#!/usr/bin/env python3
"""EEA 2022 measured N/P FW outlets versus independent 2022 Norway BMS Lucas grids.

NO butterfly species/occurrence/abundance fields accessed; no EQR scores.
Source-locked official EEA 2022 annual reported positive measured N and P
mass; unique original full reporter+plant ID and one freshwater discharge point.
Strict physical proximity uses grid extent + 5m coordinate uncertainty.
Proximity is NOT direct wastewater exposure or butterfly puddling.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np,pandas as pd
import gesui_2022_sebms_events_to_eea_wwtp_source_only as eea

NORWAY="NO"
BELGIUM="BE"
EXPECTED_SOURCE_PLANTS=5492
EXPECTED_SOURCE_DCP=24416
PRIMARY_NEAR_KM=2.0
SECONDARY_NEAR_KM=5.0
REMOTE_KM=10.0
MIN_NEAR=10
MIN_REMOTE=10
MIN_MATCHED_COUNTIES=2
SOURCE_QC_SCHEMA="gesui_2022_Norway_NINA_source_only_locality_Lucas_grid_independence_QC_v01"

BBOX={
 NORWAY:(4.0,55.0,33.0,72.0),
 BELGIUM:(2.0,49.0,7.0,52.0),
}

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(1048576),b""):h.update(chunk)
    return h.hexdigest()

def dkm(la,lo,lb,lon):
    r=math.pi/180
    dlat=(la-lb)*r;dlon=(lo-lon)*r
    q=math.sin(dlat/2)**2+math.cos(la*r)*math.cos(lb*r)*math.sin(dlon/2)**2
    return 6371*2*math.atan2(math.sqrt(q),math.sqrt(max(0,1-q)))

def source_eea():
    w,wp=eea.download_layer(eea.PLANT,eea.BOTH_WHERE,eea.PLANT_FIELDS,EXPECTED_SOURCE_PLANTS)
    d,dp=eea.download_layer(eea.OUTFALL,"1=1",eea.OUT_FIELDS,EXPECTED_SOURCE_DCP)
    w=eea.normalize_key(w);d=eea.normalize_key(d)
    w=w[w.compound_join_key.notna()].copy()
    d=d[d.compound_join_key.notna()].copy()
    pcounts=w.compound_join_key.value_counts()
    dcounts=d.compound_join_key.value_counts()
    w=w[w.compound_join_key.map(pcounts).eq(1)].copy()
    w["n_full_original_outfalls"]=w.compound_join_key.map(dcounts).fillna(0).astype(int)
    w=w[w.n_full_original_outfalls.eq(1)].copy()
    full=w.merge(d.drop_duplicates("compound_join_key"),on="compound_join_key",
                 suffixes=("_WWTP","_DCP"),validate="one_to_one",how="inner")
    for col in ("uwwNDischargeMeasured","uwwPDischargeMeasured",
                "dcpLatitude","dcpLongitude"):
        full[col]=pd.to_numeric(full[col],errors="coerce")
    full=full[full.dcpWaterBodyType.astype("string").str.strip().str.upper().eq("FW")]
    full=full[full.dcpLatitude.between(45,73)&full.dcpLongitude.between(-13,35)]
    reporter=(full["repCode_WWTP"] if "repCode_WWTP" in full.columns
              else full["repCode"]).astype("string").str.strip().str.upper()
    full=full.copy()
    full["original_report_code"]=reporter
    result={};info={}
    for iso in (NORWAY,BELGIUM):
        west,south,east,north=BBOX[iso]
        geographic=full[full.dcpLatitude.between(south,north) &
                        full.dcpLongitude.between(west,east)]
        # Source QA country labels from original prefix and geolocated source.
        # Do NOT invent a generic SE_ prefix; Sweden previously encoded SE369.
        code_ok=geographic.original_report_code.str.startswith(iso,na=False)
        matched=geographic[code_ok].copy()
        repcount=Counter(geographic.original_report_code.fillna("MISSING").astype(str))
        info[iso]={
           "all_qualifying_reporter_code_geographic_outlets_n":len(geographic),
           "native_country_prefix_matched_FW_outlets_n":len(matched),
           "geography_candidate_reporter_codes":dict(repcount.most_common(25)),
           "original_reporter_codes_matching_country":dict(
                  Counter(matched.original_report_code.astype(str))),
           "geographic_country_extent_bbox":BBOX[iso],
           "annual_reported_measured_NP_fields_used":True,
           "year_specific_individual_lab_date_not_verified":True}
        result[iso]=matched
    qc={
      "original_2022_EEA_both_positive_measured_mass_facility_rows":EXPECTED_SOURCE_PLANTS,
      "original_2022_EEA_discharge_point_rows":EXPECTED_SOURCE_DCP,
      "original_one_to_one_FW_outfalls_valid_in_Wider_EU_extent":len(full),
      "official_2022_plant_query_page_digest":hashlib.sha256(json.dumps(wp,sort_keys=True).encode()).hexdigest(),
      "official_2022_DCP_query_page_digest":hashlib.sha256(json.dumps(dp,sort_keys=True).encode()).hexdigest(),
      "country_source_checks":info,
      "original_source_is_not_delivered_instream_organic_N_or_sodium":True,
    }
    return result,qc

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-grid-qc",type=Path,required=True)
    p.add_argument("--source-grid-csv",type=Path,required=True)
    p.add_argument("--flanders-events-qc",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    prior=json.loads(args.source_grid_qc.read_text())
    if prior["schema"]!=SOURCE_QC_SCHEMA or not prior["grid_from_original_locality_20_segment_criterion_pass"]:
        raise ValueError("Norway independent grid not verified in frozen Event-only source")
    if sha(args.source_grid_csv)!=prior["source_only_2022_verified_grid_csv_sha256"]:
        raise ValueError("Original Lucas grids changed after source freeze")
    source_event=json.loads(args.flanders_events_qc.read_text())
    flanders=source_event["archive_reports"]["Flanders_Meetnetten"]["original_event_diagnostics"]
    if flanders["2022_repeat_visit_public_geo_500m_sites_n"]!=0:
        raise ValueError("Flanders original coordinate-uncertainty gate unexpectedly changed")
    original_grid=list(csv.DictReader(args.source_grid_csv.open()))
    if len(original_grid)!=prior["original_2022_source_only_site_count_before_grid_independence_verification"]:
        raise ValueError("Original 2022 Norway grid population mismatch")
    sources,sourceQC=source_eea()
    rows=[]
    for grid in original_grid:
        la=float(grid["grid_centroid_lat"]);lo=float(grid["grid_centroid_lon"])
        radius=float(grid["grid_max_section_distance_from_center_km"])
        gps=float(grid["coordinateUncertaintyInMeters_max"])/1000
        distances=[dkm(la,lo,float(x.dcpLatitude),float(x.dcpLongitude))
                   for x in sources[NORWAY].itertuples(index=False)]
        nearest=min(distances) if distances else None
        nearest_pessimistic=(nearest+radius+gps) if nearest is not None else None
        nearest_optimistic=max(0.,nearest-radius-gps) if nearest is not None else None
        near2=nearest_pessimistic is not None and nearest_pessimistic<=PRIMARY_NEAR_KM
        near5=nearest_pessimistic is not None and nearest_pessimistic<=SECONDARY_NEAR_KM
        far10=nearest_optimistic is None or nearest_optimistic>=REMOTE_KM
        rows.append({
           "original_source_grid_hash":grid["original_locality_grid_key_SHA256"],
           "county":grid["county"],"lat":grid["grid_centroid_lat"],
           "lon":grid["grid_centroid_lon"],
           "source_2022_independent_visit_days":grid["original_2022_unique_visit_dates"],
           "source_2022_event_section_rows":grid["original_event_rows_2022"],
           "grid_radius_km":radius,"public_uncertainty_km":gps,
           "nearest_2022_registered_measured_NP_FW_outfall_center_distance_km":nearest,
           "conservative_all_grid_sections_near_2km":near2,
           "conservative_all_grid_sections_near_5km":near5,
           "conservative_all_grid_sections_far_10km":far10,
           "distance_status":"NO_ELIGIBLE_REGISTERED_NP_FW_SOURCE" if nearest is None else "computed",
           "no_lycaenidae_species_or_count_data_read":True,
        })
    q={}
    near2=[x for x in rows if x["conservative_all_grid_sections_near_2km"]]
    near5=[x for x in rows if x["conservative_all_grid_sections_near_5km"]]
    far10=[x for x in rows if x["conservative_all_grid_sections_far_10km"]]
    county=lambda coll:set(x["county"] for x in coll)
    for code,s in (("strict_2km",near2),("secondary_5km",near5)):
        overlap=county(s)&county(far10)
        q[code]={
         "independent_Lucas_grid_near_n":len(s),
         "independent_Lucas_grid_remote_ge_10km_n":len(far10),
         "source_near_counties":len(county(s)),
         "source_remote_counties":len(county(far10)),
         "within_county_near_far_overlap_n":len(overlap),
         "overlap_counties":sorted(overlap),
         "predeclared_near_ge_10_remote_ge_10_and_ge_2_common_counties_pass":
              len(s)>=MIN_NEAR and len(far10)>=MIN_REMOTE and len(overlap)>=MIN_MATCHED_COUNTIES,
        }
    a=args.out
    a.parent.mkdir(parents=True,exist_ok=True)
    csvp=a.parent/"NORWAY_2022_INDEPENDENT_LUCAS_GRIDS_VS_REPORTING_WWTP_SOURCE_ONLY.csv"
    with csvp.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    output={
      "schema":"gesui_2022_Norway_Flanders_publisher_Event_only_to_EEA_WWTP_proximity_feasibility_v01",
      "date_UTC":__import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat(),
      "original_public_2022_Norway_Event_source_QC_SHA256":sha(args.source_grid_qc),
      "original_public_2022_Norway_independent_grid_CSV_SHA256":sha(args.source_grid_csv),
      "original_public_2022_Flanders_Event_source_QC_SHA256":sha(args.flanders_events_qc),
      "original_EEA_2022_measured_NP_FW_source_only":sourceQC,
      "source_independent_Norway_grid_units":len(rows),
      "geographic_source_distance_bands_km":{
          "strict_near":PRIMARY_NEAR_KM,"secondary_near":SECONDARY_NEAR_KM,"remote":REMOTE_KM},
      "strict_original_individual_public_location_uncertainty_le_500m":True,
      "all_segment_source_exposure_footprint_bound_conservative":True,
      "source_only_country_gates":{NORWAY:q,BELGIUM:{
        "strict_2022_event_sites_after_public_500m_accuracy_gate":0,
        "predeclared_source_distance_gate_pass":False,
        "reason":"original Belgian 2022 public transect coordinates generalized >=1km, observed uncertainty 7071m"
      }},
      "original_source_only_grid_distance_CSV":csvp.name,
      "original_source_only_grid_distance_CSV_SHA256":sha(csvp),
      "NO_BUTTERFLY_SPECIES_ABUNDANCE_EQR_OR_PUDDLING_DATA_OPENED":True,
      "close_proximity_is_not_effluent_dose":True,
      "registered_source_zero_is_not_no_sewage":True,
      "outcome_data_permission_only_if_gate_pass":any(
         z["predeclared_near_ge_10_remote_ge_10_and_ge_2_common_counties_pass"] for z in q.values()),
      "caveats":[
         "Norwegian Lucas squares, not 50m segments or visits, are independent monitoring units",
         "Annual treatment outfall N/P mass not directly measured environmental water protein/sodium",
         "Geographic proximity not directed in-stream exposure; terrestrial butterflies can disperse away from river",
         "Outfalls not reported in EEA 2022 are not proof of absence of wastewater",
         "No observation of Lycaenidae behavior or biological response read"
      ],
    }
    args.out.write_text(json.dumps(output,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_NORWAY_2022_BMS_GRID_TO_EEA_NP_SOURCE_ONLY_FEASIBILITY",
          json.dumps({"independent_grids":len(rows),
              "WWTP_NO":len(sources[NORWAY]),"WWTP_BE":len(sources[BELGIUM]),
              "norway_gates":q,"flemish_2022_eligible_sites":0},sort_keys=True),flush=True)

if __name__=="__main__":main()
