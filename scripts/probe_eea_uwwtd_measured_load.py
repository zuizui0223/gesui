#!/usr/bin/env python3
"""Probe OFFICIAL EU public UWWTD2022 facility N/P reported discharge-load fields.

Read official ArcGIS REST layer schema + aggregate COUNT queries + only three
non-sensitive source examples. NO plant, nutrient-receptor, or causal analysis.
Published 2026-03-03; primary reference year 2022 (not historic 2011-2024).
"""
from __future__ import annotations
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
import argparse,hashlib,json,requests

BASE="https://water.discomap.eea.europa.eu/arcgis/rest/services/WISE_UWWTD"
SERVICES={
 "treatment_2022":BASE+"/TreatmentPlants_2022/MapServer/0",
 "discharge_2022":BASE+"/DischargePoints_2022/MapServer/0"}
FIELDS=["uwwCode","repCode","uwwName","uwwState","uwwLatitude","uwwLongitude",
        "uwwNRemoval","uwwPRemoval","uwwNDischargeMeasured","uwwPDischargeMeasured",
        "uwwNDischargeCalculated","uwwPDischargeCalculated",
        "uwwNDischargeEstimated","uwwPDischargeEstimated",
        "uwwWasteWaterTreated","uwwBeginLife"]
COUNTRY_SAMPLES=("DE","FR","AT","CZ","DK","BE")
def fetch(url,params=None):
    r=requests.get(url,params=params,timeout=(12,55),
       headers={"Accept":"application/json",
                "User-Agent":"gesui-official-EEA-UWWTD-measured-load-source-audit/0.1"})
    r.raise_for_status()
    raw=r.content
    if len(raw)>1_500_000:raise ValueError("UWWTD source exceeded conservative metadata cap")
    data=json.loads(raw)
    if "error" in data:raise RuntimeError("EEA REST error "+str(data["error"])[:300])
    return data,{"url":r.url,"status_code":r.status_code,
                 "content_length":len(raw),"body_sha256":hashlib.sha256(raw).hexdigest()}
def count(base,where):
    d,rcpt=fetch(base+"/query",{"f":"json","where":where,"returnCountOnly":"true"})
    if "count" not in d:raise ValueError("REST count field missing "+where)
    return {"where":where,"count":int(d["count"]),"source_query_sha256":rcpt["body_sha256"]}
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    results={}
    for label,base in SERVICES.items():
        meta,metadata_receipt=fetch(base,{"f":"json"})
        fields=meta.get("fields",[])
        names=[x["name"] for x in fields]
        summary={
          "official_arcgis_rest_layer":base,
          "layer_name":meta.get("name"),
          "metadata_body_sha256":metadata_receipt["body_sha256"],
          "max_record_count":meta.get("maxRecordCount"),
          "all_field_names":names,
          "has_required_P_load": "uwwPDischargeMeasured" in names,
          "has_required_N_load": "uwwNDischargeMeasured" in names,
          "has_treatment_removal_fields":("uwwPRemoval" in names and "uwwNRemoval" in names),
          "source_response_fields":{x["name"]:{"type":x.get("type"),"alias":x.get("alias")}
                   for x in fields if x.get("name") in FIELDS}
        }
        q=[("all","1=1")]
        if label=="treatment_2022":
            q += [
             ("reported_measured_N","uwwNDischargeMeasured > 0"),
             ("reported_measured_P","uwwPDischargeMeasured > 0"),
             ("reported_measured_both","uwwNDischargeMeasured > 0 AND uwwPDischargeMeasured > 0"),
             ("N_data_any_status","uwwNDischargeMeasured IS NOT NULL"),
             ("P_data_any_status","uwwPDischargeMeasured IS NOT NULL"),
             ("NRemoval_flag_numeric","uwwNRemoval IS NOT NULL"),
             ("PRemoval_flag_numeric","uwwPRemoval IS NOT NULL"),
             ("has_coordinates","uwwLatitude IS NOT NULL AND uwwLongitude IS NOT NULL"),
             ("has_plant_ID","uwwCode IS NOT NULL"),
             ("reporter_DE","repCode = 'DE'"),
             ("reporter_AT","repCode = 'AT'")
            ]
        summary["counts"]={}
        for key,where in q:
            try:summary["counts"][key]=count(base,where)
            except Exception as e:summary["counts"][key]={"error":str(e)[:250]}
        if label=="treatment_2022":
            try:
                sample,rec=fetch(base+"/query",{
                    "f":"json","where":"uwwPDischargeMeasured > 0",
                    "outFields":",".join(FIELDS),
                    "returnGeometry":"false","resultRecordCount":"3",
                    "orderByFields":"OBJECTID ASC"})
                feats=sample.get("features",[])
                # Avoid redistributing names or exact WWTP coordinates in the audit.
                values=[x.get("attributes",{}) for x in feats]
                numeric=["uwwNDischargeMeasured","uwwPDischargeMeasured",
                         "uwwWasteWaterTreated","uwwState","uwwPRemoval","uwwNRemoval"]
                summary["positive_P_load_source_sample_redacted"]=[
                 {"country_code":a.get("repCode"),
                  "measured_N_load_reported":a.get("uwwNDischargeMeasured"),
                  "measured_P_load_reported":a.get("uwwPDischargeMeasured"),
                  "treatment_volume_reported":a.get("uwwWasteWaterTreated"),
                  "uwwState_code_not_decoded":a.get("uwwState"),
                  "N_removal_code_not_decoded":a.get("uwwNRemoval"),
                  "P_removal_code_not_decoded":a.get("uwwPRemoval")}
                 for a in values]
                summary["sample_response_sha256"]=rec["body_sha256"]
            except Exception as e:summary["sample_error"]=type(e).__name__+": "+str(e)[:200]
        results[label]=summary
        print("GESUI_EEA_UWWTD_PUBLIC_API_SOURCE_GATE",
             json.dumps({k:v for k,v in summary.items() if k!="all_field_names"},sort_keys=True),flush=True)
    out={
      "schema":"gesui_EEA_UWWTD2022_official_public_measured_discharge_load_source_gate_v01",
      "retrieved_utc":datetime.now(timezone.utc).isoformat(),
      "published_date":"2026-03-03",
      "reference_year":"2022",
      "source_dataset":"EEA UWWTD 2023 data call, FINAL PUBLIC VERSION, Oct 2025",
      "EEA_UWWTD_public_datahub":"https://www.eea.europa.eu/en/datahub/datahubitem-view/6244937d-1c2c-47f5-bdf1-33ca01ff1715",
      "data_dictionary":"https://dd.eionet.europa.eu/datasets/latest/UWWTDArt15",
      "primary_units_require_data_dictionary_validation":True,
      "source_coordinates_or_station_rows_bulk_downloaded":False,
      "source_is_not_historical_2011_2024":True,
      "underlying_measured_P_and_N_compared_to_Peters_receptors":False,
      "no_plant_ecological_outcome_queried":True,
      "source_metadata":results
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_EEA_UWWTD_REPORTED_LOAD_SOURCE_QC_COMPLETE",
          json.dumps({"both_measured":results["treatment_2022"]["counts"].get("reported_measured_both"),
           "P_reported":results["treatment_2022"]["counts"].get("reported_measured_P"),
           "N_reported":results["treatment_2022"]["counts"].get("reported_measured_N")},sort_keys=True),flush=True)
if __name__=="__main__":main()
