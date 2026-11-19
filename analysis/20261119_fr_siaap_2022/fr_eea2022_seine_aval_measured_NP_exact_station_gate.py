#!/usr/bin/env python3
"""Source-only EEA UWWTD 2022 reported measured N/P discharge for Seine Aval.

Exact station identifier join first; name candidates are NOT confirmed
crosswalks. One unambiguous direct freshwater outfall required to associate
the facility's annual reported export with a river location. Calculated or
estimated N/P values must never replace missing measured values. No biota.
"""
import argparse
import hashlib
import json
import re
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen

STATION_SANDRE="037800501000"
API="https://water.discomap.eea.europa.eu/arcgis/rest/services/WISE_UWWTD"
PLANTS=API+"/TreatmentPlants_2022/MapServer/0"
OUTFALLS=API+"/DischargePoints_2022/MapServer/0"
FIELDS_PLANT=("repCode","uwwCode","uwwName","uwwNDischargeMeasured","uwwPDischargeMeasured",
              "uwwNDischargeCalculated","uwwPDischargeCalculated",
              "uwwNDischargeEstimated","uwwPDischargeEstimated","uwwState",
              "uwwLatitude","uwwLongitude")
FIELDS_DCP=("repCode","uwwCode","dcpCode","dcpName","dcpWaterBodyType",
            "dcpLatitude","dcpLongitude","dcpSurfaceWaters")
QUERIES={
 "exact_SANDRE":"uwwCode = '037800501000'",
 "name_Seine_Aval":"UPPER(uwwName) LIKE '%SEINE AVAL%'",
 "name_Acheres":"UPPER(uwwName) LIKE '%ACHERES%'"
}
CAP_BYTES=850_000

def fetch(base,query):
    url=base+"/query?"+urlencode(query)
    if not url.startswith(API+"/") or "/query?" not in url:
        raise ValueError("Noncanonical EEA public source query")
    req=Request(url,headers={"Accept":"application/json","User-Agent":"gesui-university-2022-source-only-exact-station-audit/1.0"},method="GET")
    with urlopen(req,timeout=45) as resp:
        body=resp.read(CAP_BYTES+1)
        if len(body)>CAP_BYTES:raise ValueError("Source response too large")
        if resp.status!=200:raise ValueError("Official EEA API HTTP non-200")
    doc=json.loads(body)
    if doc.get("error"):raise ValueError("EEA query response error "+str(doc["error"])[:200])
    return doc,{"bytes":len(body),"SHA256":hashlib.sha256(body).hexdigest(),"http":200}

def attrs(doc):
    if not isinstance(doc,dict) or not isinstance(doc.get("features"),list):
        raise ValueError("Source query has no features field")
    if doc.get("exceededTransferLimit"):
        raise ValueError("Truncated source query cannot establish unique records")
    return [q["attributes"] for q in doc["features"] if isinstance(q,dict) and isinstance(q.get("attributes"),dict)]

def match_exact(plant_candidates,outlets):
    exact=[]
    for p in plant_candidates:
        if str(p.get("uwwCode") or "").strip()!=STATION_SANDRE:
            continue
        rep=str(p.get("repCode") or "")
        if not rep.upper().startswith("FR"):
            # Code alone might collide between countries.
            continue
        exact.append(p)
    if len(exact)!=1:
        return {"status":"EXACT_STEU_NOT_UNIQUELY_CONFIRMED",
                "exact_fr_station_rows":len(exact),"reported_measured_annual_N_P_t":None,
                "confirmed_one_direct_Fw_outfall":False}
    x=exact[0]
    rep=x["repCode"]
    d=[o for o in outlets if str(o.get("repCode") or "").strip()==rep and
       str(o.get("uwwCode") or "").strip()==STATION_SANDRE]
    fw=[o for o in d if str(o.get("dcpWaterBodyType") or "").upper()=="FW"]
    n,p=x.get("uwwNDischargeMeasured"),x.get("uwwPDischargeMeasured")
    def positive_num(v):
        try:
            f=float(v)
            return f>0 and f<1e12
        except (ValueError,TypeError):
            return False
    hasmeasured=positive_num(n) and positive_num(p)
    isfw=(len(d)==1 and len(fw)==1)
    confirmed=hasmeasured and isfw
    return {
      "status":"EXACT_EEA_2022_MEASURED_FRESHWATER_SOURCE_MATCH" if confirmed else "EXACT_STEU_PRESENT_INCOMPLETE_SOURCE",
      "exact_fr_station_rows":1,
      "station_identifier":STATION_SANDRE,
      "original_reporter_code":rep,
      "reported_measured_annual_N_P_t":
         {"N_t_per_2022":float(n),"P_t_per_2022":float(p)} if confirmed else None,
      "reported_N_measurement_available":positive_num(n),
      "reported_P_measurement_available":positive_num(p),
      "same_key_outfall_count":len(d),
      "direct_freshwater_outfall_count":len(fw),
      "confirmed_one_direct_Fw_outfall":isfw,
      "outlet_original_waterbody_code":d[0].get("dcpWaterBodyType") if len(d)==1 else None,
      "facility_record_count_not_annual_plant_farm_N_P_balance":1,
      "calculated_or_estimated_never_substituted_for_measured":True,
      "product_specific_2022_land_N_P_mass_missing":True,
      "paired_measured_water_N_P_and_2022_actual_field_N_P_mass_GO":False,
      "biotic_response_and_mechanism_GO":False
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    result={"schema":"gesui_FR_SeineAval_EEA2022_measured_NP_exact_station_SOURCE_only_v1",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "eea_primary_year":2022,
      "eea_units_for_original_measured_discharge":"tonnes_per_year",
      "unit_primary_definitions":{
       "N":"https://dd.eionet.europa.eu/dataelements/104885",
       "P":"https://dd.eionet.europa.eu/dataelements/104853"},
      "official_EEA_public_services":{"plants":PLANTS,"outfalls":OUTFALLS},
      "current_FR_SANDRE_station_id":STATION_SANDRE,
      "SOURCE_coverage_retrieved":False,
      "source_quality":"ACCESS_STOP",
      "EEA_NP_outlet_match":None,
      "unconfirmed_name_matches":[],
      "primary_2022_sludge_actual_spreading_product_N_P_mass_unverified":True,
      "no_EQR_NUTRI_HNS_used":True
    }
    try:
        plants,receipts=[],{}
        for title,where in QUERIES.items():
            q,rec=fetch(PLANTS,{"f":"json","where":where,"outFields":",".join(FIELDS_PLANT),
                                "returnGeometry":"false","resultRecordCount":"35"})
            matches=attrs(q)
            if len(matches)>35:raise ValueError("Unexpected response size")
            receipts[title]={"count":len(matches),**rec}
            plants.extend(matches)
            if title!="exact_SANDRE":
                result["unconfirmed_name_matches"].extend([
                    {"candidate_lookup":title,"uwwCode":str(v.get("uwwCode") or "")[:100],
                     "repCode":str(v.get("repCode") or "")[:100],
                     "name":str(v.get("uwwName") or "")[:140]}
                    for v in matches])
        # Keep canonical country+code key and avoid multiple search aliases.
        unique={}
        for x in plants:
            k=(str(x.get("repCode")),str(x.get("uwwCode")))
            if k in unique and unique[k]!=x:
                raise ValueError("Same EEA station key conflicts across queries")
            unique[k]=x
        exact_candidates=[x for x in unique.values() if str(x.get("uwwCode") or "").strip()==STATION_SANDRE and str(x.get("repCode") or "").upper().startswith("FR")]
        outfalls=[]
        if len(exact_candidates)==1:
            p=exact_candidates[0]
            rep=str(p.get("repCode") or "")
            if not re.fullmatch(r"[A-Za-z0-9_.-]{2,100}",rep):
                raise ValueError("Unsafe original reporter code format")
            d,rec=fetch(OUTFALLS,{"f":"json","where":f"repCode = '{rep}' AND uwwCode = '{STATION_SANDRE}'",
                                   "outFields":",".join(FIELDS_DCP),"returnGeometry":"false","resultRecordCount":"35"})
            outfalls=attrs(d)
            receipts["exact_original_station_outfall"]={"count":len(outfalls),**rec}
        result["EEA_NP_outlet_match"]=match_exact(list(unique.values()),outfalls)
        result["source_query_receipts"]=receipts
        result["SOURCE_coverage_retrieved"]=True
        result["source_quality"]=result["EEA_NP_outlet_match"]["status"]
    except Exception as e:
        result["source_error"]=type(e).__name__+": "+str(e)[:240]
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf8")
    print("GESUI_FR_SIAAP_SEINEAVAL_EEA2022_EXACT_SOURCE",json.dumps({
      "status":result["source_quality"],"error":result.get("source_error"),
      "match":result["EEA_NP_outlet_match"],"source_receipts":result.get("source_query_receipts"),
      "unconfirmed_name_candidates":len(result["unconfirmed_name_matches"]),
    },sort_keys=True))

if __name__=="__main__":
    main()
