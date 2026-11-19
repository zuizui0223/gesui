#!/usr/bin/env python3
"""SOURCE-first official Seine-Normandie treated-effluent autosurveillance catalog.

Does NOT download large resource CSVs, infer 2022 effluent flux, or use any
biology. The public data.gouv.fr dataset is 64833d77507be41507260862,
publisher Agence de l'eau Seine-Normandie; original monthly resources are
listed in its API. We freeze actual metadata before selecting a month.
"""
import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

DATASET_ID="64833d77507be41507260862"
DATASET_API=f"https://www.data.gouv.fr/api/1/datasets/{DATASET_ID}/"
DATASET_PAGE="https://www.data.gouv.fr/datasets/assainissement-stations-depuration-suivi-de-la-qualite-des-eaux-traitees"
MAX_META_BYTES=5_000_000
RESOURCE_YEAR_MONTH=re.compile(r"autosurveillance.donnees.(20\d{2}).(0[1-9]|1[012])",re.I)

def parse_catalog(meta):
    if not isinstance(meta,dict) or not isinstance(meta.get("resources"),list):
        raise ValueError("Original public French catalog is not a dataset with resources")
    published=meta.get("title") or ""
    resources=meta["resources"]
    if len(resources)>1000:
        raise ValueError("Unbounded official source resources")
    months={}
    candidates=[]
    for r in resources:
        title=str(r.get("title") or r.get("filename") or "")
        m=RESOURCE_YEAR_MONTH.search(title)
        if m:
            year,month=m.groups()
            key=f"{year}-{month}"
            if key in months:raise ValueError("Duplicate source year-month titles")
            months[key]={"year":int(year),"month":int(month),"title":title,
                         "resource_id":r.get("id"),"format":r.get("format"),
                         "size_bytes":r.get("filesize"),
                         "official_distribution_url":r.get("url")}
        else:
            candidates.append({"title":title,"format":r.get("format")})
    years=sorted(set(x["year"] for x in months.values()))
    for key in ("2022-01","2022-12"):
        if key not in months:pass
    exact22=sorted(x for x in months if x.startswith("2022-"))
    exact18=sorted(x for x in months if x.startswith("2018-"))
    return {
      "dataset_title":published,
      "original_resource_count":len(resources),
      "parsed_year_month_csv_resources":len(months),
      "source_years_with_monthly_resources":years,
      "resource_names_2022":exact22,
      "resource_names_2018":exact18,
      "has_12_distinct_months_for_2022":len(exact22)==12,
      "has_12_distinct_months_for_2018":len(exact18)==12,
      "2022_monthly_resource_metadata":[months[k] for k in exact22],
      "2018_monthly_resource_metadata":[months[k] for k in exact18],
      "other_resource_descriptions":candidates[:40],
      "no_original_data_rows_opened":True,
      "reported_treated_effluent_NGL_PT_not_yet_identified":True,
      "flow_weighting_and_A4_outlet_not_yet_audited":True,
      "SIAAP_Seine_Aval_2022_plant_match_not_yet_confirmed":True,
      "2022_annual_measured_water_N_P_mass_not_yet_confirmed":True,
      "actual_2022_sludge_field_N_P_mass_also_not_yet_verified":True,
      "paired_water_land_2022_mass_GO":False,
      "biology_EQR_NUTRI_HNS_values_read":False
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    payload={
      "schema":"gesui_FR_AESN_2022_2018_WWTP_autosurveillance_metadata_source_gate_v1",
      "official_catalog":DATASET_PAGE,
      "official_API":DATASET_API,
      "utc":datetime.now(timezone.utc).isoformat(),
      "source_retrieval_status":"ACCESS_STOP",
      "source_facts":None,
      "no_2022_SIAAP_N_P_numbers_imputed":True
    }
    try:
        req=Request(DATASET_API,headers={"User-Agent":"gesui-source-only-2022-effluent-catalog/1.0","Accept":"application/json"},method="GET")
        with urlopen(req,timeout=38) as resp:
            data=resp.read(MAX_META_BYTES+1)
            if len(data)>MAX_META_BYTES:raise ValueError("Official catalog exceeds bounded metadata cap")
            if resp.status!=200:raise ValueError("Official catalog HTTP not 200")
            payload["original_response_bytes"]=len(data)
            payload["original_response_SHA256"]=hashlib.sha256(data).hexdigest()
            payload["official_HTTP_status"]=resp.status
        meta=json.loads(data)
        payload["source_facts"]=parse_catalog(meta)
        payload["source_retrieval_status"]="OFFICIAL_CATALOG_PARSED"
    except Exception as e:
        payload["source_error"]=type(e).__name__+": "+str(e)[:250]
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(payload,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf8")
    concise={k:payload.get(k) for k in ("source_retrieval_status","source_error","original_response_bytes")}
    if payload["source_facts"]:
        z=payload["source_facts"]
        concise.update({"years":z["source_years_with_monthly_resources"],
                       "2022_month_count":len(z["resource_names_2022"]),
                       "2018_month_count":len(z["resource_names_2018"]),
                       "first_2022_resources":z["2022_monthly_resource_metadata"][:2]})
    print("GESUI_AESN_PUBLIC_AUTOSURVEILLANCE_SOURCE_CATALOG",json.dumps(concise,sort_keys=True))

if __name__=="__main__":
    main()
