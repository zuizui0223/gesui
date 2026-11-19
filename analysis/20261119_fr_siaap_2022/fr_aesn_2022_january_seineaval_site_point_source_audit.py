#!/usr/bin/env python3
"""AESN 2022-Jan source-only exact Seine Aval STATION/POINT/parameter audit.

Download one authentic public monthly CSV (official data.gouv.fr distribution)
with SHA256, no raw data redistribution, no measured nutrient values opened.
Count station-specific records, point and SANDRE parameter codes, sample-date
and quality-unit labels. Do not use ecological outcomes or estimate flux.
"""
import argparse,csv,hashlib,json,re,shutil,time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request,urlopen

CATALOG="https://www.data.gouv.fr/api/1/datasets/64833d77507be41507260862/"
NAME="autosurveillance-donnees-2022-01.csv"
ID="037800501000"
PREF_ID="FR"+ID
REQUIRED={"Année","Mois","Station_Code SANDRE","Station_Nom","Point_Prelevement",
          "Prelevement_Date","Parametre_Code SANDRE","Parametre_Libelle",
          "Resultat_Unite_Symbole","Qualification"}
MONITORING_PARAMETERS={"1551":"NGL_azote_global","1350":"PT_phosphore_total"}
MAX_SIZE=150_000_000
MAX_CATALOG=5_000_000

def approved_url(doc):
    vals=[r for r in doc.get("resources",[]) if r.get("title")==NAME or r.get("filename")==NAME]
    if len(vals)!=1:raise ValueError("One official 2022-01 resource expected")
    r=vals[0]
    u=str(r.get("url") or "")
    p=urlparse(u)
    if p.scheme!="https" or p.netloc!="static.data.gouv.fr" or not p.path.endswith("/"+NAME):
        raise ValueError("Official resource URL host/name unexpected")
    if str(r.get("format") or "").lower()!="csv":raise ValueError("Source resource not CSV")
    n=int(r.get("filesize") or 0)
    if not (50_000_000<n<MAX_SIZE):raise ValueError("Unexpected 2022 original CSV size")
    return u,n

def get_catalog():
    req=Request(CATALOG,headers={"Accept":"application/json","User-Agent":"gesui-2022-aesn-source-only-site-audit/1.0"},method="GET")
    with urlopen(req,timeout=40) as r:
        raw=r.read(MAX_CATALOG+1)
        if r.status!=200 or len(raw)>MAX_CATALOG:raise ValueError("Original catalog absent or too large")
    return json.loads(raw),hashlib.sha256(raw).hexdigest()

def download_exact(url,expected,target):
    # The same original public URL, with bounded retries for transient upstream
    # timeouts; NEVER use unverified mirrors or guessed endpoints.
    errors=[]
    for attempt in range(1,4):
        try:
            req=Request(url,headers={"Accept":"text/csv,application/octet-stream",
                                    "User-Agent":"gesui-2022-aesn-original-site-only/1.0"},method="GET")
            sha=hashlib.sha256()
            amount=0
            with urlopen(req,timeout=110) as r:
                if r.status!=200 or urlparse(r.url).netloc!="static.data.gouv.fr":
                    raise ValueError("Original source transfer origin or HTTP != 200")
                with target.open("wb") as out:
                    while True:
                        b=r.read(1_048_576)
                        if not b:break
                        amount+=len(b)
                        if amount>MAX_SIZE:raise ValueError("Monthly CSV download exceeded size cap")
                        sha.update(b);out.write(b)
            if amount!=expected:raise ValueError("Official resource expected bytes mismatch")
            return {"download_attempt":attempt,"original_full_file_bytes":amount,
                    "original_full_file_SHA256":sha.hexdigest(),
                    "authentic_monthly_CSV_file_not_redistributed":True}
        except Exception as e:
            errors.append(type(e).__name__+":"+str(e)[:120])
            target.unlink(missing_ok=True)
            if attempt<3:time.sleep(1+attempt)
    raise RuntimeError("Official CSV download unavailable: "+" | ".join(errors))

def normalize(s):
    return re.sub(r"\s+"," ",str(s or "")).strip().upper()

def is_station_code(v):
    return normalize(v) in {ID,PREF_ID}

def count_site(path):
    exact_rows=0
    candidates=Counter()
    points=Counter()
    key_parameters=defaultdict(Counter)
    parameter_units=defaultdict(Counter)
    quality=Counter()
    dates=Counter()
    point_types=defaultdict(Counter)
    unrelated_rows=0
    mismatched_source_year=0
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        reader=csv.DictReader(f,delimiter=";")
        fields=reader.fieldnames or []
        if len(fields)!=30 or not REQUIRED.issubset(fields):
            raise ValueError("Original AESN CSV schema changed; do not process")
        for row in reader:
            if row.get("Année","").strip()!="2022" or row.get("Mois","").strip() not in ("1","01"):
                mismatched_source_year+=1
                continue
            code=str(row.get("Station_Code SANDRE") or "").strip()
            name=normalize(row.get("Station_Nom"))
            if is_station_code(code):
                exact_rows+=1
                p=str(row.get("Point_Prelevement") or "").strip() or "MISSING"
                param=str(row.get("Parametre_Code SANDRE") or "").strip()
                label=str(row.get("Parametre_Libelle") or "").strip()[:150]
                unit=str(row.get("Resultat_Unite_Symbole") or "").strip()[:40] or "MISSING"
                points[p]+=1
                quality[str(row.get("Qualification") or "").strip()[:80] or "MISSING"]+=1
                dates[str(row.get("Prelevement_Date") or "").strip()[:12]]+=1
                if param in MONITORING_PARAMETERS:
                    key_parameters[p][param]+=1
                    parameter_units[(p,param,unit,label)]+=1
                point_types[p][str(row.get("Parametre_Type") or "").strip()[:60]]+=1
            elif "SEINE AVAL" in name or "ACHERES" in name or "ACHÈRES" in name:
                candidates[(code,name)]+=1
            else:
                unrelated_rows+=1
    return {
        "full_CSV_rows_other_STEUs_not_redistributed":unrelated_rows,
        "exact_SeineAval_station_code_rows":exact_rows,
        "other_name_candidates_not_original_ID_confirmed":[
            {"station_SANDRE":k[0],"station_name":k[1],"record_count":v}
            for k,v in candidates.most_common(20)],
        "original_source_year_month_mismatch_rows":mismatched_source_year,
        "exact_code_station_point_row_counts":dict(points),
        "exact_code_point_parameter_NGL_PT_counts":{p:dict(q) for p,q in key_parameters.items()},
        "exact_code_point_parameter_unit_and_label":[{
            "point":p,"SANDRE_parameter":param,"unit_symbol":unit,"parameter_label":label,
            "rows":count} for (p,param,unit,label),count in sorted(parameter_units.items())],
        "exact_code_monitoring_qualification_counts":dict(quality),
        "exact_code_distinct_sample_date_count":len(dates),
        "exact_code_point_parameter_types":{p:dict(z) for p,z in point_types.items()},
        "numeric_Resultat_Valeur_never_read":True,
        "no_biological_outcomes_read":True,
        "NGL_PT_quality_and_units_not_yet_sufficient_for_annual_load":True,
        "2022_annual_reported_water_N_P_load_confirmed":False,
        "paired_2022_water_land_mass_and_causal_response_GO":False,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--work",type=Path,required=True)
    a=ap.parse_args()
    a.work.mkdir(parents=True,exist_ok=True)
    result={"schema":"gesui_AESN_2022_01_full_original_CSV_SeineAval_site_param_source_QA_v1",
            "source_catalog":CATALOG,"UTC":datetime.now(timezone.utc).isoformat(),
            "source_status":"ACCESS_STOP",
            "numeric_wastewater_N_P_analysis_started":False,
            "biological_EQR_NUTRI_HNS_values_read":False}
    target=a.work/NAME
    try:
        catalog,catalogsha=get_catalog()
        u,size=approved_url(catalog)
        result["source_catalog_SHA256"]=catalogsha
        result["original_2022_01_resource_url"]=u
        result["expected_official_CSV_bytes"]=size
        result["original_transfer_receipt"]=download_exact(u,size,target)
        result["source_station_parameter_coverage"]=count_site(target)
        result["source_status"]="ORIGINAL_JAN2022_SITE_PARAMETER_COVERAGE_AUDITED"
    except Exception as e:
        result["source_error"]=type(e).__name__+": "+str(e)[:500]
    finally:
        target.unlink(missing_ok=True)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf8")
    info=result.get("source_station_parameter_coverage") or {}
    print("GESUI_AESN_2022JAN_SEINE_AVAL_ORIGINAL_CSV_SOURCE",json.dumps({
      "status":result["source_status"],"error":result.get("source_error"),
      "source_file_sha":(result.get("original_transfer_receipt") or {}).get("original_full_file_SHA256"),
      "exact_station_rows":info.get("exact_SeineAval_station_code_rows"),
      "name_candidates":info.get("other_name_candidates_not_original_ID_confirmed"),
      "points":info.get("exact_code_station_point_row_counts"),
      "NGL_PT_at_point":info.get("exact_code_point_parameter_NGL_PT_counts"),
      "source_unit_label":info.get("exact_code_point_parameter_unit_and_label")
    },ensure_ascii=False,sort_keys=True))

if __name__=="__main__":
    main()
