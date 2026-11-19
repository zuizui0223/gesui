#!/usr/bin/env python3
"""2022 AESN monthly WWTP performance: SOURCE-only station/point/unit audit.

Government public data.gouv.fr ~53MB original performance summary. Fetch at
verified official URL in bounded byte ranges; verify total length/Content-Range,
record full SHA, scan original rows but do NOT read Resultat_Valeur. Report
Seine Aval exact SANDRE code 037800501000 station/point/parameter, quality flags,
units, temporal coverage. No ecology. Flux is not deemed measured if only its
column name or units exist.
"""
import argparse,csv,hashlib,json,re,time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request,urlopen

CATALOG="https://www.data.gouv.fr/api/1/datasets/assainissement-stations-depuration-indicateurs-de-performance/"
NAME="stepcol-indicateurs-mensuels-2022.csv"
STEU="037800501000"
STEU_PREFIX="FR"+STEU
PARAMS={"1551":"NGL_azote_global","1350":"PT_phosphore_total"}
EXPECTED_COLS=(
 "Annee","Mois","Perimetre_Gestion_Code","Perimetre_Gestion_Libelle",
 "Station_Type","Station_Code SANDRE","Station_Nom","Departement_Numero",
 "Commune_Code INSEE","Commune_Libelle","XL93","YL93",
 "Station_Capacite nominale","Point_Prelevement_Code",
 "Point_Prelevement_Libelle","Parametre_Code SANDRE","Parametre_Libelle",
 "Resultat_Valeur","Resultat_Unite_Code","Resultat_Unite_Symbole","Etat_calcul"
)
MAX_CATALOG=2_000_000
CHUNK=8*1024*1024

def original_resource(catalog):
    matches=[]
    for x in catalog.get("resources",[]):
        u=str(x.get("url") or "")
        if NAME in (x.get("title"),x.get("filename")) or u.endswith("/"+NAME):
            matches.append(x)
    if len(matches)!=1:raise ValueError("Original 2022 monthly performance resource not unique")
    row=matches[0]
    url=str(row.get("url") or "")
    parse=urlparse(url)
    if parse.scheme!="https" or parse.netloc!="static.data.gouv.fr" or not parse.path.endswith("/"+NAME):
        raise ValueError("Source URL not official static.data.gouv.fr")
    size=int(row.get("filesize") or 0)
    if not 30_000_000<=size<=90_000_000:raise ValueError("Original CSV size out of expected bounds")
    return url,size

def fetch_catalog():
    r=Request(CATALOG,headers={"Accept":"application/json","User-Agent":"gesui-university-AESN-station-source-audit/1.0"},method="GET")
    with urlopen(r,timeout=40) as resp:
        data=resp.read(MAX_CATALOG+1)
        if len(data)>MAX_CATALOG or resp.status!=200:raise ValueError("Metadata source unavailable")
    return json.loads(data),hashlib.sha256(data).hexdigest()

def download_ranges(url,size,target,get=urlopen):
    sha=hashlib.sha256()
    count=0
    receipts=[]
    with target.open("wb") as out:
        for start in range(0,size,CHUNK):
            end=min(size-1,start+CHUNK-1)
            expected=end-start+1
            error=None
            for attempt in range(1,4):
                try:
                    req=Request(url,headers={"Range":f"bytes={start}-{end}",
                         "User-Agent":"gesui-university-AESN-source-range/1.0",
                         "Accept":"text/csv,application/octet-stream"},method="GET")
                    with get(req,timeout=65) as resp:
                        if urlparse(resp.url).netloc!="static.data.gouv.fr":
                            raise ValueError("Untrusted redirect")
                        if resp.status!=206:
                            raise ValueError("Range not honored; refusing full parallel downloads")
                        range_header=resp.headers.get("Content-Range","")
                        match=re.fullmatch(r"bytes\s+(\d+)-(\d+)/(\d+)",range_header)
                        if not match or tuple(map(int,match.groups()))!=(start,end,size):
                            raise ValueError("Incorrect source content range")
                        block=resp.read(expected+1)
                        if len(block)!=expected:
                            raise ValueError("Original source range truncated or too long")
                        h=hashlib.sha256(block).hexdigest()
                    sha.update(block)
                    out.write(block)
                    count+=len(block)
                    receipts.append({"start":start,"end":end,"size":len(block),"sha256":h,"attempt":attempt})
                    break
                except Exception as e:
                    error=f"{type(e).__name__}:{str(e)[:140]}"
                    if attempt<3:time.sleep(attempt)
            else:raise RuntimeError(f"Official 2022 performance source range {start}-{end} failed: {error}")
    if count!=size or target.stat().st_size!=size:
        raise ValueError("Original official performance CSV size not reproduced")
    return {"original_total_bytes":size,"original_SHA256":sha.hexdigest(),
            "source_range_receipts":receipts,
            "original_measurement_rows_not_redistributed":True}

def norm(x):
    return re.sub(r"\s+"," ",str(x or "")).strip().upper()

def audit_rows(p):
    station_rows=0
    name_candidates=Counter()
    site_names=Counter()
    parameter_point=Counter()
    parameter_unit=Counter()
    month_point_param=defaultdict(set)
    states=Counter()
    months=Counter()
    source_bad_year=0
    other_stations=0
    with Path(p).open(newline="",encoding="utf-8-sig") as source:
        reader=csv.DictReader(source,delimiter=";")
        if tuple(reader.fieldnames or ())!=EXPECTED_COLS:
            raise ValueError("Unexpected monthly source column schema")
        for row in reader:
            year=str(row["Annee"] or "").strip()
            month=str(row["Mois"] or "").strip()
            code=norm(row["Station_Code SANDRE"])
            name=norm(row["Station_Nom"])
            if code in (STEU,STEU_PREFIX):
                if year!="2022" or not month.isdigit() or not 1<=int(month)<=12:
                    source_bad_year+=1
                    continue
                station_rows+=1
                mo=f"{int(month):02d}"
                p_code=norm(row["Point_Prelevement_Code"]) or "MISSING"
                p_label=norm(row["Point_Prelevement_Libelle"])[:100]
                k=str(row["Parametre_Code SANDRE"] or "").strip()
                ku=str(row["Resultat_Unite_Code"] or "").strip()[:60] or "MISSING"
                us=str(row["Resultat_Unite_Symbole"] or "").strip()[:100] or "MISSING"
                state=norm(row["Etat_calcul"])[:100] or "MISSING"
                site_names[name]+=1
                months[mo]+=1
                states[state]+=1
                parameter_point[(p_code,p_label,k)]+=1
                parameter_unit[(p_code,k,ku,us,state)]+=1
                month_point_param[(p_code,k)].add(mo)
                # Crucial: NEVER access row['Resultat_Valeur'] here.
            elif "SEINE AVAL" in name or "ACHERES" in name or "ACHÈRES" in name:
                name_candidates[(code,name)]+=1
            else:other_stations+=1
    values={
       "source_station_code_exact_matches":[STEU,STEU_PREFIX],
       "exact_station_record_rows":station_rows,
       "exact_station_reported_names":dict(site_names),
       "exact_station_month_counts":dict(sorted(months.items())),
       "exact_station_record_year_mismatches":source_bad_year,
       "other_station_rows_not_published":other_stations,
       "nonmatching_name_candidates":[{"code":c,"name":n,"row_count":v}
          for (c,n),v in name_candidates.most_common(20)],
       "exact_station_calculation_state_counts":dict(states),
       "exact_station_point_parameter_counts":[{"point":p,"label":l,
                                                 "SANDRE_parameter":k,"rows":v}
          for (p,l,k),v in sorted(parameter_point.items())],
       "exact_station_parameter_unit_state_counts":[{"point":p,"SANDRE_parameter":k,
                        "unit_code":u,"unit_symbol":s,"calculation_state":c,"rows":v}
          for (p,k,u,s,c),v in sorted(parameter_unit.items())],
       "NGL_PT_1551_1350_month_coverage_by_point":[{"point":p,"SANDRE_parameter":k,
                    "distinct_month_count":len(mos),"months":sorted(mos)}
          for (p,k),mos in sorted(month_point_param.items()) if k in PARAMS],
       "original_Resultat_Valeur_NEVER_READ":True,
       "source_result_unit_not_necessarily_mass":True,
       "monthly_indicators_not_independent_source_STEUs":True,
       "flux_calculation_status_not_validated_for_mass_balance":True,
       "same_2022_SeineAval_actual_sludge_N_P_mass_missing":True,
       "no_EQR_NUTRI_HNS_opened":True,
       "paired_water_land_nutrient_mediation_GO":False,
    }
    return values

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--work",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)
    target=args.work/NAME
    result={"schema":"gesui_AESN2022_station_monthly_performance_exact_SeineAval_SOURCE_ONLY_v1",
            "utc":datetime.now(timezone.utc).isoformat(),"status":"ACCESS_STOP",
            "numeric_result_values_read":False,"biology_data_read":False}
    try:
        c,cs=fetch_catalog()
        url,size=original_resource(c)
        result["official_2022_monthly_source_url"]=url
        result["official_catalog_sha256"]=cs
        result["original_download_receipt"]=download_ranges(url,size,target)
        result["source_station_point_parameter_coverage"]=audit_rows(target)
        result["status"]="ORIGINAL_2022_PERFORMANCE_CSV_SOURCE_ONLY_COVERAGE"
    except Exception as exc:
        result["source_error"]=type(exc).__name__+": "+str(exc)[:600]
    finally:
        target.unlink(missing_ok=True)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    d=result.get("source_station_point_parameter_coverage",{})
    print("GESUI_AESN2022_PERFORMANCE_ORIGINAL_SOURCE_ONLY",json.dumps({
        "status":result["status"],"error":result.get("source_error"),
        "file_sha":(result.get("original_download_receipt") or {}).get("original_SHA256"),
        "station_rows":d.get("exact_station_record_rows"),
        "names":d.get("exact_station_reported_names"),
        "months":d.get("exact_station_month_counts"),
        "name_only_candidates":d.get("nonmatching_name_candidates"),
        "NGL_PT_point_months":d.get("NGL_PT_1551_1350_month_coverage_by_point"),
        "original_value_never_read":True},ensure_ascii=False))

if __name__=="__main__":main()
