#!/usr/bin/env python3
"""ROSEAU 2017-22 original public sludge CSV: true record grain & kommune QA.

Source-only. Outcome-agnostic. One official year per invocation. Public URLs,
POST keys and published JS sha are fixed. Do not save raw rows or sum duplicated
station totals repeated per material/parameter assay.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,unicodedata
from collections import Counter,defaultdict
from decimal import Decimal,InvalidOperation
from datetime import datetime,timezone
from pathlib import Path
import requests

BASES=(
  "https://www.assainissement.developpement-durable.gouv.fr",
  "https://assainissement.developpement-durable.gouv.fr",
  "https://www2.assainissement.developpement-durable.gouv.fr",
)
PAGE="/pages/data/basededonneeboue.php"
JS="/public/js/exportbasededonneeBoue_min.js"
HANDLER="/controllers/GetDataExportBoueHandler.php"
JSSHA="f4cf1e2fe84e0ccf608ece6c5d363426960e2b730e21900d48265eba4a78f0e5"
PINNED_2018_SHA="bb9660f13e3dc0b7db359cfc414f14e4510db3b51d73ce0baa4bd5cca821f6c2"
YEARS={str(y) for y in range(2017,2023)}
MAXBYTES=95_000_000
MINBYTES=3_485
FIELD={
  "station":"Code SANDRE la station",
  "station_name":"Nom de la station",
  "dept":"Numéro département",
  "region":"Nom de la région",
  "bassin":"Code du Bassin",
  "tranche":"Tranche obligation",
  "spreading":"Quantité épandage agricole (tMS/an)",
  "produced":"Prod boues sans réactif (tMS/an)",
  "commune":"Commune d'épandage",
  "material":"Référence de la matière",
  "sample_date":"Date du prélèvement",
  "sample_year":"Année du prélèvement",
  "parameter":"Référence paramètre",
  "lot":"Numéro de lot",
}
RE_NUM=re.compile(r"^[+-]?\d+(?:[.,]\d+)?$")
INSEE=re.compile(r"^(?:\d{5}|2[AB]\d{3})$")
EXPLICIT_INSEE_IN_LABEL=re.compile(r"(?:\(|\[)\s*((?:\d{5}|2[AB]\d{3}))\s*(?:\)|\])")

def sha(raw):return hashlib.sha256(raw).hexdigest()
def decode(raw):
    for enc in ("utf-8-sig","utf-8","cp1252"):
        try:return raw.decode(enc),enc
        except UnicodeDecodeError:pass
    raise ValueError("Unknown official CSV encoding")
def numeric(x):
    x=(x or "").strip().replace("\u00a0","").replace("\u202f","").replace(" ","")
    if not x or not RE_NUM.fullmatch(x):return None
    try:z=Decimal(x.replace(",","."))
    except InvalidOperation:return None
    return z if z.is_finite() and z>=0 else None
def type_commune(raw):
    s=(raw or "").strip()
    if not s:return "empty"
    if INSEE.fullmatch(s.upper()):return "exact_INSEE_code"
    if EXPLICIT_INSEE_IN_LABEL.search(s.upper()):return "label_plus_INSEE_code"
    if re.search(r"[;,/\n]",s):return "multi_or_delimited_text_UNRESOLVED"
    if s.isdigit():return "other_numeric_NOT_validated_INSEE"
    return "text_name_or_unknown"
def normalize(x):
    return unicodedata.normalize("NFKC",(x or "")).strip().strip('"').strip().strip('"')
def classify(raw):
    rawtext,encoding=decode(raw)
    if "<html" in rawtext[:500].lower() or "aucune donn" in rawtext[:500].lower():
        raise ValueError("No official CSV records returned")
    io_obj=io.StringIO(rawtext)
    reader=csv.reader(io_obj,delimiter="|")
    headers=[normalize(h) for h in next(reader)]
    if len(headers)<30 or len(set(headers))!=len(headers):
        raise ValueError("Expected official 34-column pipe-separated exporter; headers invalid")
    missing=set(FIELD.values())-set(headers)
    if missing:raise ValueError("Required original government fields missing "+str(sorted(missing)))
    ix={key:headers.index(val) for key,val in FIELD.items()}
    count=0
    combos=Counter()
    stations=set()
    station_commune=set()
    station_material=set()
    param_and_sample=set()
    region=set();depts=set();bassin=set()
    spreading=Counter()
    production=Counter()
    commune_types=Counter()
    commune_values=set()
    commune_literal_per_station=defaultdict(set)
    per_station_quantity=defaultdict(set)
    per_station_production=defaultdict(set)
    per_site_fields=defaultdict(set)
    positive_total_and_commune=0
    negative_or_unparsed_qty=0
    explicit_insee=set()
    name_only=set()
    sample_years=Counter()
    sample_date_rows=0
    duplicates=Counter()
    for fields in reader:
        if not any(v.strip() for v in fields):continue
        if len(fields)!=len(headers):raise ValueError("Raw government CSV row width changed")
        count+=1
        station=fields[ix["station"]].strip()
        commune=fields[ix["commune"]].strip()
        qty=numeric(fields[ix["spreading"]])
        prod=numeric(fields[ix["produced"]])
        region.add(fields[ix["region"]].strip())
        depts.add(fields[ix["dept"]].strip())
        bassin.add(fields[ix["bassin"]].strip())
        stations.add(station)
        commune_types[type_commune(commune)]+=1
        if commune:commune_values.add(commune)
        if commune and station:station_commune.add((station,commune))
        if station and fields[ix["material"]].strip():station_material.add((station,fields[ix["material"]].strip()))
        if station:commune_literal_per_station[station].add(commune)
        if station:per_station_quantity[station].add(str(qty) if qty is not None else "missing")
        if station:per_station_production[station].add(str(prod) if prod is not None else "missing")
        if station:
            per_site_fields[station].add((fields[ix["sample_date"]].strip(),fields[ix["parameter"]].strip(),
                                         fields[ix["lot"]].strip()))
        if fields[ix["sample_date"]].strip():sample_date_rows+=1
        if fields[ix["sample_year"]].strip():sample_years[fields[ix["sample_year"]].strip()]+=1
        if station and fields[ix["parameter"]].strip():
            param_and_sample.add((station,fields[ix["sample_date"]].strip(),fields[ix["parameter"]].strip()))
        spreading["positive" if qty is not None and qty>0 else
                  "valid_zero" if qty==0 else "invalid_or_blank"]+=1
        production["positive" if prod is not None and prod>0 else
                   "valid_zero" if prod==0 else "invalid_or_blank"]+=1
        if qty is None:negative_or_unparsed_qty+=1
        if qty is not None and qty>0 and commune:positive_total_and_commune+=1
        if type_commune(commune)=="exact_INSEE_code":explicit_insee.add(commune.upper())
        elif type_commune(commune)=="label_plus_INSEE_code":
            match=EXPLICIT_INSEE_IN_LABEL.search(commune.upper())
            if match:explicit_insee.add(match.group(1))
        elif commune:name_only.add(commune)
        duplicates[(station,commune,str(qty),fields[ix["sample_date"]].strip(),fields[ix["parameter"]].strip())]+=1
    if count<=0:raise ValueError("No government rows")
    nonblank_station={z for z in stations if z}
    station_bad={s:len(q) for s,q in per_station_quantity.items() if len(q)>1 and s}
    station_prod_bad={s:len(q) for s,q in per_station_production.items() if len(q)>1 and s}
    # 2018 sample labels are NOT the application date.
    return {
      "raw_records":count,"raw_columns":len(headers),"encoding":encoding,
      "official_unique_stations_nonempty":len(nonblank_station),
      "official_distinct_station_commune_literal_pairs":len(station_commune),
      "distinct_source_commune_literal_strings":len(commune_values),
      "distinct_explicit_INSEE_codes_only_not_validated_against_official_INSEE_geometry":len(explicit_insee),
      "distinct_non_INSEE_commune_literals":len(name_only),
      "official_unique_station_material_pairs":len(station_material),
      "official_unique_station_sampledate_parameter":len(param_and_sample),
      "original_analytical_record_duplicate_signatures":sum(n-1 for n in duplicates.values() if n>1),
      "distinct_regions_in_original_station_records":len({s for s in region if s}),
      "distinct_departments_in_original_station_records":len({s for s in depts if s}),
      "distinct_bassins_in_original_station_records":len({s for s in bassin if s}),
      "raw_rows_with_sample_date":sample_date_rows,
      "sample_year_labels":dict(sample_years),
      "raw_quantity_positive_zero_invalid_or_blank":dict(spreading),
      "raw_production_positive_zero_invalid_or_blank":dict(production),
      "raw_rows_with_positive_station_level_agricultural_qty_AND_nonempty_commune":positive_total_and_commune,
      "raw_commune_literal_types":dict(commune_types),
      "stations_where_agricultural_spreading_quantity_varies_between_duplicate_rows":len(station_bad),
      "stations_where_production_quantity_varies_between_duplicate_rows":len(station_prod_bad),
      "station_count_with_multiple_commune_literal_values":sum(len(z)>1 for z in commune_literal_per_station.values()),
      "looks_like_one_row_per_independent_spreading_event":False,
      "commune_field_is_not_an_application_parcel_coordinate":True,
      "sample_year_is_not_assumed_to_be_spreading_year":True,
      "station_total_dry_mass_not_summed_over_analysis_sample_rows":True,
    }
def retrieve(year):
    errors=[]
    for base in BASES:
        c=requests.Session()
        c.headers.update({"User-Agent":"gesui-france-routeC-original-sludge-source-QA/20261119"})
        try:
            page=c.get(base+PAGE,timeout=(12,35))
            page.raise_for_status()
            if str(year) not in page.text:raise ValueError("Official page lacks requested year")
            js=c.get(base+JS,timeout=(12,35))
            js.raise_for_status()
            if sha(js.content)!=JSSHA:raise ValueError("Official exporter JS source changed")
            if 'fetch("/controllers/GetDataExportBoueHandler.php"' not in js.text:
                raise ValueError("Published official JS no longer names the public handler")
            for key in ("bassin","region","annee","format","tranche","departement"):
                if f'formData.append("{key}"' not in js.text:raise ValueError("Official form key absent "+key)
            body={"bassin":"ALL","region":"ALL","annee":str(year),"format":"csv",
                  "tranche":"ALL","departement":"ALL"}
            res=c.post(base+HANDLER,data=body,timeout=(15,110),
                        headers={"Referer":base+PAGE,"Origin":base,"Accept":"text/csv,text/plain,*/*"})
            res.raise_for_status()
            raw=res.content
            if not MINBYTES<=len(raw)<=MAXBYTES:raise ValueError("Official export byte size outside gate")
            details=classify(raw)
            report={"year":str(year),"host":base,"official_public_handler_path":HANDLER,
                    "status":"ORIGINAL_CSV_RETRIEVED_QA","official_response_content_type":res.headers.get("Content-Type"),
                    "raw_csv_sha256":sha(raw),"raw_csv_bytes":len(raw),
                    "official_js_sha256":sha(js.content),
                    "official_page_sha256":sha(page.content),
                    "source_row_grain":details,
                    "earlier_2018_sha_matched":sha(raw)==PINNED_2018_SHA if str(year)=="2018" else None,
                    "true_independent_application_events_not_verified":True,
                    "source_missing_not_true_zero":True}
            if year=="2018" and not report["earlier_2018_sha_matched"]:
                report["warning"]="Source byte content changed relative to independently archived 2018 record: provenance drift"
            return report
        except Exception as e:
            errors.append({"government_hostname":base,"error":str(e)[:300]})
    return {"year":str(year),"status":"ACCESS_STOP","official_hosts_attempted":errors}
def main():
    a=argparse.ArgumentParser()
    a.add_argument("--year",choices=sorted(YEARS),required=True)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    result=retrieve(args.year)
    data={"schema":"gesui_FR_ROSEAU_sludge_routeC_2017to22_SOURCE_only_record_grain_20261119_v1",
          "UTC":datetime.now(timezone.utc).isoformat(),
          "year":args.year,"result":result,
          "no_raw_government_csv_row_republished":True,
          "no_biotic_EQR_NUTRI_HNS_values_read":True,
          "not_using_absence_of_commune_as_negative_control":True,
          "NUTS2_old_2x2_SOURCE_STOP_remains":True}
    (args.out/f"FR_ROSEAU_BOUES_{args.year}_SOURCE_GRAIN_QA.json").write_text(
      json.dumps(data,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_FR_SLUDGE_SOURCE_GRAIN",json.dumps({
      "year":args.year,"status":result["status"],
      "raw_rows":result.get("source_row_grain",{}).get("raw_records"),
      "stations":result.get("source_row_grain",{}).get("official_unique_stations_nonempty"),
      "commune_strings":result.get("source_row_grain",{}).get("distinct_source_commune_literal_strings"),
      "commune_format":result.get("source_row_grain",{}).get("raw_commune_literal_types"),
      "spreading_qty":result.get("source_row_grain",{}).get("raw_quantity_positive_zero_invalid_or_blank"),
      "sha_2018_match":result.get("earlier_2018_sha_matched")
    },ensure_ascii=False),flush=True)
if __name__=="__main__":main()
