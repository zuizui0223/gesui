#!/usr/bin/env python3
"""USGS Kotalik industrial effluent 2021/2022 true PFAS data unit audit.

Source-only, SHA256-pinned original 2 annual PFAS assay files, 2022 nutrients,
metabolism and decomposition. Never treat rows of different analytes as
biological replicates; true specimen UUID is NOT in the published header.
Never treat the one industrial stream as municipal WWTP replicate.
Does not read response magnitudes and does not fit an ecological effect.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,re,time,urllib.error,urllib.parse,urllib.request
from pathlib import Path

ITEM="68af0ddcd4be0253bb2c7c6b"
API="https://www.sciencebase.gov/catalog/item/"+ITEM+"?format=json"
PUBDOI="10.5066/P13BSPG4"
SCHEMA="GESUI_20261118_KOTALIK_INDUSTRIAL_REAL_2021_22_SITE_DATE_STAGE_AND_FUNCTION_SOURCE_ONLY_V1"
ORIG_FILES={
 "SoyoilFoodWeb_Table1_2022Results.txt":("05f09778feb45fad755c64216ac37b95246a80dd3c15ebcbbbf31155b78f8970",2839),
 "SoyoilFoodWeb_Table6_2021Results.txt":("2fdd7c259c0ac360253b3045d3cf6e0bd38a5fad22fefed73800560e9ca77d6e",1268),
 "SoyoilFoodWeb_Table13_Nutrients.txt":("eafe63566c2f071dc42d21eac438ee56f1887f7f10f1900f8356e0661f3e95f8",81),
 "SoyoilFoodWeb_Table14_Decomposition.txt":("058fa757dab0e2b55d41f197595d6cc8dfd0b29fd05ca82c50d68ef1411cfa8e",12),
 "SoyoilFoodWeb_Table15_Metabolism.txt":("83aac8818151ce265ce1330ee391038d6832f99a28d8787afb5ebb2343448395",36),
 "SoyoilFoodWeb_Table11_Fieldparameters.txt":("cfd0831eef6e686386d3ed938862066bc9d69854d0339c3d22204386ad3ecf02",18)
}
ORIG_COLS="Site Sample_datetime Sample_type Lifestage Parameter CAS Result_value Remark Reporting_level Unit Report_level_code Method".split()
MAX=1_000_000
def sha(b):return hashlib.sha256(b).hexdigest()
def fetch(url):
    if not url.startswith("https://www.sciencebase.gov/"):raise RuntimeError("Not the original USGS catalog file URL")
    last=None
    for attempt in range(3):
        try:
            r=urllib.request.Request(url,headers={
              "User-Agent":"gesui-industrial-2021-22-source-true-taxonomic-stage-audit/20261118",
              "Accept":"application/json,text/plain,*/*"})
            with urllib.request.urlopen(r,timeout=75) as h:
                b=h.read(MAX+1)
                if len(b)>MAX:raise RuntimeError("Unexpectedly large source file")
                return b
        except (urllib.error.URLError,urllib.error.HTTPError,TimeoutError) as e:
            last=e
            if isinstance(e,urllib.error.HTTPError) and e.code not in (429,500,502,503,504):
                break
            time.sleep(2+attempt*3)
    raise RuntimeError("Original USGS source unavailable:"+str(last)[:180])

def read_table(meta,name):
    f=meta[name]
    b=fetch(f["downloadUri"])
    expected,n=ORIG_FILES[name]
    if sha(b)!=expected:raise RuntimeError("USGS original 2021-22 source file SHA changed "+name)
    text=b.decode("utf-8-sig","replace")
    lines=list(csv.DictReader(io.StringIO(text),delimiter="\t"))
    if len(lines)!=n:raise RuntimeError("Original source file expected n differs "+name)
    if len(lines[0])<=1:raise RuntimeError("Original tab-separated header parse failed "+name)
    return lines
def date(s):
    s=str(s or "").strip()
    if not re.match(r"^20[12]\d(?:0[1-9]|1[0-2])(?:0[1-9]|[12]\d|3[01])",s):
        return None
    return s[:8]

def scan_biology(rr,year):
    fields=list(rr[0])
    if fields!=ORIG_COLS:raise RuntimeError("Original USGS biological PFAS fields drift")
    by=collections.Counter()
    analytes=collections.defaultdict(set)
    dates=collections.defaultdict(set)
    nonreport=collections.Counter()
    missing=collections.Counter()
    for r in rr:
        loc=r["Site"].strip()
        lifestage=r["Lifestage"].strip()
        typ=r["Sample_type"].strip()
        d=date(r["Sample_datetime"])
        if d is None:missing["invalid_or_missing_collection_date"]+=1
        if d and d[:4]!=str(year):raise RuntimeError("Mixed original calendar years in same named annual file")
        key=(loc,typ,lifestage)
        by[key]+=1
        dates[key].add(d)
        analytes[key].add(r["Parameter"].strip())
        nonreport[(loc,r["Remark"].strip(),r["Report_level_code"].strip())]+=1
        if not str(r["Result_value"]).strip() or r["Result_value"].strip().upper() in ("NA","N/A"):
            missing["value_missing"]+=1
    support=[]
    for (loc,typ,stage),n in sorted(by.items()):
        good=sorted(d for d in dates[loc,typ,stage] if d is not None)
        support.append({"original_file_year":year,"source_location":loc,
          "source_sample_type":typ,"source_lifestage":stage,
          "original_analyte_assay_rows_not_biological_replicates":n,
          "original_distinct_PFAS_analytes":len(analytes[loc,typ,stage]),
          "original_distinct_collection_dates":len(good),
          "source_calendar_dates":good,
          "non_unique_field_specimen_ID_NOT_recorded":True})
    return {"year":year,"n_analyte_assays":len(rr),
        "original_n_distinct_location_type_stage_triplets":len(support),
        "source_site_type_lifestage_counts":support,
        "original_nonreport_flag_counts":[{"site":a,"remark":b,"report_code":c,"n":v}
          for (a,b,c),v in sorted(nonreport.items())],
        "source_analyte_name_levels":sorted({r["Parameter"].strip() for r in rr}),
        "source_sample_unit_levels":sorted({r["Unit"].strip() for r in rr}),
        "source_missing":dict(missing),
        "no_sample_uuid_and_per_analyte_rows_not_independent_replicates":True}
def scan_function(rows,name):
    col=list(rows[0])
    loccol="Site"
    if loccol not in col:raise RuntimeError("Published source functional table lacks Site column")
    dtcol=next((x for x in ("Sample_date","Date") if x in col),None)
    if not dtcol:raise RuntimeError("Published source function table lacks collection date")
    key=collections.Counter()
    for r in rows:key[(r["Site"].strip(),str(r[dtcol]).strip())]+=1
    return {"file":name,"rows":len(rows),"original_columns":col,
      "original_source_location_date_counts":[{"location":a,"date":b,"n_rows":n}
       for (a,b),n in sorted(key.items())],
      "functional_response_numbers_not_read":True,
      "up_down_measures_are_one_industrial_river_not_multibasin":True}
def run(args):
    args.out.mkdir(parents=True,exist_ok=True)
    result={"schema":SCHEMA,"original_USGS_DOI":PUBDOI,"original_item":ITEM,
      "industrial_food_processing_plant_not_municipal_WWTP":True,
      "physical_independent_impacted_river_reach_network_n":1,
      "source_2021_2022_numerical_response_NOT_OPENED":True}
    b=fetch(API);source=json.loads(b)
    if source.get("id")!=ITEM or len(source.get("files",[]))!=16:
        raise RuntimeError("Original USGS immutable item/16-source-file contract changed")
    f={x["name"]:x for x in source["files"]}
    if not set(ORIG_FILES).issubset(f):raise RuntimeError("Original source file set changed")
    result["original_full_metadata_sha256"]=sha(b)
    result["source_file_sha256"]={k:v[0] for k,v in ORIG_FILES.items()}
    data={x:read_table(f,x) for x in ORIG_FILES}
    bio={
      "2021":scan_biology(data["SoyoilFoodWeb_Table6_2021Results.txt"],2021),
      "2022":scan_biology(data["SoyoilFoodWeb_Table1_2022Results.txt"],2022)}
    func={
      "nutrient":scan_function(data["SoyoilFoodWeb_Table13_Nutrients.txt"],"Nutrients"),
      "decomposition":scan_function(data["SoyoilFoodWeb_Table14_Decomposition.txt"],"Decomposition"),
      "metabolism":scan_function(data["SoyoilFoodWeb_Table15_Metabolism.txt"],"Metabolism"),
      "field_parameters":scan_function(data["SoyoilFoodWeb_Table11_Fieldparameters.txt"],"Fieldparameters")}
    result["status"]="PASS_TRUE_2021_2022_INDUSTRIAL_SOURCE_SITE_DATE_STAGE_QA_NO_NUMERIC_EFFECT"
    result["bio_source_stages"]=bio
    result["ecosystem_function_source_stages"]=func
    outfile=args.out/"KOTALIK_2021_2022_USGS_ORIGINAL_INDUSTRIAL_PFAS_BIOTA_FUNCTION_SOURCE_QA.json"
    outfile.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("KOTALIK_2021_22_INDUSTRIAL_SOURCE="+json.dumps({
      "status":result["status"],
      "biota_2021_n":bio["2021"]["n_analyte_assays"],
      "biota_2022_n":bio["2022"]["n_analyte_assays"],
      "2021_categories":[(x["source_location"],x["source_sample_type"],x["source_lifestage"],
           x["original_distinct_collection_dates"]) for x in bio["2021"]["source_site_type_lifestage_counts"]],
      "2022_categories":[(x["source_location"],x["source_sample_type"],x["source_lifestage"],
           x["original_distinct_collection_dates"]) for x in bio["2022"]["source_site_type_lifestage_counts"]],
      "function_2022_n":{k:len(v["original_source_location_date_counts"]) for k,v in func.items()}},
      ensure_ascii=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())
