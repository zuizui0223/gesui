#!/usr/bin/env python3
"""True same-reach, same-date industrial wastewater nutrients and terrestrial
leaf litter processing in stream: USGS SOYOIL source 2022.

Source-lock published 81 original nutrient chemistry rows and 12 litter
mass bags, both original July19/Aug02 observations. Different 2026 and
2025 papers reuse this exact ONE industrial outfall: never an independent
replication, municipal sewage, plant→butterfly subsidy or WWTP×year model.
Raw original 3 replicates per Site×Date×analyte are analytical samples;
within-date 3 litter bags per reach do not give independent rivers.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,statistics,time,urllib.error,urllib.request
from pathlib import Path

URL="https://www.sciencebase.gov/catalog/item/68af0ddcd4be0253bb2c7c6b?format=json"
META="0833b852675d08a83af304503d6c06f3cb7c8977af3d85948ed7b55364bcd409"
FILES={
 "nutrient":("SoyoilFoodWeb_Table13_Nutrients.txt",
  "eafe63566c2f071dc42d21eac438ee56f1887f7f10f1900f8356e0661f3e95f8",81),
 "leaf":("SoyoilFoodWeb_Table14_Decomposition.txt",
  "058fa757dab0e2b55d41f197595d6cc8dfd0b29fd05ca82c50d68ef1411cfa8e",12)
}
DATES=("20220705","20220719","20220802")
OVERLAP=("20220719","20220802")
SITES=("Upstream","Downstream","Outfall")
PARAMS={"NH4-N":"micrograms nitrogen per liter",
        "NO3-N":"milligrams nitrogen per liter",
        "SRP":"micrograms phosphorus per liter"}
SCHEMA="GESUI_20261118_SOYOIL_TRUE_2022_3LOC_3DATE_3NUTRIENTS_AND_2DATE_LITTER_FUNCTION_V1"

def sha(b):return hashlib.sha256(b).hexdigest()
def get(url):
    if not url.startswith("https://www.sciencebase.gov/"):raise RuntimeError("Not trusted USGS")
    e=None
    for i in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"gesui-genuine-Soyoil-2022-original-water-leaf-source-alignment/1"})
            with urllib.request.urlopen(req,timeout=90) as r:
                b=r.read(2_000_001)
                if len(b)>2_000_000:raise RuntimeError("Source bytes over hard cap")
                return b
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as ex:
            e=ex
            if isinstance(ex,urllib.error.HTTPError) and ex.code not in [429,500,502,503,504]:break
            time.sleep(2+i*3)
    raise RuntimeError(str(e)[:240])
def source():
    meta=get(URL)
    if sha(meta)!=META:raise RuntimeError("Pinned USGS original metadata source changed")
    fs={x["name"]:x for x in json.loads(meta)["files"]}
    data={}
    for k,(name,signed,n) in FILES.items():
        raw=get(fs[name]["downloadUri"])
        if sha(raw)!=signed:raise RuntimeError("Original byte SHA drift "+k)
        reader=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace")),delimiter="\t"))
        if len(reader)!=n:raise RuntimeError("Original source n drift "+k)
        data[k]=reader
    return data
def numeric(v):
    try:n=float(str(v).strip())
    except (ValueError,TypeError):return None
    return n if math.isfinite(n) else None

def audit_nutrient(rows):
    expected={(site,date,pa) for site in SITES for date in DATES for pa in PARAMS}
    groups=collections.defaultdict(list)
    for r in rows:
        site=r["Site"].strip();date=r["Sample_date"].strip();p=r["Parameter"].strip()
        if (site,date,p) not in expected or r["Unit"].strip()!=PARAMS[p]:
            raise RuntimeError("Original USGS source site-date-parameter or units changed")
        groups[(site,date,p)].append(numeric(r["Result_value"]))
    if set(groups)!=expected or any(len(v)!=3 for v in groups.values()):
        raise RuntimeError("Expected genuine 3source×3date×3 analyte×3 replicate 81-row support unavailable")
    # Prespecified source quality: >=2 original finite per critical site/date/chemistry
    if any(sum(x is not None for x in v)<2 for v in groups.values()):
        raise RuntimeError("Source finite genuine nutrient replicates under two; no effects fitted")
    data=[]
    for (site,date,p),v in sorted(groups.items()):
        q=[x for x in v if x is not None]
        data.append({"site":site,"date":date,"parameter":p,"original_unit":PARAMS[p],
          "n_original_chemistry_replicates":len(v),"n_original_numeric_replicates":len(q),
          "source_mean":statistics.mean(q),"source_sd":statistics.stdev(q) if len(q)>1 else None,
          "source_median":statistics.median(q)})
    deltas=[]
    for date in DATES:
        for p in PARAMS:
            vals={r["site"]:r for r in data if r["date"]==date and r["parameter"]==p}
            u,d,o=(*(vals[k]["source_mean"] for k in SITES),)
            deltas.append({"sample_date":date,"parameter":p,"unit":PARAMS[p],
              "upstream_original_mean":u,"downstream_original_mean":d,"outfall_original_mean":o,
              "DS_minus_US":d-u,"DS_to_US_ratio_if_original_U_positive":d/u if u>0 else None,
              "outfall_to_downstream_if_positive":o/d if d>0 else None,
              "original_upstream_downstream_data_same_year_same_calendar_date":True})
    return data,deltas

def audit_leaf(rows):
    cells=collections.defaultdict(list)
    bags=set()
    for r in rows:
        site,date,bag=r["Site"].strip(),r["Sample_date"].strip(),r["Bag_number"].strip()
        if site not in ("Upstream","Downstream") or date not in OVERLAP or r["Unit"]!="grams":
            raise RuntimeError("Original field leafbag site/date/unit changed")
        if (site,date,bag) in bags:raise RuntimeError("Duplicate true leaf bag")
        bags.add((site,date,bag))
        initial=numeric(r["Initial_mass"]);final=numeric(r["Final_mass"])
        if initial is None or final is None or not (initial>0 and 0<=final<=initial):
            raise RuntimeError("Original leaf pack final mass not physical")
        cells[(site,date)].append((initial-final)/initial)
    if len(bags)!=12 or set(cells)!={(s,d) for s in ("Upstream","Downstream") for d in OVERLAP} or any(len(v)!=3 for v in cells.values()):
        raise RuntimeError("Original real 12bag source before ecology outcome did not pass fixed support")
    out=[]
    for date in OVERLAP:
        upstream=statistics.mean(cells[("Upstream",date)])
        downstream=statistics.mean(cells[("Downstream",date)])
        out.append({"date":date,"mean_original_leaf_mass_lost_upstream":upstream,
          "mean_original_leaf_mass_lost_downstream":downstream,
          "DS_minus_US_fraction_mass_lost":downstream-upstream,
          "original_bags_each_side":3})
    return out

def run(args):
    data=source()
    nutrients,nd= audit_nutrient(data["nutrient"])
    leaves=audit_leaf(data["leaf"])
    match=[]
    for l in leaves:
        a=[r for r in nd if r["sample_date"]==l["date"]]
        if len(a)!=3:raise RuntimeError("One-to-many source date not three nutrient types")
        for r in a:match.append({**r,**l})
    if len(match)!=6:raise RuntimeError("Original same-day nutrient-leaf lineages not 2 dates ×3 parameters")
    result={"schema":SCHEMA,
      "status":"PASS_TRUE_SAME_DATE_2022_PUBLISHED_NUTRIENT_LAND_TO_WATER_LITTER_DESCRIPTIVE_ONE_OUTFALL",
      "published_source_USGS_DOI":"10.5066/P13BSPG4",
      "published_original_2026_leafprocessing_paper_DOI":"10.1016/j.aquatox.2026.107966",
      "source_units_exact":PARAMS,
      "USGS_original_source_sha256":{k:v[1] for k,v in FILES.items()},
      "n_original_chemistry_analyte_rows":81,
      "original_chemistry_source_3_sites_by_3_dates_by_3_params_by_3_replicates":True,
      "original_chemistry_true_calendar_dates":list(DATES),
      "original_litter_decomposition_measurement_dates":list(OVERLAP),
      "same_source_same_year_same_day_matched_nutrient_function_cells":6,
      "measured_original_water_nutrient_cell_values":nutrients,
      "original_same_day_site_paired_water_nutrient_contrasts":nd,
      "original_source_litter_fraction_mass_lost":leaves,
      "land_origin_leaf_litter_to_stream_is_reverse_of_aquatic_to_riparian_spider_source":True,
      "single_industrial_outfall_not_municipal_sewage":True,
      "original_replicate_n_81_does_not_mean_81_independent_WWTP":True,
      "two_date_association_not_sufficient_for_valid_wastewater_causal_mediation":True,
      "no_cross_compartment_actual_nutrient_flux_dose_estimated":True,
      "author_2026_field_result_already_published_not_independent_new_study":True}
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"SOYOIL_2022_TRUE_NP_SRPs_vs_SAME_DAY_LITTER_LOSS_SOURCE.json").write_text(
       json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    with (args.out/"SOYOIL_2022_NUTRIENT_LITTER_SAME_DAY_TWO_HARVESTS.csv").open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=list(match[0]));writer.writeheader();writer.writerows(match)
    print("SOYOIL_SAME_DAY_ORIGINAL_NUTRIENT_LEAF_PROCESSING="+json.dumps({
      "status":result["status"],"matched_date_n":2,
      "date_specific_down_minus_upstream_leaf_loss":{x["date"]:x["DS_minus_US_fraction_mass_lost"] for x in leaves},
      "nutrient_delta_by_date":[{"date":x["sample_date"],"param":x["parameter"],
         "DS_minus_US":x["DS_minus_US"],"unit":x["unit"]} for x in nd]},ensure_ascii=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())
