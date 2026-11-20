#!/usr/bin/env python3
"""Reanalyse original 2022 field leaf litter processing at the USGS SOYOIL
industrial effluent upstream/downstream reaches, source SHA-locked.

Replicates prior published Zachritz et al. (2026) field observation:
original leafpack mass loss (n=12 bags), respiration (36 technical readings
nested in the SAME 12 bags), water nutrient assay source (81 analyte rows).
It is land→water leaf litter resource processing, not water→terrestrial subsidy.
One stream, one industrial outfall, only 2 original timepoints. Do not
claim independent new experimental effect, independent WWTP replication,
or nutrient-specific PFAS causality from the observational field contrast.
"""
from __future__ import annotations
import argparse,csv,collections,hashlib,io,json,math,statistics,time,urllib.error,urllib.request
from pathlib import Path

URL="https://www.sciencebase.gov/catalog/item/68af0ddcd4be0253bb2c7c6b?format=json"
METADATA_SHA="0833b852675d08a83af304503d6c06f3cb7c8977af3d85948ed7b55364bcd409"
FILES={
 "PFAS_2022":("SoyoilFoodWeb_Table1_2022Results.txt",
  "05f09778feb45fad755c64216ac37b95246a80dd3c15ebcbbbf31155b78f8970",2839),
 "NUTRIENTS":("SoyoilFoodWeb_Table13_Nutrients.txt",
  "eafe63566c2f071dc42d21eac438ee56f1887f7f10f1900f8356e0661f3e95f8",81),
 "LITTER":("SoyoilFoodWeb_Table14_Decomposition.txt",
  "058fa757dab0e2b55d41f197595d6cc8dfd0b29fd05ca82c50d68ef1411cfa8e",12),
 "RESPIRATION":("SoyoilFoodWeb_Table15_Metabolism.txt",
  "83aac8818151ce265ce1330ee391038d6832f99a28d8787afb5ebb2343448395",36),
 "FIELD":("SoyoilFoodWeb_Table11_Fieldparameters.txt",
  "cfd0831eef6e686386d3ed938862066bc9d69854d0339c3d22204386ad3ecf02",18)
}
SCHEMA="GESUI_20261118_SOYOIL_ONE_OUTFALL_2022_ORIGINAL_LITTER_LOSS_RESPIRATION_3BAGS2DATES_EXPLORATORY_V1"
EXPECTED={"20220719":{"Downstream":3,"Upstream":3},
          "20220802":{"Downstream":3,"Upstream":3}}
RESP_DATES={"20220719":"20220714","20220802":"20220802"}

def sha(b):return hashlib.sha256(b).hexdigest()
def get(url,max_bytes=6_000_000):
    if not url.startswith("https://www.sciencebase.gov/"):
        raise RuntimeError("Only authentic USGS public source files")
    err=None
    for i in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"gesui-20261118-soyoil-authentic-source-ecoprocess-analysis",
                "Accept":"application/octet-stream,*/*"})
            with urllib.request.urlopen(req,timeout=90) as r:
                b=r.read(max_bytes+1)
                if len(b)>max_bytes:raise RuntimeError("Source original file over size cap")
                return b
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
            err=e
            if isinstance(e,urllib.error.HTTPError) and e.code not in (429,500,502,503,504):break
            time.sleep(2+i*3)
    raise RuntimeError(f"Published source retrieval failed: {str(err)[:250]}")

def fetchall():
    m=get(URL)
    if sha(m)!=METADATA_SHA:raise RuntimeError("Original USGS metadata changed after source QA, stop")
    md=json.loads(m)
    if md.get("id")!="68af0ddcd4be0253bb2c7c6b":raise RuntimeError("Wrong source dataset ID")
    fs={f["name"]:f for f in md["files"]}
    out={}
    for k,(name,signed_sha,n) in FILES.items():
        if name not in fs:raise RuntimeError("Original true file absent "+name)
        raw=get(fs[name]["downloadUri"])
        if sha(raw)!=signed_sha:raise RuntimeError("Original published source SHA mismatch "+name)
        rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace")),delimiter="\t"))
        if len(rows)!=n:raise RuntimeError(f"Original {name} true row count changed: {len(rows)} != {n}")
        if not rows or any("\t" in key for key in rows[0]):raise RuntimeError("Source delimiter/schema mismatch")
        out[k]=rows
        print("SOYOIL_SOURCE_VERIFIED="+json.dumps({"file":name,"row_n":len(rows),"sha256":signed_sha}),flush=True)
    return out

def numeric(v,key):
    text=str(v.get(key,"")).strip()
    if not text or text.lower() in ("na","nan"):raise RuntimeError("Unexpected missing numeric original "+key)
    x=float(text)
    if not math.isfinite(x):raise RuntimeError("Non-finite original numeric "+key)
    return x

def primary(data):
    leaf=data["LITTER"];res=data["RESPIRATION"]
    cols_leaf={"Site","Sample_date","Bag_number","Initial_mass","Final_mass","Unit"}
    cols_resp={"Site","Sample_date","Bag_number","Replicate","Oxygen_time_mass","Method"}
    if not cols_leaf.issubset(leaf[0]) or not cols_resp.issubset(res[0]):
        raise RuntimeError("Original true bag/oxygen columns missing")
    leaf_index={}; per_date=collections.Counter()
    for r in leaf:
        site=r["Site"];date=r["Sample_date"];bag=r["Bag_number"]
        if (site,date,bag) in leaf_index:raise RuntimeError("Duplicate source leafpack bag")
        if r["Unit"]!="grams" or (site not in ("Upstream","Downstream")):
            raise RuntimeError("Unknown original mass units/site")
        start=numeric(r,"Initial_mass");end=numeric(r,"Final_mass")
        if not (start>0 and 0<=end<=start):raise RuntimeError("Invalid original bag mass")
        per_date[(date,site)]+=1
        leaf_index[(site,date,bag)]={
           "site":site,"leaf_harvest_date":date,"bag_number":bag,
           "initial_mass_original_grams":start,"final_mass_original_grams":end,
           "fraction_mass_remaining":end/start,
           "fraction_mass_lost":(start-end)/start}
    if {d:{s:per_date[(d,s)] for s in ("Upstream","Downstream")} for d in EXPECTED}!=EXPECTED:
        raise RuntimeError("Original site/date litter bag support drift")
    if sum(per_date.values())!=12:raise RuntimeError("Original litter source has unexpected dates or sites")
    metabolic=collections.defaultdict(list)
    resp_dates=collections.Counter()
    for r in res:
        site=r["Site"];date=r["Sample_date"];bag=r["Bag_number"]
        if site not in ("Upstream","Downstream"):raise RuntimeError("Unexpected respiration source site")
        if str(r.get("Method"))!="Bott 2006":raise RuntimeError("Original respiration method changed")
        d=numeric(r,"Oxygen_time_mass")
        if d>=0:raise RuntimeError("Oxygen_time_mass was expected negative as oxygen decline")
        key=(site,bag)
        metabolic[key].append({"date":date,"technical_replicate":r["Replicate"],"oxygen_decline_magnitude":-d})
        resp_dates[(date,site)]+=1
    if len(metabolic)!=12 or any(len(a)!=3 for a in metabolic.values()):
        raise RuntimeError("Oxygen 36 readings do not give exactly 12 bags × 3 technical replicates")
    if {(site,d) for d in EXPECTED for site in ("Upstream","Downstream")} != {(site,date) for (site,date,bag) in leaf_index}:
        raise RuntimeError("Bag-date expected source discrepancy")
    output=[]
    for (_site,harvest_date,bag),orig in sorted(leaf_index.items()):
        key=(_site,bag)
        q=metabolic.get(key,[])
        if len(q)!=3:raise RuntimeError("Source bag mass/respiration ID overlap not exact")
        rd={x["date"] for x in q}
        if rd!={RESP_DATES[harvest_date]}:
            raise RuntimeError("Respiration vs bag harvest date differs from published 14/19 and 2/2 source")
        ids=[x["technical_replicate"] for x in q]
        if len(set(ids))!=3:raise RuntimeError("Oxygen technical replicate names duplicate")
        output.append({**orig,
           "respiration_sampling_date":next(iter(rd)),
           "respiration_replicate_n":len(q),
           "mean_original_oxygen_time_mass_negative":-statistics.mean(x["oxygen_decline_magnitude"] for x in q),
           "source_oxygen_decline_magnitude_per_mass":statistics.mean(x["oxygen_decline_magnitude"] for x in q),
           "oxygen_units_original_source_column_not_independently_confirmed":True})
    return output

def summary(bags):
    date_records=[]
    for harvest in sorted(EXPECTED):
        allgroups={}
        for site in ("Upstream","Downstream"):
            q=[r for r in bags if r["leaf_harvest_date"]==harvest and r["site"]==site]
            if len(q)!=3:raise RuntimeError("Source bag count changed at date")
            allgroups[site]={
              "n_original_distinct_bags":3,
              "n_original_nested_respiration_lab_replicates":9,
              "mean_fraction_mass_lost":statistics.mean(x["fraction_mass_lost"] for x in q),
              "sd_fraction_mass_lost":statistics.stdev(x["fraction_mass_lost"] for x in q),
              "mean_oxygen_decline_magnitude_per_mass":statistics.mean(x["source_oxygen_decline_magnitude_per_mass"] for x in q),
              "sd_between_bags_oxygen_decline_magnitude_per_mass":statistics.stdev(x["source_oxygen_decline_magnitude_per_mass"] for x in q)}
        u=allgroups["Upstream"];d=allgroups["Downstream"]
        date_records.append({"leaf_harvest_date":harvest,
          "respiration_date":RESP_DATES[harvest],
          "groups":allgroups,
          "DS_minus_US_fraction_mass_lost":d["mean_fraction_mass_lost"]-u["mean_fraction_mass_lost"],
          "DS_minus_US_oxygen_decline_magnitude_per_mass":d["mean_oxygen_decline_magnitude_per_mass"]-u["mean_oxygen_decline_magnitude_per_mass"],
          "DS_to_US_leaf_mass_loss_ratio_if_positive":d["mean_fraction_mass_lost"]/u["mean_fraction_mass_lost"],
          "DS_to_US_oxygen_decline_magnitude_ratio_if_positive":d["mean_oxygen_decline_magnitude_per_mass"]/u["mean_oxygen_decline_magnitude_per_mass"]})
    return date_records

def source_environment(data):
    # Only structural source metrics; no PFAS effect chosen after outcome.
    pfas=data["PFAS_2022"]; nut=data["NUTRIENTS"]; field=data["FIELD"]
    need=[("PFAS_2022",pfas,("Site","Sample_datetime","Sample_type","Lifestage","Parameter","Result_value","Remark","Unit")),
          ("NUTRIENTS",nut,("Site","Sample_date","Parameter","Result_value","Unit")),
          ("FIELD",field,("Site","Date","Water_temperature","Specific_conductance","Discharge"))]
    out={}
    for name,rs,fields in need:
        if not set(fields).issubset(rs[0]):raise RuntimeError("Source environmental record fields missing "+name)
        counts=collections.Counter(str(r.get("Site") or "") for r in rs)
        out[name]={
          "original_assay_rows_n":len(rs),
          "original_site_label_counts":dict(counts),
          "source_2022_parameter_levels_n":len(set(r.get("Parameter") for r in rs if r.get("Parameter"))) if name!="FIELD" else None,
          "source_2022_stage_type_levels":dict(collections.Counter(str(r.get("Sample_type") or "") for r in rs))
                    if name=="PFAS_2022" else None,
          "original_row_is_analyte_or_field_observation_not_independent_wwtp":True}
    return out

def run(args):
    args.out.mkdir(parents=True,exist_ok=True)
    original=fetchall()
    bags=primary(original)
    agg=summary(bags)
    env=source_environment(original)
    status="PASS_ORIGINAL_ONE_INDUSTRIAL_OUTFALL_LAND_TO_WATER_BAG_RESPIRATION_EXPLORATORY"
    r={"schema":SCHEMA,"status":status,
       "source_DOI":"10.5066/P13BSPG4",
       "published_2025_PFAS_foodweb_DOI":"10.1021/acs.est.5c04867",
       "published_2026_leaf_litter_DOI":"10.1016/j.aquatox.2026.107966",
       "year_2021_2022_same_industrial_outfall_NOT_independent_published_study_units":True,
       "original_Source_Pinned_SHA256":{k:x[1] for k,x in FILES.items()},
       "original_mass_loss_bags_n":len(bags),
       "original_respiration_technical_readings_n":len(original["RESPIRATION"]),
       "independent_river_source_areas_n":1,
       "original_independent_exposure_reaches_US_DS_n":2,
       "original_leaf_harvest_dates_n":2,
       "original_oxygen_respiration_dates_n":2,
       "respiration_and_harvest_July_differ_by_five_days":True,
       "observed_original_bag_date_and_group_descriptive_statistics":agg,
       "original_12_bags_with_mean_of_3_technical_O2_readings":bags,
       "original_environment_2022_source_structure":env,
       "no_2022_PFAS_dose_effect_or_nutrient_specific_causal_model_estimated":True,
       "not_independent_multiriver_meta_analysis_or_terrestrial_biomass_exposure":True,
       "measured_land_origin_leaf_resources_are_being_processed_IN_water":True,
       "original_paper_laboratory_28d_PFOS_mesocosm_is_NOT_represented_by_USGS_field_csv":True}
    (args.out/"SOYOIL_2022_AUTHENTIC_LEAF_TO_STREAM_DECOMPOSITION_RESPIRATION_DESCRIPTIVE.json").write_text(
        json.dumps(r,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    with (args.out/"SOYOIL_2022_ORIGINAL_12_BAG_MASS_LOSS_VS_RESPIRATION.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(bags[0]));w.writeheader();w.writerows(bags)
    print("SOYOIL_ORIGINAL_LAND_TO_WATER_FUNCTION="+json.dumps({
      "status":status,"n_bags":12,"n_oxygen_reads":36,
      "results":[{"date":x["leaf_harvest_date"],
          "DS_minus_US_fraction_mass_lost":x["DS_minus_US_fraction_mass_lost"],
          "DS_minus_US_oxygen_decline_magnitude":x["DS_minus_US_oxygen_decline_magnitude_per_mass"]}
          for x in agg]},ensure_ascii=False),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    run(ap.parse_args())
