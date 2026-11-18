#!/usr/bin/env python3
"""USGS 2022 industrial effluent: true paired *river function* descriptive audit.

Published original Kotalik et al. 2025 tables at DOI 10.5066/P13BSPG4.
Fixed endpoints BEFORE opening numeric response values:
  (a) every published nutrient Parameter, sample day July5 July19 Aug2;
  (b) original decomposition mass-loss fraction, dates July19 Aug2;
  (c) original oxygen_time_mass, dates July14 Aug2, 3 original leaf/bag
      identities per site/date and 3 technical reps per bag (not n=9).
One source stream; *no independent WWTP treatment replicates*, no causal
N/P mediation or ecological service increase claim, no confirmatory p-values.
Do NOT confuse oxygen consumption / negative signed rate with productivity.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,statistics,time,urllib.error,urllib.request
from pathlib import Path
ITEM="68af0ddcd4be0253bb2c7c6b"
URL="https://www.sciencebase.gov/catalog/item/"+ITEM+"?format=json"
FILES={
 "nutrients":("SoyoilFoodWeb_Table13_Nutrients.txt","eafe63566c2f071dc42d21eac438ee56f1887f7f10f1900f8356e0661f3e95f8",81),
 "decomposition":("SoyoilFoodWeb_Table14_Decomposition.txt","058fa757dab0e2b55d41f197595d6cc8dfd0b29fd05ca82c50d68ef1411cfa8e",12),
 "metabolism":("SoyoilFoodWeb_Table15_Metabolism.txt","83aac8818151ce265ce1330ee391038d6832f99a28d8787afb5ebb2343448395",36)
}
SCHEMA="GESUI_20261118_USGS_INDUSTRIAL_2022_ECOSYSTEM_FUNCTION_PAIRED_SOURCE_DESCRIPTIVE_V1"
def sha(b):return hashlib.sha256(b).hexdigest()
def fetch(uri):
    if not uri.startswith("https://www.sciencebase.gov/"):raise RuntimeError("Source not official USGS")
    for attempt in range(3):
        try:
            req=urllib.request.Request(uri,headers={"User-Agent":"gesui-20261118-industrial-2022-function-original-source-audit","Accept":"application/json,text/plain,*/*"})
            with urllib.request.urlopen(req,timeout=75) as r:
                b=r.read(1_000_001)
                if len(b)>1_000_000:raise RuntimeError("Original source file too big")
                return b
        except (urllib.error.URLError,urllib.error.HTTPError,TimeoutError):
            if attempt==2:raise
            time.sleep(2+attempt*2)
def tofloat(s):
    try:
        v=float(s)
        return v if math.isfinite(v) else None
    except (TypeError,ValueError):return None
def load():
    meta=json.loads(fetch(URL))
    if meta.get("id")!=ITEM or len(meta.get("files",[]))!=16:raise RuntimeError("USGS source item changed")
    index={r["name"]:r for r in meta["files"]}
    out={}
    for name,(filename,original_sha,n) in FILES.items():
        raw=fetch(index[filename]["downloadUri"])
        if sha(raw)!=original_sha:raise RuntimeError("Original full raw source changed: "+name)
        rr=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace")),delimiter="\t"))
        if len(rr)!=n:raise RuntimeError("Original source data rows drift "+name)
        if len(rr[0])<4:raise RuntimeError("Original TAB source misparsed "+name)
        out[name]=rr
    return out
def fixed_gate(source):
    N=source["nutrients"];D=source["decomposition"];M=source["metabolism"]
    # predeclared source groups/dates; do NOT change after reading signs.
    n=collections.Counter((r["Site"],r["Sample_date"]) for r in N)
    d=collections.Counter((r["Site"],r["Sample_date"]) for r in D)
    m=collections.Counter((r["Site"],r["Sample_date"]) for r in M)
    N_dates=("20220705","20220719","20220802")
    D_dates=("20220719","20220802")
    M_dates=("20220714","20220802")
    for dt in N_dates:
        if any(n[s,dt]!=9 for s in ("Upstream","Outfall","Downstream")):
            raise RuntimeError("Original nutrient site+date+replicate observation source drift")
    for dt in D_dates:
        if any(d[s,dt]!=3 for s in ("Upstream","Downstream")):
            raise RuntimeError("Original three true decomposer bags per site/date drift")
    for dt in M_dates:
        if any(m[s,dt]!=9 for s in ("Upstream","Downstream")):
            raise RuntimeError("Original three bags by three technical O2 reps source drift")
    return {"nutrient":list(N_dates),"decomposition":list(D_dates),"metabolism":list(M_dates)}
def nutr(s,source_dates):
    rr=s["nutrients"]
    params=sorted(set(r["Parameter"] for r in rr))
    output=[]
    for p in params:
      p_rows=[r for r in rr if r["Parameter"]==p]
      units=set(r["Unit"] for r in p_rows)
      if len(units)!=1:raise RuntimeError("Original nutrient parameter unit inconsistent "+p)
      for dt in source_dates["nutrient"]:
        vals={}
        for site in ("Upstream","Outfall","Downstream"):
          raw=[r for r in p_rows if r["Site"]==site and r["Sample_date"]==dt]
          numbers=[v for r in raw if (v:=tofloat(r["Result_value"])) is not None]
          vals[site]={"n_original_source_replicates":len(raw),
                      "n_original_numeric_readable":len(numbers),
                      "mean":statistics.mean(numbers) if numbers else None,
                      "raw_replicate_codes":sorted(set(r["Replicate"] for r in raw))}
        a=vals["Upstream"]["mean"];b=vals["Downstream"]["mean"]
        output.append({"analyte":p,"date":dt,"reported_unit":list(units)[0],
           "upstream":vals["Upstream"],"outfall":vals["Outfall"],"downstream":vals["Downstream"],
           "down_minus_upstream_in_original_units":b-a if a is not None and b is not None else None})
    return output
def decomp(s,source_dates):
    rr=s["decomposition"];out=[]
    for dt in source_dates["decomposition"]:
      by={}
      for site in ("Upstream","Downstream"):
        q=[r for r in rr if r["Site"]==site and r["Sample_date"]==dt]
        values=[]
        for r in q:
          a=tofloat(r["Initial_mass"]);b=tofloat(r["Final_mass"])
          if a is None or a<=0 or b is None or b<0:continue
          values.append((a-b)/a)
        by[site]={"n_original_independent_bags":len(q),"n_valid_original_bags":len(values),
                  "mean_fraction_initial_mass_lost":statistics.mean(values) if values else None,
                  "bag_ids":sorted(r["Bag_number"] for r in q),
                  "source_original_mass_unit":sorted(set(r["Unit"] for r in q))}
      a=by["Upstream"]["mean_fraction_initial_mass_lost"];b=by["Downstream"]["mean_fraction_initial_mass_lost"]
      out.append({"date":dt,"upstream":by["Upstream"],"downstream":by["Downstream"],
        "down_minus_upstream_fraction_lost":b-a if a is not None and b is not None else None})
    return out
def metab(s,source_dates):
    rr=s["metabolism"];out=[]
    for dt in source_dates["metabolism"]:
      by={}
      for site in ("Upstream","Downstream"):
        q=[r for r in rr if r["Site"]==site and r["Sample_date"]==dt]
        bag=collections.defaultdict(list)
        for r in q:
          val=tofloat(r["Oxygen_time_mass"])
          if val is not None:bag[r["Bag_number"]].append(val)
        perbag=[statistics.mean(v) for v in bag.values() if v]
        by[site]={"n_original_assay_rows":len(q),
         "original_bag_id_n":len(set(r["Bag_number"] for r in q)),
         "n_bags_with_numeric_rate":len(perbag),
         "mean_of_bag_mean_reported_Oxygen_time_mass":statistics.mean(perbag) if perbag else None,
         "individual_bag_means_signed_as_reported":dict(sorted((k,statistics.mean(v)) for k,v in bag.items() if v)),
         "reported_value_sign_not_altered":True}
      a=by["Upstream"]["mean_of_bag_mean_reported_Oxygen_time_mass"]
      b=by["Downstream"]["mean_of_bag_mean_reported_Oxygen_time_mass"]
      out.append({"date":dt,"upstream":by["Upstream"],"downstream":by["Downstream"],
        "down_minus_upstream_original_signed_Oxygen_time_mass":b-a if a is not None and b is not None else None})
    return out
def run(args):
    data=load()
    source_dates=fixed_gate(data)
    n=nutr(data,source_dates);d=decomp(data,source_dates);m=metab(data,source_dates)
    report={"schema":SCHEMA,
       "status":"SOURCE_VERIFIED_2022_INDUSTRIAL_ONE_STREAM_UPDOWN_FUNCTION_DESCRIPTIVE_ONLY",
       "public_original_USGS_item":ITEM,
       "original_source_DOI":"10.5066/P13BSPG4",
       "original_2025_publication_DOI":"10.1021/acs.est.5c04867",
       "original_data_source_SHA256":{k:v[1] for k,v in FILES.items()},
       "original_sampled_physical_receiving_streams_n":1,
       "original_distinct_outfalls_in_this_study_n":1,
       "source_original_visit_dates":source_dates,
       "fixed_all_nutrient_analytes_not_direction_selected":True,
       "nutrient_source_concentration_by_3_locations_3_dates":n,
       "decomposition_fraction_initial_mass_lost_by_2_dates":d,
       "metabolism_signed_original_oxygen_rate_by_2_dates":m,
       "apparent_tech_replicates_not_treatment_independent":True,
       "signed_oxygen_measurement_not_claimed_to_be_GPP_or_total_NPP":True,
       "industrial_treatment_not_municipal_wastewater_exposure":True,
       "this_is_a_secondary_reanalysis_of_published_data_not_independent_reproduction":True,
       "only_two_temporal_river_function_pairs_not_confirmatory_ci_or_pvalues":True,
       "no_PFAS_or_spider_biological_outcome_numbers_used_in_function_exploration":True,
       "no_cross_domain_transfer_or_trophic_mediation_fit":True}
    args.out.mkdir(parents=True,exist_ok=True)
    path=args.out/"KOTALIK_2022_USGS_ONE_STREAM_REAL_NUTRIENT_METABOLISM_DECOMPOSITION_DESCRIPTIVE.json"
    path.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("KOTALIK_ORIGINAL_2022_FUNCTION="+json.dumps({
      "status":report["status"],
      "nutrient_analytes":sorted(set(q["analyte"] for q in n)),
      "nutrient_down_up_by_analyte":[{"name":p,"per_date":[q["down_minus_upstream_in_original_units"] for q in n if q["analyte"]==p]} for p in sorted(set(q["analyte"] for q in n))],
      "decomposition_down_up_fraction":[x["down_minus_upstream_fraction_lost"] for x in d],
      "reported_signed_oxygen_rate_down_up":[x["down_minus_upstream_original_signed_Oxygen_time_mass"] for x in m]},ensure_ascii=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())
