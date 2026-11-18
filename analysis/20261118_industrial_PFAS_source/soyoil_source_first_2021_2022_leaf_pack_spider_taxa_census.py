#!/usr/bin/env python3
"""Source-only 2021/2022 USGS SOYOIL physical outfall to insect and riparian
spider PFAS measurement opportunity. No analyte effects, no pseudo-replicate
"sample" inference from long-format chemical assay rows.

Treat 2025 foodweb and 2026 microbiome/leaf litter articles using this
same outfall as related studies, not three independent natural basins.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,re,time,urllib.request,urllib.error
from pathlib import Path

METADATA="https://www.sciencebase.gov/catalog/item/68af0ddcd4be0253bb2c7c6b?format=json"
META_SHA="0833b852675d08a83af304503d6c06f3cb7c8977af3d85948ed7b55364bcd409"
FILES={"2021":("SoyoilFoodWeb_Table6_2021Results.txt","2fdd7c259c0ac360253b3045d3cf6e0bd38a5fad22fefed73800560e9ca77d6e",1268),
 "2022":("SoyoilFoodWeb_Table1_2022Results.txt","05f09778feb45fad755c64216ac37b95246a80dd3c15ebcbbbf31155b78f8970",2839),
 "NUTRI":("SoyoilFoodWeb_Table13_Nutrients.txt","eafe63566c2f071dc42d21eac438ee56f1887f7f10f1900f8356e0661f3e95f8",81)}
SCHEMA="GESUI_20261118_SOYOIL_ORIGINAL_OUTFALL_UPDOWN_2021_22_PFAS_INSECT_SPIDER_LEAF_PACK_SOURCE_GATE_V1"
def sha(x):return hashlib.sha256(x).hexdigest()
def get(uri):
    if not uri.startswith("https://www.sciencebase.gov/"):raise RuntimeError("Non-USGS source")
    for i in range(3):
        try:
            req=urllib.request.Request(uri,headers={"User-Agent":"gesui-Soyoil-original-2021-22-insect-spider-source-gate/20261118"})
            with urllib.request.urlopen(req,timeout=90) as f:
                raw=f.read(7_000_000)
                if len(raw)>=7_000_000:raise RuntimeError("Unexpected original source size")
                return raw
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError):
            if i==2:raise
            time.sleep(1+i*2)

def source():
    meta=get(METADATA)
    if sha(meta)!=META_SHA:raise RuntimeError("USGS metadata SHA drift")
    f={v["name"]:v for v in json.loads(meta)["files"]}
    out={}
    for key,(name,signed,n) in FILES.items():
        raw=get(f[name]["downloadUri"])
        if sha(raw)!=signed:raise RuntimeError("Original published "+name+" sha drift")
        r=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace")),delimiter="\t"))
        if len(r)!=n:raise RuntimeError("Original source row count changed")
        out[key]=r
    return out

def run(args):
    tables=source()
    audit={}
    for year in ("2021","2022"):
        tab=tables[year]
        if set(["Site","Sample_datetime","Sample_type","Lifestage","Parameter","Unit","Remark","Result_value"]).difference(tab[0]):
            raise RuntimeError("Original assay column contract changed")
        cells=collections.Counter()
        full=collections.Counter()
        units=collections.defaultdict(set)
        analytes=collections.defaultdict(set)
        codes=collections.Counter()
        for r in tab:
            site=str(r["Site"]).strip()
            stage=str(r["Lifestage"]).strip()
            sample=str(r["Sample_type"]).strip()
            date=str(r["Sample_datetime"]).strip()[:8]
            k=(site,date,sample,stage)
            cells[k]+=1
            full[(site,sample,stage)]+=1
            units[(site,sample,stage)].add(str(r["Unit"]).strip())
            analytes[(site,sample,stage)].add(str(r["Parameter"]).strip())
            codes[(site,sample,stage,str(r["Remark"]).strip())]+=1
        audit[year]={
          "original_chemical_assay_rows_n":len(tab),
          "n_original_site_date_sample_type_lifestage_cells":len(cells),
          "original_assay_site_type_lifestage_support":[{"site":k[0],"sample_type":k[1],"lifestage":k[2],
              "original_analyte_assay_rows_not_independent_biological_replicates":v,
              "original_distinct_analyte_names_n":len(analytes[k]),
              "original_reported_units":sorted(units[k]),
              "original_sample_result_remarks":{q[3]:n for q,n in codes.items() if q[:3]==k}}
              for k,v in sorted(full.items())],
          "original_site_date_type_stage_observation_census":[
            {"site":s,"original_date":d,"sample_type":t,"lifestage":l,"n_analyte_rows":n}
            for (s,d,t,l),n in sorted(cells.items())],
          "original_chemical_concentration_values_not_read_or_tested":True}
    nutrient=collections.Counter((q["Site"],q["Sample_date"],q["Parameter"],q["Unit"]) for q in tables["NUTRI"])
    audit["nutrients"]={"original_n_assay_rows":81,
        "original_source_site_date_parameters":[{"site":k[0],"date":k[1],"parameter":k[2],"unit":k[3],"original_replicate_or_chemistry_rows":v}
          for k,v in sorted(nutrient.items())],
        "original_n_numeric_magnitudes_opened":0}
    source_sites={x["site"] for x in audit["2022"]["original_assay_site_type_lifestage_support"]}
    if source_sites!={"Upstream","Downstream","Outfall"}:
        raise RuntimeError("Unexpected original source PFAS site levels")
    nearstream=[x for x in audit["2022"]["original_assay_site_type_lifestage_support"]
      if x["sample_type"] in ("Leaf Pack","Tetragnathidae","Trichoptera","Chironomidae","Lepidoptera","Surface Water","Effluent")]
    result={"schema":SCHEMA,
      "status":"PASS_SOURCE_ONLY_PFAS_2021_2022_OUTFALL_TO_LEAF_PACK_TO_INSECTS_SPIDERS_SUPPORT",
      "USGS_original_data_DOI":"10.5066/P13BSPG4",
      "one_industrial_effluent_outfall_one_stream":True,
      "original_source_file_sha256":{v[0]:v[1] for v in FILES.values()},
      "source_years":["2021","2022"],
      "original_all_dates_stage_cells":audit,
      "selected_interesting_source_stage_support_but_without_outcome_numeric_read":nearstream,
      "source_2022_Lepidoptera_is_order_not_butterfly_or_larval_host_assay":True,
      "source_leaf_pack_chemical_assay_has_no_bag_number_for_join_to_original_12_decomposition_bags":True,
      "source_no_independent_multibasin_material_transfer_or_nutrient_NP_effect":True,
      "source_same_outfall_2025_foodweb_2026_microbiome_2026_litter_papers_not_three_independent_treatments":True}
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"SOYOIL_REAL_2021_2022_USGS_PFAS_LEAF_PACK_INSECT_SPIDER_SOURCE_SUPPORT.json").write_text(
       json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("SOYOIL_SOURCE_DUAL_DIRECTION="+json.dumps({
       "status":result["status"],"2022_assay_n":len(tables["2022"]),
       "2021_assay_n":len(tables["2021"]),
       "source_2022_selected_stage_rows":nearstream[:25],
       "source_2022_chemistry_nutrient_rows":81},ensure_ascii=False)[:7500],flush=True)
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    run(ap.parse_args())
