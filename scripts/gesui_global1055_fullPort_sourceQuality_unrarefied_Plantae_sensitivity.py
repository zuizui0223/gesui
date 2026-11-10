#!/usr/bin/env python3
"""Posthoc reporting-eligibility selection bias audit, ALL original 1055 ports.

Old all-Plantae raw positive was reported on 752 GBIF 50-record eligible
sites; this audits raw indexed richness and counts on FULL 1055 successfully
queried WPI original source ports, including 303 low-coverage sites.
No new GBIF requests or new WWp source sites. No rarefaction at sparse ports.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

D=Path("results/global_all_Plantae_WWp2015_20261110")
SOURCE=D/"ORIGINAL_ALL1055_42COUNTRY_WWp2015_BEFORE_NEW808_GBIF_SOURCE.json"
PLANT=D/"ORIGINAL_ALL1055_42COUNTRY_GBIF_PLANTAE_2011_24_50RECORD_RAREFACTION_SOURCE.json"
PREV=D/"GLOBAL_ALL1055_WWp2015_ALL_PLANTAE_50RAREFIED_42COUNTRY_MACRO_AND_COUNTRY_HOLDOUT.json"
WWP="log_WWp15"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def fit(df,y):
    form=y+"~log_WWTP+WC_built+WC_crop+abs_lat+sinlon+coslon+"+WWP+"+C(country)"
    model=smf.ols(form,data=df).fit(cov_type="cluster",cov_kwds={"groups":df.country})
    return {
        "n_source_ports":len(df),
        "n_original_countries":int(df.country.nunique()),
        "beta_WWp":float(model.params[WWP]),
        "country_cluster95":list(map(float,model.conf_int().loc[WWP])),
        "nominal_p":float(model.pvalues[WWP]),
        "R2":float(model.rsquared),
        "formula":form,
        "not_true_plant_abundance_or_causal_sewage_effect":True
    }

def run(root):
    x=json.loads((root/SOURCE).read_text())
    gb=json.loads((root/PLANT).read_text())
    prev=json.loads((root/PREV).read_text())
    if gb["original_all1055_WWp_WPI_country_source_frame_SHA256"]!=sha(root/SOURCE):
        raise ValueError("All Plantae raw recording source changed")
    if prev["original_247old_plus_808new_allPlantae_GBIF_1055_result_SHA256"]!=sha(root/PLANT):
        raise ValueError("Original positive 752 models not on identical GBIF source")
    s={z["original_port_id"]:z for z in x["source_all_1055_public_country_port_site_metadata_no_new_GBIF_response"]}
    p={z["original_port_id"]:z for z in gb["original_1055_port_all_Plantae_source_GBIF_facet_summaries"]}
    if len(s)!=1055 or len(p)!=1055 or set(s)!=set(p):
        raise ValueError("Original WPI population not complete")
    rows=[]
    for key in sorted(s):
        a=s[key];q=p[key]
        if q["GBIF_status"]!="COMPLETE" or q["source_country"]!=a["original_source_country"]:
            raise ValueError("All1055 original source API or country incomplete")
        if q["indexed_GBIF_all_Plantae_2011_2024_records"]<q["faceted_original_speciesKey_records"] or q["faceted_original_speciesKey_records"]<q["faceted_original_distinct_Plantae_speciesKeys"]:
            raise ValueError("GBIF data impossible: mapped species keys more than indexed source records")
        lat=float(a["original_port_lat"]);lon=float(a["original_port_lon"])
        rows.append({
          "port_id":int(key),
          "country":a["original_source_country"],
          "source_plant50_quality":bool(q["source_quality_min50_faceted_records_min10_species_coverage70pct"]),
          WWP:float(a["source_WWp_log1p"]),
          "log_WWTP":math.log1p(float(a["source_WWTP_count_25km"])),
          "WC_built":float(a["source_ESA2021_built_fraction"]),
          "WC_crop":float(a["source_ESA2021_crop_fraction"]),
          "abs_lat":abs(lat),"sinlon":math.sin(math.radians(lon)),"coslon":math.cos(math.radians(lon)),
          "log_plant_speciesKeys":math.log1p(int(q["faceted_original_distinct_Plantae_speciesKeys"])),
          "log_plant_records":math.log1p(int(q["indexed_GBIF_all_Plantae_2011_2024_records"])),
          "GBIF_speciesKey_coverage":q["facet_record_coverage"],
        })
    df=pd.DataFrame(rows)
    assert len(df)==1055 and df.country.nunique()==42
    ok=df[df.source_plant50_quality].copy()
    low=df[~df.source_plant50_quality].copy()
    if len(ok)!=752 or len(low)!=303 or ok.country.nunique()!=38:
        raise ValueError("Original source plant-quality split no longer 752/303")
    allmodels={k:{} for k in ("unfiltered_1055","prior_GBIF50_supported_752","below_original_GBIF50_quality_303")}
    for k,sub in (("unfiltered_1055",df),("prior_GBIF50_supported_752",ok),("below_original_GBIF50_quality_303",low)):
        for y in ("log_plant_speciesKeys","log_plant_records"):
            allmodels[k][y]=fit(sub,y)
    for y,prevkey in (("log_plant_speciesKeys","raw_Plantae_speciesKey"),
                      ("log_plant_records","all_Plantae_GBIF_record_count")):
        a=allmodels["prior_GBIF50_supported_752"][y]
        original=prev["all1055_2015WWp_and_2011_24_all_Plantae_countryFE_models"][prevkey]
        if abs(a["beta_WWp"]-original["beta_on_log1p_modeled_WWp2015"])>1e-8 or any(abs(u-v)>1e-8 for u,v in zip(a["country_cluster95"],original["country_cluster95"])):
            raise ValueError("Original 752 site's raw positive beta and 95CI not exactly reproduced")
    influence={}
    for y in ("log_plant_speciesKeys","log_plant_records"):
        values=[]
        for country in sorted(df.country.unique()):
            rest=df[df.country!=country]
            m=smf.ols(y+"~log_WWTP+WC_built+WC_crop+abs_lat+sinlon+coslon+"+WWP+"+C(country)",data=rest).fit()
            values.append({"country_omitted":str(country),"WWp_beta":float(m.params[WWP])})
        influence[y]={
            "omitted_original_country_n":len(values),
            "positive_n":sum(v["WWp_beta"]>0 for v in values),
            "beta_range":[min(v["WWp_beta"] for v in values),max(v["WWp_beta"] for v in values)],
            "all_original_country_omission_results":values,
            "not_out_of_country_prediction":True
        }
    return {
        "schema":"gesui_global1055_42country_WWp_allPlantae_raw_GBIF_reporting_gate_fullsource_sensitivity_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "source_original_all1055_sites_SHA256":sha(root/SOURCE),
        "original_GBIF_all1055_source_SHA256":sha(root/PLANT),
        "original_previous_full752_macro_result_SHA256":sha(root/PREV),
        "all_1055_WPI_ports_in_source_42_countries_included_for_raw_record_analysis":True,
        "all_1055_original_GBIF_count_API_successful_without_requery":True,
        "quality_supported_752_old_beta_exactly_reproduced":True,
        "no_zero_imputation_of_failed_GBIF_queries":True,
        "GBIF_record_quality_752_subset_and_303_unsupported_is_posthoc_outcome_based_selection":True,
        "two_recorded_Plantae_endpoints_by_full_and_2_quality_strata":allmodels,
        "all_source_42country_leave_country_influence":influence,
        "no_50record_rarefaction_on_303_low_recording_sites":True,
        "2015_WWp_is_production_model_not_plant_nutrient_dose":True
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    result=run(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "results":{k:{name:{"beta":v["beta_WWp"],"ci":v["country_cluster95"],"n":v["n_source_ports"]}
                     for name,v in x.items()}
                   for k,x in result["two_recorded_Plantae_endpoints_by_full_and_2_quality_strata"].items()},
        "LOCO":{k:{n:v[n] for n in ("positive_n","beta_range")} for k,v in result["all_source_42country_leave_country_influence"].items()}
    },indent=2))
if __name__=="__main__":
    main()
