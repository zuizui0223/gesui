#!/usr/bin/env python3
"""EXPLORATORY country-matched Lycaenidae photo REPORTING fraction vs WWTP sources.

Crucial gate: read original source+2021 Plantae-only selection committed
before using any 2022 butterfly records, and STOP unless source eligibility
was sufficient. This is an opportunistic iNaturalist observer distribution,
NOT butterfly abundance, larval performance, puddling or sewage causality.
"""
from __future__ import annotations
import argparse,hashlib,json,math,time,urllib.parse,urllib.request,urllib.error
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

API="https://api.inaturalist.org/v1/observations"
FAMILIES={"Lycaenidae":47923,"Pieridae":48508,"Papilionidae":47223}
MIN_TOTAL_PHOTOS_PER_TILE=20
MIN_SUPPORTED_TILES=25
MIN_SUPPORTED_PER_COUNTRY=10
MIN_SUPPORTED_BLOCKS=10
YEAR="2022"
SEASON_START="2022-04-01"
SEASON_END="2022-10-31"
N_BOOTSTRAP=999
RNG_SEED=20261108

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for ch in iter(lambda:f.read(1048576),b""):h.update(ch)
    return h.hexdigest()

def api_total(tile,place_id,taxonid):
    params={
      "taxon_id":int(taxonid),"place_id":int(place_id),
      "swlat":float(tile["lat_floor"]),"swlng":float(tile["lon_floor"]),
      "nelat":float(tile["lat_ceiling"]),"nelng":float(tile["lon_ceiling"]),
      "d1":SEASON_START,"d2":SEASON_END,"photos":"true","geo":"true",
      "geoprivacy":"open","taxon_geoprivacy":"open",
      "quality_grade":"research","captive":"false","acc_below":500,
      "per_page":1}
    url=API+"?"+urllib.parse.urlencode(params)
    receipt={"source_url":url,"http_status":None,"total_results":None}
    for i in range(4):
        try:
            req=urllib.request.Request(url,headers={
              "User-Agent":"gesui-2022-DE-AT-lycaenid-photo-proportion-source-locked/0.1",
              "Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=35) as r:
                data=r.read(600000)
                j=json.loads(data)
                if not isinstance(j.get("total_results"),int):
                    raise ValueError("Missing public iNaturalist total_results")
                receipt.update({"http_status":r.status,
                     "total_results":int(j["total_results"]),
                     "source_JSON_SHA256":hashlib.sha256(data).hexdigest()})
                return receipt
        except urllib.error.HTTPError as e:
            receipt.update({"http_status":e.code,"error":"HTTP "+str(e.code)})
            if e.code not in (429,500,502,503) or i==3:break
        except Exception as e:
            receipt["error"]=type(e).__name__+": "+str(e)
            if i==3:break
        time.sleep(1.6*(i+1))
    return receipt

def calc_source(x,field):
    return np.log1p(np.asarray([float(v[field]) for v in x]))

def design(rows,source_field,need_Y=True):
    country=np.asarray([x["countryCode"] for x in rows])
    if not all(x in ("DE","AT") for x in country):raise ValueError("Unexpected country")
    src=calc_source(rows,source_field)
    plant=np.log1p(np.asarray([float(x["2021_Plantae_photo_count"]) for x in rows]))
    latitude=np.asarray([float(x["lat_floor"])+0.5 for x in rows])
    longitude=np.asarray([float(x["lon_floor"])+0.5 for x in rows])
    def std(arr):
        v=float(np.std(arr))
        return (arr-float(np.mean(arr)))/v if v>0 else np.zeros_like(arr)
    # Standardize source *within country* to ensure country-FE slope
    # reflects intra-country variation, not DE vs AT different monitoring.
    src_z=np.empty_like(src)
    for c in ("DE","AT"):
        mask=country==c
        s=float(np.std(src[mask]))
        if s<1e-10:raise ValueError("Within-country source predictor has no variation: "+source_field+" "+c)
        src_z[mask]=(src[mask]-np.mean(src[mask]))/s
    X=np.column_stack([np.ones(len(rows)),src_z,std(plant),std(latitude),std(longitude),(country=="AT").astype(float)])
    if need_Y:
        y=np.asarray([math.log((x["Lycaenidae_2022_photo_n"]+.5) /
                 (x["three_family_total_2022"]-x["Lycaenidae_2022_photo_n"]+.5)) for x in rows])
    else:y=None
    return X,y,{"predictor_columns":["intercept","within_country_z_source_log1p","z_log1p_2021_Plantae_photos","z_latitude","z_longitude","AT_indicator"],
                "source_predictor":source_field}

def fit(rows,field):
    X,Y,names=design(rows,field)
    rank=int(np.linalg.matrix_rank(X))
    if rank<X.shape[1]:
        return {"error":"Rank-deficient fixed model, no coefficient interpreted","rank":rank,"n":len(rows)}
    b=np.linalg.lstsq(X,Y,rcond=None)[0]
    err=Y-X@b
    return {"n":len(rows),"rank":rank,"source_coefficient":float(b[1]),
            "coefficients":{c:float(v) for c,v in zip(names["predictor_columns"],b)},
            "root_mean_sq_residual_descriptive":float(np.sqrt(np.mean(err**2))),
            "source_units":"per within-country one SD of log1p original EEA 2022 source predictor",
            "no_causal_meaning":True}

def result_summary(rows,source_fields):
    country=Counter(x["countryCode"] for x in rows)
    blocksets=set(x["countryCode"]+"|"+x["geographic_2degree_block"] for x in rows)
    primary=source_fields[0]
    output={
      "primary_support_n":len(rows),
      "country_n":dict(country),
      "independent_2degree_blocks_n":len(blocksets),
      "raw_total_2022_butterfly_family_photo_records":int(sum(x["three_family_total_2022"] for x in rows)),
      "raw_family_photo_counts":{f:sum(x[f+"_2022_photo_n"] for x in rows) for f in FAMILIES},
      "unweighted_lycaenidae_photo_reporting_share_mean":float(np.mean([x["Lycaenidae_2022_photo_n"]/x["three_family_total_2022"] for x in rows])),
      "country_summary":{},
      "sewage_causality_not_inferred":True,
    }
    for c in ("DE","AT"):
        rr=[x for x in rows if x["countryCode"]==c]
        if not rr:continue
        share=np.asarray([x["Lycaenidae_2022_photo_n"]/x["three_family_total_2022"] for x in rr])
        src=calc_source(rr,primary)
        rho,p=spearmanr(src,share)
        effort_rho,effort_p=spearmanr(src,[float(x["2021_Plantae_photo_count"]) for x in rr])
        output["country_summary"][c]={
          "n_cells":len(rr),"photo_share_median":float(np.median(share)),
          "photo_share_mean_unweighted":float(np.mean(share)),
          "source_logcount_to_photo_share_rho":float(rho) if np.isfinite(rho) else None,
          "source_logcount_to_photo_share_nominal_p_exploratory":float(p) if np.isfinite(p) else None,
          "source_logcount_to_2021_plant_effort_rho":float(effort_rho) if np.isfinite(effort_rho) else None,
          "sum_2022_three_family_photos":int(sum(x["three_family_total_2022"] for x in rr)),
        }
    output["three_predeclared_source_predictors"]={f:fit(rows,f) for f in source_fields}
    # 2 degree region bootstrap refits fixed model. 999 cluster resamples.
    group=defaultdict(list)
    for r in rows:group[r["countryCode"]+"|"+r["geographic_2degree_block"]].append(r)
    blocklist=list(group)
    rng=np.random.default_rng(RNG_SEED)
    boot=[]
    if len(blocklist)>=MIN_SUPPORTED_BLOCKS:
        for k in range(N_BOOTSTRAP):
            ids=rng.choice(blocklist,size=len(blocklist),replace=True)
            sample=[x for key in ids for x in group[key]]
            counts=Counter(r["countryCode"] for r in sample)
            if min(counts.get("AT",0),counts.get("DE",0))<5:continue
            try:
                b=fit(sample,primary).get("source_coefficient")
                if b is not None and np.isfinite(b):boot.append(b)
            except (ValueError,np.linalg.LinAlgError):continue
    output["primary_source_count_2degree_cluster_bootstrap"]={
       "n_successful_draws":len(boot),"seed":RNG_SEED,"n_requested":N_BOOTSTRAP,
       "95pct_descriptive_interval":[float(x) for x in np.quantile(boot,[.025,.975])] if len(boot)>=300 else None,
       "cluster_independence_is_imperfect":True}
    outlier_removed=[]
    for c in ("DE","AT"):
        rr=[x for x in rows if x["countryCode"]==c]
        if rr:
            rm=max(rr,key=lambda x:(int(x[primary]),float(x["lat_floor"]),float(x["lon_floor"])))
            outlier_removed.append(rm)
    pruned=[x for x in rows if x not in outlier_removed]
    output["remove_highest_source_count_tile_each_country_fixed_sensitivity"]={
       "n_removed":len(outlier_removed),"n_remaining":len(pruned),
       "source_count_model":fit(pruned,primary)}
    srcP=calc_source(rows,source_fields[1]);srcN=calc_source(rows,source_fields[2])
    corr,_=spearmanr(srcP,srcN)
    output["separate_P_N_2022_annual_WWTP_reported_mass_log1p_rank_correlation"]=float(corr)
    output["do_not_infer_independent_P_N_nutrient_mechanisms"]=True
    return output

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--2021-plant-source-gate",type=Path,required=True,dest="gate")
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    prior=json.loads(a.gate.read_text())
    if prior["schema"]!="gesui_DE_AT_original_EEA2022_1degree_tile_2021_Plant_only_effort_freeze_v01":
        raise ValueError("Wrong independently committed prior-year non-butterfly source gate")
    if not prior["no_2022_butterfly_photo_records_or_2022_butterfly_taxon_used"]:
        raise ValueError("2021 prior-year/source-only gate contaminated by 2022 butterfly responses")
    out={
      "schema":"gesui_2022_DE_AT_EEA_FW_WWTP_macro_Lycaenidae_vs_two_other_families_photo_reporting_v01",
      "source_only_2021_Plantae_gate_sha256":sha(a.gate),
      "2021_Plant_effort_gate_success":prior["gate_2021_effort_eligibility_before_2022_butterfly_outcome"],
      "original_EEA2022_source_manifest_sha256":prior["source_only_original_2022_EEA_manifest_sha256"],
      "original_source_positive_1degree_cells":prior["original_preselected_source_positive_DE_AT_grid_n"],
      "year":YEAR,"butterfly_photo_window":[SEASON_START,SEASON_END],
      "butterfly_taxon_IDs":FAMILIES,
      "min_2022_targetgroup_photos_per_eligible_cell":MIN_TOTAL_PHOTOS_PER_TILE,
      "min_2022_photo_eligible_total":MIN_SUPPORTED_TILES,
      "min_2022_photo_eligible_per_country":MIN_SUPPORTED_PER_COUNTRY,
      "min_country_geographic_2degree_blocks":MIN_SUPPORTED_BLOCKS,
      "no_butterfly_photo_pixels_or_species_userIDs_read":True,
      "no_actual_butterfly_abundance_or_puddling_measured":True,
      "year_specific_EEA_source_not_direct_instream_substrate_nutrient_dose":True,
      "no_wastewater_butterfly_causal_effect_identified":True,
      "utc":datetime.now(timezone.utc).isoformat()
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    if not out["2021_Plant_effort_gate_success"]:
        out.update({"stage":"STOP_BEFORE_2022_BUTTERFLY_API","butterfly_2022_API_queries":0,
          "conclusion":"Independent preoutcome 2021 plant observational-effort gate did not pass. No 2022 butterfly results read."})
        a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
        print("GESUI_MACRO_2021_PLANT_EFFORT_STOP_NO_BUTTERFLY_OBSERVATIONS",flush=True)
        return
    # API source uses fixed original country-place IDs verified by 2021 non-butterfly effort phase.
    countries=prior["official_iNaturalist_country_places"]
    ids={k:int(v["official_place"]["id"]) for k,v in countries.items()}
    rows=[x for x in prior["source_and_2021_plant_only_full_cell_records"] if x["2021_plant_effort_eligible"]]
    obs=[]
    for i,x in enumerate(rows):
        traces={f:api_total(x,ids[x["countryCode"]],id) for f,id in FAMILIES.items()}
        counts={k:v["total_results"] for k,v in traces.items()}
        allok=all(isinstance(v,int) for v in counts.values())
        rec={
           **{k:x[k] for k in (
              "countryCode","lat_floor","lon_floor","lat_ceiling","lon_ceiling",
              "geographic_2degree_block","2021_Plantae_photo_count",
              "public_EEA_both_measured_positive_unique_FW_outfalls",
              "public_EEA_total_reported_P_t_year","public_EEA_total_reported_N_t_year")},
           "original_sources_2022_qualifying_n":int(x["public_EEA_both_measured_positive_unique_FW_outfalls"]),
           "Lycaenidae_2022_photo_n":counts.get("Lycaenidae"),
           "Pieridae_2022_photo_n":counts.get("Pieridae"),
           "Papilionidae_2022_photo_n":counts.get("Papilionidae"),
           "three_family_total_2022":sum(counts.values()) if allok else None,
           "all_family_API_success":allok,
           "source_receipts_2022":traces,
           "not_individual_abundance":True}
        obs.append(rec)
        if i%7==0 or i==len(rows)-1:
            print("GESUI_DE_AT_2022_3FAMILY_PUBLIC_PHOTO_COUNTS",
                json.dumps({"queried_cells":i+1,"total_cells":len(rows),"country":x["countryCode"],
                  "original_EEA_WWTP_source_count":rec["original_sources_2022_qualifying_n"],
                  "photo_counts":counts},sort_keys=True),flush=True)
        time.sleep(.45)
    supported=[r for r in obs if r["three_family_total_2022"] is not None and
                r["three_family_total_2022"]>=MIN_TOTAL_PHOTOS_PER_TILE]
    c=Counter(r["countryCode"] for r in supported)
    blocks=set(r["countryCode"]+"|"+r["geographic_2degree_block"] for r in supported)
    gate=len(supported)>=MIN_SUPPORTED_TILES and c.get("DE",0)>=MIN_SUPPORTED_PER_COUNTRY and c.get("AT",0)>=MIN_SUPPORTED_PER_COUNTRY and len(blocks)>=MIN_SUPPORTED_BLOCKS
    out.update({
       "source_only_2021_plant_effort_eligible_cells":len(rows),
       "2022_family_API_queries":len(rows)*3,
       "2022_cells_with_all_three_API_valid":sum(r["all_family_API_success"] for r in obs),
       "2022_cells_all_three_API_missing":sum(not r["all_family_API_success"] for r in obs),
       "2022_three_family_photos_ge20_support_cells":len(supported),
       "2022_three_family_photos_ge20_support_by_country":dict(c),
       "2022_supported_distinct_country_2degree_blocks":len(blocks),
       "2022_photo_share_model_gate":gate,
       "all_fixed_cells_2022_family_photo_source_receipts":obs,
    })
    if gate:
        fields=("public_EEA_both_measured_positive_unique_FW_outfalls",
                "public_EEA_total_reported_P_t_year",
                "public_EEA_total_reported_N_t_year")
        out["exploratory_wastewater_source_vs_public_photo_reporting_descriptive"]=result_summary(supported,fields)
        out["interpretation"]="Exploratory country- and botanical-effort-adjusted photo reporting association across original source-positive 1degree cells ONLY. Not sewage cause or butterfly population response."
    else:
        out["stage"]="STOP_BEFORE_2022_PHOTO_SHARE_REGRESSION"
        out["interpretation"]="2022 family photo data support insufficient for predeclared country and geographic block model. No wastewater effect inferred."
    a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_DE_AT_2022_LYCAENIDAE_MACRO_REPORTING_SOURCE_GATED_RESULT",
      json.dumps({"source_cells":len(obs),"supported":len(supported),
        "supported_country":dict(c),"2deg_blocks":len(blocks),"model_gate":gate,
        "desc":out.get("exploratory_wastewater_source_vs_public_photo_reporting_descriptive",{}).get("three_predeclared_source_predictors")},sort_keys=True),flush=True)

if __name__=="__main__":main()
