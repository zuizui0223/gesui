#!/usr/bin/env python3
"""Explicitly POSTHOC: does original 2022 DE/AT WWTP reported source P still predict
Lycaenidae *photo reporting composition* in heldout 2° blocks after independent
ESA WorldCover 2021 built-up/cropland/country land area and Plantae photo effort?

2021 WorldCover categorical COG source all 53 regions must be Git-frozen
before this script is run. Exact original 2022 response and 2021 plant effort
must match archived SHA. No new observations, photographs or species data queried.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

PREV_PHOTO_SHA=None
PREV_EEA_SHA="926ab310904be048d5d36bbfebd5eabfe0b318279800142f2e2a86bee924f4ce"
SEED=20261108
SHUFFLE_N=199
BOOT_N=999
MODELS={
 "outfall_count":"public_EEA_both_measured_positive_unique_FW_outfalls",
 "reported_P_t_per_year":"public_EEA_total_reported_P_t_year",
 "reported_N_t_per_year":"public_EEA_total_reported_N_t_year",
}
MIN_SOURCE_VALID=36
MIN_DE=20
MIN_AT=10
MIN_SPATIAL=12

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def rkey(c):
    return (c["countryCode"],int(float(c["lat_floor"])),int(float(c["lon_floor"])))

def set_block(x):return x["countryCode"]+"|"+x["geographic_2degree_block"]

def objective(rows):
    return np.asarray([math.log((float(r["Lycaenidae_2022_photo_n"])+0.5)/
       (float(r["Pieridae_2022_photo_n"])+float(r["Papilionidae_2022_photo_n"])+0.5)) for r in rows])

def feature_tr_te(train,test,field,include_land,lowres,override=None):
    def vals(sample,name,op=lambda a:float(a)):
        return np.asarray([op(r[name]) for r in sample],dtype=float)
    cc_t=np.asarray([1.0 if r["countryCode"]=="AT" else 0. for r in train])
    cc_v=np.asarray([1.0 if r["countryCode"]=="AT" else 0. for r in test])
    base_t=[np.ones(len(train)),cc_t]
    base_v=[np.ones(len(test)),cc_v]
    feats=[
       ("2021_Plantae_photo_count",lambda v:math.log1p(float(v))),
       ("lat_floor",lambda v:float(v)+.5),
       ("lon_floor",lambda v:float(v)+.5),
    ]
    suff= "country_terrestrial_land_area_km2_approx"
    if include_land:
        if lowres:
            feats.extend([
              ("wc_built_up_frac_256",float),("wc_crop_frac_256",float),
              ("wc_terrestrial_km2_256",lambda v:math.log1p(float(v))),
            ])
        else:
            feats.extend([
              ("wc_built_up_frac_512",float),("wc_crop_frac_512",float),
              ("wc_terrestrial_km2_512",lambda v:math.log1p(float(v))),
            ])
    for fieldname,conv in feats:
        a=vals(train,fieldname,conv);b=vals(test,fieldname,conv)
        m=float(np.mean(a));s=float(np.std(a))
        if s<=1e-13:raise ValueError("Zero training predictor variance "+fieldname)
        base_t.append((a-m)/s);base_v.append((b-m)/s)
    X0=np.column_stack(base_t);T0=np.column_stack(base_v)
    orig=MODELS[field]
    def source(sample):
        if override is None:return vals(sample,orig)
        return np.asarray([override[rkey(r)] for r in sample],float)
    tr=np.log1p(source(train));te=np.log1p(source(test))
    zt=np.zeros(len(train),float);zv=np.zeros(len(test),float)
    for country in ("DE","AT"):
        mask=np.asarray([r["countryCode"]==country for r in train])
        mtest=np.asarray([r["countryCode"]==country for r in test])
        if mask.sum()<3:raise ValueError("Insufficient trained regional source per country")
        x=tr[mask];m=float(np.mean(x));s=float(np.std(x))
        if s<=1e-13:raise ValueError("No distinct source intensity within country "+country)
        zt[mask]=(x-m)/s
        zv[mtest]=(te[mtest]-m)/s
    X1=np.column_stack([X0,zt]);T1=np.column_stack([T0,zv])
    return X0,T0,X1,T1

def OLS_predict(X,T,y):
    if np.linalg.matrix_rank(X)<X.shape[1]:
        raise ValueError("Rank deficient 2deg spatial block-heldout training")
    return T @ np.linalg.lstsq(X,y,rcond=None)[0]

def cv(rows,field,include_land=True,lowres=False,override=None,full=False):
    original_groups=sorted(set(set_block(r) for r in rows))
    preds=[]
    for group in original_groups:
        tr=[r for r in rows if set_block(r)!=group]
        te=[r for r in rows if set_block(r)==group]
        X0,T0,X1,T1=feature_tr_te(tr,te,field,include_land,lowres,override)
        yy=objective(tr);yt=objective(te)
        baseline=OLS_predict(X0,T0,yy)
        plus=OLS_predict(X1,T1,yy)
        for r,y,p0,p1 in zip(te,yt,baseline,plus):
            preds.append({
              "key":rkey(r),"country":r["countryCode"],"spatial_2degree_group":group,
              "observed_y":float(y),"baseline_prediction":float(p0),
              "WWTP_source_added_prediction":float(p1),
              "baseline_squared_error":float((y-p0)**2),
              "source_added_squared_error":float((y-p1)**2),
              "baseline_absolute_error":float(abs(y-p0)),
              "source_added_absolute_error":float(abs(y-p1)),
            })
    base=float(np.mean([p["baseline_squared_error"] for p in preds]))
    added=float(np.mean([p["source_added_squared_error"] for p in preds]))
    result={
       "n_original_1degree_country_cells":len(rows),
       "n_source_2degree_geographic_test_blocks":len(original_groups),
       "baseline_heldout_MSE":base,
       "source_added_heldout_MSE":added,
       "relative_MSE_improvement_percent":100*(base-added)/base if base>0 else None,
       "baseline_heldout_MAE":float(np.mean([p["baseline_absolute_error"] for p in preds])),
       "source_added_heldout_MAE":float(np.mean([p["source_added_absolute_error"] for p in preds])),
       "source_field":MODELS[field],
       "with_WorldCover2021_built_crop_country_land_area":include_land,
       "WorldCover_resampled_nearest_pixel_side":256 if lowres else 512 if include_land else None,
       "country_and_2deg_block_holdout_no_region_leakage":True,
    }
    if full:
        result["heldout_predictions_every_cell"]=preds
        result["country_specific_heldout_MSE"]={}
        for c in ("DE","AT"):
            pp=[p for p in preds if p["country"]==c]
            a=float(np.mean([p["baseline_squared_error"] for p in pp]))
            b=float(np.mean([p["source_added_squared_error"] for p in pp]))
            result["country_specific_heldout_MSE"][c]={
               "n":len(pp),"baseline_MSE":a,"added_source_MSE":b,
               "relative_MSE_improvement_percent":100*(a-b)/a if a>0 else None}
    return result

def make_shuffled_source(rows,field,rng,within_block):
    original=MODELS[field]
    by=defaultdict(list)
    for r in rows:
        group=(r["countryCode"],r["geographic_2degree_block"] if within_block else "ALL")
        by[group].append(r)
    shuffled={}
    for group,sub in by.items():
        v=[float(r[original]) for r in sub]
        for r,x in zip(sub,rng.permutation(v)):shuffled[rkey(r)]=float(x)
    return shuffled

def block_boot_ci(rows,full_cv):
    if not "heldout_predictions_every_cell" in full_cv:return None
    pred=full_cv["heldout_predictions_every_cell"]
    by=defaultdict(list)
    for p in pred:by[p["spatial_2degree_group"]].append(p)
    keys=list(by)
    rng=np.random.default_rng(SEED+91)
    draws=[]
    for _ in range(BOOT_N):
        newkeys=rng.choice(keys,size=len(keys),replace=True)
        selection=[p for key in newkeys for p in by[key]]
        a=float(np.mean([p["baseline_squared_error"] for p in selection]))
        b=float(np.mean([p["source_added_squared_error"] for p in selection]))
        if a>0:draws.append(100*(a-b)/a)
    return {"resamples":len(draws),"seed":SEED+91,
       "interval_95pct_descriptive":[float(q) for q in np.quantile(draws,[.025,.975])] if len(draws)>0 else None,
       "spatial_2deg_prediction_error_bootstrap_not_independent_experiment":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--land-source-only",type=Path,required=True)
    p.add_argument("--original-2022-family-photo",type=Path,required=True)
    p.add_argument("--original-2021-plant-gate",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    land=json.loads(a.land_source_only.read_text())
    photo=json.loads(a.original_2022_family_photo.read_text())
    plant=json.loads(a.original_2021_plant_gate.read_text())
    if land["schema"]!="gesui_worldcover2021_original_categorical_COG_country_masked_EEA_DE_AT_1deg_source_only_v01":
        raise ValueError("Wrong original ESA categorical WorldCover source")
    if land["NO_2022_BUTTERFLY_PHOTO_RESULTS_OPENED"] is not True or land["source_only_processed_cells"]!=53:
        raise ValueError("WorldCover source did not verify all 53 cells BLIND to butterfly")
    if land["original_EEA_source_manifest_sha256"]!="926ab310904be048d5d36bbfebd5eabfe0b318279800142f2e2a86bee924f4ce":
        raise ValueError("WorldCover target not source-matched to original EEA 1156 WWTP sources")
    if photo["source_only_2021_Plantae_gate_sha256"]!=sha(a.original_2021_plant_gate):
        raise ValueError("Original 2022 family photo data not tied to earlier source/effort gate SHA")
    if not plant["gate_2021_effort_eligibility_before_2022_butterfly_outcome"] or not photo["2022_photo_share_model_gate"]:
        raise ValueError("Original posthoc source study did not pass original gate")
    by_land={rkey(r):r for r in land["source_cells"]}
    if len(by_land)!=53:raise ValueError("Invalid original independent WorldCover geographic source keys")
    original=[r for r in photo["all_fixed_cells_2022_family_photo_source_receipts"]
              if isinstance(r.get("three_family_total_2022"),int) and
                   r["three_family_total_2022"]>=photo["min_2022_targetgroup_photos_per_eligible_cell"]]
    if len(original)!=44 or Counter(r["countryCode"] for r in original)!=Counter({"DE":29,"AT":15}):
        raise ValueError("Outcome-dependent sample not immutable original 44 cells")
    selected=[]
    rejected=[]
    for r in original:
        k=rkey(r)
        if k not in by_land:raise ValueError("Unlisted EEA 2022 WorldCover land source region "+str(k))
        wc=by_land[k]
        if not wc["original_landcountry_geometry_and_NODATA_gate"] or not wc["WorldCover_resolution_reliability_gate"]:
            rejected.append({"cell":k,"reason":"pre-declared country land fraction or original categorical source lowres instability","original_land_frac":wc["orig_WorldCover_categorical_sample_main"]["country_terrestrial_land_fraction_cell"]})
            continue
        main=wc["orig_WorldCover_categorical_sample_main"]
        sens=wc["orig_WorldCover_categorical_sample_secondary_256"]
        row=dict(r)
        for kstr,src in (("512",main),("256",sens)):
            row["wc_built_up_frac_"+kstr]=float(src["built_up_share_of_country_terrestrial"])
            row["wc_crop_frac_"+kstr]=float(src["cropland_share_of_country_terrestrial"])
            row["wc_tree_frac_"+kstr]=float(src["tree_share_of_country_terrestrial"])
            row["wc_grass_frac_"+kstr]=float(src["grassland_share_of_country_terrestrial"])
            row["wc_terrestrial_km2_"+kstr]=float(src["country_terrestrial_land_area_km2_approx"])
            row["wc_country_terrestrial_landfraction_"+kstr]=float(src["country_terrestrial_land_fraction_cell"])
        if row["wc_terrestrial_km2_512"]<=0:raise ValueError("Accepted WorldCover country-land area zero")
        selected.append(row)
    c=Counter(r["countryCode"] for r in selected)
    blocks=set(set_block(r) for r in selected)
    srcdev=[max(abs(r["wc_built_up_frac_512"]-r["wc_built_up_frac_256"]),
                abs(r["wc_crop_frac_512"]-r["wc_crop_frac_256"])) for r in selected]
    stable_n=sum(v<=.1 for v in srcdev)
    enough=(len(selected)>=MIN_SOURCE_VALID and c.get("DE",0)>=MIN_DE and
            c.get("AT",0)>=MIN_AT and len(blocks)>=MIN_SPATIAL and
            stable_n/max(1,len(selected))>=.9)
    root={
     "schema":"gesui_DE_AT_2022_original_WWTP_source_P_N_Lycaenidae_2021_WorldCover_urban_crop_posthoc_2deg_CV_v01",
     "timestamp":datetime.now(timezone.utc).isoformat(),
     "POSTHOC_after_prior_positive_2022_photo_model_and_P_holdout_results":True,
     "original_2022_family_photo_sha256":sha(a.original_2022_family_photo),
     "original_2021_nonbutterfly_effort_sha256":sha(a.original_2021_plant_gate),
     "independent_2021_ESA_WorldCover_53cell_source_SHA256":sha(a.land_source_only),
     "original_2021_NaturalEarth_country50m_geojson_sha256":land["original_natural_earth_admin0_source"]["original_geojson_SHA256"],
     "original_worldcover_categorical_class_codes":land["original_WorldCover_categorical_raster_codes"],
     "original_2022_photo_supported_countries_1deg_cells":len(original),
     "landcover_quality_eligible_original_photo_cells":len(selected),
     "landcover_quality_by_country":dict(c),
     "landcover_quality_2degree_block_count":len(blocks),
     "landcover_quality_excluded_original_2022_photo_cells":rejected,
     "original_worldcover_country_land_area_and_resolution_gate_pass":enough,
     "original_WorldCover_512_vs256_built_crop_max_diff_eligible":max(srcdev) if srcdev else None,
     "original_WorldCover_512_vs256_stable_eligible_count":stable_n,
     "original_WorldCover_full_53cell_source_coverage":land["source_land_quality_pass_count"],
     "all_original_same_sample_cells_2022_photo_counts_and_WorldCover_categorical_2021_covariates":[{
         "countryCode":r["countryCode"],"lat_floor":r["lat_floor"],"lon_floor":r["lon_floor"],
         "country_2deg_block":r["geographic_2degree_block"],
         "2021_Plantae_photo_n":r["2021_Plantae_photo_count"],
         "2022_Lycaenidae_photo_n":r["Lycaenidae_2022_photo_n"],
         "2022_total_3family_photos":r["three_family_total_2022"],
         "registered_2022_EEA_both_measured_FW_outfalls":r["public_EEA_both_measured_positive_unique_FW_outfalls"],
         "registered_P_t_year":r["public_EEA_total_reported_P_t_year"],
         "registered_N_t_year":r["public_EEA_total_reported_N_t_year"],
         "WC2021_512_built_land_fraction":r["wc_built_up_frac_512"],
         "WC2021_512_crop_land_fraction":r["wc_crop_frac_512"],
         "WC2021_512_country_terrestrial_area_km2":r["wc_terrestrial_km2_512"],
         "WC2021_256_built_land_fraction":r["wc_built_up_frac_256"],
         "WC2021_256_crop_land_fraction":r["wc_crop_frac_256"]} for r in selected],
     "remaining_caveats":["One-degree source-positive tile observations from iNaturalist are nonprobability photographed submissions, not numbers of butterflies",
        "Original 2022 WWTP annual discharge mass is not nutrient chemistry measured at terrestrial butterfly feeding substrate",
        "ESA WorldCover 2021 categorical raster class fraction is an approximation resampled from original 10m 3degree COG over country50m mask",
        "2021 previous-year plant photo reporting is imperfect proxy for independent butterfly target-group effort",
        "Source-only 1degree region geometry is not directly hydrologically connected to source WWTP","Year/individual lab sample date of EEA reported mass unverified",
        "Posthoc model and environmental covariates were selected after earlier positive association was seen"],
     "no_new_2022_butterfly_API_or_photo_IDs_queried":True,
     "no_sewage_causality_or_puddling_claim":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    if not enough:
        root.update({"status":"STOP_BEFORE_LANDCOVER_ADJUSTED_MODEL",
            "reason":"predeclared original landcover source and 2022 photograph compatible coverage gate failed"})
        a.out.write_text(json.dumps(root,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
        print("GESUI_DE_AT_ESAWC2021_LANDCOVER_GATE_STOP_NO_WWTP_EFFECT",flush=True)
        return
    # Comparison always on exact same 43 land-QC-valid 2022 photo cells; source mass
    # and photo outcomes are never used for dynamic replacement or tuning.
    results={}
    for field in MODELS:
        base=cv(selected,field,include_land=False,lowres=False,full=True)
        land_adjusted=cv(selected,field,include_land=True,lowres=False,full=True)
        land_256=cv(selected,field,include_land=True,lowres=True,full=False)
        results[field]={
         "same_43cell_original_without_landcover":base,
         "new_ESAWorldCover2021_adjusted_original512":land_adjusted,
         "fixed_ESAWC2021_resolution_256_sensitivity":land_256,
         "added_landcover_reduction_in_baseline_MSE_percent":100*(base["baseline_heldout_MSE"]-land_adjusted["baseline_heldout_MSE"])/base["baseline_heldout_MSE"],
        }
        if field=="reported_P_t_per_year":
            results[field]["descriptive_2deg_block_resampling_CI_original512"]=block_boot_ci(selected,land_adjusted)
        print("GESUI_LYCAENIDAE_ESAWC2021_2DEG_OUTOFREGION_CONTROLLED",
          json.dumps({"source":field,"n":len(selected),"source_gain_old_same43_pct":base["relative_MSE_improvement_percent"],
            "source_gain_with_2021_built_crop_land512_pct":land_adjusted["relative_MSE_improvement_percent"],
            "source_gain_with_land256_pct":land_256["relative_MSE_improvement_percent"],
            "land_baseline_mse":land_adjusted["baseline_heldout_MSE"],
            "source_mse":land_adjusted["source_added_heldout_MSE"]},sort_keys=True),flush=True)
    # Additional fixed permutation & leverage only for P, and the source-only land quality sample.
    actual=results["reported_P_t_per_year"]["new_ESAWorldCover2021_adjusted_original512"]["relative_MSE_improvement_percent"]
    prng=np.random.default_rng(SEED+100)
    p_perm={}
    for blockwise in (False,True):
        sims=[];failed=Counter()
        for _ in range(SHUFFLE_N):
            perm=make_shuffled_source(selected,"reported_P_t_per_year",prng,blockwise)
            try:
                val=cv(selected,"reported_P_t_per_year",True,False,override=perm)["relative_MSE_improvement_percent"]
                if val is not None and math.isfinite(val):sims.append(float(val))
                else:failed["nonfinite"]+=1
            except (ValueError,np.linalg.LinAlgError):failed["rank_or_training"]+=1
        k="within_original_2deg_block" if blockwise else "within_country"
        p_perm[k]={
          "n_requested":SHUFFLE_N,
          "valid":len(sims),"failed":dict(failed),
          "null_median_gain_percent":float(np.median(sims)) if sims else None,
          "source_mass_shuffled_95pct_gain_range":[float(v) for v in np.quantile(sims,[.025,.975])] if len(sims)>10 else None,
          "pseudo_tail_fraction_equal_or_better_plus1":float((1+sum(v>=actual for v in sims))/(1+len(sims))) if sims else None,
          "not_randomized_WWTP_intervention":True}
    maxP=[]
    for country in ("DE","AT"):
        rr=[r for r in selected if r["countryCode"]==country]
        maxP.append(max(rr,key=lambda r:(float(r["public_EEA_total_reported_P_t_year"]),
              float(r["lat_floor"]),float(r["lon_floor"]))))
    pruned=[r for r in selected if r not in maxP]
    root["landcover_adjusted_source_axis_results"]=results
    root["posthoc_P_mass_spatial_label_permutations"]=p_perm
    root["posthoc_P_source_mass_top_tile_country_only_removal"]={
       "original_mass_only_removed":[{"country":r["countryCode"],"lat_floor":r["lat_floor"],
               "lon_floor":r["lon_floor"],"source_P_t_per_year":r["public_EEA_total_reported_P_t_year"],
               "no_outcome_used_for_selection":True} for r in maxP],
       "same41cell_ESA2021_landcover_adjusted_MSE_comparison":cv(pruned,"reported_P_t_per_year",True,False,full=True),
    }
    root["source_landcover_spearman_between_original_2022_EEA_Pmass_and_2021_built_crop_fractions"]={}
    for country in ("DE","AT"):
        rr=[r for r in selected if r["countryCode"]==country]
        s=np.log1p([float(r["public_EEA_total_reported_P_t_year"]) for r in rr])
        rho_b,_=spearmanr(s,[r["wc_built_up_frac_512"] for r in rr])
        rho_c,_=spearmanr(s,[r["wc_crop_frac_512"] for r in rr])
        root["source_landcover_spearman_between_original_2022_EEA_Pmass_and_2021_built_crop_fractions"][country]={
            "n":len(rr),"rho_Psource_vs_ESA2021_built":float(rho_b),"rho_Psource_vs_ESA2021_crop":float(rho_c),
            "not_pollution_landcover_causality":True}
    root["status"]="POSTHOC_LANDCOVER_OUTOFREGION_PREDICTION_AND_PLACEBO_COMPLETED"
    a.out.write_text(json.dumps(root,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_DE_AT_ESAWC2021_INDEPENDENT_LANDCOVER_POSTHOC_2022_LYCAENIDAE_SOURCE_PREDICTION_COMPLETE",
      json.dumps({"n":len(selected),"bycountry":dict(c),"landcover_adjusted_P_gain_pct":actual,
        "P_country_perm_tail":p_perm["within_country"]["pseudo_tail_fraction_equal_or_better_plus1"],
        "P_2deg_perm_tail":p_perm["within_original_2deg_block"]["pseudo_tail_fraction_equal_or_better_plus1"],
        "top_P_mass_each_country_removed_P_gain_pct":root["posthoc_P_source_mass_top_tile_country_only_removal"]["same41cell_ESA2021_landcover_adjusted_MSE_comparison"]["relative_MSE_improvement_percent"]},sort_keys=True),flush=True)

if __name__=="__main__":main()
