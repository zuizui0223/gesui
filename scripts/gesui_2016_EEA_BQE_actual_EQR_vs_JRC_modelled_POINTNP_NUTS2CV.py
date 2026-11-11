#!/usr/bin/env python3
"""Actual 2016 EEA benthic/phytobenthos normalized EQR × modelled 2014–18 JRC N/P.

Important independent source: original EEA 2016 Biology monitoring site IDs and
GISCO NUTS 2016 geography were frozen BEFORE reading 2016 numeric normalized EQR.
JRC PsN/PsP are MODELLED 2014–18 regional domestic+INDUSTRIAL point source
tonnes per NUTS2 km², not measured UWWTD 2016 site-specific nutrients.
NUTS2-disjoint fivefold CV, country×method and lat/lon/region area controls.
No causal wastewater effect, no pooled 2018/22 directly measured WWTP claims.
"""
from __future__ import annotations
import argparse,hashlib,json,math,zipfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
import statsmodels.formula.api as smf

SRC=Path("results/eea2016_thirdyear_source_20261111")
SOURCE=SRC/"2016_EEA_BENTHIC_GISCO_NUTS2_GREEN_POINT_NP_SOURCE_GATE_BEFORE_EQR.json"
FOLDS=SRC/"2016_EEA_BENTHIC_JRC_NUTS2_FROZEN_5FOLDS_NO_EQR.json"
META=SRC/"2016_BQE_PREOUTCOME_SOURCE_ONLY/2016_FIXED_BQE_MONITORING_SITE_METADATA_NO_EQR.csv"
BQE_SOURCE=SRC/"2016_BQE_PREOUTCOME_SOURCE_ONLY/2016_BQE_PREOUTCOME_SOURCE_QC.json"
ZIP_MEMBER="eea_t_waterbase-biology-2026_p_1990-2024_v01_r00/Waterbase_v2025_1_T_WISE2_BiologyEQRData.csv"
EXPECTED="f11001a57172b618bee95f58fd0ec51c8e71be7b9c01a70927dde21b7ca6e2a9"
BIO={"EEA_13-01-4":"benthic_invertebrates","EEA_124-04-9":"phytobenthos","EEA_14-02-8":"fish"}
YEAR="2016"
ALPHA=10.
SEED=20161111
VALID_MIN={"benthic_invertebrates":{"n":800,"countries":8,"NUTS2":60,"blocks":15},
           "phytobenthos":{"n":600,"countries":8,"NUTS2":40,"blocks":12}}
SOURCE_COLS=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
             "phenomenonTimeReferenceYear","parameterWaterBodyCategory",
             "observedPropertyDeterminandBiologyEQRCode","procedureClassificationSystem",
             "resultNormalisedEQRValue","resultObservationStatus"]

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(2*1024*1024),b""):h.update(b)
    return h.hexdigest()
def sha_member(z,p):
    h=hashlib.sha256()
    with z.open(p) as f:
        for b in iter(lambda:f.read(2*1024*1024),b""):h.update(b)
    return h.hexdigest()
def read_original(root,eeazip):
    src=json.loads((root/SOURCE).read_text())
    folds=json.loads((root/FOLDS).read_text())
    bqe=json.loads((root/BQE_SOURCE).read_text())
    if folds["original_GISCO2016_JRC2014to18_official_source_full_site_manifest_SHA256"]!=sha(root/SOURCE):
        raise ValueError("Original source-only 2016 study cohort changed after fold freezing")
    if folds["original_2016_EEA_numeric_EQR_outcome_VALUES_not_read"] is not True or bqe["pre_2016_Biology_EQR_values_read"] is not False:
        raise ValueError("Original source gate failed EQR pre-outcome separation")
    if not src["precommitted_2016_benthic_NP_modelled_NUTS2_SUPPORT_GATE"]["can_open_2016_EEA_numeric_normalised_EQR_and_test_POSTHOC_correlation"]:
        raise ValueError("Original source-first country/region/join support gate failed")
    with zipfile.ZipFile(eeazip) as z:
        if sha_member(z,ZIP_MEMBER)!=EXPECTED:
            raise ValueError("Original official EEA 2016 Biology EQR member source SHA has changed")
        frame=pd.read_csv(z.open(ZIP_MEMBER),usecols=SOURCE_COLS,dtype="string",low_memory=False)
    frame=frame.loc[
        frame.phenomenonTimeReferenceYear.eq(YEAR)&
        frame.parameterWaterBodyCategory.eq("RW")&
        frame.observedPropertyDeterminandBiologyEQRCode.isin(BIO)
    ].copy()
    keys=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme"]
    for k in keys:
        frame[k]=frame[k].astype("string").str.strip()
    frame["sitekey"]=frame[keys].fillna("").agg("|".join,axis=1)
    frame["siteSHA"]=frame.sitekey.map(lambda x:hashlib.sha256(x.encode()).hexdigest())
    frame["bio"]=frame.observedPropertyDeterminandBiologyEQRCode.map(BIO)
    if frame.duplicated(["siteSHA","bio"]).any():
        raise ValueError("Original official 2016 BQE has duplicate site and element records")
    # Source geographic metadata only. EQR numeric measured in original EEA CSV.
    coords=pd.read_csv(root/META,usecols=[
        "site_key","observedPropertyDeterminandBiologyEQRCode",
        "lat","lon","region_3deg_block","countryCode","procedureClassificationSystem"],
        dtype={"site_key":"string"},low_memory=False)
    coords["siteSHA"]=coords.site_key.map(lambda x:hashlib.sha256(x.encode()).hexdigest())
    coords["bio"]=coords.observedPropertyDeterminandBiologyEQRCode.map(BIO)
    if coords.duplicated(["siteSHA","bio"]).any():
        raise ValueError("Frozen 2016 EEA source site/element metadata contains duplicates")
    origins=src["original_2016_source_site_and_modelled_NUTS2_PN_values_HASHED_IDs_NO_EQR"]
    positions={}
    for entry in origins:
        k=(entry["original_2016_BQE_site_key_SHA256"],entry["original_2016_BQE_group"])
        if k in positions:raise ValueError("Original 2016 geo/JRC source site duplicated")
        positions[k]=entry
    if len(positions)!=len(origins) or len(positions)<4700:
        raise ValueError("Original 2016 source map incomplete")
    numeric=pd.to_numeric(frame.resultNormalisedEQRValue,errors="coerce")
    frame["EQRest"]=numeric
    frame["EQRest_valid"]=numeric.between(0,1,inclusive="both")&np.isfinite(numeric)
    frame["source_2016_country_method"]=frame.countryCode.astype(str)+"|"+frame.procedureClassificationSystem.fillna("UNKNOWN").astype(str)
    # Double check SOURCE-ONLY key/code/country/method before biological values.
    matched=frame.loc[[
        (a,b) in positions for a,b in zip(frame.siteSHA,frame.bio)
    ]].copy()
    matched=matched.merge(coords[["siteSHA","bio","lat","lon","region_3deg_block",
                                  "countryCode","procedureClassificationSystem"]],
                          on=["siteSHA","bio"],how="inner",
                          validate="one_to_one",suffixes=("","_FROZEN"))
    if (matched.countryCode!=matched.countryCode_FROZEN).any() or (
        matched.procedureClassificationSystem.fillna("UNKNOWN")!=matched.procedureClassificationSystem_FROZEN.fillna("UNKNOWN")).any():
        raise ValueError("Original EEA 2016 station classification/source original method changed")
    samples=[]
    for m in matched.itertuples(index=False):
        x=positions[(m.siteSHA,m.bio)]
        originalfold=folds["source_biological_classes"][m.bio]["fixed_region_to_fivefold_source_only"][x["original_GISCO2016_NUTS2"]]
        if str(m.countryCode)!=x["original_2016_EEA_country"]:
            raise ValueError("Original JRC NUTS2 source country differs from live EQR record")
        samples.append({
           "site_sha256":m.siteSHA,
           "bio":m.bio,
           "country":str(m.countryCode),
           "country_method":str(m.source_2016_country_method),
           "nuts2":x["original_GISCO2016_NUTS2"],
           "geo3deg":str(m.region_3deg_block),
           "original_fixed_nuts2_fold":int(originalfold),
           "lat":float(m.lat),"lon":float(m.lon),
           "logArea":math.log1p(float(x["original_GISCO2016_NUTS2_km2"])),
           "logPointN":float(x["log1p_modelled_POINT_N_tonnes_per_km2"]),
           "logPointP":float(x["log1p_modelled_POINT_P_tonnes_per_km2"]),
           "EQRest":float(m.EQRest) if pd.notna(m.EQRest) and math.isfinite(float(m.EQRest)) else None,
           "EQRest_valid":bool(m.EQRest_valid),
           "status":str(m.resultObservationStatus) if pd.notna(m.resultObservationStatus) else "MISSING"
        })
    df=pd.DataFrame(samples)
    if df.empty:raise ValueError("No original 2016 EQR source station matched frozen JRC sites")
    if df.duplicated(["site_sha256","bio"]).any():raise ValueError("Original source BQE site count duplicated after numeric join")
    for bio in BIO.values():
        goal=folds["source_biological_classes"][bio]["n_sites_source_fixed_before_EQR"]
        if len(df[df.bio==bio])!=goal:
            raise ValueError(f"Source-only 2016 BQE site population changed at numeric outcome stage {bio}")
    return df,src,folds,bqe

def source_qc(df,bio,source):
    q=df.loc[df.bio==bio]
    valid=q.loc[q.EQRest_valid].copy()
    invalid=q.loc[~q.EQRest_valid]
    stats={
      "original_source_sites_pre_EQR":len(q),
      "original_EEA2016_numeric_EQR_valid_sites":len(valid),
      "missing_or_invalid_numeric_EQR":len(invalid),
      "number_valid_original_EQR_countries":valid.country.nunique(),
      "number_valid_original_EQR_NUTS2_regions":valid.nuts2.nunique(),
      "number_valid_original_EQR_3deg_blocks":valid.geo3deg.nunique(),
      "numeric_EQR_min":float(valid.EQRest.min()) if len(valid) else None,
      "numeric_EQR_median":float(valid.EQRest.median()) if len(valid) else None,
      "numeric_EQR_max":float(valid.EQRest.max()) if len(valid) else None,
      "valid_site_original_country_counts":dict(Counter(valid.country)),
      "missing_or_invalid_by_country":dict(Counter(invalid.country)),
      "original_status_of_missing_or_invalid":dict(Counter(invalid.status)),
      "sites_above_normalized_EQR_1_not_silently_clamped":sum((q.EQRest.fillna(0)>1)),
      "sites_negative_normalized_EQR_not_silently_zeroed":sum((q.EQRest.fillna(0)<0)),
      "numeric_EEA_BQE_not_species_richness_or_population_abundance":True
    }
    return valid,stats

def build_predictor(n):
    cont=["lat","lon","logArea"]
    if n!="BASE":cont.append("logPoint"+n)
    transformer=ColumnTransformer([
      ("numeric",StandardScaler(),cont),
      ("country_method",OneHotEncoder(handle_unknown="ignore",sparse_output=False),["country_method"])
    ],remainder="drop",sparse_threshold=0)
    return Pipeline([("transform",transformer),("model",Ridge(alpha=ALPHA))])

def means_by_group(g,cols,group):
    return {m:float(np.mean([np.mean(z[m]) for _,z in g.groupby(group)])) for m in cols}

def compute_prediction(g):
    g=g.sort_values(["original_fixed_nuts2_fold","nuts2","site_sha256"]).reset_index(drop=True).copy()
    models={"baseline":"BASE","upstream_modelPointN":"N","upstream_modelPointP":"P"}
    pred={k:np.full(len(g),np.nan) for k in models}
    for f in range(1,6):
        train=g.original_fixed_nuts2_fold.ne(f)
        test=g.original_fixed_nuts2_fold.eq(f)
        if test.sum()==0 or train.sum()<100:raise ValueError("Frozen 2016 NUTS2 spatial CV fold empty")
        intersect=set(g.loc[train,"nuts2"])&set(g.loc[test,"nuts2"])
        if intersect:raise ValueError("One original NUTS2 predictor region leaked across source-frozen CV folds")
        for k,flag in models.items():
            model=build_predictor(flag)
            model.fit(g.loc[train],g.loc[train,"EQRest"])
            pred[k][test.to_numpy()]=model.predict(g.loc[test])
    for k in models:
        if not np.isfinite(pred[k]).all():raise ValueError("Out-of-region EEA BQE EQR model predictions incomplete")
        g["pred_"+k]=pred[k]
        g["sqerr_"+k]=(g.EQRest-pred[k])**2
        g["abserr_"+k]=np.abs(g.EQRest-pred[k])
    loss={key:float(g["sqerr_"+key].mean()) for key in models}
    abs_loss={key:float(g["abserr_"+key].mean()) for key in models}
    rel={k:100*(loss["baseline"]-loss[k])/loss["baseline"] for k in ("upstream_modelPointN","upstream_modelPointP")}
    # Honest explanatory OOS comparisons at independent region and country grain.
    grp_nuts=g.groupby("nuts2")[["sqerr_"+key for key in models]].mean()
    grp_country=g.groupby("country")[["sqerr_"+key for key in models]].mean()
    equal_nuts={k:float(grp_nuts["sqerr_"+k].mean()) for k in models}
    equal_country={k:float(grp_country["sqerr_"+k].mean()) for k in models}
    def gain(o,k):return 100*(o["baseline"]-o[k])/o["baseline"]
    rng=np.random.default_rng(SEED+len(g))
    # Independent NUTS2-region resampling is *descriptive*, not a prospectively
    # registered causal confidence interval; geographic autocorrelation remains.
    boots={}
    for label,groupframe in (("NUTS2",grp_nuts),("country",grp_country)):
        n=len(groupframe)
        raw=groupframe.to_numpy(float)
        sims={k:[] for k in ("upstream_modelPointN","upstream_modelPointP")}
        for _ in range(999):
            chosen=raw[rng.integers(0,n,n)]
            val=chosen.mean(axis=0)
            if val[0]<=0:continue
            for k,ix in (("upstream_modelPointN",1),("upstream_modelPointP",2)):
                sims[k].append(float((val[0]-val[ix])/val[0]*100))
        boots[label]={k:{"descriptive_percentile_95":list(map(float,np.quantile(a,[.025,.975]))),
                        "n_nonparametric_resamples":len(a)}
                       for k,a in sims.items()}
    by_region=[]
    for nuts,d in g.groupby("nuts2"):
        by_region.append({"nuts2":nuts,"n_original_2016_BQE_sites":len(d),
            "country":d.country.iloc[0],
            "source_frozen_test_fold":int(d.original_fixed_nuts2_fold.iloc[0]),
            "mean_EQR_original_numeric":float(d.EQRest.mean()),
            "mean_logPointN_modelled2014to18_density":float(d.logPointN.iloc[0]),
            "mean_logPointP_modelled2014to18_density":float(d.logPointP.iloc[0]),
            "CV_region_MSE":{k:float(d["sqerr_"+k].mean()) for k in models},
            "BQE_quality_mean_is_aggregate_not_individual_EEA_score":True
        })
    return {
      "year_of_original_EEA_BQE_observation":2016,
      "model_source_year":"2014–2018 averaged JRC GREEN point-source modeled",
      "model_grid":"NUTS2 2016 regional allocation",
      "forecast_validation":"5-fold source-frozen geographically disjoint NUTS2 regions, no source same-region train-test leakage",
      "ridge_original_alpha":ALPHA,
      "baseline_predictors":"training-standardized lat+lon+log(NUTS2 area) and country|procedureClassificationSystem",
      "added_modelled_source_predictors":"log1p(source tonnes per km2) of N or P, independently",
      "source_test_EEA_BQE_site_count":len(g),
      "source_test_independent_NUTS2_regions":g.nuts2.nunique(),
      "source_test_countries":g.country.nunique(),
      "source_test_original_3deg_blocks":g.geo3deg.nunique(),
      "MSE_site_equal":loss,
      "MAE_site_equal":abs_loss,
      "MSE_gain_pct_site_equal":rel,
      "MSE_NUTS2_equal":equal_nuts,
      "MSE_country_equal":equal_country,
      "MSE_gain_pct_NUTS2_equal":{k:gain(equal_nuts,k) for k in rel},
      "MSE_gain_pct_country_equal":{k:gain(equal_country,k) for k in rel},
      "descriptive_bootstrap_95gain_pct":boots,
      "frozen_fold_case_counts":{str(i):int(g.original_fixed_nuts2_fold.eq(i).sum()) for i in range(1,6)},
      "source_JRC_green_region_sampled_EEA_EQR_summary":by_region,
      "same_source_magnitude_across_sites_in_NUTS2_NOT_site_independent_treatments":True
    }

def adjusted_OLS(g):
    g=g.copy()
    g["country_method"]=g.country_method.astype(str)
    forms={"pointN":"EQRest~logPointN+logArea+lat+lon+C(country_method)",
           "pointP":"EQRest~logPointP+logArea+lat+lon+C(country_method)",
           "both_N_and_P":"EQRest~logPointN+logPointP+logArea+lat+lon+C(country_method)"}
    out={}
    for label,formula in forms.items():
        model=smf.ols(formula,data=g).fit(cov_type="cluster",cov_kwds={"groups":g.nuts2},
                      use_t=False)
        coef={}
        for nutrient in ("logPointN","logPointP"):
            if nutrient in model.params:
                coef[nutrient]={
                    "beta_per_source_log1p_tonnes_km2":float(model.params[nutrient]),
                    "NUTS2_cluster_95":list(map(float,model.conf_int().loc[nutrient])),
                    "exploratory_nominal_p":float(model.pvalues[nutrient]),
                }
        out[label]={"formula":formula,"source_site_n":len(g),"distinct_source_region_n":g.nuts2.nunique(),
             "R2":float(model.rsquared),
             "original_BQE_country_method_adjusted_source_log_magnitude_coefficient":coef,
             "not_randomized_or_confirmatory":True}
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--EEA-biology-zip",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    frame,source,folds,bqe=read_original(root,a.EEA_biology_zip)
    reports={}
    for bio in ("benthic_invertebrates","phytobenthos","fish"):
        group,qc=source_qc(frame,bio,source)
        d={"precommitted_source_country_region_outcome_QC":qc}
        gate=VALID_MIN.get(bio)
        permitted=bool(gate and len(group)>=gate["n"] and group.country.nunique()>=gate["countries"] and
            group.nuts2.nunique()>=gate["NUTS2"] and group.geo3deg.nunique()>=gate["blocks"])
        d["precommitted_numeric_EQR_support_gate_PASS"]=permitted
        if permitted:
            d["raw_group_EEA_BQE_EQR_vs_JRC_model_source_spearman"]={
                "rho_NP_2014to18_modelledN":float(spearmanr(group.EQRest,group.logPointN).statistic),
                "rho_NP_2014to18_modelledP":float(spearmanr(group.EQRest,group.logPointP).statistic),
                "not_independent_NUTS2_regions":True
            }
            d["observational_country_method_adjusted_OLS_NUTS2cluster"]=adjusted_OLS(group)
            d["original_NUTS2_disjoint_holdout_Ridge10"]=compute_prediction(group)
        else:
            d["STOP_reason"]="Source-only EQR numeric outcome or country/region gate not reached; no ecological model"
        reports[bio]=d
        print("EEA2016_GREEN_2014_18_BQE_RESULT",bio,json.dumps({
             "qc":qc,"gate":permitted,
             "raw":d.get("raw_group_EEA_BQE_EQR_vs_JRC_model_source_spearman"),
             "CV":{k:v for k,v in d.get("original_NUTS2_disjoint_holdout_Ridge10",{}).items()
                   if k in ("MSE_gain_pct_site_equal","MSE_gain_pct_NUTS2_equal",
                            "MSE_gain_pct_country_equal","source_test_EEA_BQE_site_count")},
             "reg":d.get("observational_country_method_adjusted_OLS_NUTS2cluster")
        },ensure_ascii=False),flush=True)
    result={
       "schema":"gesui_EEA2016_official_benthic_EQR_JRC_GREEN_modelled_2014to18_POINT_NP_NUTS2016_source_locked_results_v01",
       "UTC":datetime.now(timezone.utc).isoformat(),
       "official_EEA_Biology2026_original_full_zip_SHA256":sha(a.EEA_biology_zip),
       "official_EEA2016_original_member_SHA256":EXPECTED,
       "original_2016_sites_GISCO2016_JRC_modelled_SOURCE_before_EQR_SHA256":sha(root/SOURCE),
       "original_original5fold_SOURCE_before_EQR_SHA256":sha(root/FOLDS),
       "no_hyperparameter_alpha_search_Ridge10":True,
       "no_2016_measured_UWWTD_plant_WWTP_NP_load_claimed":True,
       "source_JRC_green2014to18_domestic_AND_industrial_modelled_region_point_load_tonnes":True,
       "no_EEA2016_real_plant_species_richness_assessed":True,
       "actual_2016_EEA_BQE_EQR_is_normalized_biological_condition_not_species_abundance":True,
       "new_results_are_POSTHOC_after_2018_2022_trueWWTP_positive_predictive_gains":True,
       "no_SEWAGE_CAUSAL_ECOLOGY_CONCLUSION":True,
       "results":reports,
       "full_individual_2016_EEA_BQE_EQR_observation_values_NOT_REPUBLISHED":True
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
if __name__=="__main__":main()
