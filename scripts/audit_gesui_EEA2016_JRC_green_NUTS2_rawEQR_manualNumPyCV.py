#!/usr/bin/env python3
"""Independent RAW ORIGINAL EEA 2016 EQR × JRC NUTS2 GREEN model-source audit.

Separate EEA original ZIP parser, original source hash matching, manual
NumPy ridge regression normal equations with train-fold-only normalization,
manual 2016 region-disjoint predictions and model comparison. Not calling
the primary sklearn/statsmodels implementation. Reproduces 2016 benthic
and phytobenthos site/regional MSE and ecological raw Spearman estimates.
NO 2016 measured WWTP emission or causal wastewater effect.
"""
from __future__ import annotations
import argparse,hashlib,json,zipfile,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr

ROOT=Path("results/eea2016_thirdyear_source_20261111")
A=ROOT/"2016_EEA_BENTHIC_GISCO_NUTS2_GREEN_POINT_NP_SOURCE_GATE_BEFORE_EQR.json"
F=ROOT/"2016_EEA_BENTHIC_JRC_NUTS2_FROZEN_5FOLDS_NO_EQR.json"
SITES=ROOT/"2016_BQE_PREOUTCOME_SOURCE_ONLY/2016_FIXED_BQE_MONITORING_SITE_METADATA_NO_EQR.csv"
RESULT=ROOT/"2016_EEA_RIVER_BENTHIC_EQR_JRC_GREEN_MODEL_POINTNP_NUTS2_HELDOUT_RESULTS.json"
EQRFILE="eea_t_waterbase-biology-2026_p_1990-2024_v01_r00/Waterbase_v2025_1_T_WISE2_BiologyEQRData.csv"
EXPECTED="f11001a57172b618bee95f58fd0ec51c8e71be7b9c01a70927dde21b7ca6e2a9"
Y="2016"
CODES={"EEA_13-01-4":"benthic_invertebrates","EEA_124-04-9":"phytobenthos","EEA_14-02-8":"fish"}
ALPHA=10.
NEED=("countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
      "parameterWaterBodyCategory","observedPropertyDeterminandBiologyEQRCode",
      "phenomenonTimeReferenceYear","procedureClassificationSystem","resultNormalisedEQRValue")

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(2*1024*1024),b""):h.update(b)
    return h.hexdigest()
def member_sha(z,name):
    h=hashlib.sha256()
    with z.open(name) as v:
        for chunk in iter(lambda:v.read(2*1024*1024),b""):h.update(chunk)
    return h.hexdigest()
def sample(root,archive):
    orig=json.loads((root/A).read_text())
    folds=json.loads((root/F).read_text())
    first=json.loads((root/RESULT).read_text())
    if first["original_2016_sites_GISCO2016_JRC_modelled_SOURCE_before_EQR_SHA256"]!=sha(root/A) or (
        first["original_original5fold_SOURCE_before_EQR_SHA256"]!=sha(root/F)):
        raise ValueError("2016 original SOURCE selected after ecological EQR: inconsistent hashes")
    if folds["original_GISCO2016_JRC2014to18_official_source_full_site_manifest_SHA256"]!=sha(root/A) or not (
        folds["original_2016_EEA_numeric_EQR_outcome_VALUES_not_read"]):
        raise ValueError("Original prior source-only sample/folds not immutable before biology")
    positions={}
    for x in orig["original_2016_source_site_and_modelled_NUTS2_PN_values_HASHED_IDs_NO_EQR"]:
        k=(x["original_2016_BQE_site_key_SHA256"],x["original_2016_BQE_group"])
        if k in positions:raise ValueError("Original GISCO NUTS2016 site×element duplicated")
        positions[k]=x
    # Strictly extract raw original source normalized EQRs, not any newer/synthetic site file.
    with zipfile.ZipFile(archive) as z:
        if member_sha(z,EQRFILE)!=EXPECTED:
            raise ValueError("Raw EEA EQR scientific member digest mismatch")
        raw=pd.read_csv(z.open(EQRFILE),usecols=list(NEED),dtype="string",low_memory=False)
    raw=raw.loc[
        raw.phenomenonTimeReferenceYear.eq(Y)&
        raw.parameterWaterBodyCategory.eq("RW")&
        raw.observedPropertyDeterminandBiologyEQRCode.isin(CODES)
    ].copy()
    raw["key"]=raw[["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme"]].fillna("").astype(str).agg("|".join,axis=1)
    raw["source_id"]=raw.key.map(lambda v:hashlib.sha256(v.encode()).hexdigest())
    raw["bio"]=raw.observedPropertyDeterminandBiologyEQRCode.map(CODES)
    raw["sample_id"]=list(zip(raw.source_id,raw.bio))
    if raw.sample_id.duplicated().any():
        raise ValueError("Raw original 2016 BQE has duplicate source site×BQE")
    raw=raw.loc[raw.sample_id.isin(positions)].copy()
    # Frozen year2016 geographic coords inspected WITHOUT outcome at earlier SOURCE stage.
    sites=pd.read_csv(root/SITES,usecols=["site_key","lat","lon","observedPropertyDeterminandBiologyEQRCode","procedureClassificationSystem"],dtype="string")
    sites["source_id"]=sites.site_key.map(lambda v:hashlib.sha256(str(v).encode()).hexdigest())
    sites["bio"]=sites.observedPropertyDeterminandBiologyEQRCode.map(CODES)
    points={(x.source_id,x.bio):(float(x.lat),float(x.lon),str(x.procedureClassificationSystem))
            for x in sites.itertuples(index=False)}
    xrows=[]
    for x in raw.itertuples(index=False):
        key=(x.source_id,x.bio)
        meta=positions[key]
        if key not in points:raise ValueError("Original source station ID missing geographic data")
        lat,lon,original_method=points[key]
        method=str(x.procedureClassificationSystem) if pd.notna(x.procedureClassificationSystem) else "UNKNOWN"
        if method!=original_method and not (
            method=="UNKNOWN" and original_method in ("<NA>","nan","None")
        ):
            raise ValueError("Original EEA2026 BQE assessment procedure changed after 2016 source")
        if str(x.countryCode)!=meta["original_2016_EEA_country"]:
            raise ValueError("Original EEA2016 country does not match JRC/GISCO model source")
        val=pd.to_numeric(x.resultNormalisedEQRValue,errors="coerce")
        score=float(val) if pd.notna(val) else math.nan
        xrows.append({
          "key":x.source_id,"bio":x.bio,"country":str(x.countryCode),
          "country_method":str(x.countryCode)+"|"+method,
          "nuts2":meta["original_GISCO2016_NUTS2"],
          "fold":int(folds["source_biological_classes"][x.bio]["fixed_region_to_fivefold_source_only"][meta["original_GISCO2016_NUTS2"]]),
          "lat":lat,"lon":lon,
          "logArea":math.log1p(float(meta["original_GISCO2016_NUTS2_km2"])),
          "pointN":float(meta["log1p_modelled_POINT_N_tonnes_per_km2"]),
          "pointP":float(meta["log1p_modelled_POINT_P_tonnes_per_km2"]),
          "EQRest":score
        })
    df=pd.DataFrame(xrows)
    if df.duplicated(["key","bio"]).any():raise ValueError("Original EQR source SHA key duplicate")
    return df,first,positions

def design(train,test,predictor):
    keys=["lat","lon","logArea"]
    if predictor!="baseline":keys.append("pointN" if predictor=="N" else "pointP")
    A=train[keys].to_numpy(float)
    B=test[keys].to_numpy(float)
    mu=A.mean(axis=0);sigma=A.std(axis=0,ddof=0)
    sigma=np.where(sigma>0,sigma,1.0)
    A=(A-mu)/sigma
    B=(B-mu)/sigma
    cats=sorted(train.country_method.unique())
    onehot_train=np.zeros((len(train),len(cats)))
    onehot_test=np.zeros((len(test),len(cats)))
    category={cat:i for i,cat in enumerate(cats)}
    for row,cat in enumerate(train.country_method):onehot_train[row,category[cat]]=1
    for row,cat in enumerate(test.country_method):
        if cat in category:onehot_test[row,category[cat]]=1
    return np.column_stack((A,onehot_train)),np.column_stack((B,onehot_test))
def manual_ridge_predict(train,test,predictor):
    a,b=design(train,test,predictor)
    y=train.EQRest.to_numpy(float)
    mx=a.mean(axis=0)
    my=y.mean()
    cent=a-mx
    coefficients=np.linalg.solve(cent.T@cent+ALPHA*np.eye(cent.shape[1]),cent.T@(y-my))
    return my+(b-mx)@coefficients

def check_models(df,original,bio):
    d=df[df.bio==bio].copy()
    expected=original["results"][bio]
    nonnumeric=d.EQRest.isna()|~np.isfinite(d.EQRest)|~d.EQRest.between(0,1)
    if int(nonnumeric.sum())!=expected["precommitted_source_country_region_outcome_QC"]["missing_or_invalid_numeric_EQR"]:
        raise ValueError("Original 2016 numeric EQR omission counts do not reproduce: "+bio)
    d=d.loc[~nonnumeric].copy().reset_index(drop=True)
    cv=expected["original_NUTS2_disjoint_holdout_Ridge10"]
    if len(d)!=cv["source_test_EEA_BQE_site_count"] or d.nuts2.nunique()!=cv["source_test_independent_NUTS2_regions"]:
        raise ValueError("Original 2016 EQR source outcome case size or region count not reproduced")
    y=d.EQRest.to_numpy(float)
    preds={k:np.full(len(d),np.nan) for k in ("baseline","N","P")}
    for fold in range(1,6):
        train=d.loc[d.fold!=fold]
        test=d.loc[d.fold==fold]
        if set(train.nuts2)&set(test.nuts2):
            raise ValueError("Original REGION2026 JRC posthoc spatial holdout leaking same NUTS2")
        for key in preds:
            preds[key][(d.fold==fold).to_numpy()]=manual_ridge_predict(train,test,key)
    errors={k:(y-v)**2 for k,v in preds.items()}
    if not all(np.isfinite(a).all() for a in errors.values()):
        raise ValueError("Numerical independent site-level EQR ridge audit contains missing predictions")
    result={}
    aliases={"baseline":"baseline","N":"upstream_modelPointN","P":"upstream_modelPointP"}
    for k,v in errors.items():
        orig=cv["MSE_site_equal"][aliases[k]]
        mse=float(v.mean())
        if abs(mse-orig)>1e-8:
            raise ValueError(f"Independent numpy Ridge10 original 2016 EQR MSE mismatched {bio}/{k}: {mse} expected {orig}")
        region=pd.DataFrame({"nuts2":d.nuts2.values,"sq":v}).groupby("nuts2").sq.mean()
        orig_region=cv["MSE_NUTS2_equal"][aliases[k]]
        if abs(float(region.mean())-orig_region)>1e-8:
            raise ValueError(f"Independent EEA2016 NUTS2-equal MSE mismatch {bio}/{k}")
        result[k]={"manual_numpy_MSE_site_equal":mse,"manual_numpy_MSE_equal_NUTS2":float(region.mean())}
    nspearman=float(spearmanr(y,d.pointN).statistic)
    pspearman=float(spearmanr(y,d.pointP).statistic)
    earlier=expected["raw_group_EEA_BQE_EQR_vs_JRC_model_source_spearman"]
    if abs(nspearman-earlier["rho_NP_2014to18_modelledN"])>1e-12 or abs(pspearman-earlier["rho_NP_2014to18_modelledP"])>1e-12:
        raise ValueError("Independent original EEA2016 EQR and JRC 2014-18 source Spearman differs")
    return {
        "independent_original_EEA2016_numeric_EQR_n":len(d),
        "independent_2016_numeric_EQR_country_count":d.country.nunique(),
        "independent_original_EEA2016_biology_NUTS2_count":d.nuts2.nunique(),
        "strict_NUTS2_disjoint_SOURCE_ONLY_5fold_predicted":True,
        "manual_separate_2016_year_NUTS2_source_predictor_Ridge10_reproductions":result,
        "independent_spearman_original_2016_EQR_vs_modelled_JRC_N":nspearman,
        "independent_spearman_original_2016_EQR_vs_modelled_JRC_P":pspearman,
        "same_GISCO2016_site_metadata_and_original_NUTS2_source_IDs":True,
        "no_2016_site_specific_WWTP_measured_LOAD_claimed":True
    }
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--EEA-zip",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    df,first,source=sample(root,a.EEA_zip)
    results={}
    for k in ("benthic_invertebrates","phytobenthos"):
        results[k]=check_models(df,first,k)
        print("INDEPENDENT_EEA2016_JRC_MODEL_SOURCE_QA",k,json.dumps(results[k]),flush=True)
    receipt={
       "schema":"gesui_original_EEA2016_benthic_phytobenthos_EQR_manualNumPy_2014_18_GREEN_pointNP_79NUTS2_independent_CV_v01",
       "official_full_EEA_Biology2026_ZIP_SHA256":sha(a.EEA_zip),
       "official_uncompressed_EEA2016_EQR_member_SHA256":EXPECTED,
       "source_geo_2016_GISCO_JRC2014to18_SHA256":sha(root/A),
       "source_frozen2016_NUTS2_DISJOINT_CV_SHA256":sha(root/F),
       "source_2016_original_BQE_site_geography_SHA256":sha(root/SITES),
       "source_first_2016_EQR_observational_results_SHA256":sha(root/RESULT),
       "independent_sklearn_and_statsmodels_free_NumPy_Ridge10_original_source_reproduction":results,
       "not_2016_measured_WWTP_specific_loads":True,
       "not_causal_or_confirmatory_claim":True
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
if __name__=="__main__":main()
