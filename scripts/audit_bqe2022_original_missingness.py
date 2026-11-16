#!/usr/bin/env python3
"""2022 EEA Biology EQR source attrition audit, with explicit nullable-value handling.

Diagnoses why original official EEA 2022 river biological site records did not
enter the 2022 WWTP load → EQR prediction models. This DOES NOT refit, tune,
impute, reclassify, or change existing positive predictive results.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import zipfile
import numpy as np
import pandas as pd

ORIGINAL_EQR_MEMBER=(
    "eea_t_waterbase-biology-2026_p_1990-2024_v01_r00/"
    "Waterbase_v2025_1_T_WISE2_BiologyEQRData.csv"
)
EXPECTED_EQR_SHA256="f11001a57172b618bee95f58fd0ec51c8e71be7b9c01a70927dde21b7ca6e2a9"
EXPECTED_SOURCE_CSV_SHA256="e12df511d64c4c580ef2422dccdbbd11a1309c84d29d467fd4527206e41fa636"
EXPECTED_N={"EEA_13-01-4":(3282,2475),"EEA_124-04-9":(3663,2784),
            "EEA_14-02-8":(1300,1243)}
BQE={"EEA_13-01-4":"benthic_invertebrates",
     "EEA_124-04-9":"phytobenthos","EEA_14-02-8":"fish"}
KEY_FIELDS=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme"]
FIELDS=KEY_FIELDS+[
    "parameterWaterBodyCategory",
    "observedPropertyDeterminandBiologyEQRCode",
    "phenomenonTimeReferenceYear",
    "resultEQRValue","resultNormalisedEQRValue",
    "resultObservationStatus",
    "procedureClassificationSystem",
    "resultEcologicalStatusClassValue",
]
def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for part in iter(lambda:f.read(1024*1024),b""):
            h.update(part)
    return h.hexdigest()

def original_2022(path):
    with zipfile.ZipFile(path) as z:
        digest=hashlib.sha256()
        with z.open(ORIGINAL_EQR_MEMBER) as f:
            for data in iter(lambda:f.read(1024*1024),b""):
                digest.update(data)
        if digest.hexdigest()!=EXPECTED_EQR_SHA256:
            raise ValueError("Official uncompressed EEA 2026 EQR data changed")
        data=pd.read_csv(z.open(ORIGINAL_EQR_MEMBER),usecols=FIELDS,
                         dtype="string",low_memory=False)
    data=data[
        data.phenomenonTimeReferenceYear.eq("2022") &
        data.parameterWaterBodyCategory.eq("RW") &
        data.observedPropertyDeterminandBiologyEQRCode.isin(BQE)
    ].copy()
    for f in KEY_FIELDS:data[f]=data[f].str.strip()
    data["site_key"]=data[KEY_FIELDS].fillna("").agg("|".join,axis=1)
    if data.duplicated(["site_key","observedPropertyDeterminandBiologyEQRCode"]).any():
        raise ValueError("Official EQR source has duplicate year/site/biological element keys")
    return data

def normalized_eqr_valid(s):
    """True/False ndarray; a nullable numeric NA is ALWAYS False, not NA."""
    numeric=pd.to_numeric(s,errors="coerce")
    a=numeric.to_numpy(dtype=float,na_value=np.nan)
    return np.isfinite(a)&(a>=0),a

def _rate(x,denom):
    return float(x/denom) if denom else None

def summarize(source,raw,group,code,expected_observed):
    # Source-only registry includes exact geography/exposure but deliberately
    # excluded any biological score; all original values are read here.
    df=source[source.observedPropertyDeterminandBiologyEQRCode.eq(code)].copy()
    right=raw[raw.observedPropertyDeterminandBiologyEQRCode.eq(code)].copy()
    if df.duplicated("site_key").any() or right.duplicated("site_key").any():
        raise ValueError("Source site keys cannot be mapped unambiguously")
    merged=df.merge(right[[
        "site_key","resultEQRValue","resultNormalisedEQRValue",
        "resultObservationStatus","procedureClassificationSystem",
        "resultEcologicalStatusClassValue"
    ]],on="site_key",how="left",indicator=True,
           suffixes=("_source","_original"),validate="one_to_one")
    if len(merged)!=EXPECTED_N[code][0]:
        raise ValueError("Frozen source-only cohort was altered")
    found=merged._merge.eq("both").to_numpy()
    if not found.all():
        raise ValueError("Source-only 2022 EQR record lacks an original source row")
    valid,values=normalized_eqr_valid(merged.resultNormalisedEQRValue)
    raw_valid,raw_values=normalized_eqr_valid(merged.resultEQRValue)
    retained=valid & found
    if int(retained.sum())!=expected_observed:
        raise ValueError("Original model sample is not recovered by nonmissing EQR rule")
    merged["has_original_EQR"]=found
    merged["normalized_EQR_valid"]=retained
    merged["original_EQR_numeric"]=raw_valid
    merged["original_EQR_greater_than_one"]=raw_valid&(raw_values>1)
    merged["source_exposed"]=pd.to_numeric(
        merged.upstream_2022_unique_EEA_both_load_WWTP_count,errors="coerce").fillna(0).gt(0)
    missed=merged.loc[~merged.normalized_EQR_valid].copy()
    missing_source_value=int(np.isnan(values).sum())
    nonfinite_not_na=int((~np.isfinite(values)&~np.isnan(values)).sum())
    negative_numeric=int((np.isfinite(values)&(values<0)).sum())
    bycountry=[]
    for country,g in merged.groupby("countryCode",dropna=False):
        positive=g.source_exposed.to_numpy(dtype=bool)
        good=g.normalized_EQR_valid.to_numpy(dtype=bool)
        dropped=~good
        bycountry.append({
            "country":str(country),
            "source_2022_BQE_candidate_sites":len(g),
            "source_upstream_WWTP_positive_candidate_sites":int(positive.sum()),
            "numeric_2022_EQR_model_sites":int(good.sum()),
            "source_upstream_WWTP_positive_model_sites":int((positive&good).sum()),
            "not_modelled_due_to_original_normalized_EQR":int(dropped.sum()),
            "not_modelled_upstream_WWTP_positive":int((positive&dropped).sum()),
            "missing_fraction_upstream_positive":_rate((positive&dropped).sum(),positive.sum()),
            "missing_fraction_source_zero":_rate((~positive&dropped).sum(),(~positive).sum()),
        })
    bycountry.sort(key=lambda x:(-x["not_modelled_due_to_original_normalized_EQR"],x["country"]))
    srcpos=merged.source_exposed.to_numpy(dtype=bool)
    origstat=missed.resultObservationStatus_original.fillna("MISSING")
    upstream_pos=int(srcpos.sum())
    upstream_pos_loss=int((srcpos&~retained).sum())
    result={
        "element":group,"code":code,
        "source_site_n_before_numeric_EQR_validation":len(merged),
        "original_2022_EQR_rows_matched_by_site_key":int(found.sum()),
        "numeric_EQR_usable_sites":int(retained.sum()),
        "source_sites_excluded_due_to_non_numeric_or_missing_EQR":len(missed),
        "missing_original_normalized_EQR_values":missing_source_value,
        "nonfinite_nonmissing_values":nonfinite_not_na,
        "negative_original_numeric_EQR":negative_numeric,
        "source_total_upstream_reported_load_positive":upstream_pos,
        "model_upstream_reported_load_positive":int((srcpos&retained).sum()),
        "missing_upstream_reported_load_positive":upstream_pos_loss,
        "fraction_of_upstream_source_sites_removed":
             _rate(upstream_pos_loss,upstream_pos),
        "fraction_of_source_zero_sites_removed":
             _rate((~srcpos&~retained).sum(),(~srcpos).sum()),
        "raw_EQR_numeric_among_excluded":int(missed.original_EQR_numeric.sum()),
        "raw_EQR_greater_than_one_among_excluded":
             int(missed.original_EQR_greater_than_one.sum()),
        "original_status_codes_among_excluded":{
            str(k):int(v) for k,v in origstat.value_counts().items()},
        "original_classification_method_by_missingness":{
            "retained_count":{str(k):int(v) for k,v in
                merged.loc[merged.normalized_EQR_valid,"procedureClassificationSystem_original"].
                fillna("MISSING").value_counts().items()},
            "excluded_count":{str(k):int(v) for k,v in
                missed.procedureClassificationSystem_original.
                fillna("MISSING").value_counts().items()}
        },
        "by_country":bycountry,
        "not_an_effect_of_sewage_or_model_fit":True,
    }
    return result,merged

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--official-zip",type=Path,required=True)
    p.add_argument("--frozen-registry",type=Path,required=True)
    p.add_argument("--original-outcome-json",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if sha(a.frozen_registry)!=EXPECTED_SOURCE_CSV_SHA256:
        raise ValueError("Source-only EEA2022 physical registry SHA mismatch")
    source=pd.read_csv(a.frozen_registry,dtype="string",low_memory=False)
    if len(source)!=8245:
        raise ValueError("Source site registry original denominator changed")
    if source.duplicated(["site_key","observedPropertyDeterminandBiologyEQRCode"]).any():
        raise ValueError("Source registry contains duplicated biology site")
    prior=json.loads(a.original_outcome_json.read_text())
    original=original_2022(a.official_zip)
    result={}
    for code,group in BQE.items():
        _,n_model=EXPECTED_N[code]
        if prior["results"][group]["support"]["numeric_original_EQR_sites"]!=n_model:
            raise ValueError("Original biological fitted sample differs")
        result[group],joined=summarize(source,original,group,code,n_model)
        s=result[group]
        print("GESUI_BQE_2022_ORIGINAL_MISSINGNESS_AUDIT",json.dumps({
            "group":group,
            "source_candidates":s["source_site_n_before_numeric_EQR_validation"],
            "measured":s["numeric_EQR_usable_sites"],
            "missing":s["source_sites_excluded_due_to_non_numeric_or_missing_EQR"],
            "original_statuses_missing":s["original_status_codes_among_excluded"],
            "exposed_source_total":s["source_total_upstream_reported_load_positive"],
            "exposed_removed":s["missing_upstream_reported_load_positive"],
            "raw_EQR_present_among_normalized_missing":s["raw_EQR_numeric_among_excluded"],
            "country_missing":[(x["country"],x["not_modelled_due_to_original_normalized_EQR"])
                               for x in s["by_country"] if x["not_modelled_due_to_original_normalized_EQR"]]
        },sort_keys=True),flush=True)
    report={
       "schema":"gesui_EEA2026_Biology_EQR_original_2022_completeness_and_source_exposure_audit_v01",
       "study_year":2022,"audit_utc":datetime.now(timezone.utc).isoformat(),
       "original_EEA_Biology_uncompressed_EQR_SHA256":EXPECTED_EQR_SHA256,
       "original_source_only_2022_registry_SHA256":EXPECTED_SOURCE_CSV_SHA256,
       "original_EEA_Biology_model_run":42432384196,
       "source_region_folds_and_station_selection_unchanged":True,
       "pandas_nullable_numeric_NA_counted_explicitly":True,
       "old_biological_prediction_coefficients_untouched":True,
       "no_propagation_or_reconstruction_of_missing_normalized_EQR":True,
       "source_reported_status_definitions_need_official_dictionary_interpretation":True,
       "no_causal_WWTP_effect_inference":True,
       "findings":result
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_EEA_BQE_SOURCE_ATTRITION_FINAL",
          json.dumps({k:{"missing":v["source_sites_excluded_due_to_non_numeric_or_missing_EQR"],
                         "exposed_missing":v["missing_upstream_reported_load_positive"],
                         "fraction_positive_removed":v["fraction_of_upstream_source_sites_removed"]}
                      for k,v in result.items()},sort_keys=True),flush=True)
if __name__=="__main__":main()
