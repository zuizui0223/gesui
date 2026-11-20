#!/usr/bin/env python3
"""Source-gated ChemPop v2 orthophosphate+nitrate+ammoniacalN adjustment.

Main question: Does geographic modelled sewage-effluent percentage EDF add
prediction to original observed-family BMWP beyond fixed geographical/
agricultural/riverbed monitoring controls and same-biological-date matched
river WATER-QUALITY WINDOW SUMMARY nutrient concentrations?

This is an exploratory conditional association, NOT source-resolved N/P
input or true within-site sewage-dose manipulation. All four models use the
same source-eligible site-year samples and region-disjoint test partitions.
"""
from __future__ import annotations

import argparse, collections, hashlib, json, math, re, sys, tempfile, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.linear_model import Ridge

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis"/"20261119_europe_benthic_taxa"))
import chempop_v2_bmwp_filter as original
import chempop_v2_source_census as source

NUTRIENTS=("orthophosphate","nitrate","ammoniacalnitrogen")
LOQ_PRIMARY="WQ_MEAN_LOQ2"
SOURCE_GATES={
    "min_matching_original_samples":500,
    "min_true_original_sites":150,
    "min_site_year":200,
    "min_native_regions":4,
    "min_sites_per_region":10,
    "max_full_design_VIF_EDF_or_one_nutrient":15.0,
}
EXPECTED_SCORE_SHA="4f9be755ce6d24e645b46705dbec9c0367bc71a9ce6b335d3bf2657d8c10ac5a"
EXPECTED_ORIGINAL_N=7836
EXPECTED_ORIGINAL_SITES=707
N_BOOT=2000
SEED=20261120


def find_source(files,region,chemical):
    matches=[p for p in files if re.fullmatch(
        r"data/cp_"+re.escape(region)+r"_minv_"+re.escape(chemical)+r"stats_[0-9]{8}-[0-9]{8}[.]csv",
        p,re.I)]
    if len(matches)!=1:
        raise RuntimeError("Original nutrient source filename ambiguous/absent: "+region+" "+chemical)
    return matches[0]


def chemical_summaries(z,files,original_keys):
    chemicals={}
    summary={}
    for region in source.PREFIXES:
        for chemical in NUTRIENTS:
            path=find_source(files,region,chemical)
            stat=collections.Counter()
            recognized={}
            for row in source.csv_rows(z,path,["BIO_SITE_ID","BIO_SAMPLE_DATE",
                                              "WQ_COUNT","WQ_MEAN_LOQ1",
                                              "WQ_MEAN_LOQ2","WQ_MEAN_LOQ3",
                                              "WQ_COUNT_L_LOQ","WQ_COUNT_G_LOQ"]):
                stat["raw_original_window_rows"]+=1
                d=original.parsed_date(row["BIO_SAMPLE_DATE"])
                if d is None:
                    stat["unparseable_date"]+=1
                    continue
                key=(region,str(row["BIO_SITE_ID"]).strip(),d.date().isoformat())
                if key not in original_keys:continue
                stat["source_key_joined_to_eligible_biology"]+=1
                wq_count=source.finite_numeric(row["WQ_COUNT"])
                mean2=source.finite_numeric(row[LOQ_PRIMARY])
                if wq_count is None or wq_count<=0 or mean2 is None or mean2<0:
                    stat["missing_invalid_count_or_primary_mean"]+=1
                    continue
                mean1=source.finite_numeric(row["WQ_MEAN_LOQ1"])
                mean3=source.finite_numeric(row["WQ_MEAN_LOQ3"])
                below=source.finite_numeric(row["WQ_COUNT_L_LOQ"])
                above=source.finite_numeric(row["WQ_COUNT_G_LOQ"])
                value={"primary":float(mean2),
                       "lower_LOQ_version":float(mean1) if mean1 is not None else None,
                       "upper_LOQ_version":float(mean3) if mean3 is not None else None,
                       "n_original_chem_measurements":float(wq_count),
                       "n_below_LOQ":float(below) if below is not None else None,
                       "n_above_LOQ":float(above) if above is not None else None}
                if key in recognized:
                    if recognized[key]!=value:
                        raise RuntimeError("Conflicting original chemistry summary for identical BIO site/date")
                    stat["identical_duplicate_site_date_rows"]+=1
                    continue
                recognized[key]=value
                stat["valid_unique_bio_site_date"]+=1
                if below is not None and below>0:
                    stat["chem_summary_contains_below_LOQ"]+=1
                if mean1 is not None and mean3 is not None:
                    stat["all_three_LOQ_mean_variants_exist"]+=1
            summary[region+"|"+chemical]={
                "original_file":path,
                "source_member_SHA256":hashlib.sha256(z.read(path)).hexdigest(),
                "counters":dict(stat),
            }
            for key,v in recognized.items():
                chemicals[(key,chemical)]=v
    return chemicals,summary


def source_join_raw(raw,weights,chemicals):
    join=[]
    for row,w in zip(raw,weights):
        key=(row["region"],row["site"],row["sample_date"])
        c={name:chemicals.get((key,name)) for name in NUTRIENTS}
        if not all(v is not None for v in c.values()):continue
        join.append((row,w,c))
    candidate=[q[0] for q in join]
    common_sites={(r["region"],r["site"]) for r in candidate}
    common_sy={(r["site_region_id"],r["year"]) for r in candidate}
    by_region=collections.defaultdict(set)
    for r in candidate:
        by_region[r["region"]].add(r["site"])
    n_region=sum(len(a)>=SOURCE_GATES["min_sites_per_region"] for a in by_region.values())
    pass_gate=(len(join)>=SOURCE_GATES["min_matching_original_samples"]
        and len(common_sites)>=SOURCE_GATES["min_true_original_sites"]
        and len(common_sy)>=SOURCE_GATES["min_site_year"]
        and n_region>=SOURCE_GATES["min_native_regions"])
    original_summary={
        "matched_original_samples":len(join),
        "matched_distinct_sites":len(common_sites),
        "matched_distinct_site_year":len(common_sy),
        "regions_with_at_least10_sites":n_region,
        "region_site_count":{r:len(s) for r,s in sorted(by_region.items())},
        "source_gate_PASS":bool(pass_gate),
        "all_three_nutrients_same_BIO_SITE_ID_and_BIO_SAMPLE_DATE":True,
        "no_missing_chemical_values_imputed":True,
    }
    return join,original_summary


def site_year_matched(join,scores):
    rows=[x[0] for x in join]
    weights=[x[1] for x in join]
    df,_,_,_,_,_=original.site_year_response(rows,weights,scores)
    df["edf_log1p"]=np.log1p(df.edf_pct.to_numpy(float))
    # Water quality is a source-provided time-window summary attached to the
    # biological sample-date. Count one unique site×date even if multiple
    # original biological analysis IDs share the same source water summary.
    dates={}
    for original_row,_,q in join:
        key=(original_row["site_region_id"],original_row["year"],original_row["sample_date"])
        if key in dates and dates[key]!=q:
            raise RuntimeError("Same original site/date has inconsistent water chemistry")
        dates[key]=q
    compiled=collections.defaultdict(list)
    for (sid,year,date),vals in dates.items():
        compiled[(sid,year)].append(vals)
    chem_rows=[]
    for (sid,year),items in sorted(compiled.items()):
        record={"site_region_id":sid,"year":year,"n_original_chem_sample_dates":len(items)}
        for k in NUTRIENTS:
            value=[v[k]["primary"] for v in items]
            if not value or any(not np.isfinite(x) or x<0 for x in value):
                raise RuntimeError("Invalid original chemistry summary in paired site/year")
            record["chem_"+k]=float(np.log1p(np.mean(value)))
            record["raw_mean_"+k]=float(np.mean(value))
        chem_rows.append(record)
    q=pd.DataFrame(chem_rows)
    d=df.merge(q,on=["site_region_id","year"],how="inner",validate="one_to_one")
    if len(d)!=len(df) or len(df)!=len(q):
        raise RuntimeError("Source matched chemical and biological site-year grain mismatches")
    if d.duplicated(["site_region_id","year"]).any():
        raise RuntimeError("Duplicate matched biological site-year")
    return d, {
        "paired_original_site_year":len(d),
        "paired_original_site":int(d.site_region_id.nunique()),
        "paired_original_region":int(d.region.nunique()),
        "n_nutrient_matched_distinct_dates":len(dates),
        "source_chemistry_timewindow_not_instantaneous":True,
        "source_chemistry_siteyear_mean_of_distinct_original_sample_dates":True,
        "n_sample_dates_per_siteyear_min":int(d.n_original_chem_sample_dates.min()),
        "n_sample_dates_per_siteyear_max":int(d.n_original_chem_sample_dates.max()),
    }


def full_design(df,with_chem=True):
    X=original.model_design(df,True)
    for k in NUTRIENTS if with_chem else ():
        name="chem_"+k
        v=df[name].to_numpy(float)
        sd=float(v.std(ddof=0))
        if sd<=1e-10:raise RuntimeError("Matched analyte variance too small: "+name)
        X[name]=(v-v.mean())/sd
    if np.linalg.matrix_rank(X.to_numpy(float)) != X.shape[1]:
        raise RuntimeError("Matched original background+EDF+3 nutrients is rank deficient")
    return X


def VIF_for_original_controls(X):
    values=X.to_numpy(float)
    result={}
    targets=("EDF_log1p_std",)+tuple("chem_"+k for k in NUTRIENTS)
    for name in targets:
        idx=X.columns.get_loc(name)
        x=values[:,idx]
        C=np.delete(values,idx,axis=1)
        residual=x-C@np.linalg.lstsq(C,x,rcond=None)[0]
        num=float(np.dot(x-x.mean(),x-x.mean()))
        den=float(np.dot(residual,residual))
        if den<=1e-12:raise RuntimeError("Exact source collinearity for "+name)
        result[name]=num/den
    return result


def cluster_EDF_coeff(df,with_chem):
    X=full_design(df,with_chem)
    Y=df.y_incidence.to_numpy(float)
    w=original.site_weights(df)
    fit=sm.WLS(Y,X.to_numpy(float),weights=w).fit(
        cov_type="cluster",cov_kwds={"groups":df.site_region_id.astype(str).to_numpy()},
        use_t=False)
    focus={}
    names=("EDF_log1p_std",)+tuple("chem_"+k for k in NUTRIENTS if with_chem)
    for name in names:
        idx=X.columns.get_loc(name)
        focus[name]={"beta":float(fit.params[idx]),
                     "site_cluster_normal95":[float(x) for x in fit.conf_int()[idx]],
                     "p_nominal":float(fit.pvalues[idx])}
    return {"n_siteyear":len(df),"n_site":df.site_region_id.nunique(),
            "n_regional_blocks":df.region.nunique(),
            "coefficient":focus,"design_cols":list(X.columns),
            "not_a_causal_mediation_or_experimental_effect":True}


def region_transfer_four_models(df):
    """True region-disjoint original site held-out predictions, four nested B."""
    numeric=tuple(original.COVARIATES)
    columns={
        "B":numeric,
        "B_EDF":numeric+("edf_log1p",),
        "B_NUTRIENTS":numeric+tuple("chem_"+x for x in NUTRIENTS),
        "B_NUTRIENTS_EDF":numeric+tuple("chem_"+x for x in NUTRIENTS)+("edf_log1p",)
    }
    out={}
    for region in sorted(df.region.unique()):
        test=df[df.region==region].copy().reset_index(drop=True)
        train=df[df.region!=region].copy()
        overlap=set(train.site.astype(str)) & set(test.site.astype(str))
        train=train[~train.site.astype(str).isin(overlap)].copy().reset_index(drop=True)
        if (test.site_region_id.nunique()<10 or train.site_region_id.nunique()<200):
            out[region]={"status":"SOURCE_STOP_REGION_TRANSFER",
                         "target_n_sites":int(test.site_region_id.nunique()),
                         "train_n_sites":int(train.site_region_id.nunique())}
            continue
        preds={}
        for label,cols in columns.items():
            pre=ColumnTransformer([
                ("num",StandardScaler(),list(cols)),
                ("cat",OneHotEncoder(handle_unknown="ignore",sparse_output=False),
                 list(original.CATS))
            ],remainder="drop",sparse_threshold=0)
            A=pre.fit_transform(train)
            Z=pre.transform(test)
            weights=original.site_weights(train)
            preds[label]=Ridge(alpha=10,fit_intercept=True,solver="svd").fit(
                A,train.y_incidence.to_numpy(float),sample_weight=weights
            ).predict(Z)
        values={}
        y=test.y_incidence.to_numpy(float)
        for k,pred in preds.items():
            e=(y-pred)**2
            loss=[]
            for site in sorted(test.site_region_id.astype(str).unique()):
                idx=test.site_region_id==site
                loss.append(float(e[idx].mean()))
            values[k]=float(np.mean(loss))
        out[region]={"status":"REGION_AND_SITE_DISJOINT",
                     "train_site":train.site_region_id.nunique(),
                     "target_site":test.site_region_id.nunique(),
                     "MSE":values}
    valid=[x for x in out if out[x]["status"]=="REGION_AND_SITE_DISJOINT"]
    if len(valid)<5:
        return {"status":"SOURCE_STOP_LESS_THAN_FIVE_REGION_FOLDS","by_region":out}
    results={}
    losses={m:np.asarray([out[region]["MSE"][m] for region in valid],float)
            for m in columns}
    contrasts={
        "EDF_after_geography":("B","B_EDF"),
        "NUTRIENTS_after_geography":("B","B_NUTRIENTS"),
        "EDF_after_nutrients":("B_NUTRIENTS","B_NUTRIENTS_EDF"),
        "NUTRIENTS_after_EDF":("B_EDF","B_NUTRIENTS_EDF")
    }
    rng=np.random.default_rng(SEED)
    ix=rng.integers(0,len(valid),size=(N_BOOT,len(valid)))
    for label,(base,added) in contrasts.items():
        a=losses[base]; b=losses[added]
        bs=100*(a[ix].mean(axis=1)-b[ix].mean(axis=1))/a[ix].mean(axis=1)
        results[label]={
            "MSE_additional_gain_percent":float(100*(a.mean()-b.mean())/a.mean()),
            "fixed_region_loss_2000_bootstrap95":[float(x) for x in np.quantile(bs,[.025,.975])],
            "improved_regions":int((a>b).sum()),
            "n_valid_regions":len(valid),
            "description":"fixed FOLDS: no model refitting uncertainty"
        }
    return {"status":"REGION_DISJOINT_FOUR_NESTED_MODELS",
            "valid_regions":valid,
            "contrasts":results,
            "by_region":out,
            "not_independent_nations_or_causal_mediation":True}


def make(a):
    out={"schema":"GESUI_CHEMPOP_V2_EDF_AFTER_THREE_NP_WQ_SUMMARIES_V1",
         "date":"2026-11-20",
         "status":"SOURCE_STOP",
         "posthoc_after_EDF_BMWP_negative":True,
         "source_coverage_precondition_always_before_outcome":True,
         "named_original_nutrients":list(NUTRIENTS),
         "original_WQ_mean_definition":"WQ_MEAN_LOQ2",
         "original_EDF_static_modelled_wastewater_percentage":True,
         "not_TP_not_TN_not_WWTP_delivered_mass":True,
         "not_direct_same_day_chemistry_or_causal_path":True,
         "method":"region/site-disjoint source-only gates before fitted WLS and Ridge"}
    try:
        with tempfile.TemporaryDirectory(prefix="gesui_chem3_source_") as tmp:
            d=Path(tmp)
            file=source.retrieve(source.EIDC_URL,d/"original.zip",source.SOURCE_DOWNLOAD_LIMIT)
            scores,scoreinfo=source.score_reference(d)
            if scoreinfo["SHA256"] != EXPECTED_SCORE_SHA or len(scores)!=92:
                raise RuntimeError("Original BMWP fixed score source changed")
            out["source"]={"EIDC_ZIP_SHA256":file["SHA256"],
                           "BMWP_RData_SHA256":scoreinfo["SHA256"]}
            with zipfile.ZipFile(d/"original.zip") as z:
                files=z.namelist()
                if len(files)!=326:raise RuntimeError("Original 326 EIDC file source changed")
                raw,weights,orders,old=original.construct_eligible(z,files,scores)
                if (len(raw)!=EXPECTED_ORIGINAL_N
                    or old["source_gate"]["eligible_site_n"]!=EXPECTED_ORIGINAL_SITES):
                    raise RuntimeError("Previous source-scored 7836/707 frame drift")
                keys={(r["region"],r["site"],r["sample_date"]) for r in raw}
                chemicals,chem_src=chemical_summaries(z,files,keys)
            out["original_wq_source_files"]=chem_src
            joint,gate=source_join_raw(raw,weights,chemicals)
            out["matched_source_gate"]=gate
            if not gate["source_gate_PASS"]:
                out["status"]="SOURCE_STOP_JOINT_NP_COHORT_TOO_SMALL"
            else:
                df,joined=site_year_matched(joint,scores)
                out["matched_siteyear_source"]=joined
                if (joined["paired_original_site"]<SOURCE_GATES["min_true_original_sites"]
                    or joined["paired_original_site_year"]<SOURCE_GATES["min_site_year"]):
                    out["status"]="SOURCE_STOP_COLLAPSED_SITEYEAR_COHORT"
                else:
                    X=full_design(df,True)
                    VIF=VIF_for_original_controls(X)
                    out["joint_covariate_source_vif"]=VIF
                    if max(VIF.values())>SOURCE_GATES["max_full_design_VIF_EDF_or_one_nutrient"]:
                        out["status"]="SOURCE_MODEL_STOP_COLLINEAR_EDF_OR_NP"
                    else:
                        out["source_model_cohort"]={
                            "n_original_siteyear":len(df),
                            "n_sites":df.site_region_id.nunique(),
                            "n_regions":df.region.nunique(),
                            "n_original_matched_bio_samples":gate["matched_original_samples"],
                        }
                        out["EDF_after_all_three_chemistry_WLS"]=cluster_EDF_coeff(df,True)
                        out["EDF_without_chemistry_WLS"]=cluster_EDF_coeff(df,False)
                        out["regional_transfer"]=region_transfer_four_models(df)
                        out["status"]="SOURCE_GATE_PASS_POSTHOC_CHEM_ADJUSTED_COMPLETE"
    except Exception as ex:
        out.update({"status":"SOURCE_OR_METHOD_STOP",
                    "error":type(ex).__name__+": "+str(ex)[:750]})
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_EDF_AFTER_ORIGINAL_WQ_NP="+json.dumps({
        "status":out["status"],
        "error":out.get("error"),
        "matched_source_gate":out.get("matched_source_gate"),
        "matched_siteyear_source":out.get("matched_siteyear_source"),
        "source_VIF":out.get("joint_covariate_source_vif"),
        "cohort":out.get("source_model_cohort"),
        "EDF_nutri_WLS":out.get("EDF_after_all_three_chemistry_WLS",{}).get("coefficient"),
        "EDF_alone_WLS":out.get("EDF_without_chemistry_WLS",{}).get("coefficient"),
        "regional_transfer":out.get("regional_transfer",{}).get("contrasts"),
    },sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    make(p.parse_args())
