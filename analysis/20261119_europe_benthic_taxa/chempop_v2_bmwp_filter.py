#!/usr/bin/env python3
"""ChemPop v2 modeled effluent exposure vs BMWP family composition, source-gated.

POSTHOC observational. Primary: family-incidence equal-weight published
BMWP score per original site-year. The source-only 2010–18 cohort must
pass all fixed n_sites, sample count, regional variation and data-quality
checks BEFORE fitting any biological slope or label permutation.

Original EIDC archive 326 files, 7 English regions. BMWP scores from a
separately published 2016 biotic RData family lookup (not EA2026 official).
EDF_MEAN is a modeled site-level effluent percentage proxy, NOT measured
discharged P/N mass or time-varying sewage in each collection year.
"""
from __future__ import annotations
import argparse, collections, datetime, hashlib, json, math, sys, tempfile, zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.linear_model import Ridge

import chempop_v2_source_census as source

MIN_SAMPLES=500
MIN_SITES=200
MIN_REGIONS=5
MIN_REGION_SITES=10
MIN_REGION_EDF_DISTINCT=5
MIN_SCORED_FAMILIES=3
MIN_SCORED_FRACTION=.70
MIN_COUNT_POSITIVE_SCORED_FRACTION=.80
N_PERM=999
N_BOOT=2000
SEED=20261119
Y_MIN,Y_MAX=2010,2018
NUMERIC_SOURCE=("LC_URBAN_PER","LC_ARABLE_PER","BIO_WIDTH",
                "BIO_ALTITUDE","BIO_SILT_CLAY")
COVARIATES=("lc_urban","lc_arable","width","altitude","siltclay",
            "year","month_sin","month_cos","log_n_families","coverage")
CATS=("region","method")
MAX_EDF_PERCENT=100.0

def parsed_date(raw):
    x=str(raw or "").strip()
    for fmt in ("%Y-%m-%d","%d/%m/%Y","%Y/%m/%d","%Y-%m-%d %H:%M:%S",
                "%d-%b-%Y","%d/%m/%y"):
        try:return datetime.datetime.strptime(x,fmt)
        except ValueError:pass
    try:return datetime.datetime.fromisoformat(x.replace("Z","+00:00"))
    except ValueError:return None

def construct_eligible(z,filenames,scores):
    """Source identifiers and score coverage only; no association fit."""
    site_data={}
    site_dup=0
    all_ids_by_region=collections.defaultdict(set)
    samples={}
    scored_by_order=collections.defaultdict(collections.Counter)
    qual_counts=collections.Counter()
    year_counts=collections.Counter()
    method_counts=collections.Counter()
    source_all_families=set()
    raw_year_total=collections.Counter()
    for region in source.PREFIXES:
        name=source.find_one_by_region(filenames,region,"sitevariables\\.csv")
        for row in source.csv_rows(z,name,["BIO_SITE_ID","EDF_MEAN",*NUMERIC_SOURCE]):
            site=str(row["BIO_SITE_ID"]).strip()
            if not site:continue
            all_ids_by_region[region].add(site)
            ident=(region,site)
            if ident in site_data:
                site_dup+=1;continue
            fields={k:source.finite_numeric(row.get(k)) for k in ("EDF_MEAN",*NUMERIC_SOURCE)}
            site_data[ident]=fields
        taxa=source.find_taxa(filenames,region)
        for row in source.csv_rows(z,taxa,["BIO_SITE_ID","BIO_SAMPLE_ID","BIO_ANALYSIS_ID",
                                           "BIO_SAMPLE_DATE","BIO_ANALYSIS_METHOD",
                                           "BIO_ORD_NAME","BIO_FAM_NAME","BIO_ABUNDANCE"]):
            dt=parsed_date(row["BIO_SAMPLE_DATE"])
            if dt is None:
                qual_counts["unparseable_original_sample_date"]+=1;continue
            raw_year_total[dt.year]+=1
            if not Y_MIN<=dt.year<=Y_MAX:continue
            sid=str(row["BIO_SITE_ID"]).strip()
            aid=str(row["BIO_ANALYSIS_ID"]).strip()
            smp=str(row["BIO_SAMPLE_ID"]).strip()
            method=str(row["BIO_ANALYSIS_METHOD"]).strip() or "__unknown__"
            if not sid or not aid:
                qual_counts["missing_site_or_analysis_source_key"]+=1;continue
            key=(region,sid,aid)
            original=(smp,dt.date().isoformat(),method)
            if key not in samples:
                samples[key]={"meta":original,"date":dt,"families":set(),
                              "family_order":{}, "counts":collections.defaultdict(float),
                              "n_scored_num_count_rows":0,"n_unscored_count_rows":0,
                              "numeric_family":set(),"non_numeric_family":set()}
            v=samples[key]
            if v["meta"]!=original:
                raise RuntimeError("One original ANALYSIS_ID has incompatible sample/collection date")
            fam=str(row["BIO_FAM_NAME"] or "").strip()
            order=str(row["BIO_ORD_NAME"] or "").strip()
            if not fam:continue
            v["families"].add(fam)
            source_all_families.add(fam)
            if fam in scores:
                if order:
                    scored_by_order[fam][order]+=1
                n=source.finite_numeric(row["BIO_ABUNDANCE"])
                if n is not None and n>0:
                    v["counts"][fam]+=n
                    v["numeric_family"].add(fam)
                    v["n_scored_num_count_rows"]+=1
                else:
                    v["non_numeric_family"].add(fam)
                    v["n_unscored_count_rows"]+=1
    if site_dup:raise RuntimeError("Duplicate original sitevariable keys")
    result={
        "official_site_variable_n":len(site_data),
        "2010_18_original_analysis_n":len(samples),
        "2010_18_source_all_family_names":len(source_all_families),
        "2010_18_date_and_source_observations":dict(qual_counts),
        "n_family_source_scores":len(scores),
        "year_original_taxon_row_totals":{str(k):int(v) for k,v in sorted(raw_year_total.items())},
    }
    region_sites=collections.defaultdict(set)
    region_samples=collections.Counter()
    region_edf=collections.defaultdict(set)
    sample_rows=[]
    sample_family_weights=[]
    samples_seen=0
    counts_valid=0
    for (region,site,aid),v in samples.items():
        year=v["date"].year
        month=v["date"].month
        fam=v["families"]
        scored=sorted(fam&set(scores))
        missing=sorted(fam-set(scores))
        cov=len(scored)/len(fam) if fam else 0.0
        site_fields=site_data.get((region,site))
        if site_fields is None:
            qual_counts["no_original_sitevariable_match"]+=1
            continue
        edf=site_fields["EDF_MEAN"]
        if edf is None:
            qual_counts["missing_site_edf"]+=1
            continue
        if edf<0 or edf>MAX_EDF_PERCENT:
            raise RuntimeError("EDF_MEAN outside verified source bounded percentage scale")
        if len(scored)<MIN_SCORED_FAMILIES or cov<MIN_SCORED_FRACTION:
            qual_counts["insufficient_scored_family_coverage"]+=1
            continue
        if any(site_fields[k] is None for k in NUMERIC_SOURCE):
            qual_counts["missing_fixed_environmental_control"]+=1
            continue
        if any(not math.isfinite(site_fields[k]) for k in NUMERIC_SOURCE):
            raise RuntimeError("Nonfinite original source environmental covariate")
        # Origin of abundance units not independently guaranteed. Count CWM is
        # only a sensitivity based on *recorded positive numeric BIO_ABUNDANCE*.
        numeric=len(v["numeric_family"] & set(scored))
        all_numeric=numeric==len(scored) and not (v["non_numeric_family"]&set(scored))
        count_ok=numeric/max(len(scored),1)>=MIN_COUNT_POSITIVE_SCORED_FRACTION
        counts_valid+=int(all_numeric)
        row={
            "region":region,"site":site,"site_region_id":region+"||"+site,
            "sample_date":v["meta"][1],
            "year":year,"month":month,"method":v["meta"][2],
            "edf_pct":edf,"lc_urban":site_fields["LC_URBAN_PER"],
            "lc_arable":site_fields["LC_ARABLE_PER"],
            "width":site_fields["BIO_WIDTH"],"altitude":site_fields["BIO_ALTITUDE"],
            "siltclay":site_fields["BIO_SILT_CLAY"],"n_families":len(fam),
            "n_scored_families":len(scored),"coverage":cov,
            "positive_numeric_family_fraction":numeric/max(len(scored),1),
            "all_scored_family_recorded_numeric_positive":all_numeric,
            "count_secondary_eligible":count_ok}
        sample_rows.append(row)
        sample_family_weights.append({
            "scored":scored,
            "numeric_count":{f:v["counts"][f] for f in scored if v["counts"][f]>0}
        })
        region_samples[region]+=1
        region_sites[region].add(site)
        region_edf[region].add(round(float(edf),6))
        year_counts[year]+=1
        method_counts[row["method"]]+=1
        samples_seen+=1
    region_gates={}
    for r in source.PREFIXES:
        region_gates[r]={
            "eligible_samples":region_samples[r],
            "eligible_distinct_sites":len(region_sites[r]),
            "eligible_edf_distinct":len(region_edf[r]),
            "support_pass":len(region_sites[r])>=MIN_REGION_SITES and
                           len(region_edf[r])>=MIN_REGION_EDF_DISTINCT}
    eligible_sites=len({row["site_region_id"] for row in sample_rows})
    eligible_regions=sum(v["support_pass"] for v in region_gates.values())
    result.update({
        "source_gate":{"target_years":[Y_MIN,Y_MAX],
            "min_samples":MIN_SAMPLES,"min_sites":MIN_SITES,
            "min_regions":MIN_REGIONS,
            "min_region_sites":MIN_REGION_SITES,
            "min_region_EDF_distinct":MIN_REGION_EDF_DISTINCT,
            "min_score_coverage":MIN_SCORED_FRACTION,
            "min_scored_families":MIN_SCORED_FAMILIES,
            "eligible_sample_n":len(sample_rows),
            "eligible_site_n":eligible_sites,
            "eligible_region_n":eligible_regions,
            "region_support":region_gates,
            "eligible_year_sample_n":{str(k):int(v) for k,v in sorted(year_counts.items())},
            "eligible_methods_top10":method_counts.most_common(10),
            "all_scored_family_numeric_sample_n":counts_valid,
            "cohort_ready_for_incidence_CWM_fit":
                len(sample_rows)>=MIN_SAMPLES and eligible_sites>=MIN_SITES and
                eligible_regions>=MIN_REGIONS,
            "source_rejections":dict(qual_counts),
            "require_verified_score_identity":len(scores)==92,
            "source_is_not_true_WWTP_mass_or_P_N":True,
            "source_is_not_time_variable_EDF":True
        }
    })
    return sample_rows,sample_family_weights,scored_by_order,result

def site_year_response(sample_rows,weights,score_map):
    families=sorted(score_map)
    names={name:i for i,name in enumerate(families)}
    fvalues=np.asarray([score_map[k] for k in families],float)
    matrix_count=[]
    matrix_inc=[]
    rows=[]
    for row,raw in zip(sample_rows,weights):
        inc=np.zeros(len(fvalues),float)
        count=np.zeros(len(fvalues),float)
        for f in raw["scored"]:inc[names[f]]+=1.0/len(raw["scored"])
        csum=sum(raw["numeric_count"].values())
        if csum>0:
            for f,n in raw["numeric_count"].items():count[names[f]]+=n/csum
        matrix_inc.append(inc)
        matrix_count.append(count)
        rows.append(row)
    data=pd.DataFrame(rows)
    original_inc=np.asarray(matrix_inc,float)
    original_count=np.asarray(matrix_count,float)
    if not np.allclose(original_inc.sum(axis=1),1.0,atol=1e-12,rtol=0):
        raise RuntimeError("Incidence community weight normalization invalid")
    data["month_sin"]=np.sin(2*np.pi*(data.month.astype(float)-1)/12)
    data["month_cos"]=np.cos(2*np.pi*(data.month.astype(float)-1)/12)
    data["log_n_families"]=np.log1p(data.n_families.astype(float))
    # site×year is inference grain: mean equally over its actual analyses.
    # All high-level W predictor values are fixed at site level.
    grouped=data.groupby(["site_region_id","region","site","year"],sort=True)
    out=[]
    inc_mat=[]
    cnt_mat=[]
    count_supported=0
    for _,g in grouped:
        idx=g.index.to_numpy(int)
        first=g.iloc[0]
        edf_values=g.edf_pct.to_numpy(float)
        if np.max(edf_values)-np.min(edf_values)>1e-12:
            raise RuntimeError("Source site/year EDF not constant")
        new={k:float(first[k]) for k in
             ("year","edf_pct","lc_urban","lc_arable","width","altitude","siltclay")}
        new["region"]=first.region
        new["site"]=first.site
        new["site_region_id"]=first.site_region_id
        new["method"]=g.method.astype(str).value_counts().index[0]
        new["month_sin"]=float(g.month_sin.mean())
        new["month_cos"]=float(g.month_cos.mean())
        new["log_n_families"]=float(g.log_n_families.mean())
        new["coverage"]=float(g.coverage.mean())
        new["original_sample_n"]=len(g)
        A=original_inc[idx].mean(axis=0)
        B=original_count[idx].mean(axis=0)
        new["y_incidence"]=float(A@fvalues)
        new["count_supported"]=bool(g.count_secondary_eligible.all())
        if new["count_supported"]:
            count_supported+=1
            new["y_recorded_abundance"]=float(B@fvalues)
        else:
            new["y_recorded_abundance"]=np.nan
        out.append(new);inc_mat.append(A);cnt_mat.append(B)
    df=pd.DataFrame(out)
    A=np.asarray(inc_mat,float)
    B=np.asarray(cnt_mat,float)
    # The label-permutation universe must consist of species actually
    # observable in the frozen 2010-18 scored source cohort; unused
    # dictionary entries must NOT act as phantom exchangeability slots.
    observed=A.sum(axis=0)>0
    A=A[:,observed]
    B=B[:,observed]
    families=[fam for fam,seen in zip(families,observed) if seen]
    fvalues=fvalues[observed]
    if len(df)!=len(A):raise RuntimeError("site-year-source grain mismatch")
    if not np.allclose(A@fvalues,df.y_incidence,atol=1e-10):
        raise RuntimeError("First-principles CWM incidence reconstruction mismatch")
    if df.duplicated(["site_region_id","year"]).any():
        raise RuntimeError("Duplicate source site-year after group")
    return df,A,B,families,fvalues,count_supported

def model_design(df,include_focal=True):
    df=df.reset_index(drop=True)
    controls=list(COVARIATES)
    X=pd.DataFrame(index=df.index)
    for c in controls:
        v=df[c].to_numpy(float)
        if not np.isfinite(v).all():
            raise RuntimeError("Model original source background contains missing covariate "+c)
        sd=np.std(v)
        X[c]=(v-v.mean())/(sd if sd>1e-9 else 1.0)
    if include_focal:
        raw=np.log1p(df.edf_pct.to_numpy(float))
        sd=raw.std(ddof=0)
        if sd<1e-9:raise RuntimeError("No EDF across-source variance")
        X["EDF_log1p_std"]=(raw-raw.mean())/sd
    for cat in CATS:
        one=pd.get_dummies(df[cat].astype(str),prefix=cat,drop_first=True,dtype=float)
        for name in one:X[name]=one[name].to_numpy(float)
    X=sm.add_constant(X,has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy(float))!=len(X.columns):
        raise RuntimeError("Source controls+EDF collinear/rank deficient")
    return X

def site_weights(df):
    ids=df.site_region_id.astype(str)
    n=ids.value_counts()
    return (1./ids.map(n)).to_numpy(float)

def fit_cluster(df,response):
    # Observation is site-year. Every site gets weight one.
    X=model_design(df,True)
    y=df[response].to_numpy(float)
    if not np.isfinite(y).all():raise RuntimeError("Nonfinite scored CWM source response")
    w=site_weights(df)
    model=sm.WLS(y,X.to_numpy(float),weights=w).fit(
        cov_type="cluster",cov_kwds={"groups":df.site_region_id.astype(str).to_numpy()},
        use_t=False
    )
    name="EDF_log1p_std"
    i=list(X.columns).index(name)
    return {"beta":float(model.params[i]),
            "cluster_site_normal95":[float(v) for v in model.conf_int()[i]],
            "nominal_p":float(model.pvalues[i]),
            "n_site_year":len(df),
            "n_distinct_site":df.site_region_id.nunique(),
            "n_region":df.region.nunique(),
            "design_cols":len(X.columns),
            "adjustment_cols":list(X.columns),
            "mean_original_BMWP":float(y.mean()),
            "sitecluster_not_annual_causality":True}

def null_score_shuffle(df,A,families,values,order_by_family,B=999):
    X=model_design(df,True)
    cols=list(X.columns)
    i=cols.index("EDF_log1p_std")
    x=X.iloc[:,i].to_numpy(float)
    C=X.drop(columns=["EDF_log1p_std"]).to_numpy(float)
    w=site_weights(df);sw=np.sqrt(w)
    residual=x-C@np.linalg.lstsq(sw[:,None]*C,sw*x,rcond=None)[0]
    den=float((w*residual)@residual)
    coeff=(w*residual)@A/den
    observed=float(coeff@values)
    recovered=fit_cluster(df,"y_incidence")["beta"]
    if abs(observed-recovered)>1e-9:
        raise RuntimeError("Original source EDF WLS coefficient not independently reproduced from site×family")
    families_by_order=collections.defaultdict(list)
    unknown=[]
    conflicting=[]
    for i,fam in enumerate(families):
        orders=order_by_family.get(fam,collections.Counter())
        if not orders:
            unknown.append(fam);continue
        if len(orders)>1:
            conflicting.append(fam)
            continue
        order=next(iter(orders))
        families_by_order[order].append(i)
    if len(unknown)+len(conflicting)>len(families)*.1:
        return {"status":"SOURCE_STOP_AMBIGUOUS_FAMILY_ORDER","n_unknown_order":len(unknown),
                "n_conflicting_order":len(conflicting)}
    rng=np.random.default_rng(SEED)
    vgroups=[np.asarray(x,int) for x in families_by_order.values() if len(x)>1]
    draws=np.zeros(B,float)
    for b in range(B):
        v=values.copy()
        for ids in vgroups:
            v[ids]=rng.permutation(values[ids])
        draws[b]=float(coeff@v)
    center=float(draws.mean())
    p=float((np.count_nonzero(np.abs(draws-center)>=np.abs(observed-center)-1e-12)+1)/(B+1))
    return {"status":"WITHIN_ORDER_PERMUTATION",
            "observed_beta":observed,
            "null_center":center,
            "two_sided_centered_tail_p":p,
            "null_percentile95":[float(a) for a in np.quantile(draws,[.025,.975])],
            "B":B,"n_scored_families":len(families),
            "n_shuffle_orders":len(vgroups),
            "n_unresolved_order_families":len(unknown),
            "n_multiorder_families":len(conflicting),
            "not_phylogenetic_or_causal":True}

def region_transfer(df,response):
    preds={}
    valid_regions=[]
    by=[]
    for region in sorted(df.region.unique()):
        test=df[df.region==region].copy()
        train=df[df.region!=region].copy()
        overlapping=set(train.site)&set(test.site)
        if overlapping:
            train=train[~train.site.isin(overlapping)].copy()
        if test.site.nunique()<10 or train.site.nunique()<200:
            preds[region]={"status":"SOURCE_STOP_REGION","test_site":int(test.site.nunique()),
                           "train_site":int(train.site.nunique())}
            continue
        def prediction(cols):
            column=ColumnTransformer([
                ("num",StandardScaler(),list(cols)),
                ("cat",OneHotEncoder(handle_unknown="ignore",sparse_output=False),
                 list(CATS))],remainder="drop",sparse_threshold=0)
            X=column.fit_transform(train);Z=column.transform(test)
            w=site_weights(train)
            return Ridge(alpha=10,fit_intercept=True,solver="svd").fit(
                X,train[response].to_numpy(float),sample_weight=w
            ).predict(Z)
        base=prediction(COVARIATES)
        extra=prediction(list(COVARIATES)+["edf_log1p"])
        y=test[response].to_numpy(float)
        unique=test.site_region_id.astype(str).unique()
        error0=[];error1=[]
        for sid in unique:
            mask=test.site_region_id==sid
            error0.append(float(np.mean((y[mask]-base[mask])**2)))
            error1.append(float(np.mean((y[mask]-extra[mask])**2)))
        a=float(np.mean(error0));b=float(np.mean(error1))
        preds[region]={"status":"REGION_FULL_HOLDOUT",
                       "test_site":len(unique),"train_site":train.site_region_id.nunique(),
                       "baseline_MSE":a,"add_EDF_MSE":b,
                       "EDF_added_MSE_gain_pct":100*(a-b)/a}
        valid_regions.append(region)
        by.append((a,b))
    if len(valid_regions)<MIN_REGIONS:
        return {"status":"SOURCE_STOP_FEWER_5_VALID_REGIONS",
                "valid_region_names":valid_regions,"by_region":preds}
    a=np.asarray([q[0] for q in by]);b=np.asarray([q[1] for q in by])
    rng=np.random.default_rng(SEED)
    ix=rng.integers(0,len(a),size=(N_BOOT,len(a)))
    boot=100*(a[ix].mean(axis=1)-b[ix].mean(axis=1))/a[ix].mean(axis=1)
    return {"status":"REGION_DISJOINT_MODEL_TRANSFER",
            "valid_regions":valid_regions,"n_regions":len(valid_regions),
            "n_improved_regions":int(np.sum(a>b)),
            "region_equal_W_added_MSE_gain_pct":float(100*(a.mean()-b.mean())/a.mean()),
            "fixed_error_2000_region_boot95":[float(v) for v in np.quantile(boot,[.025,.975])],
            "by_region":preds,
            "interval_is_fixed_errors_only":True}

def main(args):
    output={"schema":"GESUI_CHEMPOP_V2_BMWP_EDF_TRAIT_SELECTION_V1",
            "source":{"original_EIDC_UUID":source.EIDC_UUID,
                      "separate_BMWP_score_commit":source.SCORE_COMMIT},
            "exploratory_posthoc":True,
            "not_independent_contemporaneous_WWTP_NP_mass":True,
            "does_not_assert_within_site_EDF_temporal_change":True,
            "status":"SOURCE_STOP",
            "protocol_before_outcome":"docs/PROTOCOL_20261119_CHEMPOP_V2_BMWP_INC_COUNT_TRANSFER.md"}
    try:
        with tempfile.TemporaryDirectory(prefix="gesui_ChemPop_BMWP_") as tmp:
            temp=Path(tmp)
            arch=source.retrieve(source.EIDC_URL,temp/"original.zip",source.SOURCE_DOWNLOAD_LIMIT)
            scores,sourcelog=source.score_reference(temp)
            output["source"].update({"EIDC_ZIP_SHA256":arch["SHA256"],
                                     "EIDC_original_bytes":arch["source_archive_bytes"],
                                     "score_rdata_SHA256":sourcelog["SHA256"]})
            with zipfile.ZipFile(temp/"original.zip") as z:
                files=z.namelist()
                if len(files)!=326 or len(set(files))!=326:raise RuntimeError("Original publisher ZIP contents drift")
                raw,weight,orders,source_state=construct_eligible(z,files,scores)
                output["source_gate"]=source_state
        if not source_state["source_gate"]["cohort_ready_for_incidence_CWM_fit"]:
            output.update({"status":"SOURCE_STOP_SOURCE_GATES",
                           "reason":"Existing pre-outcome site/sample/region source threshold not passed"})
        else:
            df,A,B,families,values,nsupported=site_year_response(raw,weight,scores)
            # First inspect geometry/variance/covariate rank BEFORE opening regression.
            by_region=df.groupby("region").site_region_id.nunique()
            if len(df)<250 or df.site_region_id.nunique()<MIN_SITES or sum(by_region>=10)<MIN_REGIONS:
                output.update({"status":"SOURCE_STOP_SITEYEAR_SUPPORT",
                               "reason":"Site-year collapse reduced source support",
                               "n_original_siteyear":len(df),"n_site":df.site_region_id.nunique(),
                               "region_site_n":by_region.to_dict()})
            else:
                df["edf_log1p"]=np.log1p(df.edf_pct)
                # All predictors available. Build matrix rank/gates pre outcome.
                design=model_design(df,True)
                if np.linalg.matrix_rank(design.to_numpy(float))<design.shape[1]:
                    raise RuntimeError("Model source environmental design not full rank")
                res={"status":"POSTHOC_OBSERVATIONAL_COEFFICIENTS",
                     "n_site_year":len(df),
                     "n_original_2010_18_samples":len(raw),
                     "n_distinct_sites":df.site_region_id.nunique(),
                     "n_regions":df.region.nunique(),
                     "n_count_secondary_supported_siteyears":nsupported,
                     "source_FAMILY_scoring_values_n":len(values),
                     "source_method_counts":df.method.value_counts().to_dict(),
                     "recorded_abundance_not_performance_or_density":True,
                     "EDF_is_modelled_site_level_not_time_varying":True,
                     "analysis":{}}
                incfit=fit_cluster(df,"y_incidence")
                perm=null_score_shuffle(df,A,families,values,orders,B=N_PERM)
                reg=region_transfer(df,"y_incidence")
                res["analysis"]["incidence_primary"]={"fit":incfit,"within_order_label_null":perm,
                                                      "out_region_predictive_added_EDF":reg}
                if nsupported>=200 and df.loc[df.count_supported].site_region_id.nunique()>=100:
                    other=df[df.count_supported].copy()
                    res["analysis"]["recorded_abundance_secondary"]={
                        "fit":fit_cluster(other,"y_recorded_abundance"),
                        "out_region_predictive_added_EDF":region_transfer(other,"y_recorded_abundance"),
                        "n_source_eligible_siteyear":len(other)}
                else:
                    res["analysis"]["recorded_abundance_secondary"]={
                        "status":"SOURCE_STOP_ABUNDANCE_COVERAGE","n_count_siteyear":nsupported}
                output.update(res)
    except Exception as e:
        output.update({"status":"SOURCE_OR_METHOD_STOP",
                       "error":type(e).__name__+": "+str(e)[:700]})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    digest={"status":output["status"],"source":output["source"],
            "source_gate":output.get("source_gate",{}).get("source_gate",{}),
            "n_site_year":output.get("n_site_year"),
            "n_sites":output.get("n_distinct_sites"),
            "analysis":output.get("analysis"),
            "error":output.get("error")}
    print("GESUI_CHEMPOP_BMWP_EDF_TRAIT_TEST="+json.dumps(
        digest,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",required=True,type=Path)
    main(parser.parse_args())
