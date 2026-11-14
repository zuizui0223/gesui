#!/usr/bin/env python3
"""Source-authentic five-family butterfly comparison in 454 EU H7×period rows.

Posthoc *observational* family composition only. No measured wastewater
chemistry, direct puddling/feeding or plant/butterfly interactions. The five
family source shards are restored from original pre-existing Actions artifacts,
not re-fetched or re-selected from fitted signs.
"""
import collections
import csv
import hashlib
import json
import math
import os
import pathlib
import numpy as np

FAMILIES=("Lycaenidae","Nymphalidae","Pieridae","Papilionidae","Hesperiidae")
SHA_H6="372b800b85abe2f395241987cea2000ff901cf5aff31bf1d1fb825304f2495fc"
SHA_ATLAS="ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686"
SOURCE_SHARD_SCHEMA="HYBAS7_234_EU_BOTH_EEA_SOURCE_PLANT_S50_AND_FIVE_BUTTERFLY_FAMILIES_original_v1"
OUT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_20261114_5BUTTERFLY_FAMILY_MACRO"
OUT.mkdir(parents=True,exist_ok=True)

def read_json(key,sha=None):
    p=pathlib.Path(os.environ[key])
    b=p.read_bytes()
    got=hashlib.sha256(b).hexdigest()
    if sha and got!=sha:
        raise RuntimeError("Original source SHA mismatch "+key+": "+got)
    return json.loads(b),got

def read_csv(key):
    p=pathlib.Path(os.environ[key])
    b=p.read_bytes()
    with p.open(encoding="utf-8-sig",newline="") as handle:
        rows=list(csv.DictReader(handle))
    return rows,hashlib.sha256(b).hexdigest()

def count_check(key,val):
    if not isinstance(val,int) or val<0:
        raise RuntimeError("Invalid original source count "+key)
    return val

def feature_numbers(row):
    fields=("log_local_H7_area_km2","log1p_original_EQR_site_n",
        "log1p_GBIF_Animalia_source_records","log1p_GBIF_Plantae_source_records",
        "log1p_GBIF_5family_butterfly_source_records",
        "log1p_original_EEA_site_DIS_AV_CMS","crop_up","urban_up","naturalout")
    vals=[float(row[x]) for x in fields]
    if not all(math.isfinite(z) for z in vals):
        raise RuntimeError("Background source covariate missing / invalid")
    return vals

def source():
    h6,sha_h6=read_json("ORIGINAL_TRUE_H6_GEOGRAPHY",SHA_H6)
    atlas0,sha_atlas=read_json("ORIGINAL_H7_HYDROATLAS_SOURCE",SHA_ATLAS)
    geom=h6["original_454_H7_study_basin_geographic_parent_mapping"]
    h6_parent={str(q["HYBAS_ID_L07"]):str(q["official_parent_HYBAS_ID_L06"]) for q in geom}
    if len(h6_parent)!=232 or len(set(h6_parent.values()))!=90:
        raise RuntimeError("Original genuine 232 H7 / 90 H6 shape changed")
    atlas=atlas0["original_1219_H7_covariate_by_real_ID"]
    all_families={}
    shards={}
    for i in range(4):
        j,sh=read_json("ORIGINAL_FAMILY_SHARD_"+str(i))
        qa,_=read_json("ORIGINAL_FAMILY_QA_"+str(i))
        if qa.get("mode")!="full" or qa.get("shard")!=i or \
            qa.get("original_expected_field_queries")!=len(j):
            raise RuntimeError("Frozen original 5-family shard QA mismatch "+str(i))
        for q in j:
            field=q["source_field"]
            if field not in FAMILIES:continue
            k=(str(q["HYBAS_ID_L07"]),q["period"],field)
            if k in all_families:raise RuntimeError("Duplicated original family-period field")
            all_families[k]=q
        shards[str(i)]=sh
    if len(all_families)!=234*2*5:
        raise RuntimeError("Original 234×2×5 source support incomplete")

    rows,raw_sha=read_csv("ORIGINAL_454_EEA_GBIF_SOURCE_ROWS")
    if len(rows)!=454 or len({(str(r["HYBAS_ID_L07"]),r["PERIOD"]) for r in rows})!=454:
        raise RuntimeError("Original 454-row same-unit panel shifted")
    if len({str(r["HYBAS_ID_L07"]) for r in rows})!=232:
        raise RuntimeError("Original 232-basin panel shifted")

    for r in rows:
        h=str(r["HYBAS_ID_L07"]);p=r["PERIOD"]
        if h not in h6_parent or h not in atlas:
            raise RuntimeError("Missing genuine H6 or HydroATLAS geography")
        r["H6_parent"]=h6_parent[h]
        x=atlas[h]
        r["crop_up"]=float(x["crp_pc_use"])
        r["urban_up"]=float(x["urb_pc_use"])
        r["naturalout"]=math.log1p(float(x["dis_m3_pyr"]))
        family_counts={}
        for f in FAMILIES:
            s=all_families[(h,p,f)]
            if s["status"]!="ok":
                raise RuntimeError("A source butterfly family GBIF query was NOT complete "+str((h,p,f)))
            family_counts[f]=count_check(f,s["source_record_index_count"])
        n=sum(family_counts.values())
        # The original signed 454 CSV saved log1p(total) and the Lycaenidae
        # share/logodds, NOT the untransformed family counts. Recover their
        # exact intended comparison from original source shards only.
        original_total=round(math.expm1(float(r["log1p_GBIF_5family_butterfly_source_records"])))
        if n!=original_total:
            raise RuntimeError("Original five-family denominator differs from prior signed cohort")
        if n<20:
            raise RuntimeError("Prior source eligibility n>=20 source count failed")
        lyr=family_counts["Lycaenidae"]
        if abs(lyr/n-float(r["lycaenidae_fivefamily_index_share"]))>1e-10:
            raise RuntimeError("Family share fidelity failed")
        published_logodds=math.log((lyr+0.5)/(n-lyr+0.5))
        if abs(published_logodds-float(r["lycaenidae_fivefamily_index_share_logodds"]))>1e-10:
            raise RuntimeError("Original recorded Lycaenidae logodds could not be reconstructed")
        r["family_original_counts"]=family_counts
        r["native_obs_total"]=n
        r["W"]=float(r["W_static_upstream_MODEL_flow_per_UpArea"])
        if not (math.isfinite(r["W"]) and r["W"]>=0):
            raise RuntimeError("Original wastewater OUTLET modeled source invalid")
        if not (0<=r["crop_up"]<=100 and 0<=r["urban_up"]<=100):
            raise RuntimeError("HydroATLAS source percent out of range")
        feature_numbers(r)
    return rows,{"original_454_panel_sha256":raw_sha,"original_H6_source_sha256":sha_h6,
                 "original_HydroATLAS_sha256":sha_atlas,"original_family_shard_shas":shards}

def stable_country(rows):
    ct=collections.defaultdict(collections.Counter)
    for r in rows:ct[r["HYBAS_ID_L07"]][r["country"]]+=1
    return {h:sorted(counter,key=lambda v:(-counter[v],v))[0] for h,counter in ct.items()}

def outcome_values(rows):
    vals={}
    for label,a,b in [("Lycaenidae_vs_Pieridae","Lycaenidae","Pieridae"),
                      ("Lycaenidae_vs_Nymphalidae","Lycaenidae","Nymphalidae")]:
        vals[label]=np.array([math.log((r["family_original_counts"][a]+.5)/
                      (r["family_original_counts"][b]+.5)) for r in rows])
    for f in FAMILIES:
        vals[f+"_versus_other_four"]=np.array([
            math.log((r["family_original_counts"][f]+.5)/
                 (r["native_obs_total"]-r["family_original_counts"][f]+.5))
            for r in rows])
    return vals

def create_design(rows, country_effect):
    n=len(rows)
    B=np.array([feature_numbers(z) for z in rows],dtype=float)
    countries=sorted({r["country"] for r in rows})
    methods=sorted({r["original_majority_BQE_method"] for r in rows})
    period=np.array([float(z["PERIOD"]=="2020-2022") for z in rows])[:,None]
    c=np.array([[float(r["country"]==t) for t in countries[1:]] for r in rows],
               dtype=float) if country_effect else np.empty((n,0))
    m=np.array([[float(r["original_majority_BQE_method"]==t) for t in methods[1:]]
                for r in rows],dtype=float)
    baseline=np.column_stack([np.ones(n),B,c,m,period])
    W=np.array([r["W"] for r in rows],dtype=float)
    if np.std(W)<1e-10:raise RuntimeError("No WWTP geographic pressure variation")
    return baseline,np.column_stack([baseline,W]),W

def ridge_predictions(X,Y,heldout_groups,source_h6):
    n=len(Y)
    predicted=np.full(n,np.nan)
    for g in sorted(set(heldout_groups)):
        train=heldout_groups!=g
        test=~train
        if sum(train)<150 or not np.any(test):
            raise RuntimeError("True complete holdout insufficient source support")
        Z=X.copy()
        mu=np.mean(Z[train,1:],axis=0)
        sd=np.std(Z[train,1:],axis=0)
        sd=np.where(sd>1e-12,sd,1)
        Z[:,1:]=(Z[:,1:]-mu)/sd
        A=Z[train];t=Z[test]
        h=source_h6[train]
        freq=collections.Counter(h)
        weights=np.array([1.0/freq[x] for x in h])
        weights/=weights.mean()
        gram=A.T@(weights[:,None]*A)
        rhs=A.T@(weights*Y[train])
        penalty=np.eye(A.shape[1])*10.
        penalty[0,0]=0
        beta=np.linalg.solve(gram+penalty,rhs)
        predicted[test]=t@beta
    if not np.all(np.isfinite(predicted)):
        raise RuntimeError("Holdout prediction incomplete")
    return predicted

def cv_gain(Y,X0,X1,heldout,source_h6,grain,seed=20261114):
    a=(Y-ridge_predictions(X0,Y,heldout,source_h6))**2
    b=(Y-ridge_predictions(X1,Y,heldout,source_h6))**2
    if grain=="H6":
        groups=sorted(set(source_h6))
        old=np.array([np.mean(a[source_h6==g]) for g in groups])
        new=np.array([np.mean(b[source_h6==g]) for g in groups])
        rng=np.random.default_rng(seed)
        ids=rng.integers(0,len(groups),size=(2000,len(groups)))
        bootstrap=100*(old[ids].mean(axis=1)-new[ids].mean(axis=1))/old[ids].mean(axis=1)
        return {"group_equal_out_of_sample_MSE_baseline":float(old.mean()),
                "group_equal_out_of_sample_MSE_with_WWTP":float(new.mean()),
                "group_equal_MSE_gain_percent":float(100*(old.mean()-new.mean())/old.mean()),
                "H6_2000_resample_descriptive95_gain_percent":np.percentile(bootstrap,[2.5,97.5]).tolist(),
                "row_weighted_MSE_gain_percent":float(100*(a.mean()-b.mean())/a.mean())}
    if grain=="country":
        countries=sorted(set(heldout))
        per={}
        for g in countries:
            mask=heldout==g
            e0,e1=float(a[mask].mean()),float(b[mask].mean())
            per[str(g)]={"n_original_rows":int(sum(mask)),"MSE_gain_percent":float(100*(e0-e1)/e0),
                         "baseline_MSE":e0,"source_MSE":e1}
        supported=[g for g in countries if per[str(g)]["n_original_rows"]>=10 and
                    len({r for r in source_h6[heldout==g]})>=3]
        major_old=np.mean([per[str(g)]["baseline_MSE"] for g in supported])
        major_new=np.mean([per[str(g)]["source_MSE"] for g in supported])
        return {"country_holdout_row_weighted_gain_percent":float(100*(a.mean()-b.mean())/a.mean()),
                "supported_countries_rows_at_least10_H6_at_least3":list(supported),
                "supported_country_equal_MSE_gain_percent":float(100*(major_old-major_new)/major_old),
                "per_country":per}
    raise RuntimeError("Unexpected validation type")

def descriptive_beta(Y,X1,source_h6):
    Z=X1.copy()
    mu=Z[:,1:].mean(0)
    sd=Z[:,1:].std(0)
    sd=np.where(sd>1e-12,sd,1)
    Z[:,1:]=(Z[:,1:]-mu)/sd
    b=np.linalg.lstsq(Z,Y,rcond=None)[0]
    e=Y-Z@b
    groups=sorted(set(source_h6))
    bread=np.linalg.pinv(Z.T@Z)
    meat=np.zeros((Z.shape[1],Z.shape[1]))
    for g in groups:
        z=Z[source_h6==g]
        u=e[source_h6==g]
        v=z.T@u
        meat+=np.outer(v,v)
    rank=np.linalg.matrix_rank(Z)
    corr=len(groups)/(len(groups)-1)*(len(Y)-1)/(len(Y)-rank)
    var=bread@meat@bread*corr
    se=float(math.sqrt(max(0,var[-1,-1])))
    return {"descriptive_beta_per_plus1SD_original_modelled_WWTP":float(b[-1]),
            "descriptive_H6_cluster_95":[float(b[-1]-1.96*se),float(b[-1]+1.96*se)],
            "not_a_sewage_chemical_or_butterfly_abundance_causal_effect":True}

def main():
    rows,source_refs=source()
    Y=outcome_values(rows)
    h6=np.array([r["H6_parent"] for r in rows])
    country_map=stable_country(rows)
    country=np.array([country_map[r["HYBAS_ID_L07"]] for r in rows])
    fold=np.array([int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5 for x in h6])
    if len(set(h6))!=90 or min(np.bincount(fold,minlength=5))<15:
        raise RuntimeError("Authentic H6 blocked folds changed")
    X0,X1,W=create_design(rows,True)
    L0,L1,_=create_design(rows,False)
    output={"schema":"GESUI_EU454_BUTTERFLY5_ORIGINAL_FAMILY_CONTRASTS_POSTHOC_v1",
            "source_provenance":source_refs,
            "n_original_basin_period_rows":len(rows),
            "n_H7":232,"n_authentic_H6_parents":90,
            "country_holdout_group_counts":dict(collections.Counter(country)),
            "all_families":list(FAMILIES),
            "two_primary_comparisons":["Lycaenidae_vs_Pieridae","Lycaenidae_vs_Nymphalidae"],
            "exposure_type":"MODELLED HydroWASTE static outlet to H7 basin OUTLET, not measured sewage nutrient dose",
            "butterfly_response_type":"GBIF five-family observation-index count composition, not observed feeding or biomass",
            "all_new_results_are_posthoc_exploratory":True,
            "model_training":"Ridge10, group(H6)-equal train weights, train-only standardization, H6 5fold or stable-country leaveout",
            "results":{}}
    for i,(name,y) in enumerate(Y.items()):
        output["results"][name]={
            "original_observation_summary":{"mean_logratio":float(y.mean()),"SD_logratio":float(y.std())},
            "same_H6_parent_hydrological_geoblock_prediction":cv_gain(y,X0,X1,fold,h6,"H6",20261114+i),
            "country_left_out_no_country_FE_prediction":cv_gain(y,L0,L1,country,h6,"country"),
            "within_sample_H6_cluster_observational":descriptive_beta(y,X1,h6)}
    out=OUT/"RESULTS_EU454_5BUTTERFLY_FAMILY_WWTP_H6_COUNTRY_COMPARISON.json"
    out.write_text(json.dumps(output,indent=2,sort_keys=True,allow_nan=False),encoding="utf-8")
    small={k:{"H6":v["same_H6_parent_hydrological_geoblock_prediction"],
              "country_heldout":{x:y for x,y in v["country_left_out_no_country_FE_prediction"].items() if x!="per_country"},
              "beta":v["within_sample_H6_cluster_observational"]}
           for k,v in output["results"].items()}
    print("GESUI_BUTTERFLY_FAMILY_COMPARISON_START",flush=True)
    print(json.dumps({"source_sample":(len(rows),len(set(h6))),
                     "results":small},sort_keys=True,allow_nan=False),flush=True)
    print("GESUI_BUTTERFLY_FAMILY_COMPARISON_END",flush=True)

if __name__=="__main__":
    main()
