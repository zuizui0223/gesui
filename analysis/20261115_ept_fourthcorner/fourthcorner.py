#!/usr/bin/env python3
"""EPA NRSA three-survey-cycle taxon x trait x river environment ("fourth-corner-style").

Original EPA site×taxon counted individuals, source-verified taxon FFG/HABIT/PTV,
and native UID-matched river TN,TP,PHab. Tests whether the embeddedness gradient
represents loss of trait-bearing taxon incidence or changes in counted dominance.
NO phylogenetic adaptation, wastewater-specific causality, ecological flux or
formal fully identified causal path is inferred.

2026-11-15 exploratory/posthoc relative to known weighted trait results PR83.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,os
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.linear_model import Ridge

TEMP=Path(os.environ.get("RUNNER_TEMP","/tmp"))
ROOT=TEMP/"gesui_threecycle_traits"
COUNT=TEMP/"gesui_threecycle_rarefied_taxa"/"AUDIT_ONLY_ORIGINAL_EPA_3CYCLE_TAXON_COUNTS.csv"
TRAIT=ROOT/"AUDIT_ONLY_EPA_ORIGINAL_3CYCLE_TAXA_FFG_HABIT_PTV.csv"
SITE=ROOT/"AUDIT_ONLY_3CYCLE_ORIGINAL_EPT_FFG_HABIT_PTV_SAMESITE.csv"
SOURCE=TEMP/"gesui_threecycle_ept_transfer"/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv"
QA=ROOT/"NRSA_3CYCLE_ORIGINAL_EPT_FFG_HABIT_PTV_SOURCE_ONLY_SITES.json"
CYCLES=("2008_09","2013_14","2018_19")
FIELD={"clinger":("HABIT_norm","EPT_clinger_HABIT_fraction"),
       "tolerance":("PTV_norm","EPT_mean_pollution_tolerance_PTV"),
       "scraper":("FFG_norm","EPT_scraper_FFG_fraction")}
FFG=set("CF CG PA PI PR SC SH".split())
HAB=set("BU CB CN DV PK SK SP SW".split())
GEO=["lat","lon","stream_order","season_sin","season_cos","log_cond"]
NP=["log_tn","log_tp"]
PHAB=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]
CATEGORICAL=["state","eco9"]
PRIMARY=["clinger","tolerance"]

def score(v,kind):
    if kind=="tolerance":
        try:z=float(v)
        except (ValueError,TypeError):return np.nan
        return z if np.isfinite(z) and 0<=z<=10 else np.nan
    code=str(v).upper().strip()
    if not code:return np.nan
    parts=list(dict.fromkeys(x.strip() for x in code.split(",") if x.strip()))
    allowed,desired=(HAB,"CN") if kind=="clinger" else (FFG,"SC")
    if not parts or any(x not in allowed for x in parts):return np.nan
    return float(desired in parts)/len(parts)

def read():
    original=pd.read_csv(SOURCE,dtype={"source_cycle":str,"source_UID":str,"huc8":str})
    if len(original)!=2648 or original.duplicated(["source_cycle","source_UID"]).any():
        raise ValueError("PR81 source original 2648 UID population changed")
    count=pd.read_csv(COUNT,dtype={"cycle":str,"UID":str,"key":str,"order":str},keep_default_na=False)
    count=count[count.order.isin(("EPHEMEROPTERA","PLECOPTERA","TRICHOPTERA"))].copy()
    count["count"]=pd.to_numeric(count["count"],errors="raise")
    if count.duplicated(["cycle","UID","key"]).any() or not count["count"].gt(0).all():
        raise ValueError("Original EPT source count grain/count failed")
    trait=pd.read_csv(TRAIT,dtype=str,keep_default_na=False)
    if trait.duplicated(["source_cycle","key"]).any():
        raise ValueError("Original EPA source trait lookup is not unique")
    count=count.merge(trait[["source_cycle","key","order"]+[x[0] for x in FIELD.values()]],
                      left_on=["cycle","key","order"],right_on=["source_cycle","key","order"],
                      how="left",validate="many_to_one")
    if len(count)==0:raise ValueError("Original EPT trait/count merge empty")
    sites=pd.read_csv(SITE,dtype={"source_cycle":str,"source_UID":str,"huc8":str,"state":str,"eco9":str})
    if len(sites)!=2209 or sites.duplicated(["source_cycle","source_UID"]).any():
        raise ValueError("Original 2209 rare200 same-site source drift")
    qa=json.loads(QA.read_text())
    if qa["schema"]!="gesui_EPA3cycle_SOURCEONLY_ORIGINAL_EPT_SC_CN_PTV_nativeUID_true_coverage_v1":
        raise ValueError("Original trait-site QA schema changed")
    return original,count,sites,qa

def compute(count,sites,qa):
    out={}
    matrices={}
    for cycle in CYCLES:
        dsite=sites[sites.source_cycle==cycle].sort_values("source_UID").copy()
        c=count[(count.cycle==cycle)&count.UID.isin(set(dsite.source_UID))].copy()
        for kind,(original_field,original_site_field) in FIELD.items():
            e=c[["UID","key","order","count",original_field]].copy()
            e["trait"]=e[original_field].map(lambda v:score(v,kind))
            denom=e.groupby("UID")["count"].sum()
            known=e[np.isfinite(e.trait)].copy()
            numerator=(known["count"]*known.trait).groupby(known.UID).sum()
            knownmass=known.groupby("UID")["count"].sum()
            ntaxa=known.groupby("UID")["key"].nunique()
            all_taxa=e.groupby("UID")["key"].nunique()
            incidence=known.groupby("UID")["trait"].mean()
            calc_count=(numerator/knownmass)
            cover=(knownmass/denom).reindex(dsite.source_UID).fillna(0).to_numpy(float)
            published=dsite[original_site_field].to_numpy(float)
            nw=calc_count.reindex(dsite.source_UID).to_numpy(float)
            inc=incidence.reindex(dsite.source_UID).to_numpy(float)
            nknown=ntaxa.reindex(dsite.source_UID).fillna(0).to_numpy(int)
            nall=all_taxa.reindex(dsite.source_UID).fillna(0).to_numpy(int)
            taxon_coverage=nknown/np.maximum(nall,1)
            mask=cover>=0.8-1e-10
            if not np.array_equal(mask,np.isfinite(published)):
                raise ValueError(f"Original pre-outcome 80pct source cohort changed: {cycle} {kind}")
            if not np.allclose(nw[mask],published[mask],rtol=0,atol=1e-10):
                raise ValueError(f"Previous source's count-CWM not independently reproduced: {cycle} {kind}")
            q=qa["source_original_sites_by_cycle_and_trait"][cycle]["original_explicit_trait_80pct_or_more_site_eligibility"][original_site_field]
            sel=dsite.loc[mask].copy().reset_index(drop=True)
            if (len(sel),sel.huc8.nunique(),sel.state.nunique()) != (q["n_native_UID"],q["n_HUC8"],q["states"]):
                raise ValueError("Original source locked trait cohort changed "+cycle+kind)
            if not len(sel)>=250 or not sel.huc8.nunique()>=150 or not sel.state.nunique()>=20:
                raise ValueError("Original source support minimum failed "+cycle+kind)
            sel["CWM_count"]=nw[mask]
            sel["CWM_incidence"]=inc[mask]
            sel["count_minus_incidence"]=nw[mask]-inc[mask]
            sel["known_taxa_n"]=nknown[mask]
            sel["all_EPT_taxa_n"]=nall[mask]
            sel["trait_taxon_coverage"]=taxon_coverage[mask]
            sel["trait_record_coverage"]=cover[mask]
            if not np.isfinite(sel[["CWM_count","CWM_incidence"]]).all().all():
                raise ValueError("Unknown incidence trait treated as 0")
            if not np.isfinite(sel[GEO+NP+PHAB].to_numpy(float)).all():
                raise ValueError("Missing source covariate")
            # Fixed site/taxon trait-matrix for within-order trait-label null.
            ec=known[known.UID.isin(set(sel.source_UID))].copy()
            ec["id"]=ec.key+"||"+ec.order
            ids=sorted(ec.id.unique())
            si={u:i for i,u in enumerate(sel.source_UID)}
            ti={k:i for i,k in enumerate(ids)}
            cnts=np.zeros((len(sel),len(ids)),dtype=float)
            raw_trait=np.full(len(ids),np.nan)
            taxa_order={}
            for row in ec.itertuples(index=False):
                i=si[row.UID];j=ti[row.id]
                cnts[i,j]+=float(row.count)
                if np.isfinite(raw_trait[j]) and abs(raw_trait[j]-float(row.trait))>1e-12:
                    raise ValueError("Taxon trait inconsistent within original source")
                raw_trait[j]=float(row.trait)
                taxa_order[row.id]=row.order
            if not np.isfinite(raw_trait).all():raise ValueError("Taxon score missing after known restriction")
            xcount=cnts/cnts.sum(axis=1,keepdims=True)
            xinc=(cnts>0).astype(float)
            xinc/=xinc.sum(axis=1,keepdims=True)
            if not np.allclose(xcount@raw_trait,sel.CWM_count,atol=1e-10):
                raise ValueError("Site x taxon matrix fails count CWM numeric equality")
            if not np.allclose(xinc@raw_trait,sel.CWM_incidence,atol=1e-10):
                raise ValueError("Site x taxon matrix fails incidence CWM numeric equality")
            out[(cycle,kind)]=sel
            matrices[(cycle,kind)]=(xcount,xinc,raw_trait,np.array([taxa_order[k] for k in ids]),ids)
    return out,matrices

def model(sel,ykey):
    # Other GEO/NP/PHAB remain in their genuine source units; focal embedding and
    # NP scaled by within-cycle SD, as in the prior #83 within-era trait regression.
    xx=pd.concat([sel[GEO+NP+PHAB].astype(float),
                  pd.get_dummies(sel[CATEGORICAL],drop_first=True,dtype=float)],axis=1)
    for v in ["XEMBED","log_tn","log_tp"]:
        a=xx[v].to_numpy(float);sd=float(a.std(ddof=0))
        if sd<1e-10:raise ValueError("Near-constant original source env variable")
        xx[v]=(a-a.mean())/sd
    X=sm.add_constant(xx.to_numpy(float),has_constant="add")
    if np.linalg.matrix_rank(X)!=X.shape[1]:
        raise ValueError("Site×trait regression rank deficient")
    ids=sel.huc8.astype(str).to_numpy()
    counts=collections.Counter(ids)
    w=np.array([1./counts[g] for g in ids])
    fit=sm.WLS(sel[ykey].to_numpy(float),X,weights=w).fit(
           cov_type="cluster",cov_kwds={"groups":ids},use_t=False)
    res={}
    for v in ["XEMBED","log_tn","log_tp"]:
        idx=xx.columns.get_loc(v)+1
        beta=float(fit.params[idx]);interval=list(map(float,fit.conf_int()[idx]))
        res[v]={"beta":beta,"HUC8_cluster_normal95":interval,
                "p_nominal":float(fit.pvalues[idx])}
    return res

def taxon_label_permutation(sel,matrix,trait,target="XEMBED",B=399,seed=20261115):
    """Fourth-corner-style within-ORDER taxon trait assignment null.

    Holds site×taxon counts, all site geography/hydrology and trait missingness
    fixed, permutes scored TAXA_ID trait labels WITHIN each EPT order.
    The null is not randomizing geographic exposure and cannot identify causality.
    """
    abundance,incidence,values,orders,keys=matrix
    base=GEO+NP+[v for v in PHAB if v!=target]
    control=pd.concat([sel[base].astype(float),
                       pd.get_dummies(sel[CATEGORICAL],drop_first=True,dtype=float)],axis=1)
    x=sel[target].to_numpy(float)
    x=(x-x.mean())/x.std(ddof=0)
    C=sm.add_constant(control.to_numpy(float),has_constant="add")
    if np.linalg.matrix_rank(C)!=C.shape[1]:raise ValueError("Null control matrix rank-deficient")
    group=sel.huc8.astype(str).to_numpy()
    n=collections.Counter(group);w=np.array([1./n[i] for i in group])
    sw=np.sqrt(w)
    coef=np.linalg.lstsq(sw[:,None]*C,sw*x,rcond=None)[0]
    residual=x-C@coef
    den=float(np.dot(w*residual,residual))
    if den<1e-9:raise ValueError("Target environmental exposure fully predicted")
    # For each site-composition matrix L: beta_X=(xres' W L t)/(xres' W xres)
    contribution={k:(w*residual)@L/den
                  for k,L in (("CWM_count",abundance),("CWM_incidence",incidence))}
    assert len(values)==len(keys)==len(orders)
    rng=np.random.default_rng(seed)
    labelgroups=[np.where(orders==order)[0] for order in sorted(set(orders))]
    obs={k:float(np.dot(v,values)) for k,v in contribution.items()}
    simulated={k:np.empty(B) for k in obs}
    for i in range(B):
        p=values.copy()
        for ix in labelgroups:
            if len(ix)>1:p[ix]=rng.permutation(values[ix])
        for k,v in contribution.items():simulated[k][i]=float(np.dot(v,p))
    result={}
    for k in obs:
        null=simulated[k]
        center=float(null.mean())
        p=(1+int(np.sum(np.abs(null-center)>=abs(obs[k]-center)-1e-12)))/(B+1)
        result[k]={"observed_beta":obs[k],"null_center":center,
                   "within_order_trait_shuffle_two_sided_tail":float(p),
                   "null_percentile_95":list(map(float,np.percentile(null,[2.5,97.5]))),
                   "taxon_n":len(values),"taxon_order_n":len(labelgroups)}
    return result

def external_predict(dataset,source_huc,kind,metric):
    """Hold out an entire EPA survey cycle and all of its HUC8 from training."""
    GEO_N=GEO+NP;FULL=GEO+NP+PHAB
    out={}
    for era in CYCLES:
        test=dataset[(era,kind)]
        pretrain=pd.concat([dataset[(c,kind)] for c in CYCLES if c!=era],ignore_index=True)
        train=pretrain.loc[~pretrain.huc8.isin(source_huc[era])].copy()
        if (set(test.huc8)&set(train.huc8)) or len(train)<300 or train.huc8.nunique()<200:
            out[era]={"status":"source_or_training_gate_FAIL","test_n":len(test),
                       "train_n":len(train)};continue
        def pred(cols):
            pipe=ColumnTransformer([
              ("numeric",StandardScaler(),cols),
              ("cat",OneHotEncoder(handle_unknown="ignore",sparse_output=False),CATEGORICAL)
            ],remainder="drop",sparse_threshold=0)
            X=pipe.fit_transform(train);V=pipe.transform(test)
            group_sizes=train.huc8.value_counts()
            weights=(len(train)/len(group_sizes))/train.huc8.map(group_sizes)
            return Ridge(alpha=10,fit_intercept=True,solver="svd").fit(
              X,train[metric].to_numpy(float),sample_weight=weights.to_numpy(float)
            ).predict(V)
        y=test[metric].to_numpy(float)
        pa,pb=pred(GEO_N),pred(FULL)
        g,inv=np.unique(test.huc8.astype(str),return_inverse=True)
        cnt=np.bincount(inv)
        ea=np.bincount(inv,weights=(y-pa)**2)/cnt
        eb=np.bincount(inv,weights=(y-pb)**2)/cnt
        ix=np.random.default_rng(20261115).integers(0,len(g),size=(2000,len(g)))
        b=100*(ea[ix].mean(axis=1)-eb[ix].mean(axis=1))/ea[ix].mean(axis=1)
        out[era]={"status":"HUC8_and_era_independent_prediction",
                  "test_n":len(test),"test_HUC8":len(g),"train_n":len(train),
                  "train_HUC8":train.huc8.nunique(),
                  "PHab_added_MSE_gain_percent":float(100*(ea.mean()-eb.mean())/ea.mean()),
                  "fixed_OOF_HUC8_2000_bootstrap95":list(map(float,np.quantile(b,[.025,.975])))}
    return out

def bh(p):
    arr=np.asarray(p,float);n=len(arr)
    ix=np.argsort(arr)
    values=arr[ix]
    q=np.minimum.accumulate((values*n/(np.arange(n)+1))[::-1])[::-1]
    ans=np.zeros(n);ans[ix]=np.clip(q,0,1)
    return ans.tolist()

def run(outdir,permutations=399):
    original,count,sites,qa=read()
    data,matrices=compute(count,sites,qa)
    outputs={"schema":"GESUI_3ERA_EPT_SITETAXON_FOURTHCORNER_DOMINANCE_V1",
             "exploratory_posthoc_after_prior_weighted_trait_results":True,
             "EPA_original_EPT_taxonID_and_measured_counts":True,
             "trait_identity_EPA_name_ID_ORDER_verified":True,
             "n_rarefied_source_sites":2209,
             "not_true_species_richness_or_direct_functional_feeding":True,
             "not_WWTP_specific_or_causal_community_selection":True,
             "multiplicity_primary_BH12":True,
             "within_order_trait_shuffle_B":permutations,
             "results":{},"external":{}, "posthoc_strict_incidence_coverage":{}}
    tested=[]
    for era in CYCLES:
        outputs["results"][era]={}
        for kind in FIELD:
            sel=data[(era,kind)]
            matrix=matrices[(era,kind)]
            reg={k:model(sel,k) for k in ("CWM_count","CWM_incidence","count_minus_incidence")}
            null=taxon_label_permutation(sel,matrix,None,B=permutations,
                                         seed=20261115+100*CYCLES.index(era)+list(FIELD).index(kind))
            outputs["results"][era][kind]={
                "n_sites":len(sel),"HUC8":sel.huc8.nunique(),"states":sel.state.nunique(),
                "original_source_CWM_count_mean":float(sel.CWM_count.mean()),
                "original_incidence_CWM_mean":float(sel.CWM_incidence.mean()),
                "source_taxon_mean_known_per_site":float(sel.known_taxa_n.mean()),
                "known_taxa_n_le1":int(sel.known_taxa_n.le(1).sum()),
                "regression":reg,"trait_label_null":null}
            if kind in PRIMARY:
                for weighted in ("CWM_count","CWM_incidence"):
                    tested.append((era,kind,weighted,reg[weighted]["XEMBED"]["p_nominal"]))
    adj=bh([x[3] for x in tested])
    for (era,kind,weighted,p),q in zip(tested,adj):
        outputs["results"][era][kind]["regression"][weighted]["XEMBED"]["BH12_q"]=q
    # POST-RESULT stress: stricter presence quality and non-singleton support.
    # The original core coefficients above are never reselected or replaced.
    for era in CYCLES:
        outputs["posthoc_strict_incidence_coverage"][era]={}
        for kind in PRIMARY:
            site=data[(era,kind)]
            strict=site[(site.known_taxa_n>=3)&
                        (site.trait_taxon_coverage>=.80)&
                        (site.original_EPT_counted_individuals>=10)].copy()
            diag={"strict_n":len(strict),"strict_HUC8":strict.huc8.nunique(),
                  "strict_states":strict.state.nunique(),
                  "cutoffs_not_pre_registered":True,
                  "min_3_scored_taxa":True,
                  "min_80pct_incidence_taxon_scored":True,
                  "min_10_counted_EPT_individuals":True}
            if len(strict)>=250 and strict.huc8.nunique()>=150 and strict.state.nunique()>=20:
                diag["status"]="exploratory_posthoc_numeric_sensitivity"
                diag["slopes"]={k:model(strict,k) for k in
                                ("CWM_count","CWM_incidence","count_minus_incidence")}
            else:
                diag["status"]="SOURCE_GATE_FAIL_NO_SLOPES"
            outputs["posthoc_strict_incidence_coverage"][era][kind]=diag
    source_huc={c:set(original.loc[original.source_cycle==c,"huc8"]) for c in CYCLES}
    for kind in PRIMARY:
        outputs["external"][kind]={
            metric:external_predict(data,source_huc,kind,metric)
            for metric in ("CWM_count","CWM_incidence")}
    outdir.mkdir(parents=True,exist_ok=True)
    out=outdir/"FOURTHCORNER_EPT_TRAIT_FILTER_RESULTS.json"
    out.write_text(json.dumps(outputs,indent=2,sort_keys=True,allow_nan=False))
    print("FOURTHCORNER_EPT_FILTER_SUMMARY",json.dumps({
        "original_n":2209,"primary":{era:{k:{
           "n":outputs["results"][era][k]["n_sites"],
           "weighted_XEMBED":outputs["results"][era][k]["regression"]["CWM_count"]["XEMBED"],
           "incidence_XEMBED":outputs["results"][era][k]["regression"]["CWM_incidence"]["XEMBED"],
           "difference_XEMBED":outputs["results"][era][k]["regression"]["count_minus_incidence"]["XEMBED"],
           "taxon_label_null":outputs["results"][era][k]["trait_label_null"],
         } for k in PRIMARY} for era in CYCLES},
        "external":outputs["external"],
        "strict_source_sensitivity":outputs["posthoc_strict_incidence_coverage"]},sort_keys=True,allow_nan=False),flush=True)
    return outputs

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--outdir",type=Path,default=TEMP/"gesui_fourthcorner")
    a.add_argument("--permutations",type=int,default=399)
    arg=a.parse_args()
    run(arg.outdir,arg.permutations)
