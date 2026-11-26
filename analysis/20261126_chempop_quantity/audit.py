#!/usr/bin/env python3
"""Fixed ChemPop matched-chemistry abundance/richness/concentration audit.
Actual recorded counts only; sample totals are not river density or biomass.
Source replay before response models. No raw source rows exported."""
import argparse, collections, hashlib, inspect, json, pathlib, sys, tempfile, zipfile
import numpy as np, pandas as pd
from scipy.stats import hypergeom
import statsmodels.api as sm
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.linear_model import Ridge
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261125_chempop_resolution"))
import chemistry_order as prior
f,src,wq=prior.f,prior.src,prior.wq
METHOD="Actual abundance of taxa in the sample"
CHEM=["chem_"+k for k in wq.NUTRIENTS]
BASE=["lc_urban","lc_arable","width","altitude","siltclay","year","month_sin","month_cos","edf_log1p"]
SEED=20261126
def augment_parser():
    # Exact anchored changes to the pinned prior parser: retain all native
    # reported integer counts before the original score/family filters.
    # Original keys, order, eligibility and score computations are unchanged.
    s=inspect.getsource(f.construct_eligible)
    replacements=[
      ('"family_order":{}, "counts":collections.defaultdict(float),',
       '"family_order":{}, "counts":collections.defaultdict(float), "all_counts":collections.Counter(), "unknown_count":0, "invalid_count_rows":0, "native_count_rows":0,'),
      ('            if not fam:continue\n            v["families"].add(fam)',
       '''            if method == "Actual abundance of taxa in the sample":
                v["native_count_rows"]+=1
                n0=source.finite_numeric(row["BIO_ABUNDANCE"])
                if n0 is None or n0<0 or n0!=int(n0):
                    v["invalid_count_rows"]+=1
                elif fam:
                    v["all_counts"][fam]+=int(n0)
                else:
                    v["unknown_count"]+=int(n0)
            if not fam:continue
            v["families"].add(fam)'''),
      ('        region_samples[region]+=1',
       '''        sample_family_weights[-1].update(all_counts=dict(v["all_counts"]),unknown_count=v["unknown_count"],invalid_count_rows=v["invalid_count_rows"],native_count_rows=v["native_count_rows"])
        region_samples[region]+=1''')]
    for a,b in replacements:
        assert s.count(a)==1,("prior parser anchor drift",a)
        s=s.replace(a,b)
    ns=dict(f.__dict__);exec(compile(s,"<source-locked augmented parser>","exec"),ns)
    return ns["construct_eligible"]
def rare(n,N,m=50):
    n=np.asarray(n,dtype=int);assert N>=m and np.all(n>0) and n.sum()<=N
    a=hypergeom.sf(0,N,n,m)
    # Direct without-replacement absence recurrence, independent implementation.
    miss=np.ones(len(n))
    for i in range(m):miss*=np.maximum(0,1-n/(N-i))
    np.testing.assert_allclose(a,1-miss,atol=5e-12,rtol=0)
    return float(a.sum())
def metrics(w):
    counts=np.array([n for n in w["all_counts"].values() if n>0],int)
    unknown=int(w["unknown_count"]);N=int(counts.sum())+unknown
    if w["invalid_count_rows"] or not len(counts) or N<=0:return None,"invalid_or_empty"
    if counts.sum()/N<.95:return None,"family_count_coverage_below95"
    if N<50:return None,"fewer_than50_reported_individuals"
    p=counts/counts.sum()
    return dict(total=float(N),richness=float(len(counts)),rare50=rare(counts,N),simpson=float(p@p),top_share=float(p.max()),known_count_fraction=float(counts.sum()/N)),None
def ci(loss):
    a=np.asarray(loss);g=100*(1-a[:,1].mean()/a[:,0].mean())
    rng=np.random.default_rng(SEED);b=a[rng.integers(0,len(a),(5000,len(a)))].mean(1)
    return float(g),*map(float,np.quantile(100*(1-b[:,1]/b[:,0]),[.025,.975]))
def predict(tr,te,y,cols):
    prep=ColumnTransformer([("numeric",StandardScaler(),cols),("region",OneHotEncoder(handle_unknown="ignore",sparse_output=False),["region"])],sparse_threshold=0)
    X=prep.fit_transform(tr);Z=prep.transform(te);w=1/tr.groupby("site_region_id").site_region_id.transform("size").to_numpy(float);w*=len(tr)/w.sum()
    p=Ridge(alpha=10,solver="svd").fit(X,tr[y],sample_weight=w).predict(Z)
    # Independent augmented least squares; same scaling, weights, penalty.
    a=np.column_stack((np.ones(len(X)),X));b=np.column_stack((np.ones(len(Z)),Z))
    pen=np.diag([0.]+[np.sqrt(10.)]*X.shape[1]);sw=np.sqrt(w)
    q=np.linalg.lstsq(np.vstack((a*sw[:,None],pen)),np.r_[tr[y].to_numpy()*sw,np.zeros(pen.shape[0])],rcond=None)[0]
    np.testing.assert_allclose(p,b@q,atol=1e-8,rtol=0)
    return p
def regression(d,y):
    # Directional index fixed as equal mean of three standardized log1p
    # concentrations. Not a dose, not total N or P.
    x=d[BASE+CHEM].astype(float);sd=x.std(ddof=0);assert (sd>0).all()
    x=(x-x.mean())/sd
    x["nutrient_index"]=x[CHEM].mean(axis=1);x["nutrient_index"]/=x["nutrient_index"].std(ddof=0)
    reg=pd.get_dummies(d.region,drop_first=True,dtype=float).set_axis(d.index)
    X=sm.add_constant(pd.concat((x[BASE+["nutrient_index"]],reg),axis=1)).to_numpy(float)
    assert np.linalg.matrix_rank(X)==X.shape[1]
    w=1/d.groupby("site_region_id").site_region_id.transform("size").to_numpy(float)
    fit=sm.WLS(d[y],X,weights=w).fit(cov_type="cluster",cov_kwds={"groups":d.site_region_id},use_t=False)
    j=1+len(BASE);q=np.linalg.lstsq(X*np.sqrt(w)[:,None],d[y].to_numpy()*np.sqrt(w),rcond=None)[0]
    np.testing.assert_allclose(q,np.asarray(fit.params),atol=1e-9,rtol=0)
    lo,hi=np.asarray(fit.conf_int())[j]
    return dict(beta=float(q[j]),lo=float(lo),hi=float(hi),p=float(fit.pvalues.iloc[j]))
def main(work,out):
    out.mkdir(parents=True,exist_ok=True)
    status={"status":"SOURCE_STOP","posthoc":True,"no_sewage_causation":True}
    try:
        prev=prior.receipt(work/"chemistry.zip","18c764318c2143a40565b8d58be350d91c7b8e3b371c9705038401db8730ff43")
        order=prior.receipt(work/"order.zip","785d977b8eb81a80d51f13b870d5c8ce0a399c20debc02af5283a5824922a03a")
        with tempfile.TemporaryDirectory() as tmp:
            p=pathlib.Path(tmp)
            retrieved=src.retrieve(src.EIDC_URL,p/"original.zip",src.SOURCE_DOWNLOAD_LIMIT)
            scores,srec=src.score_reference(p)
            assert srec["SHA256"]==wq.EXPECTED_SCORE_SHA and len(scores)==92
            with zipfile.ZipFile(p/"original.zip") as z:
                names=z.namelist();assert len(names)==len(set(names))==326
                expected=dict(order["original_14_critical_member_SHA256"])
                for v in prev["original_wq_source_files"].values():expected[v["original_file"]]=v["source_member_SHA256"]
                assert len(expected)==35
                for k,v in expected.items():assert prior.digest(z.read(k))==v,k
                raw,ww,orders,qa=augment_parser()(z,names,scores)
                assert len(raw)==7836 and qa["source_gate"]["eligible_site_n"]==707
                chem,_=wq.chemical_summaries(z,names,{(r["region"],r["site"],r["sample_date"]) for r in raw})
            joint,gate=wq.source_join_raw(raw,ww,chem)
            assert gate==prev["matched_source_gate"]
            old,joined=wq.site_year_matched(joint,scores)
            assert joined==prev["matched_siteyear_source"]
            assert (len(joint),len(old),old.site_region_id.nunique(),old.region.nunique())==(5000,2624,591,7)
            # Verify original incidence calculation and exact scored-count reconstruction.
            actual=[(r,w,c) for r,w,c in joint if r["method"]==METHOD]
            status.update(source_35_members=expected,source_wrapper=retrieved["SHA256"],source_cohort=gate,method_counts=dict(collections.Counter(r["method"] for r,w,c in joint)),actual_matched_samples=len(actual),actual_matched_sites=len({r["site_region_id"] for r,w,c in actual}))
            rows=[];rej=collections.Counter()
            for r,w,c in actual:
                v,reason=metrics(w)
                if reason:rej[reason]+=1;continue
                # Same original scored counts, independent extra all-family parser.
                for k,n in w["numeric_count"].items():assert w["all_counts"][k]==n
                rows.append({**r,**v,**{k:float(c[k[5:]]["primary"]) for k in CHEM}})
            status["count_quality_rejections"]=dict(rej)
            d=pd.DataFrame(rows)
            if d.empty:status["status"]="STOP_NO_VALID_ACTUAL_COUNTS";return
            rc=d.groupby("region").site_region_id.nunique().to_dict()
            status["count_support"]=dict(samples=len(d),sites=int(d.site_region_id.nunique()),siteyears=int(d[["site_region_id","year"]].drop_duplicates().shape[0]),regions=rc,median_total=float(d.total.median()),median_family_count_coverage=float(d.known_count_fraction.median()))
            if len(d)<300 or d.site_region_id.nunique()<100 or status["count_support"]["siteyears"]<200 or sum(v>=10 for v in rc.values())<4:
                status["status"]="STOP_COUNT_SUPPORT";return
            # Every eligible sample retained; average sample outcomes by site-year.
            # Chemistry date duplicates collapse before site-year means.
            d["month_sin"]=np.sin(2*np.pi*(d.month-1)/12);d["month_cos"]=np.cos(2*np.pi*(d.month-1)/12);d["edf_log1p"]=np.log1p(d.edf_pct)
            d["log_total"]=np.log1p(d.total);d["log_rare50"]=np.log1p(d.rare50);d["log_richness"]=np.log1p(d.richness)
            Y=["log_total","log_rare50","simpson","top_share","log_richness"]
            key=["site_region_id","site","region","year"]
            bc=[k for k in BASE if k!="year"]
            ag=d.groupby(key,as_index=False)[bc+Y+["total","rare50","richness"]].mean()
            cd=d.drop_duplicates(key+["sample_date"]).groupby(key,as_index=False)[CHEM].mean()
            cd[CHEM]=np.log1p(cd[CHEM]);ag=ag.merge(cd,on=key,validate="one_to_one")
            assert len(ag)==status["count_support"]["siteyears"]
            status["summary_means"]={k:float(ag.groupby("site_region_id")[k].mean().mean()) for k in ["total","rare50","richness","simpson","top_share"]}
            res=[];loss=[]
            for y in Y:
                rr=regression(ag,y);ll=[]
                for reg in sorted(ag.region.unique()):
                    te=ag[ag.region==reg];tr=ag[(ag.region!=reg)&(~ag.site.isin(te.site))]
                    if te.site.nunique()<10 or tr.site.nunique()<70:continue
                    a=predict(tr,te,y,BASE);b=predict(tr,te,y,BASE+CHEM)
                    v=te.assign(e0=(te[y]-a)**2,e1=(te[y]-b)**2).groupby("site_region_id")[["e0","e1"]].mean().mean()
                    ll.append([float(v.e0),float(v.e1)]);loss.append(dict(response=y,region=reg,sites=te.site.nunique(),base=float(v.e0),full=float(v.e1)))
                if len(ll)>=4:
                    gain,lo,hi=ci(ll);rr.update(gain=gain,gain_lo=lo,gain_hi=hi,regions=len(ll),improved_regions=sum(a>b for a,b in ll))
                else:rr["prediction_status"]="STOP_REGION_SUPPORT"
                res.append(dict(response=y,**rr))
            from statsmodels.stats.multitest import multipletests
            for r,p in zip(res[:3],multipletests([r["p"] for r in res[:3]],method="holm")[1]):r["holm_primary_p"]=float(p)
            pd.DataFrame(res).to_csv(out/"results.csv",index=False);pd.DataFrame(loss).to_csv(out/"regional_losses.csv",index=False)
            status.update(status="COMPLETE",results=res)
    finally:
        (out/"receipt.json").write_text(json.dumps(status,indent=2,allow_nan=False)+"\n")
        print("QUANTITY DIVERSITY "+json.dumps(status,allow_nan=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--work",type=pathlib.Path,required=True);p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args();main(a.work,a.out)
