#!/usr/bin/env python3
"""Same-source, same-cohort chemistry-associated BMWP taxonomic resolution audit."""
import argparse, collections, hashlib, json, pathlib, sys, tempfile, zipfile
import numpy as np
import pandas as pd
import statsmodels.api as sm
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261120_chempop_bmwp"))
import wq_adjusted_edf as wq
import adversarial_order_actual as old
f=wq.original
src=wq.source
SEED=20261125
B=4999
def digest(b):return hashlib.sha256(b).hexdigest()
def receipt(p,expected):
    assert digest(p.read_bytes())==expected,p
    with zipfile.ZipFile(p) as z:
        names=[n for n in z.namelist() if n.endswith(".json")]
        assert len(names)==1
        return json.loads(z.read(names[0]))
def group_order(families,orders):
    out=collections.defaultdict(list)
    for i,k in enumerate(families):
        vals=orders.get(k,{})
        key=next(iter(vals)) if len(vals)==1 else "__singleton_"+k
        out[key].append(i)
    return [np.asarray(x,int) for x in out.values()]
def label_test(A,within,groups,X0,X1,w):
    sw=np.sqrt(w)
    Q0=np.linalg.qr(X0*sw[:,None],mode="reduced")[0]
    M=A*sw[:,None];M-=Q0@(Q0.T@M)
    Z=X1*sw[:,None];Z-=Q0@(Q0.T@Z)
    assert np.linalg.matrix_rank(Z)==Z.shape[1]
    Qz=np.linalg.qr(Z,mode="reduced")[0]
    R=Qz.T@M
    gram=M.T@M;proj=R.T@R
    def stat(S):
        num=np.einsum("bi,ij,bj->b",S,proj,S)
        den=np.einsum("bi,ij,bj->b",S,gram,S)
        assert np.all(den>1e-14)
        return num/den
    observed=stat(within[None,:])[0]
    Y=(A@within)*sw
    design=np.column_stack((X0,X1))*sw[:,None]
    r0=Y-X0*sw[:,None]@np.linalg.lstsq(X0*sw[:,None],Y,rcond=None)[0]
    r1=Y-design@np.linalg.lstsq(design,Y,rcond=None)[0]
    direct=1-(r1@r1)/(r0@r0)
    np.testing.assert_allclose(observed,direct,atol=1e-11,rtol=0)
    rng=np.random.default_rng(SEED);S=np.tile(within,(B,1))
    for k in range(B):
        for idx in groups:S[k,idx]=within[rng.permutation(idx)]
    null=stat(S)
    assert np.isfinite(null).all()
    return dict(observed_partial_R2=float(observed),null_quantiles=np.quantile(null,[.025,.5,.975]).tolist(),p_upper=float((1+np.count_nonzero(null>=observed-1e-14))/(B+1)),permutations=B,non_singleton_groups=sum(len(x)>1 for x in groups),exchangeable_families=sum(len(x) for x in groups if len(x)>1),independent_direct_partial_R2=float(direct)),null
def main(work,out):
    prior=receipt(work/"chemistry.zip","18c764318c2143a40565b8d58be350d91c7b8e3b371c9705038401db8730ff43")
    order=receipt(work/"order.zip","785d977b8eb81a80d51f13b870d5c8ce0a399c20debc02af5283a5824922a03a")
    out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        d=pathlib.Path(tmp)
        retrieved=src.retrieve(src.EIDC_URL,d/"original.zip",src.SOURCE_DOWNLOAD_LIMIT)
        scores,si=src.score_reference(d)
        assert si["SHA256"]==wq.EXPECTED_SCORE_SHA and len(scores)==92
        with zipfile.ZipFile(d/"original.zip") as z:
            names=z.namelist();assert len(names)==len(set(names))==326
            expected=dict(order["original_14_critical_member_SHA256"])
            for v in prior["original_wq_source_files"].values():expected[v["original_file"]]=v["source_member_SHA256"]
            assert len(expected)==35
            for k,v in expected.items():assert digest(z.read(k))==v,k
            raw,weights,orders,qa=f.construct_eligible(z,names,scores)
            assert len(raw)==7836 and qa["source_gate"]["eligible_site_n"]==707
            chem,cqa=wq.chemical_summaries(z,names,{(r["region"],r["site"],r["sample_date"]) for r in raw})
        joint,gate=wq.source_join_raw(raw,weights,chem)
        assert gate==prior["matched_source_gate"]
        df,joined=wq.site_year_matched(joint,scores)
        assert joined==prior["matched_siteyear_source"]
        df0,A,_,families,s,unused=f.site_year_response([v[0] for v in joint],[v[1] for v in joint],scores)
        assert df[["site_region_id","year"]].equals(df0[["site_region_id","year"]])
        assert (len(df),df.site_region_id.nunique(),df.region.nunique())==(2624,591,7)
        np.testing.assert_allclose(df.y_incidence,A@s,atol=1e-12,rtol=0)
        ref,oq=old.order_ref(families,s,orders);groups=group_order(families,orders)
        for idx in groups:np.testing.assert_allclose(ref[idx],np.mean(s[idx]),atol=1e-12,rtol=0)
        Y={"total":A@s,"between_order":A@ref,"within_order":A@(s-ref)}
        np.testing.assert_allclose(Y["total"],Y["between_order"]+Y["within_order"],atol=1e-12,rtol=0)
        X=wq.full_design(df,True);X0=f.model_design(df,True)
        vif=wq.VIF_for_original_controls(X)
        assert max(vif.values())<=wq.SOURCE_GATES["max_full_design_VIF_EDF_or_one_nutrient"]
        for k,v in prior["joint_covariate_source_vif"].items():np.testing.assert_allclose(vif[k],v,atol=1e-8,rtol=0)
        cols=["chem_"+k for k in wq.NUTRIENTS]
        assert set(X.columns)==set(X0.columns)|set(cols)
        weights=f.site_weights(df);xa=X.to_numpy(float);beta={};results={};regressions=[]
        for name,y in Y.items():
            d1=df.copy();d1["y_incidence"]=y
            fit=sm.WLS(y,xa,weights=weights).fit(cov_type="cluster",cov_kwds={"groups":df.site_region_id.astype(str).to_numpy()},use_t=False)
            independent=np.linalg.lstsq(xa*np.sqrt(weights)[:,None],y*np.sqrt(weights),rcond=None)[0]
            np.testing.assert_allclose(fit.params,independent,atol=1e-10,rtol=0)
            beta[name]=independent
            for c in ["EDF_log1p_std"]+cols:
                j=X.columns.get_loc(c);lo,hi=fit.conf_int()[j]
                regressions.append(dict(response=name,predictor=c,beta=float(fit.params[j]),lo=float(lo),hi=float(hi)))
            results[name]=wq.region_transfer_four_models(d1)
        np.testing.assert_allclose(beta["total"],beta["between_order"]+beta["within_order"],atol=1e-10,rtol=0)
        for k,v in prior["EDF_after_all_three_chemistry_WLS"]["coefficient"].items():
            np.testing.assert_allclose(beta["total"][X.columns.get_loc(k)],v["beta"],atol=1e-10,rtol=0)
        for k,v in prior["regional_transfer"]["contrasts"].items():
            np.testing.assert_allclose(results["total"]["contrasts"][k]["MSE_additional_gain_percent"],v["MSE_additional_gain_percent"],atol=1e-7,rtol=0)
        test,null=label_test(A,s-ref,groups,X0.to_numpy(float),X[cols].to_numpy(float),weights)
        cv=[]
        for name,t in results.items():
            assert len(t["valid_regions"])==7
            for key,v in t["contrasts"].items():
                cv.append(dict(response=name,contrast=key,gain=v["MSE_additional_gain_percent"],lo=v["fixed_region_loss_2000_bootstrap95"][0],hi=v["fixed_region_loss_2000_bootstrap95"][1],improved_regions=v["improved_regions"]))
        pd.DataFrame(regressions).to_csv(out/"coefficients.csv",index=False)
        pd.DataFrame(cv).to_csv(out/"regional_gains.csv",index=False)
        pd.DataFrame({"null_partial_R2":null}).to_csv(out/"within_order_label_null.csv",index=False)
        passed=bool(test["p_upper"]<.05 and results["within_order"]["contrasts"]["NUTRIENTS_after_EDF"]["fixed_region_loss_2000_bootstrap95"][0]>0)
        info=dict(status="COMPLETE",source_35_members=expected,source_wrapper=retrieved["SHA256"],score_sha=si["SHA256"],cohort=gate,order_source=oq,coefficients=regressions,region_transfer=results,within_order_label_null=test,primary_joint_criterion_pass=passed,posthoc=True,not_causal=True)
        (out/"results.json").write_text(json.dumps(info,indent=2,allow_nan=False)+"\n")
        print("CHEMISTRY ORDER AUDIT",json.dumps(dict(cohort=gate,order_source=oq,within_order_label_null=test,primary_joint_criterion_pass=passed),allow_nan=False))
        print(pd.DataFrame(regressions).to_string(index=False))
        print(pd.DataFrame(cv).to_string(index=False))
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--work",required=True,type=pathlib.Path);p.add_argument("--out",required=True,type=pathlib.Path);a=p.parse_args();main(a.work,a.out)
