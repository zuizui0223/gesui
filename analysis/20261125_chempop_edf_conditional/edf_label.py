#!/usr/bin/env python3
"""Targeted posthoc EDF score-specificity test on the identical PR162 data."""
import argparse, json, pathlib, sys, tempfile, zipfile
import numpy as np
import pandas as pd
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261125_chempop_resolution"))
import chemistry_order as prev
wq,old,f,src=prev.wq,prev.old,prev.f,prev.src
def main(work,out):
    prior=prev.receipt(work/"resolution.zip","9c04699d4fab6ba27b3e5ba99e6c17670b0c69bb912584dbe9a697391185edfc")
    out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        d=pathlib.Path(tmp)
        ret=src.retrieve(src.EIDC_URL,d/"original.zip",src.SOURCE_DOWNLOAD_LIMIT)
        scores,si=src.score_reference(d)
        assert si["SHA256"]==prior["score_sha"]==wq.EXPECTED_SCORE_SHA and len(scores)==92
        with zipfile.ZipFile(d/"original.zip") as z:
            names=z.namelist();assert len(names)==len(set(names))==326
            expected=prior["source_35_members"];assert len(expected)==35
            for k,v in expected.items():assert prev.digest(z.read(k))==v,k
            raw,weights,orders,qa=f.construct_eligible(z,names,scores)
            assert len(raw)==7836 and qa["source_gate"]["eligible_site_n"]==707
            chem,_=wq.chemical_summaries(z,names,{(r["region"],r["site"],r["sample_date"]) for r in raw})
        joint,gate=wq.source_join_raw(raw,weights,chem)
        assert gate==prior["cohort"]
        df,_=wq.site_year_matched(joint,scores)
        df0,A,_,families,s,_=f.site_year_response([v[0] for v in joint],[v[1] for v in joint],scores)
        assert df[["site_region_id","year"]].equals(df0[["site_region_id","year"]])
        assert (len(df),df.site_region_id.nunique(),df.region.nunique())==(2624,591,7)
        np.testing.assert_allclose(df.y_incidence,A@s,atol=1e-12,rtol=0)
        ref,oq=old.order_ref(families,s,orders);assert oq==prior["order_source"]
        groups=prev.group_order(families,orders)
        for idx in groups:np.testing.assert_allclose(ref[idx],np.mean(s[idx]),atol=1e-12,rtol=0)
        X=wq.full_design(df,True);X0=f.model_design(df,True)
        weights=f.site_weights(df);sw=np.sqrt(weights)
        y=A@(s-ref)
        beta=np.linalg.lstsq(X.to_numpy(float)*sw[:,None],y*sw,rcond=None)[0]
        edf="EDF_log1p_std"
        target=next(v for v in prior["coefficients"] if v["response"]=="within_order" and v["predictor"]==edf)
        np.testing.assert_allclose(beta[X.columns.get_loc(edf)],target["beta"],atol=1e-10,rtol=0)
        cols=["chem_"+k for k in wq.NUTRIENTS]
        original,original_null=prev.label_test(A,s-ref,groups,X0.to_numpy(float),X[cols].to_numpy(float),weights)
        np.testing.assert_allclose(original["observed_partial_R2"],prior["within_order_label_null"]["observed_partial_R2"],atol=1e-11,rtol=0)
        assert original["p_upper"]==prior["within_order_label_null"]["p_upper"]
        test,null=prev.label_test(A,s-ref,groups,X.drop(columns=[edf]).to_numpy(float),X[[edf]].to_numpy(float),weights)
        pd.DataFrame({"conditional_EDF_null_partial_R2":null,"replayed_chemistry_null_partial_R2":original_null}).to_csv(out/"null_draws.csv",index=False)
        result=dict(status="COMPLETE",source_members=expected,wrapper_SHA=ret["SHA256"],cohort=gate,order_source=oq,prior_EDF_coefficient_reproduced=float(beta[X.columns.get_loc(edf)]),prior_chemistry_null_reproduced=original,conditional_EDF_label_null=test,prior_EDF_fixed_region_gain=prior["region_transfer"]["within_order"]["contrasts"]["EDF_after_nutrients"],score_specificity_nominal_pass=test["p_upper"]<.05,posthoc=True,not_causal=True)
        (out/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
        print("CONDITIONAL EDF LABEL NULL",json.dumps({k:result[k] for k in ("status","prior_EDF_coefficient_reproduced","conditional_EDF_label_null","prior_EDF_fixed_region_gain","score_specificity_nominal_pass")},allow_nan=False))
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--work",required=True,type=pathlib.Path);p.add_argument("--out",required=True,type=pathlib.Path);a=p.parse_args();main(a.work,a.out)
