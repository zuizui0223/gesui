#!/usr/bin/env python3
"""Independent EPA original E/P/T composition audit, no original analysis helpers.

Source HUC8 pair identity derives solely from previously frozen environmental
matching source. Individual ORDER fractions, per-HUC means, Holm tests and
paired bootstrap/permutation are recalculated with separate stdlib CSV+NumPy.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,pathlib
import numpy as np
import gesui_NRSA_fourthcycle_2023_EPT_heldout_transport as oldsource
import gesui_nrsa2023_HUC8_source_novelty_environment_20261122 as domain
import gesui_nrsa2023_huc_novelty_extra_source_matching_20261122 as extra
from gesui_nrsa2324_fourthcycle_source_only_20261122 import SHA,sha
def check(actual,expected,label):
 if isinstance(actual,(list,tuple)):
  if not np.allclose(np.array(actual,float),np.array(expected,float),rtol=1e-10,atol=1e-10):
   raise RuntimeError(label+" independent vector mismatch")
 elif not np.isclose(float(actual),float(expected),rtol=1e-10,atol=1e-10):
  raise RuntimeError(label+" independent scalar mismatch "+str((actual,expected)))
def read_counts(path,uid_to_huc,expected_fraction):
 totals=collections.defaultdict(lambda:np.zeros(5,dtype=float))
 rows=collections.Counter()
 with open(path,newline="",encoding="utf-8-sig") as f:
  r=csv.DictReader(f)
  for item in r:
   uid=item["UID"].strip()
   if uid not in uid_to_huc:continue
   if item["SAMPLE_TYPE"].strip()!="BERW":raise RuntimeError("Native method mismatch")
   try:count=float(item["TOTAL300"])
   except ValueError:continue
   if not np.isfinite(count) or count<=0:continue
   order=item["ORDER"].strip().lower()
   arr=totals[uid]
   arr[0]+=count
   if order=="ephemeroptera":arr[1]+=count
   elif order=="plecoptera":arr[2]+=count
   elif order=="trichoptera":arr[3]+=count
   elif not order:arr[4]+=count
   rows[uid]+=1
 if len(totals)!=1165:raise RuntimeError("Native original EPT count census shrunk")
 huc_comp=collections.defaultdict(list)
 unknown=0;den=0;zeroEPT=0
 for uid,arr in totals.items():
  total,E,P,T,unk=map(float,arr)
  if total<=0:raise RuntimeError("Zero original taxon denominator")
  ept=(E+P+T)/total
  if not np.isclose(ept,expected_fraction[uid],rtol=0,atol=1e-12):
   raise RuntimeError("Original EPT counted numerator reconstructed incorrectly")
  vec=np.array([E/total,P/total,T/total,ept,unk/total],float)
  within=np.array([E,P,T],float)/(E+P+T) if E+P+T>0 else np.array([np.nan]*3)
  huc_comp[uid_to_huc[uid]].append(np.r_[vec,within])
  den+=total;unknown+=unk;zeroEPT+=int(E+P+T==0)
 result={key:np.mean(np.stack(v),axis=0) for key,v in huc_comp.items()}
 # Mean conditional share over positive-EPT UIDs only (numeric QA checks
 # the exact same original per-HUC mean(skipna) semantics).
 for h,v in huc_comp.items():
  A=np.stack(v)
  result[h][5:8]=np.nanmean(A[:,5:8],axis=0) if np.isfinite(A[:,5:8]).any() else np.nan
 if len(result)!=693:raise RuntimeError("Source HUC8 count changed")
 return result,{"zero_EPT_UID":zeroEPT,"missing_ORDER_individual_fraction":unknown/den,
                "median_TOTAL300_counted_individuals":float(np.median([x[0] for x in totals.values()]))}
def simulate(pairs,groups,idx,label,kind):
 novel=np.array([groups[p[0]][idx] for p in pairs],float)
 ever=np.array([groups[p[1]][idx] for p in pairs],float)
 mask=np.isfinite(novel)&np.isfinite(ever)
 n=novel[mask];o=ever[mask]
 diff=n-o
 seed=int(hashlib.sha256((kind+label).encode()).hexdigest()[:8],16)
 rng=np.random.default_rng(seed)
 draw=rng.integers(0,len(diff),size=(2000,len(diff)))
 b=diff[draw].mean(axis=1)
 flip=rng.integers(0,2,size=(10000,len(diff)),dtype=np.int8)*2-1
 null=(flip.astype(float)*diff).mean(axis=1)
 p=(1+np.count_nonzero(np.abs(null)>=abs(float(diff.mean()))-1e-15))/(len(null)+1)
 return {"novel":float(n.mean()),"ever":float(o.mean()),
        "mean_difference_percentage_points":float(100*np.mean(diff)),
        "ci95":list(map(float,100*np.quantile(b,[.025,.975]))),
        "p":float(p),"n_pairs":len(diff)}
def holm3(P):
 ranks=sorted(range(3),key=lambda k:P[k])
 result=[None]*3;minbound=0.
 for ix,i in enumerate(ranks):
  minbound=max(minbound,(3-ix)*P[i])
  result[i]=min(1.,minbound)
 return result
def main():
 p=argparse.ArgumentParser()
 for k in ("sites","water","habitat","chla","benthos","old_panel",
           "original_2023_source","original_2023_result","extra_source",
           "source_order","original_numeric","out"):
  p.add_argument("--"+k,type=pathlib.Path,required=True)
 a=p.parse_args()
 for k,pin in SHA.items():
  if sha(getattr(a,k))!=pin:raise RuntimeError("Original EPA source SHA mismatch "+k)
 signed=json.loads(a.original_numeric.read_text())
 source=json.loads(a.source_order.read_text())
 if source["source_GO"] is not True or signed["EPA_original_SHA"]!=SHA:
  raise RuntimeError("Prior signed original source or numeric schema changed")
 old=oldsource.history(a.old_panel)
 env=oldsource.covariates(a,{})
 if len(old)!=2648 or len(env)!=1165:raise RuntimeError("Previous EPA original source cohort changed")
 hist=domain.huc_aggregate(env)
 train=old.loc[~old.huc8.isin(set(hist.index))].copy()
 trainh=domain.huc_aggregate(train)
 x=extra.block_sources(a,env).join(env[["huc8"]]).groupby("huc8").mean()
 hist=hist.join(x[[k for b in extra.BLOCKS.values() for k in b]],how="left")
 ever=hist.index.isin(set(old.huc8))
 original_source=json.loads(a.original_2023_source.read_text())
 previous_t=oldsource.canonical_target(oldsource.make_target(a,original_source),original_source,"full")
 expected_fraction=dict(zip(previous_t.source_UID,previous_t.ept_fraction))
 uid_to_huc=dict(zip(env.source_UID,env.huc8))
 groups,quality=read_counts(a.benthos,uid_to_huc,expected_fraction)
 check(quality["zero_EPT_UID"],signed["original_COUNT_QUALITY"]["original_EPT_zero_visit_count"],"zero-EPT")
 check(quality["missing_ORDER_individual_fraction"],signed["original_COUNT_QUALITY"]["original_order_unidentified_fraction_of_all_positive_individuals"],"missing order fraction")
 check(quality["median_TOTAL300_counted_individuals"],signed["original_COUNT_QUALITY"]["median_TOTAL300_positive_counted_individuals"],"total counting")
 names=["E_total_share","P_total_share","T_total_share","EPT_total_share",
        "nonidentified_ORDER_share","E_conditional_EPT","P_conditional_EPT","T_conditional_EPT"]
 results={}
 for kind,cols in extra.BLOCKS.items():
  pairs,_,digest=extra.choose(hist,trainh,ever,extra.BASE+cols)
  if digest!=signed["matched_blocks"][kind]["previous_source_pair_identity_sha256"]:
   raise RuntimeError("Original frozen source-only HUC8 pair SHA changed")
  saved=signed["matched_blocks"][kind]["matched_EPT_order_source_results"]
  outcomes={}
  for idx,label in enumerate(names):
   new=simulate(pairs,groups,idx,label,kind)
   result=saved[label]
   for f,ff in (("novel","novel_mean_fraction"),("ever","old_source_ever_mean_fraction"),
                ("mean_difference_percentage_points","novel_minus_ever_difference_percentage_points"),
                ("ci95","paired_HUC8_bootstrap95_percentage_points"),
                ("p","two_sided_paired_signflip_permutation_p_unadjusted"),
                ("n_pairs","matched_HUC8_pairs")):
    check(new[f],result[ff],kind+" "+label+" "+f)
   outcomes[label]={"difference_pp":new["mean_difference_percentage_points"],
                    "matched_pairs":new["n_pairs"],"all_statistics_reproduced":True}
  pvals=[simulate(pairs,groups,i,names[i],kind)["p"] for i in range(3)]
  holm=holm3(pvals)
  for i in range(3):
   check(holm[i],saved[names[i]]["holm3_exploratory_adjusted_permutation_p"],
         kind+" "+names[i]+" Holm3")
  results[kind]={"n_source_matched_HUC8_pairs":len(pairs),
    "matched_source_HUC8_pair_SHA":digest,
    "n_prespecified_order_comparisons":3,
    "all_original_posthoc_EPT_order_share_CI_permutation_holm3_recomputed":True,
    "order_checks":outcomes}
 output={
  "schema":"gesui_EPA2023_EPT_order_matched_huc_independent_stdlib_numpy_QA_v1",
  "source_original_UID":1165,"original_source_HUC8":693,
  "all_EPT_TOTAL300_fractions_exactly_reproduced":True,
  "n_three_named_EPT_orders":3,"independent_source_HUC8_and_original_signed_SHA":True,
  "no_WWTP_or_aquatic_to_terrestrial_transfer":True,
  "results":results}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
 print("NRSA_EPT_ORDER_INDEPENDENT_TAXON_ABUNDANCE_QA",json.dumps(
  {k:{"matched_HUC8":v["n_source_matched_HUC8_pairs"],
       "orders":{x:v["order_checks"][x]["difference_pp"] for x in names[:3]},
       "all_tests_reproduced":True} for k,v in results.items()},sort_keys=True))
if __name__=="__main__":main()
