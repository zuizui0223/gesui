#!/usr/bin/env python3
"""Source-frozen EPA NRSA2023 E/P/T count shares in historically sampled vs novel HUC8.

Primary diagnostic: environmental-matched paired differences in original
TOTAL300-based Ephemeroptera, Plecoptera and Trichoptera fractions.
Source pair membership was independently frozen before this numerical test.
All tests exploratory/posthoc; no sewage, causal sediment, or insect emergence.
"""
from __future__ import annotations
import argparse,hashlib,json,pathlib
import numpy as np,pandas as pd
import gesui_NRSA_fourthcycle_2023_EPT_heldout_transport as past
import gesui_nrsa2023_HUC8_source_novelty_environment_20261122 as domain
import gesui_nrsa2023_huc_novelty_extra_source_matching_20261122 as extra
import gesui_nrsa2023_novel_HUC8_source_matched_EPT_prediction_20261122 as oldoutcome
from gesui_nrsa2324_fourthcycle_source_only_20261122 import SHA,sha

ORDER={"E":"ephemeroptera","P":"plecoptera","T":"trichoptera"}
PAIRS_SHA={
 "grain":"8c78c9060e0d512666062760f0467baaa54f6cfabeeaf76867ad609e3f1ea8e2",
 "hydrology_landuse":"4d21361aa86707039d4afd62cfa11c70f629c50ea01e1a0663e4c00930f254b2"}
BOOT=2000;FLIPS=10000
def get_original_ep_orders(args,rawtarget):
 columns=["UID","SITE_ID","VISIT_NO","DATE_COL","SAMPLE_TYPE","ORDER","TOTAL300"]
 b=pd.read_csv(args.benthos,usecols=columns,dtype=str,keep_default_na=False,low_memory=False)
 alluid=set(rawtarget.source_UID)
 b=b[b.UID.isin(alluid)].copy()
 if set(b.UID)!=alluid:raise RuntimeError("Original 2023 1165 selected UID taxa coverage changed")
 b["count"]=pd.to_numeric(b.TOTAL300,errors="coerce")
 b=b[b.SAMPLE_TYPE.eq("BERW")&b["count"].gt(0)].copy()
 b["order"]=b.ORDER.str.strip().str.lower()
 den=b.groupby("UID",sort=True)["count"].sum().rename("denominator")
 order_counts={}
 for label,order in ORDER.items():
  order_counts[label]=b.loc[b.order.eq(order)].groupby("UID")["count"].sum().reindex(den.index,fill_value=0)
 g=pd.DataFrame({"denominator":den,**{k:v for k,v in order_counts.items()}})
 if len(g)!=1165 or not(g.denominator>0).all():raise RuntimeError("Selected 2023 UID count denominator changed")
 g["ept_total"]=g[list(ORDER)].sum(axis=1)
 g["EPT_total_share"]=g.ept_total/g.denominator
 g["nonidentified_ORDER_count"]=b.loc[b.order.eq("")].groupby("UID")["count"].sum().reindex(g.index,fill_value=0)
 for label in ORDER:
  g[label+"_total_share"]=g[label]/g.denominator
  g[label+"_conditional_EPT"]=np.where(g.ept_total>0,g[label]/g.ept_total,np.nan)
 g["nonidentified_ORDER_share"]=g.nonidentified_ORDER_count/g.denominator
 r=rawtarget.set_index("source_UID")["ept_fraction"].reindex(g.index)
 if r.isna().any():raise RuntimeError("Original EPA target EPT source UID misalignment")
 if not np.allclose(g.EPT_total_share.to_numpy(float),r.to_numpy(float),rtol=0,atol=1e-12):
  raise RuntimeError("Source original counted EPT share cannot be reconciled with prior 4-era EPT fraction")
 if not np.allclose(g[[k+"_total_share" for k in ORDER]].sum(axis=1).to_numpy(),
                    r.to_numpy(float),rtol=0,atol=1e-12):
  raise RuntimeError("The E/P/T three shares did not exactly reconstruct original EPT fraction")
 visit_quality={
  "original_2023_UID_positive_total":len(g),
  "original_EPT_zero_visit_count":int((g.ept_total==0).sum()),
  "median_TOTAL300_positive_counted_individuals":float(g.denominator.median()),
  "q10_90_TOTAL300_positive_counted_individuals":[float(q) for q in g.denominator.quantile([.1,.9])],
  "original_order_unidentified_positive_individuals":float(g.nonidentified_ORDER_count.sum()),
  "original_order_unidentified_fraction_of_all_positive_individuals":float(g.nonidentified_ORDER_count.sum()/g.denominator.sum()),
  "all_2023_UID_EPT_fractions_exactly_reconciled_to_prior_model":True}
 cols=["EPT_total_share","nonidentified_ORDER_share"]+[k+"_total_share" for k in ORDER]+[k+"_conditional_EPT" for k in ORDER]
 c=g[cols].join(rawtarget.set_index("source_UID")[["huc8","eco9","site_id"]],how="left")
 if c.huc8.isna().any():raise RuntimeError("Missing original target HUC8 count join")
 h=c.groupby("huc8")[cols].mean()
 h["EPT_present_UID"]=c.assign(has_ept=g.ept_total.gt(0)).groupby("huc8").has_ept.sum()
 h["uid_count"]=c.groupby("huc8").size()
 if len(h)!=693:raise RuntimeError("Original target 693 HUC8 order composition changed")
 visit_quality["HUC8_with_no_positive_EPT_UID"]=int((h.EPT_present_UID==0).sum())
 visit_quality["HUC8_with_some_EPT_uid"]=int((h.EPT_present_UID>0).sum())
 return h,visit_quality
def holm(raw):
 indexed=sorted(enumerate(raw),key=lambda x:x[1])
 adj=[None]*len(raw)
 high=0
 for i,(idx,p) in enumerate(indexed):
  high=max(high,(len(raw)-i)*p)
  adj[idx]=float(min(1.,high))
 return adj
def paired_aggregate(huc,pairs,kind):
 outcomes=["E_total_share","P_total_share","T_total_share"]
 group={}
 for endpoint in outcomes+["EPT_total_share","E_conditional_EPT","P_conditional_EPT","T_conditional_EPT","nonidentified_ORDER_share"]:
  novel=np.asarray([huc.at[p[0],endpoint] for p in pairs],float)
  ever=np.asarray([huc.at[p[1],endpoint] for p in pairs],float)
  valid=np.isfinite(novel)&np.isfinite(ever)
  if endpoint in outcomes+["EPT_total_share","nonidentified_ORDER_share"] and not valid.all():
   raise RuntimeError("Missing original fully defined UID fraction in paired HUC8 "+endpoint)
  nn=novel[valid];ee=ever[valid];diff=nn-ee
  if len(nn)<40:group[endpoint]={"status":"SOURCE_STOP_INADEQUATE_POSITIVE_EPT_PAIRS","pairs":len(nn)};continue
  seed=int(hashlib.sha256((kind+endpoint).encode()).hexdigest()[:8],16)
  rng=np.random.default_rng(seed)
  boots=rng.integers(0,len(diff),size=(BOOT,len(diff)))
  ci=np.quantile(diff[boots].mean(axis=1),[.025,.975])
  flips=rng.integers(0,2,size=(FLIPS,len(diff)),dtype=np.uint8)*2-1
  null=(flips.astype(float)*diff).mean(axis=1)
  p=(1+np.count_nonzero(abs(null)>=abs(diff.mean())-1e-15))/(FLIPS+1)
  group[endpoint]={
   "status":"POSTHOC_MATCHED_HUC8_COMPLETE",
   "matched_HUC8_pairs":len(diff),
   "novel_mean_fraction":float(nn.mean()),"old_source_ever_mean_fraction":float(ee.mean()),
   "novel_minus_ever_difference_fraction":float(diff.mean()),
   "novel_minus_ever_difference_percentage_points":float(100*diff.mean()),
   "paired_HUC8_bootstrap95_percentage_points":[float(100*c) for c in ci],
   "two_sided_paired_signflip_permutation_p_unadjusted":float(p),
   "no_randomized_or_causal_interpretation":True}
 vals=[group[k]["two_sided_paired_signflip_permutation_p_unadjusted"] for k in outcomes]
 adjusted=holm(vals)
 for j,k in enumerate(outcomes):group[k]["holm3_exploratory_adjusted_permutation_p"]=adjusted[j]
 return group
def assemble(args):
 source=json.loads(args.source_order.read_text())
 if source["schema"]!="gesui_NRSA2023_EPT_taxonomic_ORDER_availability_SOURCE_ONLY_v1" or not source["source_GO"]:
  raise RuntimeError("Original independently signed pre-outcome EPT ORDER source gate STOP")
 if source["original_EPA_source_SHA256"]!=SHA:raise RuntimeError("Source EPA original SHA changed")
 ext_source=json.loads(args.extra_source.read_text())
 if ext_source["schema"]!="NRSA2023_original_four_cycle_HUC_novelty_extra_source_domains_v1":
  raise RuntimeError("Frozen original source-only extra HUC8 matching schema invalid")
 for k in SHA:
  if sha(getattr(args,k))!=SHA[k]:raise RuntimeError("Original EPA source file SHA changed "+k)
 src2023=json.loads(args.original_2023_source.read_text())
 old=past.history(args.old_panel)
 env=past.covariates(args,src2023)
 if (len(env),env.huc8.nunique(),len(old),old.huc8.nunique())!=(1165,693,2648,1030):
  raise RuntimeError("EPA original four-cycle source population changed")
 ht=domain.huc_aggregate(env)
 train=old.loc[~old.huc8.isin(set(ht.index))].copy()
 trainh=domain.huc_aggregate(train)
 if len(train)!=1185 or trainh.shape[0]!=551 or set(train.huc8)&set(ht.index):
  raise RuntimeError("Original no-target-HUC8 train source changed")
 extra_h=extra.block_sources(args,env).join(env[["huc8"]]).groupby("huc8").mean()
 ht=ht.join(extra_h[[x for q in extra.BLOCKS.values() for x in q]],how="left")
 ever=ht.index.isin(set(old.huc8))
 if (ever.sum(),(~ever).sum())!=(479,214):raise RuntimeError("Original source-novelty groups changed")
 frozen={}
 for kind,cols in extra.BLOCKS.items():
  signed=ext_source["source_only_blocks"][kind]
  pairs,regions,digest=extra.choose(ht,trainh,ever,extra.BASE+cols)
  if not signed["source_go"] or digest!=signed["source_match_pair_membership_SHA256"] or digest!=PAIRS_SHA[kind] or len(pairs)!=signed["source_only_match_pairs"]:
   raise RuntimeError("Original EPA pre-response source-only HUC8 match SHA changed "+kind)
  frozen[kind]=(pairs,digest)
 # Only after both source pair SHA exact matches, open biological order counts.
 original_t=past.canonical_target(past.make_target(args,src2023),src2023,"full")
 if len(original_t)!=1165 or original_t.huc8.nunique()!=693:
  raise RuntimeError("Original 2023 respondent target changed")
 h,quality=get_original_ep_orders(args,original_t)
 # The prior original EPT trained model must be numerically identical.
 original_ref=json.loads(args.original_2023_result.read_text())
 preds=past.all_predictions(train,original_t)
 oldoutcome.fixed_prediction_reference(preds,original_t,original_ref["results"]["full"])
 reports={}
 for kind,(pairs,digest) in frozen.items():
  reports[kind]={"previous_source_pair_identity_sha256":digest,
    "source_preoutcome_matched_HUC8_pairs":len(pairs),
    "matched_EPT_order_source_results":paired_aggregate(h,pairs,kind)}
 report={"schema":"gesui_EPA_NRSA2023_4cycle_EPT_order_source_novelty_composition_v1",
  "EPA_original_SHA":SHA,"original_source_2023_wadeable_UID":1165,
  "original_2023_huc8":693,"prior_3era_UID":2648,
  "past_years_train_HUC8":551,"historical_old_source_ever_HUC8":479,
  "historical_old_source_never_HUC8":214,
  "original_all_EPT_counted_share_exactly_reproduced":True,
  "original_3era_EPT_prediction_full_target_exactly_reproduced":True,
  "E_P_T_ORDERS_predeclared_before_composition_result":ORDER,
  "original_COUNT_QUALITY":quality,
  "matched_blocks":reports,
  "all_tests_posthoc_and_matched_not_randomized":True,
  "not_independent_ecological_replication_or_sewage_effect":True}
 return report
def main():
 ap=argparse.ArgumentParser()
 for k in ("sites","water","habitat","chla","benthos","old_panel",
           "original_2023_source","original_2023_result","extra_source","source_order","out"):
  ap.add_argument("--"+k,type=pathlib.Path,required=True)
 args=ap.parse_args();res=assemble(args)
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(res,indent=2,sort_keys=True)+"\n")
 print("NRSA_EPT_ORDER_MATCHED_HUC8_POSTHOC_RESULTS",json.dumps({
   "source":res["original_COUNT_QUALITY"],
   "matched":{kind:{"HUC8_pairs":report["source_preoutcome_matched_HUC8_pairs"],
     "EPT_total":report["matched_EPT_order_source_results"]["EPT_total_share"],
     "E_P_T":{k:report["matched_EPT_order_source_results"][k+"_total_share"]
       for k in ORDER},
     "within_EPT":{k:report["matched_EPT_order_source_results"][k+"_conditional_EPT"]
       for k in ORDER}}
      for kind,report in res["matched_blocks"].items()}},sort_keys=True),flush=True)
if __name__=="__main__":main()
