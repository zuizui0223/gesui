#!/usr/bin/env python3
"""SOURCE-ONLY original EPA 2008/2013/2018 taxon-count denominator census.

Actual same 2648 native UID panel and canonical HUC8 established by prior
PR81 source workflow. Re-download original publisher source-SHA-pinned full
Benthic Count tables (and 2008 Taxa_ID -> ORDER source). Counted-organism
opportunities, count integer, taxonomic coverage and original EPT counted
fraction concordance. NO rarefaction *richness values* or ecological fits
are calculated by this source-only file.
"""
from __future__ import annotations
import collections,hashlib,json,os,pathlib,time,urllib.request
import numpy as np,pandas as pd

TEMP=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))
ROOT=TEMP/"gesui_threecycle_rarefied_taxa"
ROOT.mkdir(parents=True,exist_ok=True)
COHORT=TEMP/"gesui_threecycle_ept_transfer"/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv"
df=pd.read_csv(COHORT,dtype={"source_UID":str,"source_cycle":str,"huc8":str,"state":str,"eco9":str})
if len(df)!=2648 or df.duplicated(["source_cycle","source_UID"]).any():
 raise RuntimeError("Original PR81 exact source n2648 per era must be verified before taxon count analysis")
expected={"2008_09":(586,437,48),"2013_14":(1116,706,48),"2018_19":(946,603,47)}
for cycle,(n,h,s) in expected.items():
 d=df[df.source_cycle.eq(cycle)]
 if len(d)!=n or d.huc8.nunique()!=h or d.state.nunique()!=s:
  raise RuntimeError("Original authentic EPA cycle cohort changed before new count response "+cycle)
q13=json.loads(pathlib.Path(os.environ["NRSA1314_SEMANTIC_SOURCE_QA"]).read_text())
if q13["schema"]!="NRSA1314_ORIGINAL_EPA_2019_SIX_TABLE_PRENUMERIC_SEMANTIC_CONTRACT_v1":
 raise RuntimeError("Prior original EPA 2013 count SHA evidence not loaded")
old13=q13["original_six_source_receipts_sha256"]["benthic"]
SHAS={
 "2008_09":"8251fbe2d2bfc699225873ec43231feda0b110bd112badb61d6e8de486870d31",
 "2013_14":old13,
 "2018_19":"4e9b3a0798b25f67c7a96d4800b94f86544bbd623abaf119358e067b158324a7",
 "2008_taxa":"6165d5a7ec8d4ee3e85149255fd967f6268d870d49c57de6302ca4e124673b26"}
URL={
 "2008_09":"https://www.epa.gov/sites/default/files/2016-11/nrsa0809bentcts.csv",
 "2013_14":"https://www.epa.gov/sites/default/files/2019-04/nrsa1314_bentcnts_04232019.csv",
 "2018_19":"https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_benthic_macroinvertebrate_count_-_data.csv",
 "2008_taxa":"https://www.epa.gov/sites/default/files/2016-06/nrsa_0809_benttaxa.csv"}
def fetch(key):
 p=ROOT/(key+".csv")
 for i in range(4):
  try:
   req=urllib.request.Request(URL[key],headers={"User-Agent":"gesui-EPA-3ERA-original-BERWW-BERW-rarefaction-source-20261114"})
   with urllib.request.urlopen(req,timeout=120) as resp,p.open("wb") as out:
    while buf:=resp.read(1048576):out.write(buf)
   h=hashlib.sha256(p.read_bytes()).hexdigest()
   if h!=SHAS[key]:raise RuntimeError("Original EPA publisher source SHA changed "+key)
   break
  except Exception:
   p.unlink(missing_ok=True)
   if i==3:raise
   time.sleep(2**i)
 d=pd.read_csv(p,dtype=str,keep_default_na=False,encoding="latin1",low_memory=False)
 p.unlink();d.columns=d.columns.str.strip().str.upper()
 return d
taxa=fetch("2008_taxa")
taxa["key"]=taxa.TAXA_ID.str.strip().str.upper()
taxa["order"]=taxa.ORDER.str.strip().str.upper()
if (taxa[taxa.key.ne("")].groupby("key").order.nunique()>1).any():
 raise RuntimeError("Original 2008 TAXA_ID -> ORDER contradicted its source")
mapping=taxa[taxa.key.ne("")].drop_duplicates("key").set_index("key").order
summary={};rows=[]
for cycle in expected:
 raw=fetch(cycle)
 d=df[df.source_cycle.eq(cycle)]
 ids=set(d.source_UID)
 met="BERWW" if cycle=="2008_09" else "BERW"
 B=raw.loc[raw.UID.isin(ids) &
           raw.SAMPLE_TYPE.str.strip().str.upper().eq(met)].copy()
 if cycle=="2008_09":
  B=B[B.SAMPLE_CAT.str.strip().str.upper().eq("P")].copy()
 B["count"]=pd.to_numeric(B.TOTAL300,errors="coerce")
 if (B["count"].dropna()<0).any():raise RuntimeError("Negative original counted organisms")
 B=B[B["count"].gt(0)].copy()
 B["key"]=B.TAXA_ID.str.strip().str.upper()
 if cycle=="2008_09":B["order"]=B.key.map(mapping).fillna("")
 else:
  if "ORDER" not in B:raise RuntimeError("Original later EPA benthic ORDER field missing")
  B["order"]=B.ORDER.str.strip().str.upper()
 B["is_ept"]=B.order.isin(("EPHEMEROPTERA","PLECOPTERA","TRICHOPTERA"))
 B["known_order"]=B.order.ne("")
 B["known_taxon"]=B.key.ne("")
 B["known_both"]=B.known_order & B.known_taxon
 # Taxon records may recur with the same UID+taxon code; count is summed
 # within that original taxon only, never counted as an extra taxon.
 if (B["count"]-B["count"].round()).abs().gt(1e-6).any():
  raise RuntimeError("Original TOTAL300 count not an integer, cannot hypergeom rarefy: "+cycle)
 group=B.groupby("UID")
 counts=group["count"].sum()
 counts_ept=(B["count"]*B["is_ept"]).groupby(B.UID).sum()
 old=d.set_index("source_UID")
 if set(counts.index)!=ids:raise RuntimeError("Some original source native UID lacks positive Benthic counts "+cycle)
 share=(counts_ept/counts).reindex(old.index)
 diff=(share-old.ept_fraction).abs()
 if diff.max()>1e-9:
  raise RuntimeError("Original count-derived EPT share disagrees with PR81 locked original source for "+cycle+":"+str(diff.max()))
 valid=(B["count"]*B["known_both"]).groupby(B.UID).sum().reindex(old.index)
 quality=valid/counts.reindex(old.index)
 raw_uid_taxon_dup=int(B.duplicated(["UID","key"]).sum())
 groupbykey=B.groupby(["UID","key"],sort=False)
 if raw_uid_taxon_dup:
  # EPT original ORDER information must agree for duplicate UID-TAXA_ID.
  same=B.groupby(["UID","key"]).order.nunique()
  if (same>1).any():raise RuntimeError("Repeated original UID TAXA_ID has conflicting ORDER "+cycle)
 by_taxon=groupbykey.agg(count=("count","sum"),order=("order","first"),known=("known_both","min")).reset_index()
 by_taxon["cycle"]=cycle
 # Preserve the exact measured specimen counts only in ephemeral runner,
 # never upload raw native source taxonomic tables or coordinates.
 rows.append(by_taxon[["cycle","UID","key","order","count","known"]])
 assert len(counts)==len(d)
 audit={
  "original_source_UIDs":len(d),"original_HUC8":d.huc8.nunique(),
  "source_positive_BERW_BERWW_taxon_rows":len(B),
  "duplicate_same_UID_TAXA_ID_rows":raw_uid_taxon_dup,
  "distinct_positive_original_taxonID_per_cycle":B.key.nunique(),
  "original_total_source_count_min":int(counts.min()),
  "original_total_source_count_quantiles":{str(i):float(v) for i,v in counts.quantile([0,.05,.25,.5,.75,.95,1]).items()},
  "n_source_original_count_at_least100":int((counts>=100).sum()),
  "n_source_original_count_at_least200":int((counts>=200).sum()),
  "n_source_original_count_at_least300":int((counts>=300).sum()),
  "source_nonempty_taxonID_and_order_coverage_mean":float(quality.mean()),
  "n_source_UID_known_order_taxonID_at_least95pct":int((quality>=.95-1e-10).sum()),
  "n_at_least200_and_95pct_known_order_ID":int(((counts.reindex(old.index)>=200)&(quality>=.95-1e-10)).sum()),
  "previously_frozen_original_EPT_share_max_absolute_discrepancy":float(diff.max()),
  "source_sample_type":met,"source_original_SHA256":SHAS[cycle],
  "non_integer_positive_original_counts_n":0}
 eligible=old[(counts.reindex(old.index)>=200)&(quality>=.95-1e-10)].copy()
 audit["200individual_eligible_n_HUC8_state"]=[
  len(eligible),eligible.huc8.nunique(),eligible.state.nunique()]
 audit["rare200_source_minimum_n250_HUC150_states20_pass"]=(len(eligible)>=250 and eligible.huc8.nunique()>=150 and eligible.state.nunique()>=20)
 print("ORIGINAL_EPA_RAREFACTION_SOURCE_CENSUS",cycle,json.dumps(audit,sort_keys=True),flush=True)
 summary[cycle]=audit
counts=pd.concat(rows,ignore_index=True)
if counts.duplicated(["cycle","UID","key"]).any():
 raise RuntimeError("Original raw EPA grouped UID+TAXA_ID not unique after original count summation")
counts.to_csv(ROOT/"AUDIT_ONLY_ORIGINAL_EPA_3CYCLE_TAXON_COUNTS.csv",index=False)
out={
 "schema":"gesui_EPA3cycle_SOURCEONLY_BERW_TAXAID_count200_rarefaction_population_v1",
 "original_frozen_2648_EPT_relative_fraction_native_UID":True,
 "source_original_publisher_count_SHA256":SHAS,
 "original_source_all_3cycle_BERW_BERWW_total300":True,
 "quality_gates_by_cycle":summary,
 "precommitted_m200_gate_all_3cycle_pass":all(x["rare200_source_minimum_n250_HUC150_states20_pass"] for x in summary.values()),
 "source_original_fraction_reconstructed_directly_from_counts_and_matched":True,
 "no_rarefied_response_values_evaluated_before_source_gate":True,
 "no_source_outcome_effect_model_fitted":True,
 "not_species_richness_or_true_naturally_occurring_abundance":True,
 "not_real_wwtp_effluent_or_emergence_service":True}
p=ROOT/"NRSA_THREECYCLE_EPT_TAXA_SOURCE_CENSUS_PRE_EFFECT.json"
p.write_text(json.dumps(out,indent=2,sort_keys=True))
print("ORIGINAL_EPA_RAREFACTION_SOURCE_GATE",out["precommitted_m200_gate_all_3cycle_pass"])
