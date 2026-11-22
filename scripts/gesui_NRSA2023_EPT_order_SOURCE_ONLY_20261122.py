#!/usr/bin/env python3
"""EPA NRSA2023/24 original EPT order taxonomy *availability* source-only gate.

Only native taxonomic ORDER labels, UID/site/visit/date/method, source QA and
environment availability; DOES NOT READ TOTAL300 counts, benthic MMI,
per-group taxonomic abundances or biological response coefficients.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,pathlib
import pandas as pd
import gesui_NRSA_fourthcycle_2023_EPT_heldout_transport as prior
from gesui_nrsa2324_fourthcycle_source_only_20261122 import SHA,sha

EPT_ORDERS=("ephemeroptera","plecoptera","trichoptera")
def main():
 ap=argparse.ArgumentParser()
 for k in ("sites","water","habitat","chla","benthos","source_qa","out"):
  ap.add_argument("--"+k,type=pathlib.Path,required=True)
 a=ap.parse_args()
 for k,v in SHA.items():
  if sha(getattr(a,k))!=v:raise RuntimeError("Original EPA published source SHA changed: "+k)
 signed=json.loads(a.source_qa.read_text())
 if signed["original_EPA_sha256"]!=SHA:
  raise RuntimeError("Original 2023 SOURCE ONLY census SHA disagrees")
 d=prior.covariates(a,signed)
 assert len(d)==1165 and d.huc8.nunique()==693 and d.site_id.nunique()==1046
 raw_header=pd.read_csv(a.benthos,nrows=0).columns.tolist()
 required=["UID","SITE_ID","VISIT_NO","DATE_COL","SAMPLE_TYPE","ORDER","TOTAL300"]
 if not set(required).issubset(raw_header):
  raise RuntimeError("Original EPA benthiccount field schema changed")
 raw=pd.read_csv(a.benthos,usecols=required[:6],dtype=str,keep_default_na=False,low_memory=False)
 raw=raw[raw.UID.isin(d.index)].copy()
 if raw.empty:raise RuntimeError("Original selected EPA benthic rows unavailable")
 original_uid=set(d.index);raw_uids=set(raw.UID)
 if raw_uids!=original_uid:
  raise RuntimeError("Original source 2023 wadeable survey UIDs without benthos taxon rows")
 metadata_check=raw.groupby("UID",sort=False).agg({
  "SITE_ID":"nunique","VISIT_NO":"nunique","DATE_COL":"nunique","SAMPLE_TYPE":"nunique"})
 if not metadata_check.eq(1).all().all():
  raise RuntimeError("Original EPA within-UID sample ID/date/type conflicts")
 s=pd.read_csv(a.sites,usecols=["UID","SITE_ID","VISIT_NO","DATE_COL"],
               dtype=str,keep_default_na=False).set_index("UID")
 b=raw.drop_duplicates("UID").set_index("UID")
 s=s.loc[sorted(raw_uids)];b=b.loc[s.index]
 checks={
  "site_mismatched":int((s.SITE_ID!=b.SITE_ID).sum()),
  "visit_mismatched":int((s.VISIT_NO!=b.VISIT_NO).sum()),
  "date_mismatched":int((pd.to_datetime(s.DATE_COL,format="mixed")!=pd.to_datetime(b.DATE_COL,format="mixed")).sum()),
  "not_original_BERW":int(b.SAMPLE_TYPE.ne("BERW").sum())}
 if any(checks.values()):raise RuntimeError("EPA taxon metadata does not match original 2023 source cohort: "+str(checks))
 names=raw.ORDER.astype(str).str.lower().str.strip()
 labels=names.value_counts().to_dict()
 for n in EPT_ORDERS:
  if n not in labels:raise RuntimeError("Original benthic taxon source missing predefined E/P/T order: "+n)
 family_field="FAMILY" in raw_header
 # Source taxonomy coverage only: taxon row labels ≠ quantitative individuals.
 data={
  "schema":"gesui_NRSA2023_EPT_taxonomic_ORDER_availability_SOURCE_ONLY_v1",
  "original_EPA_source_SHA256":SHA,
  "original_source_huc8_sample_count":693,
  "original_source_2023_UID":1165,
  "original_source_site_count":1046,
  "original_sample_method":"BERW/WADEABLE",
  "sample_ID_and_survey_date_conflicts":checks,
  "source_original_BERW_taxon_rows":len(raw),
  "median_taxon_rows_per_UID":float(raw.groupby("UID").size().median()),
  "taxon_order_label_row_profiles":dict(sorted(labels.items())),
  "taxon_order_UniqueLabels":sorted(labels),
  "original_FAMILY_column_present":family_field,
  "original_TOTAL300_field_present":True,
  "TOTAL300_benthic_individual_counts_read":False,
  "order_fraction_by_source_history_read":False,
  "all_three_prespecified_EPT_ORDERs_in_original_source":True,
  "source_GO":True,
  "note":"Taxon-row count and ORDER labeling are source taxonomy availability, NOT biological abundance. Numeric E/P/T fractions only in subsequent source-frozen workflow."}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")
 print("NRSA_EPT_ORDER_SOURCE_ONLY_GATE",json.dumps({
  "UID":len(raw_uids),"rows":len(raw),"orders":{k:labels.get(k,0) for k in EPT_ORDERS},
  "family_column":family_field,"empty_order_rows":labels.get("",0),"GO":True
  },sort_keys=True),flush=True)
if __name__=="__main__":main()
