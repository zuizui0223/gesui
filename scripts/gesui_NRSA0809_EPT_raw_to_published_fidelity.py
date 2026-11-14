#!/usr/bin/env python3
"""Actual EPA2008 BERWW TOTAL300 + original EPA TAXA_ID→ORDER, biological
measurement fidelity ONLY, not N/P or PHab effects. Signed source ancestry
and 2008 published EPT_PIND same UID cross-validation before ecology.
"""
from __future__ import annotations
import hashlib,json,os,pathlib,urllib.request,time
import numpy as np,pandas as pd

R=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui0809_count_fidelity"
R.mkdir(parents=True,exist_ok=True)
source=json.loads(pathlib.Path(os.environ["ORIGINAL_TAXON_SOURCE_QA"]).read_text())
if source["schema"]!="gesui_2008_original_source_TAXA_ID_ORDER_BERWW_TOTAL300_bridge_v1":
 raise RuntimeError("Expected original independently audited taxon-ID source")
Q=source["original_EPA_ID_order_mapping_QA"]
if not Q["eligible_to_derive_raw_EPT_fraction"] or Q["original_UIDs_ge95percent_mapped_count"]<300:
 raise RuntimeError("Source-first taxonomy eligible gate did not pass")
valid=set(Q["original_UIDs_mapped_ge95percent"])
LOCK=json.loads(pathlib.Path(os.environ["ORIGINAL_2008_FIVE_SOURCE_QA"]).read_text())
allvalid=set(LOCK["original_wadeable_primary_source_gate"]["original_UIDs"])
if len(allvalid)!=1299 or not valid<=allvalid:raise RuntimeError("Original 2008 primary wadeable eligible UID frame changed")
LINK={
 "taxa":("https://www.epa.gov/sites/default/files/2016-06/nrsa_0809_benttaxa.csv",Q["taxonomy_csv_sha256"]),
 "bentcount":("https://www.epa.gov/sites/default/files/2016-11/nrsa0809bentcts.csv",Q["original_benthic_source_sha256"]),
 "bentmmi":("https://www.epa.gov/sites/default/files/2015-09/bentcond.csv",LOCK["original_EPA_CSVs_SHA256"]["bentmmi"]),
 "siteinfo":("https://www.epa.gov/sites/default/files/2015-09/siteinfo_0.csv",Q["original_siteinfo_sha256"]),
 "phab":("https://www.epa.gov/sites/default/files/2015-09/phabmed.csv",Q["original_phab_sha256"])}
def get(k):
 url,known=LINK[k];p=R/(k+".csv")
 for i in range(3):
  try:
   req=urllib.request.Request(url,headers={"User-Agent":"gesui-EPA2008-ept-native-taxonomy-count-fidelity"})
   with urllib.request.urlopen(req,timeout=120) as resp,p.open("wb") as o:
    while x:=resp.read(1024*1024):o.write(x)
   if hashlib.sha256(p.read_bytes()).hexdigest()!=known:raise RuntimeError("Original EPA raw file sha changed "+k)
   break
  except Exception:
   p.unlink(missing_ok=True)
   if i==2:raise
   time.sleep(2**i)
 d=pd.read_csv(p,dtype=str,keep_default_na=False,encoding="latin1",low_memory=False)
 p.unlink();d.columns=d.columns.str.upper().str.strip()
 return d
T=get("taxa");C=get("bentcount");M=get("bentmmi");S=get("siteinfo");P=get("phab")
T["taxkey"]=T.TAXA_ID.str.strip().str.upper();T["ord"]=T.ORDER.str.strip().str.upper()
u=T[T.taxkey.ne("")].groupby("taxkey").ord.nunique()
if (u>1).any():raise RuntimeError("Original taxon has multiple different ORDER codes")
lookup=T[T.taxkey.ne("")].drop_duplicates("taxkey").set_index("taxkey").ord
C=C[C.SAMPLE_TYPE.str.upper().str.strip().eq("BERWW")&
    C.SAMPLE_CAT.str.upper().str.strip().eq("P")&C.UID.isin(allvalid)].copy()
C["N"]=pd.to_numeric(C.TOTAL300,errors="coerce")
C=C[C.N.gt(0)].copy()
C["taxkey"]=C.TAXA_ID.str.strip().str.upper()
C["Order"]=C.taxkey.map(lookup).fillna("")
C["EPT"]=C.Order.isin(["EPHEMEROPTERA","PLECOPTERA","TRICHOPTERA"])
C["mapped"]=C.Order.ne("")
G=C.groupby("UID")
count=G.N.sum()
known=(C.N*C.mapped).groupby(C.UID).sum().reindex(count.index,fill_value=0)
ept=(C.N*C.EPT).groupby(C.UID).sum().reindex(count.index,fill_value=0)
available=set(count.index[count.gt(0)&known.div(count).ge(.95)])
if available!=valid:raise RuntimeError("Taxon source-only eligible UID mapping derived from real field data differs")
pct_lower=ept.div(count)*100
pct_upper=(ept+count-known).div(count)*100
pct_known=ept.div(known.replace(0,np.nan))*100
mm=M[M.UID.isin(available)].set_index("UID")
pub=pd.to_numeric(mm.EPT_PIND,errors="coerce")
source_compare=pd.DataFrame({"published_EPT_percent":pub,"reconstructed_lower_pct":pct_lower.reindex(pub.index),
 "reconstructed_upper_pct":pct_upper.reindex(pub.index),"among_mapped_pct":pct_known.reindex(pub.index)}).dropna()
if len(source_compare)<80:raise RuntimeError("Original EPA published EPT_PIND source overlap too sparse for independent validation")
source_compare["absolute_lower_pp"]=(source_compare.published_EPT_percent-source_compare.reconstructed_lower_pct).abs()
source_compare["absolute_known_pp"]=(source_compare.published_EPT_percent-source_compare.among_mapped_pct).abs()
source_compare["within_missing_taxon_bounds"]=(source_compare.published_EPT_percent>=source_compare.reconstructed_lower_pct-0.01)&(source_compare.published_EPT_percent<=source_compare.reconstructed_upper_pct+0.01)
compare={
 "original_published_EPT_PIND_comparable_real_UID_n":len(source_compare),
 "median_abs_diff_pp_lower_denominator":float(source_compare.absolute_lower_pp.median()),
 "pct_above_1_percentage_point_lower":float((source_compare.absolute_lower_pp>1).mean()*100),
 "median_abs_diff_pp_excluding_unclassified_count":float(source_compare.absolute_known_pp.median()),
 "pct_within_unclassified_bound":float(source_compare.within_missing_taxon_bounds.mean()*100),
 "original_percent_published_range":list(map(float,source_compare.published_EPT_percent.agg(["min","max"]))),
 "derived_share_source_UID_n":len(available),
 "derived_positive_count_total":float(count.loc[sorted(available)].sum()),
 "unknown_taxon_share_mean_pct":float(((count-known)/count).loc[sorted(available)].mean()*100),
 "unknown_taxon_share_max_pct":float(((count-known)/count).loc[sorted(available)].max()*100)}
for label,x in [("siteinfo",S),("phab",P)]:
 if "STRAHLERORDER" in x:
  z=x.loc[x.UID.isin(allvalid),"STRAHLERORDER"].str.strip()
  compare[label+"_original_STRORDER_top_codes"]=z.value_counts().head(9).to_dict()
# Fidelity is stringent; some difference may be due published condition metric
# formula not matching simple raw specimen count; this is a legitimate FAIL.
pass_fidelity=(compare["pct_above_1_percentage_point_lower"]<=10 and len(available)>=300)
OUT={
 "schema":"gesui_EPA0809_original_BERWW_TOTAL300_EPT_ORDER_with_legacy_metric_fidelity_v1",
 "pre_effect_1280_source_original_EPT_count_fidelity_report":compare,
 "true_predeclared_fidelity_gate_pass":bool(pass_fidelity),
 "original_taxonomy_SHA256":Q["taxonomy_csv_sha256"],
 "original_count_SHA256":Q["original_benthic_source_sha256"],
 "original_source_lock_SHA256":hashlib.sha256(pathlib.Path(os.environ["ORIGINAL_TAXON_SOURCE_QA"]).read_bytes()).hexdigest(),
 "eligible_native_UIDs_for_reconstructed_share":sorted(available),
 "individual_replicate_raw_EPT_fraction_not_EPA_MMI_metric":True,
 "reported_original_EPT_fraction_based_on_count_lower_bound_not_implicitly_nonEPT":True,
 "no_instream_chemistry_or_ecological_effect_modeled":True,
 "not_sewage_causal_or_biomass_service":True
}
out=R/"NRSA0809_EPT_REAL_TOTAL300_SOURCE_PUBLISHED_FIDELITY_QA.json"
out.write_text(json.dumps(OUT,sort_keys=True,indent=2))
print("NRSA0809_EPT_COUNT_FIDELITY_QA",json.dumps(compare,sort_keys=True),flush=True)
print("NRSA0809_EPT_RECONSTRUCTION_GATE",pass_fidelity,flush=True)
