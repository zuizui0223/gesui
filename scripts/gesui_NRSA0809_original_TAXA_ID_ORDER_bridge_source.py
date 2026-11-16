#!/usr/bin/env python3
"""Source-first original EPA 2008–09 native TAXA_ID->ORDER audit and BERWW
TOTAL300 availability, before testing any stream chemistry/biodiversity signs.
"""
import hashlib,json,os,pathlib,urllib.request,time
import numpy as np,pandas as pd

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui0809_taxon_source"
ROOT.mkdir(parents=True,exist_ok=True)
base=json.loads(pathlib.Path(os.environ["NRSA0809_FROZEN_SOURCE_QA"]).read_text())
origin=json.loads(pathlib.Path(os.environ["NRSA0809_FIRST_EPA_SOURCE"]).read_text())
E=base["original_wadeable_primary_source_gate"]
eligible=set(E["original_UIDs"])
if len(eligible)!=1299:raise RuntimeError("Prior original before-numeric source fixed wadeable UID count differs")
URL={
 "taxa":"https://www.epa.gov/sites/default/files/2016-06/nrsa_0809_benttaxa.csv",
 "bentcount":"https://www.epa.gov/sites/default/files/2016-11/nrsa0809bentcts.csv",
 "phab":"https://www.epa.gov/sites/default/files/2015-09/phabmed.csv",
 "siteinfo":"https://www.epa.gov/sites/default/files/2015-09/siteinfo_0.csv"}
def fetch(label):
 p=ROOT/(label+".csv")
 for k in range(3):
  try:
   req=urllib.request.Request(URL[label],headers={"User-Agent":"gesui-2008-EPA-taxonomy-source-audit"})
   with urllib.request.urlopen(req,timeout=120) as w,p.open("wb") as out:
    while b:=w.read(1024*1024):out.write(b)
   sha=hashlib.sha256(p.read_bytes()).hexdigest()
   if label in origin["original_data_sources"] and sha!=origin["original_data_sources"][label]["source_SHA256"]:
    raise RuntimeError("Previously original source pinned EPA SHA changed")
   break
  except Exception:
   p.unlink(missing_ok=True)
   if k==2:raise
   time.sleep(2**k)
 d=pd.read_csv(p,dtype=str,keep_default_na=False,encoding="latin1",low_memory=False)
 p.unlink()
 d.columns=d.columns.str.upper().str.strip()
 return d,sha
taxa,tsh=fetch("taxa")
counts,bsh=fetch("bentcount")
phab,psh=fetch("phab")
site,ssh=fetch("siteinfo")
print("NRSA0809_TAXON_SOURCE_COLUMN_SCHEMA",json.dumps({x:list(d.columns) for x,d in
     [("taxa",taxa),("bentcount",counts)]}),flush=True)
for k in ("TAXA_ID","ORDER"):
 if k not in taxa:raise RuntimeError("Original EPA 2008 Taxa List missing "+k)
if any(k not in counts for k in ("UID","SAMPLE_TYPE","SAMPLE_CAT","TAXA_ID","TOTAL300")):
 raise RuntimeError("Original EPA 2008 source benthic Count missing mandatory source fields")
taxa["order"]=taxa.ORDER.str.strip().str.upper()
taxa["key"]=taxa.TAXA_ID.str.strip().str.upper()
taxa=taxa[taxa.key.ne("")].copy()
ambig=taxa.groupby("key").order.nunique()
ambig_keys=set(ambig[ambig>1].index)
taxa_unique=taxa.drop_duplicates(["key","order"])
bycount=taxa_unique.groupby("key").order.agg(lambda x:next(iter(x)) if len(set(x))==1 else "")
if not eligible <= set(site.UID):raise RuntimeError("Source only actual EPA SiteInfo UID disappeared")
cnt=counts[counts.UID.isin(eligible)&counts.SAMPLE_TYPE.str.strip().str.upper().eq("BERWW")&
           counts.SAMPLE_CAT.str.strip().str.upper().eq("P")].copy()
cnt["key"]=cnt.TAXA_ID.str.strip().str.upper()
cnt["n"]=pd.to_numeric(cnt.TOTAL300,errors="coerce")
if (cnt.n.dropna()<0).any():raise RuntimeError("Original source counted subsample negative")
positive=cnt[cnt.n.gt(0)].copy()
positive["original_order"]=positive.key.map(bycount).fillna("")
positive["mapped"]=positive.original_order.ne("")&~positive.key.isin(ambig_keys)
q=positive.groupby("UID").agg(original_positive_count=("n","sum"),mapped_positive_count=("n",lambda v:0)) if False else None
by_site=positive.groupby("UID").n.sum()
mapped=positive[positive.mapped].groupby("UID").n.sum().reindex(by_site.index,fill_value=0)
rowmapped=positive.groupby("UID").mapped.mean().reindex(by_site.index,fill_value=0)
v=positive[positive.mapped].n.sum()/max(1,positive.n.sum())
completes=set(by_site.index[by_site.gt(0)&mapped.div(by_site).ge(.95)])
pb=phab[phab.UID.isin(eligible)].copy()
ss=site[site.UID.isin(eligible)].copy()
samples={}
for name,d in (("siteinfo",ss),("phab",pb)):
 for col in ("STRAHLERORDER","STRAHLER_ORD","STR_ORD"):
  if col in d:
   z=pd.to_numeric(d[col],errors="coerce")
   samples[name+"_"+col]={"numeric":int(z.notna().sum()),"positive":int((z>0).sum()),
        "source_raw_nonempty":int(d[col].str.strip().ne("").sum())}
original={
 "taxonomy_csv_sha256":tsh,"original_benthic_source_sha256":bsh,
 "original_phab_sha256":psh,"original_siteinfo_sha256":ssh,
 "taxon_list_rows":len(taxa),
 "taxon_list_unique_taxonID":int(taxa.key.nunique()),
 "taxonID_with_conflicting_original_ORDER":len(ambig_keys),
 "positive_original_BERWW_taxon_rows":len(positive),
 "positive_original_TOTAL300_individuals":float(positive.n.sum()),
 "matched_original_order_positive_COUNT_rows":int(positive.mapped.sum()),
 "fraction_positive_original_taxon_rows_exact_ORDER":float(positive.mapped.mean()) if len(positive) else 0,
 "fraction_positive_original_counted_individuals_exact_ORDER":float(v),
 "original_positive_1299_source_eligible_UIDs":len(by_site),
 "original_UIDs_ge95percent_mapped_count":len(completes),
 "original_UIDs_mapped_ge95percent":sorted(completes),
 "native_stream_order_source_coverage":samples,
 "old_frozen_mmi_EPT_pct_source_count_286_from_noeffect_QA":True,
 "source_only_no_new_EPT_fraction_regression":True,
 "order_conflicts_not_resolved_by_guess":True,
 "eligible_to_derive_raw_EPT_fraction":(len(completes)>=300 and v>=.95 and len(ambig_keys)==0)
}
print("NRSA0809_TAXON_ORDER_ORIGINAL_QA",
 json.dumps({k:v for k,v in original.items() if not isinstance(v,list)},sort_keys=True),flush=True)
out={
 "schema":"gesui_2008_original_source_TAXA_ID_ORDER_BERWW_TOTAL300_bridge_v1",
 "original_1299_source_eligible_UIDs":len(eligible),
 "source_original_2008_first_source_SHA256":hashlib.sha256(pathlib.Path(os.environ["NRSA0809_FIRST_EPA_SOURCE"]).read_bytes()).hexdigest(),
 "original_2008_source_gate_sha256":hashlib.sha256(pathlib.Path(os.environ["NRSA0809_FROZEN_SOURCE_QA"]).read_bytes()).hexdigest(),
 "original_EPA_ID_order_mapping_QA":original,
 "no_original_biological_or_nutrient_numeric_effect_fit":True,
 "not_periphyton_2008_or_wastewater_load":True}
p=ROOT/"NRSA0809_BERWW_TOTAL300_TAXA_ID_ORDER_ORIGINAL_CROSSWALK_SOURCE_QA.json"
p.write_text(json.dumps(out,sort_keys=True,indent=2))
