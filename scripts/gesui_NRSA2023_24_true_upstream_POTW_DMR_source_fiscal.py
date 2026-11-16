#!/usr/bin/env python3
"""Source-only historical EPA NPDES FY2023/24/25 Discharge Monitoring Reports.

The original 195 NRSA 2023–24 UID-to-657 exact up-river registered outfalls
are frozen from independent #43449175435 before DMR number inspection.

Scan entire official fiscal archive, retaining only exact permit+external-outfall
feature+effluent location+total-N/total-P observed reporting opportunity.
Record VALUE_TYPE_CODE/DMR_UNIT_DESC separately (quantity vs concentration).
No pollutant amounts, annual mass conversion, ecosystem effects or zero controls
are produced here. Historical-source failures are NOT environmental zero.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,math,os,pathlib,zipfile
import numpy as np,pandas as pd

FY=(2023,2024,2025)
MIN_STRICT_UPSTREAM_VISITS=40
EXPECTED_TOPOLOGY_RUN=43449175435
FIELDS=[
 "EXTERNAL_PERMIT_NMBR","PERM_FEATURE_NMBR","PERM_FEATURE_TYPE_CODE",
 "MONITORING_LOCATION_CODE","MONITORING_PERIOD_END_DATE","PARAMETER_CODE",
 "DMR_VALUE_NMBR","NODI_CODE","VALUE_TYPE_CODE","DMR_UNIT_DESC",
 "DMR_VALUE_QUALIFIER_CODE","DMR_VALUE_ID","STATISTICAL_BASE_CODE"
]
CODES={"00600":"total_N","00665":"total_P"}
EFFLUENT_LOC={"1","2","EG"}

def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  while b:=f.read(1<<20):h.update(b)
 return h.hexdigest()

def normalize_feature(x):
 v=str(x).strip().upper()
 if v.isdigit() and len(v)<=3:return v.zfill(3)
 return v

def ready_sources(qa_path,geo_path):
 original=json.loads(qa_path.read_text())
 geo=json.loads(geo_path.read_text())
 if original["schema"]!="gesui_NRSA_2018_19_2023_24_UPSTREAM_POTW_SOURCE_ONLY_NO_DMR_JOIN_v1":
  raise RuntimeError("Original independently audited NLDI source-only schema changed")
 if sha(geo_path)!=original["original_geocandidate_manifest_SHA256"]:
  raise RuntimeError("Original frozen EPA ECHO geometry SHA mismatch")
 e=original["eras"]["2023_24"]
 if not e["pre_DMR_source_gate_meets_MIN40sites_25features_10HUC8_5states"]:
  raise RuntimeError("Previously fixed true-upstream source gate failed: do not scan DMR")
 if e["n_verified_different_reach_upstream_registered_outfall_site_visits"]!=195:
  raise RuntimeError("Pre-outcome source 195 UID frozen count unexpectedly changed")
 original_uids=set()
 features=set()
 sites_by_source=collections.defaultdict(set)
 original_geo={r["uid"]:r for r in geo["eras"]["2023_24"]["source_records"]}
 src_year={}
 for z in e["verified_source_only_upstream_candidate_pairs"]:
  uid=z["uid"]
  feature=(z["permit"].upper(),normalize_feature(z["feature"]))
  if uid not in original_geo:raise RuntimeError("NLDI original UID no longer in source geography")
  original_uids.add(uid);features.add(feature);sites_by_source[feature].add(uid)
 for uid in original_uids:
  date=pd.to_datetime(original_geo[uid]["date_col"],errors="coerce")
  if pd.isna(date) or date.year not in (2023,2024):
   raise RuntimeError(f"Original true-upstream native NRSA UID has invalid survey-year source: {uid}")
  src_year[uid]=int(date.year)
 if len(original_uids)!=195 or len(features)!=657:
  raise RuntimeError(f"Unexpected distinct UID/source counts {len(original_uids)}/{len(features)}")
 byyear=collections.defaultdict(set)
 for feat,ids in sites_by_source.items():
  for uid in ids:byyear[(feat,src_year[uid])].add(uid)
 return original, original_uids, features, src_year, byyear

def archive_scan(fiscal,archive,geo_source,source_qa):
 original,uids,features,uidyears,byyear=ready_sources(source_qa,geo_source)
 with zipfile.ZipFile(archive) as z:
  csvfiles=[x for x in z.namelist() if x.lower().split("/")[-1] in
            (f"npdes_dmrs_fy{fiscal}.csv",f"npdes_dmr_fy{fiscal}.csv")]
  if len(csvfiles)!=1:raise RuntimeError("No single original EPA DMR CSV in official fiscal archive: "+str(z.namelist()[:8]))
  filename=csvfiles[0]
  with z.open(filename) as f:
   hdr=pd.read_csv(f,nrows=0).columns.tolist()
  missing=set(FIELDS)-set(hdr)
  if missing:raise RuntimeError("Official EPA original DMR header changed "+repr(sorted(missing)))
  scans=collections.Counter()
  year_feature=collections.defaultdict(lambda:collections.Counter())
  type_units=collections.Counter()
  positive_feature_param_year=collections.defaultdict(set)
  reported_feature_param_year=collections.defaultdict(set)
  seen_dmr_ids=set()
  with z.open(filename) as f:
   reader=pd.read_csv(f,usecols=FIELDS,dtype=str,keep_default_na=False,
                      chunksize=150000,encoding_errors="replace",on_bad_lines="error")
   for i,df in enumerate(reader,1):
    scans["all_official_dmr_rows"]+=len(df)
    permits=df.EXTERNAL_PERMIT_NMBR.str.strip().str.upper()
    # Restrict early by original historical source-eligible permit.
    relevant_permits={p for p,feat in features}
    filt=permits.isin(relevant_permits)
    if not filt.any():continue
    d=df.loc[filt].copy()
    d["permit"]=permits.loc[filt]
    d["feature"]=d.PERM_FEATURE_NMBR.map(normalize_feature)
    eligible=np.array([(p,f) in features for p,f in zip(d.permit,d.feature)])
    d=d.loc[eligible].copy()
    if d.empty:continue
    scans["exact_permit_external_feature_raw_rows"]+=len(d)
    d=d.loc[d.PERM_FEATURE_TYPE_CODE.str.upper().str.strip().eq("EXO") &
            d.MONITORING_LOCATION_CODE.str.upper().str.strip().isin(EFFLUENT_LOC)].copy()
    if d.empty:continue
    scans["EXO_effluent_monitoring_rows"]+=len(d)
    d=d.loc[d.PARAMETER_CODE.str.strip().isin(CODES)].copy()
    if d.empty:continue
    scans["original_00600_00665_source_rows"]+=len(d)
    dates=pd.to_datetime(d.MONITORING_PERIOD_END_DATE,errors="coerce")
    years=dates.dt.year
    d=d.loc[years.isin([2023,2024])].copy()
    if d.empty:continue
    d["year"]=years.loc[d.index].astype(int)
    scans["original_calendar_2023_24_source_rows"]+=len(d)
    valid_idx=[]
    for idx,row in d.iterrows():
      pair=(row["permit"],row["feature"])
      if (pair,int(row["year"])) not in byyear:continue
      val=str(row["DMR_VALUE_NMBR"]).strip()
      try:
       numeric=float(val)
       numeric_ok=math.isfinite(numeric)
      except (TypeError,ValueError):
       numeric=math.nan;numeric_ok=False
      if str(row["NODI_CODE"]).strip() or not numeric_ok:continue
      # DMR value ID deduplicates original repeated/revised archive rows if present.
      dmr_id=str(row["DMR_VALUE_ID"]).strip()
      if dmr_id:
       if dmr_id in seen_dmr_ids:continue
       seen_dmr_ids.add(dmr_id)
      year=int(row["year"])
      param=CODES[str(row["PARAMETER_CODE"]).strip()]
      vt=str(row["VALUE_TYPE_CODE"]).strip().upper()
      units=str(row["DMR_UNIT_DESC"]).strip()
      typ="quantity" if vt.startswith("Q") else ("concentration" if vt.startswith("C") else "other")
      key=(year,pair[0],pair[1],param,typ,vt,units)
      year_feature[(year,pair[0],pair[1],param,typ)]["numeric_report_rows"]+=1
      if numeric>0:
       year_feature[(year,pair[0],pair[1],param,typ)]["strictly_positive_report_rows"]+=1
       positive_feature_param_year[(year,param,typ)].add(pair)
      reported_feature_param_year[(year,param,typ)].add(pair)
      type_units[(year,param,vt,units)]+=1
      scans["valid_original_numeric_TN_TP_monitoring_rows"]+=1
      scans["numeric_nonzero_rows"]+=int(numeric>0)
    if i%120==0:print("FY_DMR_SOURCE_PROGRESS",fiscal,"chunks",i,"raw",scans["all_official_dmr_rows"],
                      "matched",scans["valid_original_numeric_TN_TP_monitoring_rows"],flush=True)

 report={
  "schema":"gesui_NRSA2023_24_EXACT_UPSTREAM_POTW_OFFICIAL_DMR_SOURCE_ONLY_FY_v1",
  "fiscal_year":fiscal,
  "source_archive_name":archive.name,"source_archive_sha256":sha(archive),
  "official_internal_csv":filename,
  "source_NLDI_independent_6shard_manifest_SHA256":sha(source_qa),
  "source_ECHO_historical_geometry_manifest_SHA256":sha(geo_source),
  "original_frozen_upstream_NRSA_UIDs":len(uids),
  "original_frozen_upstream_distinct_external_outfall_features":len(features),
  "source_counts":dict(scans),
  "year_param_unit_value_type_summary":[{"calendar_year":y,"parameter":p,"value_type_code":v,"dmr_unit":u,"rows":n}
                                      for (y,p,v,u),n in sorted(type_units.items())],
  "exact_outfall_year_param_observed_categories":[
       {"calendar_year":y,"permit":p,"feature":f,"parameter":para,"class":typ,
        "reported_numeric_rows":d["numeric_report_rows"],
        "strictly_positive_numeric_rows":d["strictly_positive_report_rows"]}
       for (y,p,f,para,typ),d in sorted(year_feature.items())],
  "no_DMR_numeric_values_exported":True,
  "no_assumed_annual_mass_from_concentration_or_single_report":True,
  "not_instream_delivered_effluent_dose":True,
  "not_a_negative_control_or_biological_effect":True
 }
 return report

def main():
 a=argparse.ArgumentParser()
 a.add_argument("--fiscal",type=int,required=True,choices=FY)
 a.add_argument("--dmr-archive",type=pathlib.Path,required=True)
 a.add_argument("--source-qa",type=pathlib.Path,required=True)
 a.add_argument("--geo-source",type=pathlib.Path,required=True)
 a.add_argument("--out",type=pathlib.Path,required=True)
 x=a.parse_args()
 result=archive_scan(x.fiscal,x.dmr_archive,x.geo_source,x.source_qa)
 x.out.parent.mkdir(parents=True,exist_ok=True)
 x.out.write_text(json.dumps(result,indent=2,sort_keys=True))
 print("FY_DMR_SOURCE_AUDIT_RESULT",json.dumps({"FY":x.fiscal,
       "source_archive_sha256":result["source_archive_sha256"],
       "source_counts":result["source_counts"],
       "report_feature_rows":len(result["exact_outfall_year_param_observed_categories"]),
       "unit_value_type_rows":len(result["year_param_unit_value_type_summary"])},sort_keys=True))
if __name__=="__main__":main()
