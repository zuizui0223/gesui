#!/usr/bin/env python3
"""SOURCE-ONLY: real EPA FY2021/FY2022 DMR reporting on near POTW permit/features.

Same 568 original NWCA 2021 wetland source sites and 2026-11-08 official EPA
NPDES permitted-feature source ZIP locked by prior SHA. This stage never opens
2021 vegetation H_ALL/FQAI nor actual NWCA water TN/TP values. DMR reported
numeric *availability* is evaluated; numeric concentrations/loads are neither
exported nor pooled across incomparable units. No receiving wetland hydrology.
"""
from __future__ import annotations
import argparse,hashlib,json,math,zipfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from gesui_NWCA2021_568source_EPA_NPDES_registered_sewage_outfall_geoONLY import (
    original_nwca_source_sites,ndpes_outfall_source,distances_km,xyz)

B=Path("results/us_NWCA2021_plant_TNTP_20261110")
SOURCE=B/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
NPDES_SHA="f0d1310b1c07085c42527e2ede5cf97ab83d517d50548b7dadb7621a44aa6160"
FULL_CALENDAR="2021-01-01..2021-12-31"
FISCAL_COVER={2021:("2021-01-01","2021-09-30"),2022:("2021-10-01","2021-12-31")}
BASENAME={2021:"npdes_dmr_fy2021.csv",2022:"npdes_dmr_fy2022.csv"}
COLS=("EXTERNAL_PERMIT_NMBR","PERM_FEATURE_NMBR","PERM_FEATURE_TYPE_CODE",
      "MONITORING_LOCATION_CODE","MONITORING_PERIOD_END_DATE","PARAMETER_CODE",
      "DMR_VALUE_NMBR","NODI_CODE","VALUE_TYPE_CODE","DMR_UNIT_DESC",
      "DMR_VALUE_QUALIFIER_CODE")
POLLUTANTS={"00600":"Total_N_as_N","00665":"Total_P_as_P"}
EFFLUENT_LOCATION=("1","2","EG")
CHUNK=150_000

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):h.update(block)
    return h.hexdigest()

def digest(t):return hashlib.sha256(t.encode()).hexdigest()

def load_original_locked_near_feature_pairs(root,nwca_dir,outfallzip):
    source=json.loads((root/SOURCE).read_text(encoding="utf-8"))
    if source["schema"]!="gesui_US_EPA_NWCA2021_568source_NPDES2026_reported_sewage_outfall_GEO_2021permitDates_sourceONLY_v01":
        raise ValueError("Original POTW geometry source schema changed")
    if not source["no_original_NWCA_H_ALL_FQAI_or_TNTP_outcomes_accessed_in_source_stage"]:
        raise ValueError("Original outfall source not selected without seeing 2021 NWCA plant values")
    if not source["source_precommitted_GEO_eligibility_gate"]["permitted_2021_outfall_distance_association_source_support_PASS"]:
        raise ValueError("Original near10/far20 source eligibility failed")
    if source["source_EPA_ECHO_NPDES_2026_live_weekly_original_ZIP_SHA256"]!=NPDES_SHA or sha(outfallzip)!=NPDES_SHA:
        raise ValueError("2026-11-08 original ECHO NPDES national outfall ZIP SHA mismatch")
    wetlands,_=original_nwca_source_sites(root,nwca_dir)
    facilities,_=ndpes_outfall_source(outfallzip)
    facilities=facilities.loc[facilities.PERMITTED_FEATURE_LATLON&
                             facilities.PERMIT_DATES_PLAUSIBLE_2021].copy()
    if len(wetlands)!=568 or len(facilities)!=10865:
        raise ValueError("Original pre-2021 POTW/Wetland source counts differ")
    original={r["NWCA2021_source_site_visit_SHA256"]:r
        for r in source["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]}
    if len(original)!=568:raise ValueError("Original preoutcome 568 wetland geometries missing")
    coords=facilities[["LAT","LON"]].to_numpy(float)
    tree=cKDTree(xyz(coords[:,0],coords[:,1]))
    maxchord=2*math.sin((10./6371.0088)/2)
    permits_by_site={}
    pairs_by_site={}
    unique_features=set()
    for row in wetlands.itertuples(index=False):
        sid=row.source_site_visit_SHA256
        if sid not in original:raise ValueError("Original wetland source ID differs")
        point=(float(row.LAT_ANALYS),float(row.LON_ANALYS))
        candidates=tree.query_ball_point(xyz(np.array([point[0]]),np.array([point[1]]))[0],r=maxchord+1e-11)
        subset=coords[candidates] if candidates else np.zeros((0,2))
        ds=distances_km(subset,point) if candidates else np.asarray([])
        features=set()
        permitids=set()
        for index,km in zip(candidates,ds):
            if km>10:continue
            permit=str(facilities.iloc[index].PERMIT_ID).strip().upper()
            feature=str(facilities.iloc[index].FEATURE_ID).strip().upper()
            features.add((permit,feature))
            permitids.add(permit)
        if sorted(digest(pid) for pid in permitids)!=original[sid]["independent_10km_NPDES_permits_SHA256"]:
            raise ValueError("Reconstructed ORIGINAL strict <=10km 2026 POTW permit set differs for NWCA 2021 site")
        if bool(features)!=bool(original[sid]["outfall_rows_within_10km"]):
            raise ValueError("Original site proximity near10 classification modified")
        permits_by_site[sid]=permitids
        pairs_by_site[sid]=features
        unique_features.update(features)
    near={k for k,v in pairs_by_site.items() if v}
    permitunion=set().union(*permits_by_site.values())
    if len(near)!=128 or len(permitunion)!=160:
        raise ValueError("Original 2026 POTW geo source NOT 128 sites / 160 unique permits")
    return source,wetlands,permits_by_site,pairs_by_site,permitunion,unique_features

def parse_fiscal_csv(zipfile_path,fiscal,permits,pairs):
    if fiscal not in FISCAL_COVER:raise ValueError("Only FY2021 and FY2022 originally selected")
    begin,end=FISCAL_COVER[fiscal]
    records=Counter()
    permits_any=set()
    permits_numeric=set()
    permits_exo=set()
    features_exo=set()
    features_tn=set()
    features_tp=set()
    permits_tn=set()
    permits_tp=set()
    dated=set()
    category=Counter()
    pollutant=Counter()
    per_source=defaultdict(lambda:Counter())
    with zipfile.ZipFile(zipfile_path) as z:
        names=[v for v in z.namelist() if v.split("/")[-1].lower()==BASENAME[fiscal]]
        if len(names)!=1:
            raise ValueError("Original 2021/FY2022 EPA national DMR archive does not contain expected one CSV "+repr(z.namelist()[:8]))
        with z.open(names[0]) as f:
            header=pd.read_csv(f,nrows=0).columns.tolist()
        missing=set(COLS)-set(header)
        if missing:raise ValueError("EPA original DMR schema changed "+repr(missing))
        with z.open(names[0]) as f:
            reader=pd.read_csv(f,usecols=list(COLS),dtype=str,keep_default_na=False,
                              encoding="utf-8",encoding_errors="replace",
                              on_bad_lines="error",chunksize=CHUNK,low_memory=False)
            for i,chunk in enumerate(reader,1):
                records["original_all_US_DMR_rows_scanned"]+=len(chunk)
                pid=chunk.EXTERNAL_PERMIT_NMBR.astype(str).str.upper().str.strip()
                mask=pid.isin(permits)
                if not mask.any():continue
                d=chunk.loc[mask].copy()
                records["source_near160_registered_permits_raw_DMR_fiscal_rows"]+=len(d)
                d["PERMIT"]=pid.loc[mask]
                dates=pd.to_datetime(d.MONITORING_PERIOD_END_DATE,errors="coerce")
                period=dates.between(pd.Timestamp(begin),pd.Timestamp(end),inclusive="both")
                d=d.loc[period].copy()
                if d.empty:continue
                records["source_permits_calendar2021_eligible_rows"]+=len(d)
                permits_any.update(d.PERMIT)
                features=d.PERM_FEATURE_NMBR.astype(str).str.upper().str.strip()
                d["FEATURE"]=features
                for pp in d.PERMIT.unique():per_source[pp]["any_calendar2021_DMR_rows"]+=int((d.PERMIT==pp).sum())
                no_data=d.NODI_CODE.astype(str).str.strip().eq("")
                value=pd.to_numeric(d.DMR_VALUE_NMBR,errors="coerce")
                measured=no_data & np.isfinite(value)
                if measured.any():
                    rec=d.loc[measured]
                    permits_numeric.update(rec.PERMIT)
                    records["source_permits_2021_any_numeric_DMR_rows"]+=len(rec)
                    for pp in rec.PERMIT.unique():per_source[pp]["any_numeric_rows"]+=int((rec.PERMIT==pp).sum())
                exo=d.PERM_FEATURE_TYPE_CODE.astype(str).str.upper().str.strip().eq("EXO")
                location=d.MONITORING_LOCATION_CODE.astype(str).str.upper().str.strip().isin(EFFLUENT_LOCATION)
                valid=measured&exo&location
                if valid.any():
                    r=d.loc[valid]
                    permits_exo.update(r.PERMIT)
                    records["source_permits_2021_EXO_effluent_numeric_rows"]+=len(r)
                    for pr,feature in zip(r.PERMIT,r.FEATURE):
                        if (pr,feature) in pairs:features_exo.add((pr,feature))
                    records["source_permits_2021_EXO_numeric_at_registered_near_feature_rows"]+=sum((pr,fe) in pairs for pr,fe in zip(r.PERMIT,r.FEATURE))
                    for pp in r.PERMIT.unique():per_source[pp]["EXO_numeric_rows"]+=int((r.PERMIT==pp).sum())
                for code,label in POLLUTANTS.items():
                    mask2=valid & d.PARAMETER_CODE.astype(str).str.strip().eq(code)
                    if not mask2.any():continue
                    r=d.loc[mask2]
                    records["2021_original_outfall_numeric_"+label+"_rows"]+=len(r)
                    (permits_tn if code=="00600" else permits_tp).update(r.PERMIT)
                    for pr,feature in zip(r.PERMIT,r.FEATURE):
                        if (pr,feature) in pairs:
                            (features_tn if code=="00600" else features_tp).add((pr,feature))
                    for yy in r.VALUE_TYPE_CODE.astype(str).str.upper().str.strip():
                        category[(label,yy)]+=1
                    for pp in r.PERMIT.unique():per_source[pp][label+"_numeric_rows"]+=int((r.PERMIT==pp).sum())
                if i%70==0:
                    print("FY_DMR_SCAN",fiscal,i,"chunks",records["original_all_US_DMR_rows_scanned"],
                          "near_original_permits_seen",len(permits_any),"strict_EXO_measured",len(permits_exo),
                          "with_TN",len(permits_tn),"with_TP",len(permits_tp),flush=True)
    def feature_set(lst):
        return sorted(digest(p+"|"+f) for p,f in lst)
    return {
        "fiscal_archive":fiscal,
        "calendar2021_monitoring_period_subset":[begin,end],
        "original_US_DMR_rows_scanned":records["original_all_US_DMR_rows_scanned"],
        "original_target_160_source_and_record_status":dict(records),
        "source_2021_permit_ids_with_any_DMR_record_SHA256":sorted(digest(p) for p in permits_any),
        "source_2021_permit_ids_with_any_numeric_DMR_SHA256":sorted(digest(p) for p in permits_numeric),
        "source_2021_permit_ids_with_EXO_effluent_numeric_DMR_SHA256":sorted(digest(p) for p in permits_exo),
        "source_2021_near10_registered_feature_pairs_with_EXO_effluent_numeric_DMR_SHA256":feature_set(features_exo),
        "source_2021_near10_feature_pairs_with_00600_totalN_DMR_SHA256":feature_set(features_tn),
        "source_2021_near10_feature_pairs_with_00665_totalP_DMR_SHA256":feature_set(features_tp),
        "source_2021_permit_ids_with_00600_totalN_DMR_SHA256":sorted(digest(p) for p in permits_tn),
        "source_2021_permit_ids_with_00665_totalP_DMR_SHA256":sorted(digest(p) for p in permits_tp),
        "source_2021_TN_and_TP_original_DMR_value_type_codes":{k+"|"+v:n for (k,v),n in category.items()},
        "source_2021_permit_DMR_summary_original_hashes":{
            digest(p):dict(x) for p,x in sorted(per_source.items())},
        "observed_DMR_numerical_values_are_NOT_REDISRIBUTED_OR_MODELED":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--fiscal-year",type=int,required=True,choices=[2021,2022])
    p.add_argument("--nwca-source-dir",type=Path,required=True)
    p.add_argument("--npdes-outfalls-zip",type=Path,required=True)
    p.add_argument("--dmr-zip",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    src,sites,permits_by_site,pairs_by_site,permitunion,featureunion=load_original_locked_near_feature_pairs(
         root,a.nwca_source_dir,a.npdes_outfalls_zip)
    data=parse_fiscal_csv(a.dmr_zip,a.fiscal_year,permitunion,featureunion)
    data.update({
       "schema":"gesui_US_NWCA2021_568source_160permit_ECHO_FY2021FY2022_actual_DMR_record_availability_SOURCE_ONLY_v01",
       "UTC":datetime.now(timezone.utc).isoformat(),
       "original_POTW_geo_before_DMR_SHA256":sha(root/SOURCE),
       "original_NPDES_2026_11_08_outfalls_ZIP_SHA256":sha(a.npdes_outfalls_zip),
       "original_ECHO_DMR_FY_archive_SHA256":sha(a.dmr_zip),
       "original_source_568_NWCA_probability_wetlands":len(sites),
       "original_source_near10_geographic_wetlands":sum(bool(x) for x in permits_by_site.values()),
       "original_source_near10_distinct_NPDES_permit_numbers":len(permitunion),
       "original_near10_distinct_registered_source_feature_pairs":len(featureunion),
       "original_near10_feature_pairs_SHA256":sorted(digest(p+"|"+f) for p,f in featureunion),
       "original_near10_potentially_2021_sewage_permits_SHA256":sorted(digest(p) for p in permitunion),
       "original_NEAR10_NWCA_wetland_to_target_NPDES_permit_hashes":{
           sid:sorted(digest(p) for p in permits) for sid,permits in permits_by_site.items() if permits},
       "original_NEAR10_NWCA_wetland_to_target_feature_pair_hashes":{
           sid:sorted(digest(p+"|"+f) for p,f in pairs) for sid,pairs in pairs_by_site.items() if pairs},
       "original_source_not_invented_WWTP_effects_from_NWCA_HALL_or_TNTP":True,
       "actual_2021_DMR_2026_revision_survey_window_only_not_full_historic_reporting_census":True,
       "geodesic_NPDES_10km_potential_not_verified_downstream_receiver":True,
    })
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
      "fiscal":a.fiscal_year,
      "sourceNear":data["original_source_near10_geographic_wetlands"],
      "sourcePermits":data["original_source_near10_distinct_NPDES_permit_numbers"],
      "sourceFeatures":data["original_near10_distinct_registered_source_feature_pairs"],
      "originalNationwideRows":data["original_US_DMR_rows_scanned"],
      "sourceYear2021Rows":data["original_target_160_source_and_record_status"],
      "DMR_known_perm":len(data["source_2021_permit_ids_with_any_DMR_record_SHA256"]),
      "DMR_exo":len(data["source_2021_permit_ids_with_EXO_effluent_numeric_DMR_SHA256"]),
      "DMR_near_feature":len(data["source_2021_near10_registered_feature_pairs_with_EXO_effluent_numeric_DMR_SHA256"]),
      "Ntotal_permits":len(data["source_2021_permit_ids_with_00600_totalN_DMR_SHA256"]),
      "Ptotal_permits":len(data["source_2021_permit_ids_with_00665_totalP_DMR_SHA256"]),
    },indent=2))

if __name__=="__main__":
    main()
