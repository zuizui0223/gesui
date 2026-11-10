#!/usr/bin/env python3
"""US EPA official FY2021 NPDES national DMR SOURCE ONLY on pre-frozen permits.

Streams original ~576MB ZIP in place; filters to the 160 original admin NPDES
permits within 10km of 128 original NWCA wetland sites, before ANY new
candidate-butterfly/plant ecological subgroup effect model. Official 2021
permit DMR may be on another outfall under the same permit; no source-to-
wetland hydrologic proof. No unstandardized DMR values exported or added.
"""
from __future__ import annotations

import argparse,csv,hashlib,io,json,math,zipfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

GEOMETRY=Path("results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json")
DMR_URL="https://echo.epa.gov/files/echodownloads/npdes_dmrs_fy2021.zip"
EXPECTED_CSV="npdes_dmr_fy2021.csv"
PARAMS={"00600":"TOTAL_NITROGEN","00665":"TOTAL_PHOSPHORUS"}
COLUMNS=(
    "EXTERNAL_PERMIT_NMBR","PERM_FEATURE_NMBR","PARAMETER_CODE",
    "MONITORING_PERIOD_END_DATE","DMR_VALUE_NMBR","NODI_CODE",
    "STANDARD_UNIT_DESC","DMR_UNIT_DESC","DMR_VALUE_STANDARD_UNITS",
    "VALUE_TYPE_CODE","MONITORING_LOCATION_CODE",
    "DMR_VALUE_QUALIFIER_CODE"
)
REQUIRED_NEAR_SITES=40
REQUIRED_SOURCE_PERMITS=15
REQUIRED_STATES=8

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for buf in iter(lambda:f.read(1024*1024),b""):
            h.update(buf)
    return h.hexdigest()

def hash_permit(s):
    return hashlib.sha256(str(s).strip().upper().encode()).hexdigest()

def original_source_permits(root):
    f=root/GEOMETRY
    x=json.loads(f.read_text())
    if x["schema"]!="gesui_US_EPA_NWCA2021_568source_NPDES2026_reported_sewage_outfall_GEO_2021permitDates_sourceONLY_v01":
        raise ValueError("Original national NPDES point source stage not available")
    if not x["source_precommitted_GEO_eligibility_gate"]["permitted_2021_outfall_distance_association_source_support_PASS"]:
        raise ValueError("Cannot open FY2021 DMR after original source-geometric failure")
    src=x["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]
    if len(src)!=568:raise ValueError("Original 568 first-visit probability-wetland source changed")
    near=[s for s in src if s["outfall_rows_within_10km"]>0]
    if len(near)!=128 or len({s["NWCA2021_source_site_visit_SHA256"] for s in near})!=128:
        raise ValueError("Original 128 near municipal-source wetland centers changed")
    permitted={k for s in near for k in s["independent_10km_NPDES_permits_SHA256"]}
    if len(permitted)!=160 or any(len(k)!=64 for k in permitted):
        raise ValueError("Source administrative distinct 160 NPDES IDs not intact")
    return x,near,permitted

def load_national_dmr(path,allowed_sha):
    if path.stat().st_size<100_000_000:
        raise ValueError("Official EPA FY2021 DMR archive unexpectedly small")
    original={}
    hit=[]
    meta=Counter()
    param=Counter()
    with zipfile.ZipFile(path) as z:
        target=[f for f in z.namelist() if f.split("/")[-1].lower()==EXPECTED_CSV.lower()]
        if len(target)!=1:
            raise ValueError("Original FY2021 ZIP has no unique NPDES_DMR_FY2021.csv")
        csv.field_size_limit(12*1024*1024)
        with z.open(target[0]) as binary:
            source=csv.reader(io.TextIOWrapper(binary,encoding="utf-8-sig",errors="replace",newline=""))
            columns=[s.upper().strip() for s in next(source)]
            required=set(COLUMNS)-set(columns)
            if required:
                raise ValueError("EPA FY2021 original DMR source schema missing "+repr(required))
            idx={field:columns.index(field) for field in COLUMNS}
            main=idx["EXTERNAL_PERMIT_NMBR"]
            pcode=idx["PARAMETER_CODE"]
            number=idx["DMR_VALUE_NMBR"]
            nodi=idx["NODI_CODE"]
            feature=idx["PERM_FEATURE_NMBR"]
            loc=idx["MONITORING_LOCATION_CODE"]
            valtype=idx["VALUE_TYPE_CODE"]
            def rec(k):
                if k not in original:
                    original[k]={
                      "original_DMR_2021_source_rows":0,
                      "original_DMR_ANY_parameter_numeric_n":0,
                      "number_distinct_features_with_numeric_DMR":set(),
                      "NUMERIC_TN_00600":0,
                      "NUMERIC_TP_00665":0,
                      "TN_features":set(),"TP_features":set(),
                      "NODI_TN_rows":0,"NODI_TP_rows":0,
                      "unit_descriptions_by_parameter":defaultdict(Counter),
                      "value_type_by_parameter":defaultdict(Counter),
                      "monitoring_location_by_parameter":defaultdict(Counter),
                      "actual_period_end_year_by_parameter":defaultdict(Counter)
                    }
                return original[k]
            cached={}
            maxindex=max(idx.values())
            for i,row in enumerate(source,1):
                meta["original_national_DMR_csv_rows"]+=1
                if len(row)<=maxindex:
                    meta["original_national_bad_short_rows"]+=1
                    continue
                pid=row[main].strip().upper()
                if pid not in cached:
                    digest=hash_permit(pid)
                    cached[pid]=digest if digest in allowed_sha else None
                digest=cached[pid]
                if digest is None:continue
                r=rec(digest)
                r["original_DMR_2021_source_rows"]+=1
                meta["original_160_permit_source_rows"]+=1
                code=row[pcode].strip().zfill(5)
                value=row[number].strip()
                no_data=row[nodi].strip()
                if code in PARAMS:
                    param[code+"_source_DMREntries"]+=1
                    if no_data:r["NODI_"+("TN" if code=="00600" else "TP")+"_rows"]+=1
                if not value or no_data:continue
                try:
                    v=float(value)
                except ValueError:
                    meta["non_numeric_reported_DMR_original_values"]+=1
                    continue
                if not math.isfinite(v) or v<0:
                    meta["reported_DMR_negative_or_nonfinite_source_values"]+=1
                    continue
                meta["numeric_original_matched_permit_DMR_values"]+=1
                r["original_DMR_ANY_parameter_numeric_n"]+=1
                r["number_distinct_features_with_numeric_DMR"].add(row[feature].strip())
                if code in PARAMS:
                    nutrient="TN" if code=="00600" else "TP"
                    r["NUMERIC_"+nutrient+"_"+code]+=1
                    r[nutrient+"_features"].add(row[feature].strip())
                    r["unit_descriptions_by_parameter"][code].update([(row[idx["DMR_UNIT_DESC"]].strip(),row[idx["STANDARD_UNIT_DESC"]].strip())])
                    r["value_type_by_parameter"][code][row[valtype].strip()]+=1
                    r["monitoring_location_by_parameter"][code][row[loc].strip()]+=1
                    mdate=row[idx["MONITORING_PERIOD_END_DATE"]].strip()
                    r["actual_period_end_year_by_parameter"][code][mdate[:4] if mdate else "(missing)"]+=1
                    param[code+"_usable_numeric_n"]+=1
                if i%2_000_000==0:
                    print("NPDES_FY2021_NATIONAL_ROWS",i,
                          "matched160_permit_rows",meta["original_160_permit_source_rows"],
                          "permits_matched",len(original),flush=True)
    compact={}
    for k,r in original.items():
        compact[k]={
            "original_DMR_2021_source_rows":r["original_DMR_2021_source_rows"],
            "original_DMR_ANY_parameter_numeric_n":r["original_DMR_ANY_parameter_numeric_n"],
            "number_distinct_features_with_numeric_DMR":len(r["number_distinct_features_with_numeric_DMR"]),
            "NUMERIC_TN_00600":r["NUMERIC_TN_00600"],
            "NUMERIC_TP_00665":r["NUMERIC_TP_00665"],
            "unique_TN_DMR_outfall_feature_codes":len(r["TN_features"]),
            "unique_TP_DMR_outfall_feature_codes":len(r["TP_features"]),
            "NODI_TN_rows":r["NODI_TN_rows"],"NODI_TP_rows":r["NODI_TP_rows"],
            "unit_descriptions_by_parameter":{p:{str(a)+" | "+str(b):int(n) for (a,b),n in cnt.items()} for p,cnt in r["unit_descriptions_by_parameter"].items()},
            "value_type_by_parameter":{p:dict(cnt) for p,cnt in r["value_type_by_parameter"].items()},
            "monitoring_location_by_parameter":{p:dict(cnt) for p,cnt in r["monitoring_location_by_parameter"].items()},
            "actual_period_end_year_by_parameter":{p:dict(cnt) for p,cnt in r["actual_period_end_year_by_parameter"].items()},
            "FY2021_NPDES_same_PERMIT_other_outfall_may_have_been_monitored":True,
        }
    if meta["original_national_DMR_csv_rows"]<100_000:
        raise ValueError("Original national FY2021 NPDES DMR archive has unexpectedly few rows")
    return compact,dict(meta),dict(param)

def summary(root,dmr_path):
    source,near,permits=original_source_permits(root)
    data,meta,parameters=load_national_dmr(dmr_path,permits)
    rows={k:data.get(k,{}) for k in sorted(permits)}
    for k,r in rows.items():
        if not r:
            rows[k]={"original_DMR_2021_source_rows":0,
                     "original_DMR_ANY_parameter_numeric_n":0,
                     "NUMERIC_TN_00600":0,"NUMERIC_TP_00665":0}
    any_p={k for k,r in rows.items() if r["original_DMR_ANY_parameter_numeric_n"]>0}
    tn_p={k for k,r in rows.items() if r["NUMERIC_TN_00600"]>0}
    tp_p={k for k,r in rows.items() if r["NUMERIC_TP_00665"]>0}
    nutrient=tn_p|tp_p
    both=tn_p&tp_p
    groups={"any_numeric_DMR_FY2021":any_p,
            "total_N_00600_numeric_FY2021":tn_p,
            "total_P_00665_numeric_FY2021":tp_p,
            "total_N_OR_P_numeric_FY2021":nutrient,
            "total_N_AND_P_numeric_FY2021":both}
    reports={}
    for name,perm in groups.items():
        sites=[s for s in near if perm.intersection(s["independent_10km_NPDES_permits_SHA256"])]
        reports[name]={
          "original_NWCA_568_source_near10_sites_with_geographically_associated_PERMIT_reporting_n":len(sites),
          "original_reg_sewage_NPDES_permits_with_FY2021_monitoring_data":len(perm),
          "original_distinct_US_states_with_monitor_and_near_site":len({s["original_state"] for s in sites}),
          "original_ecoregions_with_monitor_and_near_site":len({s["original_ecoregion"] for s in sites}),
          "original_near_source_2021_wetland_site_hashes_only":[s["NWCA2021_source_site_visit_SHA256"] for s in sites],
          "not_verified_that_the_same_within10km_outfall_feature_or_downstream_wetland_received_the_reported_matter":True,
        }
    q=reports["total_N_OR_P_numeric_FY2021"]
    passed=bool(q["original_NWCA_568_source_near10_sites_with_geographically_associated_PERMIT_reporting_n"]>=REQUIRED_NEAR_SITES and
                q["original_reg_sewage_NPDES_permits_with_FY2021_monitoring_data"]>=REQUIRED_SOURCE_PERMITS and
                q["original_distinct_US_states_with_monitor_and_near_site"]>=REQUIRED_STATES)
    return {
        "schema":"gesui_EPA_NWCA2021_160_source_POTW_NPDES_permits_original_FY2021_DMR_TNTP_Monitoring_Source_only_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_2026_NPDES_registered_outfall_568_NWCA_source_SHA256":sha(root/GEOMETRY),
        "original_EPA_fiscalyear2021_national_NPDES_DMR_ZIP_SHA256":sha(dmr_path),
        "original_official_ECHO_national_FY2021_DMR_url":DMR_URL,
        "fiscal2021_us_EPA_monitor_period_not_same_as_calendar_2021":True,
        "original_160_administratively_related_sewage_NPDES_source_permit_HASHes_n":len(permits),
        "original_128_NWCA_2021_plant_waterjoint_wetland_sites_near10_this_permit_universe":len(near),
        "FY2021_national_original_DMR_source_row_count_QA":meta,
        "original_nutrient_monitoring_parameter_code_00600TN_00665TP_QA":parameters,
        "original_160_registered_sewage_source_permit_count_summary":reports,
        "original_160_registered_near_outfall_sewage_permits_FY2021_anyNutrient_model_support":{
           "min_2021_NWCA_near_registered_outfall_site_centers_with_any_real_FY2021_N_or_P_DMR_at_same_ADMIN_permit":REQUIRED_NEAR_SITES,
           "min_distinct_municipal_ADMIN_permits_with_nutrient_DMR":REQUIRED_SOURCE_PERMITS,
           "min_US_states_with_source_DMR_N_or_P":REQUIRED_STATES,
           "2021_FY_Nutrient_MONITORING_SOURCE_GATE_PASS":passed,
           "NO_NEW_2021_NPDES_direct_wetland_source_to_vegetation_model_PERMITTED_without_exact_outfall_ID_and_hydrological_route":True
        },
        "every_160_original_nearby_candidate_administrative_sewage_permit_DMR_availability_QA":rows,
        "this_stage_NEVER_accessed_2021_wetland_Shannon_FQAI_or_source_TNTP_values":True,
        "DMR_same_ADM_PERMIT_not_verified_near_source_PERM_FEATURE_or_wetland_receiving_hydrology":True,
        "DMR_DMR_UNIT_DESC_may_be_mass_or_concentration_not_aggregated_as_same_NutrientLoad":True,
        "no_new_sewage_ecological_causal_effect_claimed":True
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--dmr-fy2021-zip",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=summary(root,args.dmr_fy2021_zip)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
       "DMR_source":out["FY2021_national_original_DMR_source_row_count_QA"],
       "Nutrients":out["original_nutrient_monitoring_parameter_code_00600TN_00665TP_QA"],
       "source_counts":{k:{kk:vv for kk,vv in val.items() if kk!="original_near_source_2021_wetland_site_hashes_only"}
           for k,val in out["original_160_registered_sewage_source_permit_count_summary"].items()},
       "gate":out["original_160_registered_near_outfall_sewage_permits_FY2021_anyNutrient_model_support"]
    },indent=2))

if __name__=="__main__":
    main()
