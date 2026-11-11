#!/usr/bin/env python3
"""Independent offline audit of 2021-calendar NPDES TN/TP DMR SOURCE gate.

No original national FY ZIP re-download, but independently recomputes
original 160 permit/344 feature candidate population, date union, 128 near
NWCA wetland support, exact 103 same state/ecoregion near, original far 137
and preregistered monitoring-availability gate from original SOURCE receipts.
Does NOT open NWCA 2021 vegetation Shannon/TN/TP response files.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from collections import defaultdict,Counter
from pathlib import Path

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
P=D/"NWCA2021_568WETLAND_160POTW_EXACT_PERMIT_FEATURE_SOURCE_BEFORE_2021DMR.json"
G=D/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
R=D/"NWCA2021_160POTW_ORIGINAL_CAL2021_DMR_TN_TP_EXO_FEATURE_SOURCE_GATE.json"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    old=json.loads((root/P).read_text())
    geo=json.loads((root/G).read_text())
    result=json.loads((root/R).read_text())
    if result["original_568site_160permit_344feature_prior_to_DMR_SHA256"]!=sha(root/P):
        raise ValueError("Pre-DMR 2026 registered source site+outfall number SHA mismatch")
    if result["original_568site_2026ECHO_STRICT_SEWAGE_geographic_PREOUTCOME_SHA256"]!=sha(root/G):
        raise ValueError("Original source-only site/outfall geometry changed")
    candidates=old["registered_2026_site_to_permit_feature_pairs_geodesic_le10km_without_2021_biology"]
    valid={(x["original_public_NPDES_permit_id"],x["original_public_NPDES_feature_number"]) for x in candidates}
    if len(valid)!=344 or len({p for p,f in valid})!=160:
        raise ValueError("Source selection outfall feature numbers changed")
    maps=defaultdict(set);groups={}
    for rec in candidates:
        sid=rec["source_Wetland_SITE_ID_VISIT_SHA256"]
        maps[sid].add((rec["original_public_NPDES_permit_id"],rec["original_public_NPDES_feature_number"]))
        groups[sid]=(rec["source_original_state"],rec["source_original_NWCA_ecoregion"])
    if len(maps)!=128:raise ValueError("Original 128 candidate-source near wetlands changed")
    q=geo["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]
    if len(q)!=568:raise ValueError("Original 568 NWCA wetlands not in previous source QA")
    near_regions={(x["original_state"],x["original_ecoregion"]) for x in q if x["outfall_rows_within_10km"]>0}
    far_regions={(x["original_state"],x["original_ecoregion"]) for x in q
       if x["nearest_outfall_geodesic_km_within25km_else_null"] is None or
       x["nearest_outfall_geodesic_km_within25km_else_null"]>=20}
    shared=near_regions&far_regions
    matchednear={sid for sid in maps if groups[sid] in shared}
    farids={x["NWCA2021_source_site_visit_SHA256"] for x in q
       if (x["original_state"],x["original_ecoregion"]) in shared and
       (x["nearest_outfall_geodesic_km_within25km_else_null"] is None or
        x["nearest_outfall_geodesic_km_within25km_else_null"]>=20)}
    if len(shared)!=32 or len(matchednear)!=103 or len(farids)!=137:
        raise ValueError("Pre-DMR original source in 32 shared state/ecoregions no longer identical")
    records=result["full_source_permit_feature_parameter_monitoring_receipts"]
    by={}
    for row in records:
        permit=str(row["original_public_NPDES_permit_id"]).upper()
        feature=str(row["original_public_NPDES_feature_number"]).upper()
        code=str(row["original_pollutant_code"])
        if (permit,feature) not in valid or code not in ("00600","00665","50050"):
            raise ValueError("Post-DMR source feature outside original 344 feature frame")
        key=(permit,feature,code)
        if key in by:raise ValueError("Duplicate 2021 original external feature nutrient source")
        dates=row["original_2021_valid_numeric_monitor_period_end_dates_unique"]
        if dates!=sorted(set(dates)):
            raise ValueError("Calendar2021 dates duplicated across FY source files")
        if any(not re.fullmatch(r"2021-\d{2}-\d{2}",s) for s in dates):
            raise ValueError("Noncalendar2021 monitoring date slipped through EPA DMR filter")
        if row["original_2021_valid_numeric_rows_no_dedup"]<len(dates):
            raise ValueError("More reported unique dates than valid DMR numeric values")
        if row["original_2021_positive_numeric_rows_no_dedup"]>row["original_2021_valid_numeric_rows_no_dedup"]:
            raise ValueError("Positive reported value count exceeded all numeric values")
        by[key]=row
    goodN={key[:2] for key,row in by.items() if key[2]=="00600" and row["original_2021_valid_numeric_monitor_period_end_dates_unique"]}
    goodP={key[:2] for key,row in by.items() if key[2]=="00665" and row["original_2021_valid_numeric_monitor_period_end_dates_unique"]}
    both=goodN&goodP
    if len(goodN)!=result["original_2021_calendar_DMR_EXO_feature_counts"]["any_valid_numeric_TN"] or len(goodP)!=result["original_2021_calendar_DMR_EXO_feature_counts"]["any_valid_numeric_TP"]:
        raise ValueError("Original source TN/TP feature-level counts not independently reproduced")
    if len(both)!=result["original_2021_calendar_DMR_EXO_feature_counts"]["both_TN_and_TP_same_permit_external_feature"]:
        raise ValueError("Original both nutrient DMR same-feature rule changed")
    def s(ids):
        features={ft for sid in ids for ft in maps[sid]}
        return {"original_wetland_sites":len(ids),
                "original_near_2021_monitored_permit_ids":len({f[0] for f in features}),
                "n_states":len({groups[sid][0] for sid in ids})}
    anyids={sid for sid,fts in maps.items() if fts&(goodN|goodP)}
    bothids={sid for sid,fts in maps.items() if fts&both}
    a1=s(anyids);a2=s(bothids);a3=s(matchednear&bothids)
    # Only count permits with monitoring evidence, not merely all available
    # registered feature permits inside source site catchments.
    a1["original_near_2021_monitored_permit_ids"]=len({f[0] for sid in anyids for f in maps[sid]&(goodN|goodP)})
    a2["original_near_2021_monitored_permit_ids"]=len({f[0] for sid in bothids for f in maps[sid]&both})
    a3["original_near_2021_monitored_permit_ids"]=len({f[0] for sid in matchednear&bothids for f in maps[sid]&both})
    qgate=result["precommitted_DMR_SOURCE_reporting_opportunity_gate"]
    passed=(a1["original_wetland_sites"]>=20 and a1["original_near_2021_monitored_permit_ids"]>=12 and a1["n_states"]>=8 and
            a2["original_wetland_sites"]>=12 and a2["original_near_2021_monitored_permit_ids"]>=8 and a2["n_states"]>=5 and
            a3["original_wetland_sites"]>=10 and len(farids)>=60)
    if passed!=qgate["DMR_monitored_same_2021_source_feature_biological_sensitivity_allowed"]:
        raise ValueError("Source-only precommitted EPA 2021 DMR gate decision not reproduced")
    orig1=result["original_2021_DMR_sites_with_valid_TN_OR_TP_same_source_feature"]
    orig2=result["original_2021_DMR_sites_with_valid_BOTH_TN_and_TP_same_source_feature"]
    for got,ref in ((a1,orig1),(a2,orig2)):
        if (got["original_wetland_sites"]!=ref["n_original_NWCA_source_wetland_centroids"] or
            got["original_near_2021_monitored_permit_ids"]!=ref["distinct_2021_measured_external_outfall_permit_IDs"] or
            got["n_states"]!=ref["n_US_states"]):
            raise ValueError("Independent site/permit/state opportunity differs from original")
    out={
       "schema":"gesui_US_2021_NWCA_568site_160POTW_344EXO_EPA_FY2021_FY2022_DMR_source_ONLY_independent_QA_v01",
       "source_SHA256":{"preDMR_crosswalk":sha(root/P),"geo_2026original":sha(root/G),
                         "DMR_2021_2022_join":sha(root/R)},
       "original_568_site_2026_160_permit_344_feature_sample_unchanged":True,
       "original_near10_all_wetlands":len(maps),
       "matched_2021_original_near_wetlands_32state_ecoregions":len(matchednear),
       "matched_far_original_state_ecoregions":len(farids),
       "EPA_DMR_2021_original_feature_level_numeric_TN":len(goodN),
       "EPA_DMR_2021_original_feature_level_numeric_TP":len(goodP),
       "EPA_DMR_2021_original_feature_level_both_numeric_TN_TP_same_feature":len(both),
       "source_anyN_orP_NWCA_near_wetland_IDs":a1,
       "source_bothN_andP_NWCA_near_wetland_IDs":a2,
       "source_bothN_andP_previously_matched_103_NEAR_wetlands":a3,
       "original_predeclared_biology_DMR_source_availability_gate_independently_reproduced":passed,
       "original_EPA_FY2021_2022_DMR_archive_actual_2021_calendar_validation_only":True,
       "EPA_2021_numeric_NP_at_permitted_feature_not_wetland_water_effluent_provenance":True,
       "no_shannon_or_FQAI_re_read":True
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(out,indent=2))
if __name__=="__main__":
    main()
