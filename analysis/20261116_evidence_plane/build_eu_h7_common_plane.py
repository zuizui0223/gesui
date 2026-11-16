#!/usr/bin/env python3
"""EU same-unit evidence plane: 454 native H7 × source fixed 2017–19/2020–22.

This is a SOURCE JOIN (not a fresh regression or causal SEM). Native original
source files:
- Original 454 H7×period panel SHA-pinned; parent H6/period/country/method.
- 2023 EIV strict numeric accepted plant CWM 451 source-eligible periods
  from previously source-authenticated Actions #43718069015.
- Original 2018/2022 EEA annual reported WWTP mass/flow source index (454).
- Original upstream-vs-downstream reported load negative-control source (454).
- Original 2018/2022 station-level independently analyzed WWTP-reported
  positive mass and benthic EQR with TRUE 6,330 source native station-H7 map.

Different EQR indices and station subsets are retained as SEPARATE columns.
Observed WWTP effluent reported annual N/P mass is NOT measured in-stream TNTP.
No source-absent control is coded as clean water or zero actual wastewater.
No 2018/2022 survey station is independently replicated when aggregated to H7.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from pathlib import Path
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
PANEL=ROOT/"results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"
REPORTED=ROOT/"results/20261115_reported_np_hns/H7_454_REPORTED_SOURCE_INDEX_HNS_EQR_PLANT.csv"
DIRECTION=ROOT/"results/20261115_np_direction/H7_454_DIRECTIONAL_SOURCE_AUDIT_PANEL.csv"
BQE18=ROOT/"results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv"
BQE22=ROOT/"results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv"
PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
STATIONS_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
STATIONS_NAME="EEA2018_EEA2022_EPA1819_NATIVE_HYBAS7_SOURCE_ONLY_GEOGRAPHY.json"
BQE18_BLOB="0f30bfbc3dac13f6a182c8a910e9899c04778322"
BQE22_BLOB="325f199e985d9862ffee5b0edb205c61a7be3833"
REPORT_BLOB="45979613c1bc079f1c2346461eb7bdd31ac17e6a"
DIRECTION_BLOB="702131afc7d1a7bfb15e80357a506754cae7a70a"
EIV_SOURCE_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
EXPECTED_P=454
EXPECTED_NUTRI=451
PERIODS={"EEA18":"2017-2019","EEA22":"2020-2022"}
MODEL_OLD_EEA_YEAR={"EEA18":2018,"EEA22":2022}

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def gitblob(p):
    raw=p.read_bytes()
    return hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
def only(base,name):
    found=list(base.rglob(name))
    if len(found)!=1:raise RuntimeError("Source missing / ambiguous "+name)
    return found[0]
def check_keys(d,name):
    if d.duplicated(["h","period"]).any():raise RuntimeError(name+" duplicate original H7-period keys")

def locked_station_map(source_dir, fixed_h7):
    p=only(source_dir,STATIONS_NAME)
    if digest(p)!=STATIONS_SHA:raise RuntimeError("True original 6,330 station geocoded source SHA drift")
    r=json.loads(p.read_text())
    rows=r["candidate_source_site_period_rows"]
    if len(rows)!=6330:raise RuntimeError("Native EEA+EPA station original rows drift")
    seen={}
    n_with_H7=0
    for row in rows:
        panel=row.get("panel")
        if panel not in PERIODS:continue
        h=row.get("HYBAS_ID_L07")
        if h is None:continue
        h=str(h);period=PERIODS[panel]
        if (h,period) not in fixed_h7:continue
        native=str(row.get("native_id"))
        if native!=row.get("physical_id") or len(native.split("|"))!=3:
            raise RuntimeError("Native 2018/22 station key identity drift")
        k=(panel,native)
        if k in seen and seen[k]!=(h,period):
            raise RuntimeError("Native source station belongs to different true H7")
        seen[k]=(h,period)
        n_with_H7+=1
    if n_with_H7!=3525 or len(seen)!=3525:
        raise RuntimeError("True original EEA 2018+2022 source site-year 3525 mismatch")
    return seen

def site_eqr_h7(fixed_h7, station_map):
    cfg=[
      (BQE18,"EEA18","numeric_EQR_2018",
       "2018_measured_upstream_P_t_per_year",
       "2018_measured_upstream_N_t_per_year",BQE18_BLOB),
      (BQE22,"EEA22","normalised_EQR",
       "upstream_2022_EEA_P_measured_t_per_yr",
       "upstream_2022_EEA_N_measured_t_per_yr",BQE22_BLOB)
    ]
    all_data=[];status={}
    for source,panel,metric,p_mass,n_mass,sha in cfg:
        if gitblob(source)!=sha:raise RuntimeError("Authentic PR14/17 station EQR original Git blob mismatch "+str(source))
        d=pd.read_csv(source,dtype={"site_key":str},low_memory=False)
        if d.site_key.duplicated().any():raise RuntimeError("Original EQR EEA site-key population duplicated "+panel)
        if not {metric,p_mass,n_mass,"site_key"}.issubset(d.columns):
            raise RuntimeError("Original site record schema missing for "+panel)
        count_eligible=0;missing_geo=0
        for r in d[["site_key",metric,p_mass,n_mass]].itertuples(index=False,name=None):
            site,eqr,p,n=r
            pair=station_map.get((panel,site))
            if pair is None:
                missing_geo+=1;continue
            h,period=pair
            if (h,period) not in fixed_h7:raise RuntimeError("Site H7 not in true 454 panel")
            eqr=float(eqr) if pd.notna(eqr) else np.nan
            n=float(n) if pd.notna(n) else np.nan
            p=float(p) if pd.notna(p) else np.nan
            if not (np.isfinite(eqr) and 0<=eqr<=5):continue
            if not (np.isfinite(n) and np.isfinite(p) and n>=0 and p>=0):raise RuntimeError("Invalid native reported WWTP N/P site")
            all_data.append({"h":h,"period":period,"site_key":site,
                "benthic_site_EQR":eqr,
                "benthic_reported_P_tpy_index":p,"benthic_reported_N_tpy_index":n,
                "benthic_upstream_matched_report_positive":bool(n>0 or p>0),
                "benthic_source_panel":panel})
            count_eligible+=1
        status[panel]={"original_site_rows":len(d),
                       "site_rows_in_original_EU454_H7":count_eligible,
                       "original_sites_without_EU454_H7_mapping":missing_geo,
                       "station_git_blob_sha":sha,
                       "original_site_measurement_not_independent_of_H7_aggregation":True}
    v=pd.DataFrame(all_data)
    if len(v)==0:raise RuntimeError("No actual benthic observations in same original EU454")
    if v.duplicated(["h","period","site_key"]).any():
        raise RuntimeError("Site×H7-period duplicated on source site equality")
    # Site means are source descriptors, not independent observations. The
    # EEA H7 panel's own EQR assessment method may have different eligibility
    # or scaling; keep separately and never silently replace.
    g=v.groupby(["h","period"],dropna=False,sort=True).agg(
      benthic_sameH7_native_site_count=("site_key","nunique"),
      benthic_sameH7_site_equal_EQR=("benthic_site_EQR","mean"),
      benthic_sameH7_site_report_positive_count=("benthic_upstream_matched_report_positive","sum"),
      benthic_sameH7_site_P_tpy_index_mean=("benthic_reported_P_tpy_index","mean"),
      benthic_sameH7_site_N_tpy_index_mean=("benthic_reported_N_tpy_index","mean")
    ).reset_index()
    g["benthic_sameH7_any_positive_WWTP_report"]=g.benthic_sameH7_site_report_positive_count>0
    return g,status

def main(args):
    if digest(PANEL)!=PANEL_SHA:raise RuntimeError("True original 454 H7 panel SHA drift")
    for p,expected in ((REPORTED,REPORT_BLOB),(DIRECTION,DIRECTION_BLOB)):
        if gitblob(p)!=expected:raise RuntimeError("Original reported-source/direction blob changed "+str(p))
    original=pd.read_csv(PANEL,dtype={"h":str,"h6":str,"period":str},low_memory=False)
    required=["h","h6","period","country","method","W","crop","urban","riverflow",
              "EQR_original_site_equal_basin_period_mean","primary_record_mean",
              "primary_incidence_mean","plant_S50"]
    if any(v not in original.columns for v in required):
        raise RuntimeError("Original 454 source panel missing specified column "+
                           str([v for v in required if v not in original.columns]))
    if len(original)!=454 or original.h.nunique()!=232 or original.h6.nunique()!=90:
        raise RuntimeError("Original same-H7 454 × 232 H7 / 90 H6 grain changed")
    check_keys(original,"Original 454 source")
    if original.period.nunique()!=2:raise RuntimeError("Original two GBIF windows changed")
    fix=set(zip(original.h,original.period))
    trait_dir=args.eiv
    gate=json.loads(only(trait_dir,"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (not gate.get("hard_source_support_pass") or
        gate.get("schema")!="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA" or
        gate.get("source_EIV",{}).get("original_source_SHA256")!=EIV_SOURCE_SHA or
        gate.get("official_confirmed_valid_unscorable_nonbinomial_hybrids")!=175):
        raise RuntimeError("Original 4011 scored accepted plant species 451 source gate failed")
    elig=pd.read_csv(only(trait_dir,"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv"),dtype={"h":str,"period":str})
    means=pd.read_csv(only(trait_dir,"EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv"),dtype={"h":str,"period":str})
    if len(elig)!=454 or len(means)!=451:
        raise RuntimeError("Original EIV 454 opportunity/451 valid original H7-period count changed")
    check_keys(elig,"Original EIV eligibility")
    check_keys(means,"Original EIV means")
    eligible=elig[elig.full_source_trait_eligible.astype(str).str.lower()=="true"]
    if len(eligible)!=451 or set(zip(eligible.h,eligible.period))!=set(zip(means.h,means.period)):
        raise RuntimeError("Original valid EIV trait cohort not 451 exact H7-period keys")
    source_report=pd.read_csv(REPORTED,dtype={"h":str,"h6":str,"period":str})
    source_direction=pd.read_csv(DIRECTION,dtype={"h":str,"h6":str,"period":str})
    for d in (source_report,source_direction):
        if len(d)!=454 or set(zip(d.h,d.period))!=fix:raise RuntimeError("Original report/direction source does not cover full H7 panel")
        check_keys(d,"Original WWTP source")
        m=original.merge(d[["h","period","h6"]],on=["h","period"],validate="one_to_one")
        if not m.h6_x.eq(m.h6_y).all():raise RuntimeError("WWTP source true H6 parent changed")
    if int(source_report.reported_positive.sum())!=56 or int(source_direction.upstream_report_positive.sum())!=56:
        raise RuntimeError("Original 56 source-report positive H7×periods drifted")
    if not source_report.sort_values(["h","period"]).reported_positive.to_numpy().tolist()==source_direction.sort_values(["h","period"]).upstream_report_positive.to_numpy().tolist():
        raise RuntimeError("Source indicator mismatch with negative control upstream")
    station_map=locked_station_map(args.stations,fix)
    site_mean,site_qa=site_eqr_h7(fix,station_map)
    left=original[required+["nEEA","nAnimal","nPlant","nButterfly","primary_record_coverage_raw","primary_host_completeness"]].copy()
    e=eligible[["h","period","NUTRI_coverage_all_Plantae_speciesKeys",
                "exact_NUTRI_scored_vascular_accepted_species"]].copy()
    e=e.merge(means[["h","period","Plant_N_incidence","Plant_N_record"]],
              on=["h","period"],validate="one_to_one")
    m=left.merge(e,on=["h","period"],validate="one_to_one",how="left")
    if len(m)!=454:raise RuntimeError("Trait left join changed original grain")
    report=source_report[["h","period","reported_positive","reporting_class","n_source",
                           "n_mg_l","p_mg_l","log_reported_N_mgL_equivalent","log_reported_P_ugL_equivalent"]].copy()
    m=m.merge(report,on=["h","period"],validate="one_to_one",how="left")
    direction=source_direction[["h","period","upstream_report_positive","downstream_report_positive"]]
    m=m.merge(direction,on=["h","period"],validate="one_to_one",how="left")
    m=m.merge(site_mean,on=["h","period"],validate="one_to_one",how="left")
    if len(m)!=454 or m.duplicated(["h","period"]).any():
        raise RuntimeError("Evidence plane no longer exactly one row per original H7×period")
    m=m.rename(columns={"W":"modelled_static_wastewater_geography_W",
              "EQR_original_site_equal_basin_period_mean":"aquatic_EEA_original_EQR",
              "primary_incidence_mean":"butterfly_hostN_HNS_equal_species",
              "primary_record_mean":"butterfly_hostN_HNS_GBIF_record_weighted",
              "Plant_N_incidence":"plant_EIV_NUTRI_equal_species",
              "Plant_N_record":"plant_EIV_NUTRI_GBIF_record_weighted",
              "reported_positive":"EEA_upstream_positive_report_matched_flag",
              "n_mg_l":"EEA_upstream_report_N_per_model_riverflow_mgL_equivalent",
              "p_mg_l":"EEA_upstream_report_P_per_model_riverflow_mgL_equivalent",
              "plant_S50":"GBIF_plant_rarefied_S50",
              "upstream_report_positive":"EEA_upstream_positive_source_flag_independent",
              "downstream_report_positive":"EEA_downstream_positive_placebo_source_flag"})
    m["plant_EIV_original_source_eligible"]=m.plant_EIV_NUTRI_equal_species.notna()
    m["EEA_no_match_is_NOT_clean_control"]=True
    m["water_WWTP_report_N_and_P_are_NOT_inriver_TNTP"]=True
    m["water_benthic_site_report_mass_is_NOT_aquatic_exposure_dose"]=True
    m["butterfly_HNS_is_NOT_observed_larval_feeding_or_insect_count"]=True
    m["plant_EIV_NUTRI_is_NOT_measured_leaf_nitrogen"]=True
    m["native_sampling_unit"]="HydroBASINS H7 × original GBIF period"
    m["water_benthic_source_site_mean_EQR_may_differ_from_main_EEA_EQR"]=True
    # Report-matching status is kept from true source metadata. The literal
    # report-index 0 is NOT a measured zero effluent or no wastewater.
    m["water_source_report_unmatched_means_exposure_unknown"]=~m.EEA_upstream_positive_report_matched_flag.astype(bool)
    m=m.sort_values(["h","period"],kind="stable")
    allcore=m[["aquatic_EEA_original_EQR","plant_EIV_NUTRI_equal_species",
               "butterfly_hostN_HNS_equal_species"]].notna().all(axis=1)
    inc_source=m.EEA_upstream_positive_report_matched_flag.astype(bool)
    output=args.output
    output.mkdir(parents=True,exist_ok=True)
    csvpath=output/"EU_H7_PERIOD_454_WATER_PLANT_BUTTERFLY_WWTP_COMMON_PLANE.csv"
    m.to_csv(csvpath,index=False,float_format="%.12g")
    diff=m[["aquatic_EEA_original_EQR","benthic_sameH7_site_equal_EQR"]].dropna()
    empirical_qa={"schema":"GESUI_20261116_TRUE_EU_H7xPERIOD_WATER_PLANT_BUTTERFLY_PLANE_V1",
       "native_H7period_rows":len(m),"native_H7":int(m.h.nunique()),
       "native_true_H6":int(m.h6.nunique()),
       "matched_site_year_to_H7_mapping_records":len(station_map),
       "original_2022_2018_benthic_reporting_EQR_site_source":site_qa,
       "EIV_original_accepted_species_count":4011,
       "EIV_numeric_source_SHA256":EIV_SOURCE_SHA,
       "original_454_H7_panel_SHA256":PANEL_SHA,
       "original_reported_source_Git_blob_SHA1":REPORT_BLOB,
       "original_downstream_placebo_Git_blob_SHA1":DIRECTION_BLOB,
       "source_same_H7_benthic_site_groups":len(site_mean),
       "source_original_EIV_451_eligible_H7period":int(m.plant_EIV_original_source_eligible.sum()),
       "complete_same_H7_EQR_PlantN_HNS_original":int(allcore.sum()),
       "report_positive_same_H7periods":int(inc_source.sum()),
       "complete_core_and_upstream_source_positive":int((allcore&inc_source).sum()),
       "site_eqr_independent_recalculation_comparable_rows":len(diff),
       "site_eqr_vs_original_H7_panel_absolute_difference_quantile":
          diff.eval("abs(aquatic_EEA_original_EQR-benthic_sameH7_site_equal_EQR)").quantile([.05,.50,.95]).to_dict() if len(diff) else None,
       "GBIF_plant_HNS_not_field_diet_or_abundance":True,
       "no_true_local_sewage_exposure_or_causal_mediation_claim":True,
       "upstream_or_downstream_report_0_does_NOT_mean_absent_sewage":True,
       "separate_benthic_station_site_subsamples_not_counted_as_new_independent_H7":True,
       "cross_country_year_study_effects_not_pooled_into_one_coefficient":True,
       "generated_common_plane_csv_SHA256":digest(csvpath)}
    (output/"EU_H7_PERIOD_454_COMMON_PLANE_SOURCE_QA.json").write_text(
        json.dumps(empirical_qa,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_ONE_PLANE_EU_H7_PERIOD_454_SOURCE_QA="+json.dumps(
        empirical_qa,sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--eiv",required=True,type=Path)
    ap.add_argument("--stations",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    main(ap.parse_args())
