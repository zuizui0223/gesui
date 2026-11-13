#!/usr/bin/env python3
"""Original SOURCE-ONLY actual 2022 measured WWTP annual reported P/N source
coverage for exactly 273 source-frozen physically same 2022 water TP / French
EEA I2M2 benthic EQR sites in 74 original H7 basins. Outcome magnitude not
opened. Also audit 2018 EEA I2M2 same-site biological station continuity.

Both WQ and EEA source records are PRIOR derived restricted 2022 EEA WWTP
point-source inventory indices, NOT independent real-world P-flow sensors,
and 0=not matched reported source, NOT untreated/no effluent.
"""
import csv,io,json,hashlib,math,collections,os,pathlib,statistics
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_FR2022_273_TP_EQR_WWTP_SOURCE_GATE"
ROOT.mkdir(parents=True,exist_ok=True)
PAIR_SHA="448add2d6850cbc876c3e6b90bca278c6b711e64f9d363071e93bc0c18b84611"
H6_SHA="294eb69b76cd4cf40f1022a5e83278529b1d2bbceac9a38931004610a23eb052"
GITS={
 "EEA22":("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv","325f199e985d9862ffee5b0edb205c61a7be3833"),
 "EEA18":("results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv","0f30bfbc3dac13f6a182c8a910e9899c04778322"),
 "WQ22":("results/eea2022_measured_nutrient_load_river_20261107/TP_2022_reported_wwtp_load_and_measured_river_station.csv","a75a4cf140de64fd6038079a3fac6756a461f3f8")
}
def source_json(path,sha):
    raw=pathlib.Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=sha:raise RuntimeError("Original pre-outcome source SHA mismatch "+path)
    return json.loads(raw)
pairs=source_json(os.environ["PRE_ORIGINAL_273TP_EQR_SITE_MATCH"],PAIR_SHA)["original_454_H7_cohort_2022_genuine_TP_250m_source_matched_pairs"]
g=source_json(os.environ["PRE_ORIGINAL_74H7_H6_PARENT_FOLDS"],H6_SHA)["matched_original_74_H7_fold_source"]
parents={r["H7"]:r["H6"] for r in g}
if len(pairs)!=273 or len({r["EEA2022_site_key"] for r in pairs})!=273 or len({r["WQ_ID"] for r in pairs})!=273 or len(parents)!=74 or len(set(parents.values()))!=24:raise RuntimeError("Original 273/74/24 genuine geographic population changed")
def gitsha(raw):
    return hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
tables={}
for label,(path,sha) in GITS.items():
    raw=pathlib.Path(path).read_bytes()
    if gitsha(raw)!=sha:raise RuntimeError("Original publisher-derived source Git blob mismatch "+label)
    tab=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"),newline="")))
    if len(tab)!={"EEA22":2475,"EEA18":2634,"WQ22":1291}[label]:raise RuntimeError("Original 2022 source table size changed "+label)
    tables[label]=tab
old22={r["site_key"]:r for r in tables["EEA22"]}
old18={r["site_key"]:r for r in tables["EEA18"]}
WQ={r["WQ.ID"]:r for r in tables["WQ22"]}
if len(old22)!=2475 or len(old18)!=2634 or len(WQ)!=1291:
    raise RuntimeError("Original geographic site IDs not unique")
def f(v):
    try:
        n=float(v)
        return n if math.isfinite(n) else None
    except (ValueError,TypeError):return None
def dist_m(la1,lo1,la2,lo2):
    rad=math.pi/180
    dlat=(la2-la1)*rad;dlon=(lo2-lo1)*rad
    a=math.sin(dlat/2)**2+math.cos(la1*rad)*math.cos(la2*rad)*math.sin(dlon/2)**2
    return 2*6371008.8*math.asin(min(1.,math.sqrt(a)))
# Primary original actual 2022 WWTP source coverage; never use columns
# normalised_EQR or real_2022_station_median_mgL to select or calculate source.
source_rows=[];reporting=collections.Counter()
for q in pairs:
    site=old22[q["EEA2022_site_key"]]
    water=WQ[q["WQ_ID"]]
    if not(site["HYRIV_ID"]==water["HYRIV_ID"]==q["HYRIV_ID"]):
        raise RuntimeError("Original physical source native HYRIV reach mismatch")
    if site["procedureClassificationSystem"]!="I2M2" or site["countryCode"]!="FR":
        raise RuntimeError("Original FR I2M2 native source population modified")
    h=q["EEA_H7_L07"]
    d=dist_m(f(site["lat"]),f(site["lon"]),f(water["WQ.lat"]),f(water["WQ.lon"]))
    if d>250.00001 or abs(d-f(q["spherical_source_coordinate_separation_m"]))>0.3:
        raise RuntimeError("Original 273 source TP-EQR location match differs")
    p=f(site["upstream_2022_EEA_P_measured_t_per_yr"])
    n=f(site["upstream_2022_EEA_N_measured_t_per_yr"])
    p_wq=f(water["eea2022_upstream_P_measured_load_nativeunit"])
    count=f(water["eea2022_upstream_wwtp_both_positive_count"])
    if p is None or n is None or p_wq is None or count is None:
        reporting["source_index_missing_positive_reporting_unknown"]+=1
        pass
    elif p<0 or n<0 or p_wq<0 or count<0:
        raise RuntimeError("Original EEA source P/N inventory had negative load")
    if p is not None and p_wq is not None and not math.isclose(p,p_wq,rel_tol=1e-7,abs_tol=1e-6):
        reporting["two_derivation_p_mismatch_at_same_river_source"]+=1
    if (p or 0)>0:reporting["EEA22_site_positive_P"]+=1
    if (n or 0)>0:reporting["EEA22_site_positive_N"]+=1
    if (p_wq or 0)>0:reporting["WQ22_site_positive_P"]+=1
    if (count or 0)>0:reporting["WQ22_site_any_matched_WWTP_both_measured"]+=1
    source_rows.append({
      "EEA2022_site_key":q["EEA2022_site_key"],"Peters_2022_WQ_ID":q["WQ_ID"],
      "original_HYBAS_ID_L07":h,"original_H6_parent":parents[h],
      "original_2022_native_reach":q["HYRIV_ID"],
      "original_EEA_WQ_source_coordinate_offset_m":d,
      "original_EEA_2022_reported_upstream_P_inventory_tpy":p,
      "original_EEA_2022_reported_upstream_N_inventory_tpy":n,
      "original_Peters_WQ_2022_reported_upstream_P_inventory_nativeunit":p_wq,
      "original_Peters_2022_matched_positive_source_count":count,
      "original_source_P_2022_matched_positive":p is not None and p>0,
      "original_source_status_not_physical_no_WWTP_when_zero":True,
      "original_I2M2_EQR_numeric_NOT_OPENED":True,
      "original_real_2022_TP_numeric_NOT_OPENED":True
    })
if len(source_rows)!=273 or len({r["EEA2022_site_key"] for r in source_rows})!=273:
    raise RuntimeError("Source-only 2022 WWTP-report overlap one-to-one failed")

positive=[r for r in source_rows if r["original_source_P_2022_matched_positive"]]
p_h7={r["original_HYBAS_ID_L07"] for r in positive}
p_h6={r["original_H6_parent"] for r in positive}
source_2022_gate=(len(positive)>=30 and len(p_h7)>=15 and len(p_h6)>=5)
# Secondary source-only 2018 biological field-site availability, NEVER
# look at EQR scores or TP chemistry. First exact site ID, then same river
# reach & <=250m and deterministic distance tie, 1:1 EEA2018 sites.
by_reach=collections.defaultdict(list)
for r in tables["EEA18"]:
    if r["countryCode"]=="FR" and r["procedureClassificationSystem"]=="I2M2":
        by_reach[r["HYRIV_ID"]].append(r)
choices=[]
for record in source_rows:
    curr=old22[record["EEA2022_site_key"]]
    k=record["EEA2022_site_key"]
    la=f(curr["lat"]);lo=f(curr["lon"]);riv=curr["HYRIV_ID"]
    exact=old18.get(k)
    if exact is not None and exact["countryCode"]=="FR" and exact["procedureClassificationSystem"]=="I2M2" and exact["HYRIV_ID"]==riv:
        d=dist_m(la,lo,f(exact["lat"]),f(exact["lon"]))
        if d<=250:
            choices.append((0,d,k,exact["site_key"],"original_native_site_key_identical"))
    for prev in by_reach.get(riv,()):
        d=dist_m(la,lo,f(prev["lat"]),f(prev["lon"]))
        if d<=250:
            choices.append((1,d,k,prev["site_key"],"same_original_HYRIV_250m_geo_fallback"))
choices.sort(key=lambda z:(z[0],z[1],z[2],z[3]))
matched={};used18=set()
for rank,d,k,k18,mode in choices:
    if k in matched or k18 in used18:continue
    matched[k]={"original_EEA2018_site_key":k18,"original_same_physical_site_match":mode,
                "source_2018_vs_2022_station_offset_m":d,
                "original_2018_source_P_index_tpy":f(old18[k18]["2018_measured_upstream_P_t_per_year"]),
                "original_2018_source_N_index_tpy":f(old18[k18]["2018_measured_upstream_N_t_per_year"]),
                "original_2018_2022_procedure_equal_I2M2":True}
    used18.add(k18)
source_paired_2018=[]
for r in source_rows:
    if r["EEA2022_site_key"] in matched:
        source_paired_2018.append({**r,**matched[r["EEA2022_site_key"]]})
n_2018_site=len(source_paired_2018)
samekey=sum(r["original_same_physical_site_match"]=="original_native_site_key_identical" for r in source_paired_2018)
positive18=[r for r in source_paired_2018 if (r["original_2018_source_P_index_tpy"] or 0)>0]
audit={
 "schema":"GESUI_FR2022_273_TP_EQR_SAME_SITE_ACTUAL_REPORTED_2022_WWTP_P_SOURCE_COVERAGE_AND_EEA2018_CONTINUITY_v1",
 "original_source_273_physical_pairs_SHA256":PAIR_SHA,
 "original_parent_74H7_24H6_source_SHA256":H6_SHA,
 "original_genuine_EEA_GitBlob_SHA1":{k:v[1] for k,v in GITS.items()},
 "actual_2022_EEA_I2M2_Peters_WQ_source_273_original_station_pairs":len(source_rows),
 "source_2022_FR_original_H7_with_matched_water_TP_and_EQR":len({r["original_HYBAS_ID_L07"] for r in source_rows}),
 "original_h6_parents_in_FR2022_TP_EQR_source":len({r["original_H6_parent"] for r in source_rows}),
 "EEA_2022_true_reported_P_SOURCE_POSITIVE_site_pairs":len(positive),
 "EEA_2022_true_reported_P_SOURCE_POSITIVE_unique_H7":len(p_h7),
 "EEA_2022_true_reported_P_SOURCE_POSITIVE_unique_H6":len(p_h6),
 "original_2022_report_source_positive_and_two_derived_pipeline_agreement":dict(reporting),
 "original_2022_source_P_masses_of_positive_pairs_tpy":[r["original_EEA_2022_reported_upstream_P_inventory_tpy"] for r in positive],
 "prior_locked_research_gate_thresholds":{"positive_2022_EEA_WWTP_P_source_sites_min":30,
        "positive_2022_EEA_WWTP_P_source_H7_min":15,"source_positive_H6_parents_min":5},
 "causal_WWTP_P_to_river_TP_or_benthic_EQR_Model_allowed":source_2022_gate,
 "source_gate_status":"ALLOW_EXPLORATORY_ONLY" if source_2022_gate else "STOP_UNDERSUPPORTED_ACTUAL_REPORTED_WWTP_P_AT_TP_EQR_SITES",
 "source_2018_EEA_same_native_I2M2_biology_stations_matched_to_2022_TP_EQR":n_2018_site,
 "source_2018_2022_same_exact_EEA_station_key":samekey,
 "source_2018_2022_same_reach_250m_alternative_station_key":n_2018_site-samekey,
 "source_2018_2022_EEA_continuity_unique_H7":len({r["original_HYBAS_ID_L07"] for r in source_paired_2018}),
 "source_2018_2022_EEA_continuity_unique_H6":len({r["original_H6_parent"] for r in source_paired_2018}),
 "source_2018_upstream_positive_matched_WWTP_P_report_pairs":len(positive18),
 "source_2018_original_2022_EEA_biological_EQR_numeric_NOT_OPENED":True,
 "source_original_Peters_2018_TP_numeric_NOT_YET_OBTAINED":True,
 "source_2022_chemistry_and_biotic_original_scores_NOT_SELECTED_BY_MAGNITUDE":True,
 "does_not_claim_original_source_zero_means_no_treatment_plant_or_pollution":True,
 "two_previous_EEA_WWTP_source_indices_not_independent_new_effluent_measures":True,
 "no_ecological_causal_or_serial_service_claim":True
}
out={
 "schema":audit["schema"],"source_station_pair_registrations_2022":source_rows,
 "source_2018_2022_EEA_same_station_I2M2_biology_geography_only":source_paired_2018,
 "original_2022_EEA_WQ_nutrient_and_biological_response_VALUES_NOT_READ":True
}
b=(json.dumps(out,sort_keys=True,separators=(",",":"))+"\n").encode()
(ROOT/"FR2022_273_TP_EQR_SOURCE_P_EEA2018_STATION_CONTINUITY.json").write_bytes(b)
audit["original_273TP_EQR_source_nutrient_reporting_2018continuity_manifest_SHA256"]=hashlib.sha256(b).hexdigest()
(ROOT/"FR2022_273_TP_EQR_REPORTED_WWTP_P_SUPPORT_QA.json").write_text(json.dumps(audit,indent=2,sort_keys=True))
print("FR2022_273_REPORTED_P_AND_2018_SOURCE_GATE_START",flush=True)
print(json.dumps(audit,sort_keys=True),flush=True)
print("FR2022_273_REPORTED_P_AND_2018_SOURCE_GATE_END",flush=True)
# No hard failure on sparsity; a proper fail-closed STOP is retained in
# explicit original QA, and NO biology/nutrient effect model is run.
