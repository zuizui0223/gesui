#!/usr/bin/env python3
"""SOURCE ONLY USGS NLDI 200m direct river snap + downstream COMID check.

Join 4 fixed successful official USGS API source-location shards (71 original
wetland + 86 real 2021 numeric N/P-reporting permitted POTW features). Never
classify catchment-only, same COMID, absent snap or failed API as true exposure.
Any 2021 wetland plant/chemistry response remains unopened.
"""
from __future__ import annotations
import argparse,hashlib,json,time,urllib.error,urllib.parse,urllib.request
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
GEO=D/"EPA_NWCA2021_128WETLAND_160POTW_EXACT_NEAR_OUTFALL_FEATURE_SOURCE.json"
DMR=D/"EPA_NWCA2021_ORIGINAL_160POTW_2021_DMR_TNTP_SOURCE_FEASIBILITY.json"
FINAL=D/"EPA_NWCA2021_EXACT_344POTW_OUTFALL_2021_TNTP_DMR_NHDPLUS_SOURCE_GATE.json"
BASE="https://api.water.usgs.gov/nldi/linked-data"
HEAD={"Accept":"application/json","User-Agent":"gesui-NWCA2021-USGS-NLDI-downstream-gate-v1"}
MAX_NAV_KM=50
SOURCE_MIN={"wetland_sites":20,"source_features":15,"permits":10,"states":8,"original_far_controls":30}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def hashraw(raw):return hashlib.sha256(raw).hexdigest()

def api(url):
    for i in range(3):
        try:
            req=urllib.request.Request(url,headers=HEAD)
            with urllib.request.urlopen(req,timeout=58) as r:
                raw=r.read()
                data=json.loads(raw)
                if not isinstance(data,dict):raise ValueError("USGS NLDI flowline navigation non-object")
                return data,hashraw(raw),None
        except (urllib.error.HTTPError,urllib.error.URLError,ValueError,TimeoutError) as e:
            if isinstance(e,urllib.error.HTTPError) and e.code in (400,404,422):return None,None,str(e)[:140]
            if i==2:return None,None,(type(e).__name__+":"+str(e))[:140]
            time.sleep(2*(i+1))
    return None,None,"API_FAILURE"
def downstream(comid,mode):
    url=f"{BASE}/comid/{comid}/navigation/{mode}/flowlines?"+urllib.parse.urlencode(
        {"distance":MAX_NAV_KM,"f":"json","excludeGeometry":"true"})
    val,proof,err=api(url)
    if val is None:return {"comid":comid,"mode":mode,"api_error":err,"api_response_sha256":None,
                        "downstream_flowline_comids":[],"navigation_url":url}
    ids=set()
    for row in val.get("features",[]):
        v=(row.get("properties") or {}).get("nhdplus_comid",row.get("id"))
        if v is not None and str(v).isdigit():ids.add(str(v))
    return {"comid":comid,"mode":mode,"api_error":None,
            "api_response_sha256":proof,"downstream_flowline_comids":sorted(ids),
            "navigation_url":url}
def count_and_sources(shards,root):
    original=json.loads((root/GEO).read_text())
    dmr=json.loads((root/DMR).read_text())
    source=json.loads((root/FINAL).read_text())
    if source["original_pre2021_DMR_exact344geo_features_SHA256"]!=sha(root/GEO) or source["original_FY2021andFY2022_source_DMR_join_SHA256"]!=sha(root/DMR):
        raise ValueError("Original 2021 true measured POTW report/source IDs changed")
    raw={}
    shardreceipts={}
    for p in shards:
        data=json.loads(p.read_text())
        shard=data["original_source_shard_id"]
        if shard in shardreceipts:raise ValueError("Duplicate USGS NLDI official source shard")
        if data["schema"]!="gesui_NWCA2021_71wetland_86POTW_actualDMR_USGS_NLDI_comid200m_fourShard_sourceonly_v01":
            raise ValueError("Unexpected source indexing schema")
        if data["source_original_exact344_geo_SHA256"]!=sha(root/GEO) or data["source_original_2021_DMR160_FY2021FY2022_join_SHA256"]!=sha(root/DMR):
            raise ValueError("NLDI source index not from original 2021 DMR true feature cohort")
        rows=data["source_original_USGS_point_source_replies"]
        if len(rows) not in (39,40):raise ValueError("Original 157 source sample not exactly 4 fixed shards")
        for row in rows:
            key=row["source_hashed_point_id"]
            if key in raw:raise ValueError("Duplicate 2021 measured source wetland or POTW feature")
            raw[key]=row
        shardreceipts[str(shard)]=sha(p)
    if set(shardreceipts)!={"0","1","2","3"} or len(raw)!=157:
        raise ValueError("Original 71 wetland+86 true N/P source NLDI index not complete")
    if Counter(row["original_source_type"] for row in raw.values())!=Counter({"NWCA_wetland":71,"POTW_outfall":86}):
        raise ValueError("Original source 71 wetland /86 outfall not preserved")
    measured={r["NPDES_exact_permit_feature_SHA256"]
              for r in dmr["source_original_2021_exact_permit_feature_DMR_reporting_summary"]
              if r["reported_original_numeric_2021_TN00600_period_end_dates"]
                 or r["reported_original_numeric_2021_TP00665_period_end_dates"]}
    included=[x for x in original["original_source_case_site_to_exact_permit_feature_sha256_and_distance"]
              if x["NPDES_exact_permit_feature_SHA256"] in measured]
    if set(x["NWCA2021_site_visit_SHA256"] for x in included)!={k for k,r in raw.items() if r["original_source_type"]=="NWCA_wetland"}:
        raise ValueError("Original wetland DMR source eligibility and USGS indexing differ")
    if set(x["NPDES_exact_permit_feature_SHA256"] for x in included)!={k for k,r in raw.items() if r["original_source_type"]=="POTW_outfall"}:
        raise ValueError("Original POTW 2021 DMR source eligibility and USGS indexing differ")
    summary={
        "source_157_original_GBIF_free_USGS_NLDI_points":len(raw),
        "source_71_WETLAND_86_OUTFALL_original_NP2021":True,
        "USGS_position_API_failures":sum(x["NLDI_position_HTTP_status"]!=200 for x in raw.values()),
        "USGS_hydrolocation_API_failures":sum(x["NLDI_hydrolocation_HTTP_status"]!=200 for x in raw.values()),
        "USGS_catchment_found_wetlands":sum(x["original_source_type"]=="NWCA_wetland" and x["official_USGS_NLDI_catchment_comid"] is not None for x in raw.values()),
        "USGS_catchment_found_outfalls":sum(x["original_source_type"]=="POTW_outfall" and x["official_USGS_NLDI_catchment_comid"] is not None for x in raw.values()),
        "USGS_200m_direct_candidate_wetlands":sum(x["original_source_type"]=="NWCA_wetland" and x["source_strict_200m_DIRECT_stream_candidate"] for x in raw.values()),
        "USGS_200m_direct_candidate_POTW":sum(x["original_source_type"]=="POTW_outfall" and x["source_strict_200m_DIRECT_stream_candidate"] for x in raw.values()),
        "USGS_500m_model_only_wetlands":sum(x["original_source_type"]=="NWCA_wetland" and x["source_500m_modelled_sensitivity_only"] for x in raw.values()),
        "USGS_500m_model_only_POTW":sum(x["original_source_type"]=="POTW_outfall" and x["source_500m_modelled_sensitivity_only"] for x in raw.values()),
    }
    strict=[]
    same=[]
    possible_500m=[]
    for pair in included:
        r=raw[pair["NWCA2021_site_visit_SHA256"]]
        s=raw[pair["NPDES_exact_permit_feature_SHA256"]]
        if r["source_500m_modelled_sensitivity_only"] and s["source_500m_modelled_sensitivity_only"]:
            possible_500m.append(pair)
        if not(r["source_strict_200m_DIRECT_stream_candidate"] and s["source_strict_200m_DIRECT_stream_candidate"]):continue
        if r["official_USGS_NLDI_hydrolocation_comid"]==s["official_USGS_NLDI_hydrolocation_comid"]:
            same.append(pair)
        else:
            strict.append({
                "site_sha256":pair["NWCA2021_site_visit_SHA256"],
                "feature_sha256":pair["NPDES_exact_permit_feature_SHA256"],
                "permit_sha256":pair["NPDES_registered_permit_SHA256"],
                "original_state":pair["source_original_US_state"],
                "source_comid":s["official_USGS_NLDI_hydrolocation_comid"],
                "wetland_comid":r["official_USGS_NLDI_hydrolocation_comid"],
                "source_direct_snap_m":s["original_distance_to_hydrolocation_m"],
                "wetland_direct_snap_m":r["original_distance_to_hydrolocation_m"],
                "straight_geodesic_candidate_distance_km":pair["source_geodesic_outfall_to_wetland_km"],
            })
    sources={p["feature_sha256"] for p in strict}
    sites={p["site_sha256"] for p in strict}
    permits={p["permit_sha256"] for p in strict}
    states={p["original_state"] for p in strict}
    original_geo=json.loads((root/D.joinpath("EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json")).read_text())
    far=original_geo["predeclared_nearest_sewage_outfall_geodesic_source_controls"]["strict"]["matched_far_sites_from_same_state_and_ecoregion_as_a_near_site_n"]
    pass_prefilter=(len(sites)>=SOURCE_MIN["wetland_sites"] and
                    len(sources)>=SOURCE_MIN["source_features"] and
                    len(permits)>=SOURCE_MIN["permits"] and
                    len(states)>=SOURCE_MIN["states"] and
                    far>=SOURCE_MIN["original_far_controls"])
    summary.update({
        "source_exact_original_SiteFeature_2021_DMR_nutrient_candidates_n":len(included),
        "source_site_feature_both200m_same_comid_direction_unknown_n":len(same),
        "source_site_feature_both200m_DIFFERENT_comid_to_test_downstream_n":len(strict),
        "source_200m_candidate_sites_different_comids_n":len(sites),
        "source_200m_candidate_outfalls_different_comids_n":len(sources),
        "source_200m_candidate_NPDES_permits_different_comids_n":len(permits),
        "source_200m_candidate_states_different_comids_n":len(states),
        "source_site_feature_both500m_modelled_not_strict_n":len(possible_500m),
        "source_original_same_state_ecoregion_far_controls_n":far,
        "pre_navigation_direct_candidate_SUPPORT_GATE_pass":pass_prefilter
    })
    return raw,strict,summary,shardreceipts

def run(root,files):
    raw,strict,summary,receipts=count_and_sources(files,root)
    navigation=[]
    source_calls=defaultdict(dict)
    if summary["pre_navigation_direct_candidate_SUPPORT_GATE_pass"]:
        for i,comid in enumerate(sorted({p["source_comid"] for p in strict}),1):
            for mode in ("DM","DD"):
                result=downstream(comid,mode)
                source_calls[comid][mode]=result
                time.sleep(.6)
            if i%10==0:print("USGS_NLDI_2021_POTW_ROUTE",i,"/",len({p["source_comid"] for p in strict}),flush=True)
        for pair in strict:
            rec={"site_sha256":pair["site_sha256"],"feature_sha256":pair["feature_sha256"],
                 "permit_sha256":pair["permit_sha256"],"original_state":pair["original_state"],
                 "source_comid":pair["source_comid"],"wetland_comid":pair["wetland_comid"],
                 "DM_directed_different_comid":pair["wetland_comid"] in source_calls[pair["source_comid"]]["DM"]["downstream_flowline_comids"],
                 "DD_diversion_sensitivity":pair["wetland_comid"] in source_calls[pair["source_comid"]]["DD"]["downstream_flowline_comids"],
                 "original_NLDI_route_query_bounded_50networkkm_but_actual_distance_unknown":True}
            navigation.append(rec)
    confirmed=[p for p in navigation if p["DM_directed_different_comid"]]
    recsite={p["site_sha256"] for p in confirmed}
    feature={p["feature_sha256"] for p in confirmed}
    permit={p["permit_sha256"] for p in confirmed}
    state={p["original_state"] for p in confirmed}
    decision=(len(recsite)>=SOURCE_MIN["wetland_sites"] and len(feature)>=SOURCE_MIN["source_features"] and
              len(permit)>=SOURCE_MIN["permits"] and len(state)>=SOURCE_MIN["states"] and
              summary["source_original_same_state_ecoregion_far_controls_n"]>=SOURCE_MIN["original_far_controls"])
    return {
        "schema":"gesui_NWCA2021_USGS_NLDI_NHDPlusV2_2021DMR_71wetland86POTW_strict200m_DM_DD_source_GATE_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "source_original_USGS_NLDI_index_shards_SHA256":receipts,
        "source_original_344_EPA_registered_POTW_feature_geo_SHA256":sha(root/GEO),
        "source_original_2021_calendaryear_NP_DMR_status_SHA256":sha(root/DMR),
        "source_USGS_NLDI_official_API":"https://api.water.usgs.gov/nldi/linked-data",
        "source_NHDPlus_version":"NHDPlusV2 (USGS NLDI)",
        "precommitted_direct_stream_both_original_sites_200m_and_diff_comid":True,
        "precommitted_downstream_navigation_bound_km":MAX_NAV_KM,
        "original_2021_DMR_source_2026_registered_geo_reference_counts":summary,
        "official_USGS_NLDI_direct_stream_candidate_DifferentCOMID_source_pairs_only":strict,
        "USGS_source_comid_navigation_response_receipts_without_full_geography":[
          {"source_comid":comid,"mode":mode,"nav_response_sha256":v["api_response_sha256"],
           "navigation_url":v["navigation_url"],"api_error":v["api_error"],
           "flowline_reaches_count":len(v["downstream_flowline_comids"]),
           "flowlines_SHA256":hashraw(json.dumps(v["downstream_flowline_comids"],sort_keys=True).encode())}
          for comid,d in source_calls.items() for mode,v in d.items()],
        "USGS_strict200m_directed_outfall_to_wetland_potential_matches":navigation,
        "USGS_downstream_DM_navigation_performed":bool(summary["pre_navigation_direct_candidate_SUPPORT_GATE_pass"]),
        "USGS_downstream_DM_50km_outcome_source":{
          "status":"EXECUTED_NLDI_DM" if summary["pre_navigation_direct_candidate_SUPPORT_GATE_pass"] else "NOT_RUN_ORIGINAL_200M_SOURCE_SUPPORT_GATE_FAILED",
          "potential_directed_source_pairs":len(confirmed) if summary["pre_navigation_direct_candidate_SUPPORT_GATE_pass"] else None,
          "potential_directed_unique_NWCA_wetlands":len(recsite) if summary["pre_navigation_direct_candidate_SUPPORT_GATE_pass"] else None,
          "potential_directed_unique_2021_measured_POTW_features":len(feature) if summary["pre_navigation_direct_candidate_SUPPORT_GATE_pass"] else None,
          "potential_directed_unique_NPDES_administrative_permits":len(permit) if summary["pre_navigation_direct_candidate_SUPPORT_GATE_pass"] else None,
          "potential_directed_original_US_states":len(state) if summary["pre_navigation_direct_candidate_SUPPORT_GATE_pass"] else None,
          "NOT_run_is_not_a_zero_downstream_match":True,
        },
        "precommitted_original_actual_DM_FLOW_SOURCE_STUDY_GATE":{
          **SOURCE_MIN,
          "only_if_source_200m_direct_and_DifferentCOMID_and_DownstreamDM_reached":True,
          "sufficient_USGS_NLDI_source_supported_NWCA2021_wetlands_for_future_observational_vegetation_test":bool(decision)
        },
        "official_NHDPlusV2_DD_navigation_diversion_ONLY_sensitivity_not_rescue":True,
        "source_wetland_within200m_river_DOES_NOT_prove_wetland_hydrologic_inflow":True,
        "all_original_NWCA2021_H_ALL_FQAI_or_TNTP_plant_biology_NOT_ACCESSED":True,
        "NO_sewage_dose_N_P_flow_or_nutrient_mass_effect_computed":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--shard",action="append",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if len(a.shard)!=4:raise ValueError("Exactly four original USGS source point shards required")
    root=Path(__file__).resolve().parents[1]
    r=run(root,a.shard)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"position":r["original_2021_DMR_source_2026_registered_geo_reference_counts"],
                      "strict_downstream":r["USGS_downstream_DM_50km_outcome_source"],
                      "future_plant_support":r["precommitted_original_actual_DM_FLOW_SOURCE_STUDY_GATE"]},indent=2))

if __name__=="__main__":main()
