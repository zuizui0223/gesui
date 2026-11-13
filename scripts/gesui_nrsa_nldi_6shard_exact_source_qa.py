#!/usr/bin/env python3
"""Independent SOURCE-only merger for six true NLDI UT/DM × exact ECHO
external outfall source shards. Not measured effluent, in-stream N/P, biodiversity
or cause. Prevents repeated site/permit pseudo-replication and zero imputation.
"""
from __future__ import annotations
import collections,hashlib,json,os,pathlib,sys
root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_nrsa_nldi_merge"
root.mkdir(parents=True,exist_ok=True)
base=pathlib.Path(sys.argv[1])
fs=sorted(base.rglob("USGS_NLDI_NRSA_POTW_SOURCE_ROUTE_SHARD*.json"))
if len(fs)!=6:raise RuntimeError(f"Require all 6 original source shards, found {len(fs)}")
data=[json.loads(p.read_text()) for p in fs]
if {x["shard"] for x in data}!=set(range(6)):raise RuntimeError("Missing or duplicate native source shard")
if len({x["geo_source_sha256"] for x in data})!=1:raise RuntimeError("Different original geo source versions")
if len({x["original_ECHO_national_zip_sha256"] for x in data})!=1:raise RuntimeError("Different original NPDES official files")
if not all(x["no_DMR_or_biological_numeric_values_read"] for x in data):raise RuntimeError("Source gate not outcome blind")
sites=[r for x in data for r in x["site_route_records"]]
keys=[r["uid"] for r in sites]
if len(set(keys))!=len(keys):raise RuntimeError("Duplicate original visit across source shards")
out={}
for era in ("2018_19","2023_24"):
 rr=[x for x in sites if x["era"]==era]
 ok=[x for x in rr if x["status"]=="NETWORK_SOURCE_API_OK"]
 nonok=[x for x in rr if x["status"]!="NETWORK_SOURCE_API_OK"]
 bytag=collections.Counter()
 U=[]
 D=[]
 same=[]
 for x in ok:
  for src in x.get("matched_NLDI_exact_registered_pairs",[]):
   bytag[src["nldi_route_status"]]+=1
   rec={"uid":x["uid"],"state":x["state"],"huc8":x["huc8"],
        "permit":src["permit"],"feature":src["feature"],
        "nldi_id":src["nldi_exact_id"],"site_comid":x["comid"],
        "near10":src["near10"],"geodesic_km":src["geodesic_km"]}
   if src["nldi_route_status"]=="UPSTREAM_DIFFERENT_REACH":U.append(rec)
   elif src["nldi_route_status"]=="DOWNSTREAM_DIFFERENT_REACH":D.append(rec)
   elif src["nldi_route_status"]=="SAME_REACH_UNRESOLVED_ORDER":same.append(rec)
 unique_sites=set(x["uid"] for x in U)
 unique_features=set((x["permit"],x["feature"]) for x in U)
 unique_permits=set(x["permit"] for x in U)
 # Keep NPDES feature-level sources independent and preserve source register.
 min_dmr_gate= len(unique_sites)>=40 and len(unique_features)>=25 and len({x["huc8"] for x in U})>=10 and len({x["state"] for x in U})>=5
 out[era]={
  "original_eligible_geo+nativeCOMID_site_visits":len(rr),
  "n_UT_and_DM_both_API200":len(ok),
  "n_unavailable_or_error_API_site_visits":len(nonok),
  "api_status_by_kind":dict(collections.Counter(x["status"] for x in rr)),
  "n_distinct_HUC8_queried":len({x["huc8"] for x in rr}),
  "registered_external_outfall_exact_NLDI_route_type_source_pair_counts":dict(bytag),
  "n_verified_different_reach_upstream_registered_outfall_site_visits":len(unique_sites),
  "n_verified_distinct_original_outfall_features":len(unique_features),
  "n_verified_distinct_original_permits":len(unique_permits),
  "n_verified_exposure_HUC8":len({x["huc8"] for x in U}),
  "n_verified_exposure_states":len({x["state"] for x in U}),
  "n_site_visits_with_upstream_source_within10km":len({x["uid"] for x in U if x["near10"]}),
  "n_source_downstream_pairs":len(D),
  "n_source_same_reach_order_unresolved_pairs":len(same),
  "pre_DMR_source_gate_meets_MIN40sites_25features_10HUC8_5states":min_dmr_gate,
  "verified_source_only_upstream_candidate_pairs":U,
  "same_reach_and_downstream_are_not_upstream_DMR_exposure":True,
  "original_registry_zero_not_reliable_historical_no_effluent":True
 }
summ={"schema":"gesui_NRSA_2018_19_2023_24_UPSTREAM_POTW_SOURCE_ONLY_NO_DMR_JOIN_v1",
      "original_NPDES_zip_SHA256":data[0]["original_ECHO_national_zip_sha256"],
      "original_geocandidate_manifest_SHA256":data[0]["geo_source_sha256"],
      "all_shards_count":6,"original_route_sites_total":len(sites),
      "eras":out,
      "source_only_DMR_required_before_any_sewage_NP_biology":True,
      "no_assigned_measured_WWTP_mass_or_concentrations":True,
      "no_instream_biological_wwtp_effects_calculated":True,
      "NLDI_npdES_registered_station_reach_is_NOT_2026_original_ECHO_outfall_bed_measure":"exact permit+feature matched, NHDPlus2 from NLDI still secondary geospatial model",
      "conditional_next":"Only if validated source support passes and historical FY2018/19/20 or FY2023/24/25 actual N/P DMR reports meet independent nonzero, QC and monitoring-unit thresholds should ecological effects be fit."
      }
outpath=root/"NRSA_NLDI_POTW_TRUE_UPSTREAM_EXACT_SOURCE_JOIN_QA.json"
outpath.write_text(json.dumps(summ,sort_keys=True,indent=2))
for era in out:
 v=out[era]
 print("NRSA_POTW_HYDROLOGIC_SOURCE_GATE",era,json.dumps({k:v[k] for k in v if not isinstance(v[k],list)},sort_keys=True))
print("NRSA_POTW_HYDRO_SOURCE_QA_SAVED",str(outpath))
