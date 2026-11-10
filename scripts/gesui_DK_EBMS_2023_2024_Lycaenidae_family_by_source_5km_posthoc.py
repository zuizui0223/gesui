#!/usr/bin/env python3
"""POSTHOC Denmark 2023/24 standardized eBMS Lycaenidae WWTP-distance contrasts.

The 5km source sampling gate, blinded GBIF taxonomy, 107-site quality gate,
and 4 independent registered WWTP-near source groups were frozen before
any family effect fits. Not a biochemical dose, population causal test,
or a new 2023/24 treatment. Two years from same sites are NOT replications.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import defaultdict,Counter
from pathlib import Path

ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
SOURCE=ROOT+"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
SAMPLE=ROOT+"DENMARK_EBMS_2023_2024_FAMILY_SITEYEAR_COUNTS_QC_BEFORE_WWTP_CONTRAST.json"
CLUSTER=ROOT+"DENMARK_2023_24_EVENT_SITE_WWTP_CLUSTER_PSEUDOREPLICATION_PREOUTCOME_QA.json"
FAMS=("Lycaenidae","Pieridae","Papilionidae")
YEAR=(2023,2024)
NEAR="near5";REMOTE="remote10"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):
    return json.loads(path.read_text())
def density(item,fam):
    return item["family_individual_counts"][fam]/item["distinct_original_visit_days"]
def logodds(item):
    # Species recording conditional on at least one of three target families.
    if item["target_3family_count"]<=0:return None
    counts=item["family_individual_counts"]
    return math.log((counts["Lycaenidae"]+.5)/(counts["Pieridae"]+counts["Papilionidae"]+.5))
def mean(nums):
    return sum(nums)/len(nums) if nums else None
def group_class(item):
    d=item["source_2022_nearest_qualified_WWTP_km"]
    return NEAR if d<=5 else REMOTE if d>=10 else "intermediate"
def aggregate(items):
    x={"original_repeated_monitor_site_coordinate_rows":len(items),
       "siteyears_with_at_least_one_of_3_families":sum(z["target_3family_count"]>0 for z in items),
       "siteyears_without_any_3_family_counts":sum(z["target_3family_count"]==0 for z in items),
       "original_distinct_visit_days":sum(z["distinct_original_visit_days"] for z in items)}
    x["observed_individualCount_by_family"]={f:sum(z["family_individual_counts"][f] for z in items) for f in FAMS+("Nymphalidae","Hesperiidae")}
    nvis=x["original_distinct_visit_days"]
    x["butterfly_family_individuals_per_original_visit_day"]={
        f:sum(density(z,f) for z in items)/len(items) if items else None for f in FAMS
    }
    x["pooled_individualCount_per_visit_day"]={
        f:x["observed_individualCount_by_family"][f]/nvis if nvis else None for f in FAMS
    }
    logs=[v for z in items if (v:=logodds(z)) is not None]
    x["site_equal_logodds_Lycaenidae_vs_PieridaePapilionidae_given_3family_detected"]=mean(logs)
    l=x["observed_individualCount_by_family"]["Lycaenidae"]
    other=sum(x["observed_individualCount_by_family"][f] for f in ("Pieridae","Papilionidae"))
    x["pooled_logodds_given_positive_total_count"]= math.log((l+.5)/(other+.5)) if l+other>0 else None
    x["zero_target_denominator_is_not_assigned_artificial_logodds_zero"]=True
    return x

def hav(a,b):
    r=math.pi/180
    dlat=(a["lat"]-b["lat"])*r
    dlon=(a["lon"]-b["lon"])*r
    t=math.sin(dlat/2)**2+math.cos(a["lat"]*r)*math.cos(b["lat"]*r)*math.sin(dlon/2)**2
    return 12742*math.asin(min(1,math.sqrt(t)))

def source_wwtp_groups(source,by_hash,panel):
    sourceout=source["original_2022_EEA_DK_qualifying_Fw_outfalls"]
    groups=defaultdict(lambda:{"source_site_ids":set(),"year_rows":defaultdict(list)})
    for item in panel:
        if group_class(item)!=NEAR:continue
        spot=by_hash[item["source_rounded5_site_hash"]]
        matches=sorted([(hav(spot,z),z["outfall_id_SHA256"]) for z in sourceout])
        d,key=matches[0]
        if abs(d-item["source_2022_nearest_qualified_WWTP_km"])>1e-8:
            raise ValueError("New near WWTP cluster identity mismatch source")
        groups[key]["source_site_ids"].add(item["source_rounded5_site_hash"])
        groups[key]["year_rows"][item["site_year"]].append(item)
    return [
        {"outfall_id_sha256":key,
         "n_source_near_site_coordinate_rows":len(vals["source_site_ids"]),
         "year_summary":{str(y):aggregate(vals["year_rows"][y]) for y in YEAR}}
        for key,vals in sorted(groups.items())
    ]

def compare(bycell):
    result=[]
    for cell,year in sorted(bycell):
        grouped=bycell[(cell,year)]
        if not grouped[NEAR] or not grouped[REMOTE]:continue
        near=aggregate(grouped[NEAR]);far=aggregate(grouped[REMOTE])
        o=near["site_equal_logodds_Lycaenidae_vs_PieridaePapilionidae_given_3family_detected"]
        p=far["site_equal_logodds_Lycaenidae_vs_PieridaePapilionidae_given_3family_detected"]
        q=near["pooled_logodds_given_positive_total_count"]
        r=far["pooled_logodds_given_positive_total_count"]
        result.append({
            "region_1degree":cell,"year":year,
            "near_source_5km":near,"remote_10km_or_more":far,
            "effect_lycaenidae_site_equal_family_logodds_near_minus_remote":o-p if o is not None and p is not None else None,
            "effect_lycaenidae_pooled_family_logodds_near_minus_remote":q-r if q is not None and r is not None else None,
            "effect_Lycaenidae_absolute_individualCount_per_visit_near_minus_remote":
                near["pooled_individualCount_per_visit_day"]["Lycaenidae"]-far["pooled_individualCount_per_visit_day"]["Lycaenidae"],
            "within_region_contrast_not_independent_WWTP_causal_effect":True,
            "shared_region_year_composition_qualifies_only_if_target_family_recorded_on_both_sides":
                o is not None and p is not None
        })
    return result

def run(root):
    sp=root/SOURCE;sample=root/SAMPLE;cp=root/CLUSTER
    source=load(sp);counts=load(sample);cluster=load(cp)
    assert counts["pre_effect_model_family_count_observation_gate_pass"] is True
    assert counts["original_taxon_WWTP_species_count_data_quality_QC"]["original_event_to_occurrence_link_fraction"]>=.95
    assert counts["original_taxon_WWTP_species_count_data_quality_QC"]["visit_days_mismatch"]==0
    assert cluster["distinct_nearest_WWTP_outfalls_with_ge5km_near_sites"]==4
    assert source["source_distance_and_yearvisit_eligibility"]["5km"]["butterfly_original_occurrence_extension_may_be_accessed_ONLY_if_true"] is True
    assert source["source_distance_and_yearvisit_eligibility"]["2km"]["butterfly_original_occurrence_extension_may_be_accessed_ONLY_if_true"] is False
    sites=source["original_2023_2024_repeated_event_site_source_metadata"]
    sites={x["site_id_SHA256"]:x for x in sites if x["repeat_2023_2024_dates_gate"] and x["precise_500m_explicit_uncertainty_gate"] and x["rounded5_site_coord_scatter_le2m"]}
    assert len(sites)==107
    panel=counts["site_year_family_counts_original_public_CC_BY_NC_aggregate"]
    assert len(panel)==214 and len({x["source_rounded5_site_hash"] for x in panel})==107
    bycell=defaultdict(lambda:defaultdict(list))
    summary={}
    for x in panel:
        assert x["source_rounded5_site_hash"] in sites
        assert x["site_year"] in YEAR
        st=sites[x["source_rounded5_site_hash"]]
        assert x["original_geo_1deg_cell"]==st["country_1deg_lat_lon_cell"]
        assert abs(x["source_2022_nearest_qualified_WWTP_km"]-st["nearest_qualified_EEA2022_DK_FW_outfall_km"])<1e-8
        bycell[(x["original_geo_1deg_cell"],x["site_year"])][group_class(x)].append(x)
    for year in YEAR:
        within=[x for x in panel if x["site_year"]==year]
        summary[str(year)]={name:aggregate([x for x in within if group_class(x)==name]) for name in (NEAR,REMOTE,"intermediate")}
    contrasts=compare(bycell)
    common=source["source_distance_and_yearvisit_eligibility"]["5km"]["country_1deg_cells_with_both_near_and_remote"]
    assert common==3 and len(contrasts)==6
    year_summary={}
    for y in YEAR:
        all_recs=[x for x in contrasts if x["year"]==y]
        valid=[x for x in all_recs if x["shared_region_year_composition_qualifies_only_if_target_family_recorded_on_both_sides"]]
        slopes=[x["effect_lycaenidae_site_equal_family_logodds_near_minus_remote"] for x in valid]
        pooled=[x["effect_lycaenidae_pooled_family_logodds_near_minus_remote"] for x in valid]
        year_summary[str(y)]={
            "total_source_supported_matched_geographical_cells":len(all_recs),
            "composition_data_eligible_regions":len(valid),
            "zero_target_family_near_side_regions":[x["region_1degree"] for x in all_recs if not x["near_source_5km"]["siteyears_with_at_least_one_of_3_families"]],
            "regionwise_site_equal_Lycaenidae_family_logodds_near_minus_remote":slopes,
            "regionwise_pooled_family_logodds_near_minus_remote":pooled,
            "valid_region_effect_signs":{"positive":sum(z>0 for z in slopes),"negative":sum(z<0 for z in slopes)},
            "equal_region_logodds_mean_ONLY_IF_ALL3_POSITIVE_TARGET_SUPPORTED":mean(slopes) if len(valid)==3 else None,
            "equal_region_pooled_logodds_mean_ONLY_IF_ALL3_POSITIVE_TARGET_SUPPORTED":mean(pooled) if len(valid)==3 else None,
            "leave_one_region_out_range_IF_ALL3_SUPPORTED":[min([(sum(slopes)-x)/2 for x in slopes]),max([(sum(slopes)-x)/2 for x in slopes])] if len(valid)==3 else None,
            "not_a_significant_or_causal_4outfall_effect":True
        }
    groupunits=source_wwtp_groups(source,sites,panel)
    if len(groupunits)!=4:
        raise ValueError("Actual near butterflies must stay grouped under four source outfall IDs")
    out={
        "schema":"gesui_DK_eBMS_2023_2024_source_near5_vs_remote10_standardized_Lycaenidae_posthoc_region_contrast_v01",
        "analysis_posthoc_after_DE_AT_2022_iNat_positive":True,
        "original_frozen_sha256":{
            "source_event_and_2022_EEA_WWTP":sha(sp),
            "blinded_original_GBIF_taxonomy_and_licensed_event_occurrence_familycount_panel":sha(sample),
            "pre_butterfly_4outfall_spatial_pseudoreplication_source_QA":sha(cp)
        },
        "original_species_counts_provenance":"Denmark EBMS v1.4, 2023/24, CC-BY-NC 4.0",
        "original_2022_EEA_annual_NP_reported_FW_outfall_proxy_not_2023_2024_butterfly_nutrient_exposure":True,
        "2023_2024_siteyear_count_data":summary,
        "original_3_shared_1degree_region_year_contrasts":contrasts,
        "year_level_same_region_summary":year_summary,
        "near_5km_effective_independent_outfall_group_2023_2024_summary":groupunits,
        "primary_three_shared_region_x_two_year_result_NOT_inferentially_supported_if_any_cell_year_has_no_three_family_target_count":any(x["composition_data_eligible_regions"]<3 for x in year_summary.values()),
        "near5_sourceWWTP_independent_outfall_units_n":len(groupunits),
        "no_2km_source_near_butterfly_effect_estimated":True,
        "no_true_population_change_or_effluent_causality_or_adult_puddling_observed":True,
        "no_individual_record_or_personal_observer_IDs_republished":True,
        "outcome_decision":"Descriptive evidence only. Family ratios are not defined where neither of three target families is observed, years and sites within 4 source outfalls are correlated, and annual 2022 WWTP load is a static geographic proxy."
    }
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=run(root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print(json.dumps({
        "2023_2024_regional_matched_effect_support":out["year_level_same_region_summary"],
        "near_5km_effective_outfalls":out["near5_sourceWWTP_independent_outfall_units_n"],
        "positive_adult_WWTP_effect_proven":False
    },ensure_ascii=False,indent=2))
if __name__=="__main__":main()
