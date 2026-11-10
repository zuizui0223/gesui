#!/usr/bin/env python3
"""Independent offline spatial source unit audit of frozen Denmark EBMS/EEA data.

This reads only 2023-2024 original Event-site metadata and 2022 original WWTP
outfalls; does NOT read any butterfly occurrence, family or count outcomes.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import defaultdict
from pathlib import Path

SOURCE="results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"

def hav(a,b):
    r=math.pi/180
    t=(a["lat"]-b["lat"])*r
    s=(a["lon"]-b["lon"])*r
    q=math.sin(t/2)**2+math.cos(a["lat"]*r)*math.cos(b["lat"]*r)*math.sin(s/2)**2
    return 12742*math.asin(min(1.,math.sqrt(q)))

def clusters(items,radius):
    links=list(range(len(items)))
    def root(i):
        while links[i]!=i:
            links[i]=links[links[i]]
            i=links[i]
        return i
    for i,a in enumerate(items):
        for j in range(i+1,len(items)):
            if hav(a,items[j])<=radius:
                links[root(i)]=root(j)
    groups=defaultdict(list)
    for i,item in enumerate(items):
        groups[root(i)].append(item)
    return list(groups.values())

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    p=root/SOURCE
    original=json.loads(p.read_text())
    assert original["source_distance_and_yearvisit_eligibility"]["5km"]["butterfly_original_occurrence_extension_may_be_accessed_ONLY_if_true"] is True
    assert original["source_distance_and_yearvisit_eligibility"]["2km"]["butterfly_original_occurrence_extension_may_be_accessed_ONLY_if_true"] is False
    assert original["original_EEA2022_outfall_source_QC"]["unique_Danish_qualified_2022_FW_bothNP_outfalls_n"]==203
    sites=original["original_2023_2024_repeated_event_site_source_metadata"]
    good=[x for x in sites if x["repeat_2023_2024_dates_gate"] and x["precise_500m_explicit_uncertainty_gate"] and x["rounded5_site_coord_scatter_le2m"]]
    assert len(good)==107
    near=[x for x in good if x["nearest_qualified_EEA2022_DK_FW_outfall_km"]<=5]
    remote=[x for x in good if x["nearest_qualified_EEA2022_DK_FW_outfall_km"]>=10]
    assert len(near)==23 and len(remote)==38
    source=original["original_2022_EEA_DK_qualifying_Fw_outfalls"]
    byout=defaultdict(list)
    for s in near:
        ranked=sorted(((hav(s,z),z["outfall_id_SHA256"]) for z in source))
        distance,key=ranked[0]
        assert abs(distance-s["nearest_qualified_EEA2022_DK_FW_outfall_km"])<1e-8
        byout[key].append(s)
    out={
        "schema":"gesui_DK_ebms_2023_2024_5km_near_falldown_to_outfall_clusters_source_only_v01",
        "original_event_source_SHA256":hashlib.sha256(p.read_bytes()).hexdigest(),
        "public_DK_reported_2022_annual_P_source_n":len(source),
        "original_year_repeated_event_eligible_sites":len(good),
        "near5_source_proximity_sites":len(near),
        "remote_ge10km_sites":len(remote),
        "distinct_nearest_WWTP_outfalls_with_ge5km_near_sites":len(byout),
        "nearest_WWTP_outfall_near_site_cluster_n":{k:len(v) for k,v in byout.items()},
        "no_butterfly_outcomes_read":True,
        "geo_region_compare_5km":original["source_distance_and_yearvisit_eligibility"]["5km"],
        "radius_cluster_analyses":{},
        "interpretation":"Source screening 5km gate passed but 23 public site-coordinate rows lie near only four distinct WWTP sources; site rows within a transect or close sites are correlated, not independent WWTP treatment replicates. Do not infer family benefit without source-level clustering and habitat chemistry."
    }
    for rad in (0.25,0.5,1.,2.):
        d={}
        for name,items in (("near",near),("remote",remote)):
            cs=clusters(items,rad)
            d[name]={
                "original_site_coordinate_rows":len(items),
                "spatial_connected_components_n":len(cs),
                "component_size_desc":sorted((len(g) for g in cs),reverse=True),
                "independent_unit_warning":True}
        out["radius_cluster_analyses"][f"{rad:.2f}km"]=d
    assert out["radius_cluster_analyses"]["0.50km"]["near"]["spatial_connected_components_n"]==4
    assert out["radius_cluster_analyses"]["0.50km"]["remote"]["spatial_connected_components_n"]==8
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"near_independent_outfall_clusters":len(byout),
                      "near_sites":len(near),"remote_sites":len(remote),
                      "clusters500m":out["radius_cluster_analyses"]["0.50km"]},indent=2))
if __name__=="__main__":main()
