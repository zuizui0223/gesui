#!/usr/bin/env python3
"""Independent pure-Python OLS arithmetic audit of frozen ESA WorldCover and EBMS.

Different implementation from original numpy.linalg.lstsq, no API or raw
licensed occurrence data. Checks original 6-area source-balanced site sample,
denominator eligibility and each previously fitted regional habitat slope.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import defaultdict,Counter
from pathlib import Path

ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
LAND=ROOT+"DK2025_NEW399_286SITE_WORLD_COVER2021_500M_NO_FAMILY_SOURCE.json"
FAM=ROOT+"DENMARK_2025_NEW399_SITE_LYCAENIDAE_FAMILY_POST_SOURCE_GATE_COUNTS.json"
FIT=ROOT+"DK2025_NEW399_WORLD_COVER2021_LYCAENIDAE_3VS5_WITHIN6_REGION_HABITAT_QC.json"
COV=("built","crop","grass","tree")

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def solve_normal(X,y):
    """Independent pivoted standard-library normal equation, no numpy."""
    n=len(X);p=len(X[0])
    a=[[0.]*(p+1) for _ in range(p)]
    for i in range(n):
        for j in range(p):
            a[j][p]+=X[i][j]*y[i]
            for k in range(p):a[j][k]+=X[i][j]*X[i][k]
    for col in range(p):
        q=max(range(col,p),key=lambda r:abs(a[r][col]))
        if abs(a[q][col])<1e-9:return None
        a[q],a[col]=a[col],a[q]
        div=a[col][col]
        for k in range(col,p+1):a[col][k]/=div
        for j in range(p):
            if j==col:continue
            factor=a[j][col]
            for k in range(col,p+1):a[j][k]-=factor*a[col][k]
    return [row[p] for row in a]

def eligible(row,endpoint):
    c=row["count"]
    if endpoint=="three":
        return sum(c[f] for f in ("Lycaenidae","Pieridae","Papilionidae"))>0
    if endpoint=="five":
        return sum(c[f] for f in ("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae"))>0
    return True

def response(row,endpoint):
    c=row["count"]
    if endpoint=="three":
        return math.log((c["Lycaenidae"]+.5)/(c["Pieridae"]+c["Papilionidae"]+.5))
    if endpoint=="five":
        return math.log((c["Lycaenidae"]+.5)/(c["Pieridae"]+c["Papilionidae"]+c["Nymphalidae"]+c["Hesperiidae"]+.5))
    name={"Lycaenidae":"Lycaenidae","Pieridae":"Pieridae","Nymphalidae":"Nymphalidae","Hesperiidae":"Hesperiidae"}[endpoint]
    return math.log((c[name]+.5)/row["days"])

def fit(rows,endpoint,landcover):
    a=[r for r in rows if eligible(r,endpoint)]
    areas=sorted({r["region"] for r in a})
    if any(len({r["source"] for r in a if r["region"]==z})!=2 for z in areas):return None
    X=[[1.,float(r["source"]=="near5"),*[float(r["region"]==z) for z in areas[1:]],
        math.log(r["days"]),*([r["land"][k] for k in COV] if landcover else [])] for r in a]
    y=[response(r,endpoint) for r in a]
    b=solve_normal(X,y)
    return {"n_sites":len(a),"n_regions":len(areas),
            "original_slope_near_minus_remote":b[1],
            "exp_slope":math.exp(b[1])} if b else None

def check(a,b,what):
    if a is None or b is None:
        if a is not b:raise ValueError("Missing source model comparison "+what)
    elif abs(a-b)>1e-7:raise ValueError(f"Independent model mismatch {what}: {a} vs {b}")

def main():
    q=argparse.ArgumentParser();q.add_argument("--out",type=Path,required=True);args=q.parse_args()
    root=Path(__file__).resolve().parents[1]
    files={k:root/x for k,x in (("land2021",LAND),("family2025",FAM),("original_model",FIT))}
    land=json.loads(files["land2021"].read_text())
    family=json.loads(files["family2025"].read_text())
    report=json.loads(files["original_model"].read_text())
    assert report["source_original_SHA256"]["worldcover2021"]==sha(files["land2021"])
    assert report["source_original_SHA256"]["original_family_source"]==sha(files["family2025"])
    assert len(land["original_site_geo_worldcover_source_counts"])==286
    w={x["source_original_public_site_SHA256"]:x for x in land["original_site_geo_worldcover_source_counts"]}
    f=family["2025_original_399_site_accepted_family_counts"]
    a=[]
    for v in f:
        sid=v["original_public_site_SHA256"]
        if not v["source_gate_2plus_any_Occurrence_dates"] or v["source_band"] not in ("near5","remote10"):continue
        t=w[sid]
        if t["source_QC"]!="COMPLETE":continue
        if t["original_2025_any_published_Occurrence_linked_dates"]!=v["linked_visit_days"]:
            raise ValueError("Original reported-event count differs")
        a.append({"sid":sid,"source":v["source_band"],"region":v["original_geo_1deg_cell"],
                  "days":v["linked_visit_days"],"count":v["family_individualCount"],
                  "land":t["fractions"],"cluster":t["source_near_WWTP_5km_component_SHA256"]})
    groups=defaultdict(set)
    for r in a:groups[r["region"]].add(r["source"])
    common=sorted(k for k,g in groups.items() if len(g)==2)
    if common!=report["within_1degree_comparable_geographic_regions"] or len(common)!=6:
        raise ValueError("Original same six spatial comparisons changed")
    x=[r for r in a if r["region"] in common]
    assert len(x)==139 and Counter(r["source"] for r in x)==Counter({"near5":61,"remote10":78})
    endpoints={"logLycaenidaePerDate":"Lycaenidae","logitThree":"three",
               "logitFive":"five","logPieridaePerDate":"Pieridae",
               "logNymphalidaePerDate":"Nymphalidae","logHesperiidaePerDate":"Hesperiidae"}
    result={}
    for key,endpoint in endpoints.items():
        result[key]={}
        for column,with_land in (("region_plus_logSurveyDays",False),
                                 ("region_plus_logSurveyDays_plus_WorldCover2021",True)):
            value=fit(x,endpoint,with_land)
            orig=report["models"][key][column]
            if value["n_sites"]!=orig["sample_site_n"]:
                raise ValueError("Original response denominator sample changed")
            check(value["original_slope_near_minus_remote"],orig["coefficient_log_near_minus_remote"],key+"/"+column)
            result[key][column]=value
        # True spatial leave-one-area-out model; source similarity only a
        # descriptive robustness check, not independent replications.
        omitted={}
        for region in common:
            reduced=[r for r in x if r["region"]!=region]
            vv=fit(reduced,endpoint,True)
            omitted[region]=vv["exp_slope"] if vv else None
        values=[v for v in omitted.values() if v is not None]
        result[key]["leave_one_1degree_area_out_adjusted_exponentiated_slopes"]=omitted
        result[key]["LOO_adjusted_positive_ratio_count"]=sum(v>1 for v in values)
        result[key]["LOO_valid_region_drops"]=len(values)
        result[key]["LOO_adjusted_ratio_range"]=[min(values),max(values)] if values else None
    near_components=sorted({r["cluster"] for r in x if r["source"]=="near5"})
    family_component={}
    for key,endpoint in endpoints.items():
        vals=[]
        for part in near_components:
            reduced=[r for r in x if not(r["source"]=="near5" and r["cluster"]==part)]
            vv=fit(reduced,endpoint,True)
            if vv:vals.append(vv["exp_slope"])
        family_component[key]={
            "near_5km_spatial_components_attempted":len(near_components),
            "leave_one_near_component_adjusted_valid_n":len(vals),
            "source_component_adjusted_ratio_minmax":[min(vals),max(vals)] if vals else None,
            "source_component_drops_pseudoreplication_not_hydrology_proof":True
        }
    output={
        "schema":"gesui_DK2025_NEW399_WorldCover2021_original_6area_adjusted_family_independent_OLS_QA_v01",
        "source_SHA256":{k:sha(v) for k,v in files.items()},
        "independent_standard_library_normal_equation_not_numpy_lstsq":True,
        "all_six_outcomes_and_both_regions_v_land_model_slopes_exactly_reproduced":True,
        "original_six_geographies":common,
        "source_region_2025_original_public_site_n":len(x),
        "original_2025_worldcover2021_500m_quality_sample":286,
        "independent_family_model_arithmetic":result,
        "near_source_5km_component_leave_one_out":family_component,
        "2021_WorldCover_classes_not_host_species_or_butterfly_dose":True,
        "all_WWTP_proximity_inference_noncausal_POSTHOC":True
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"all_models_exact":output["all_six_outcomes_and_both_regions_v_land_model_slopes_exactly_reproduced"],
        "endpoints":{k:{"region_adj_ratio":v["region_plus_logSurveyDays"]["exp_slope"],
                        "land_adj_ratio":v["region_plus_logSurveyDays_plus_WorldCover2021"]["exp_slope"],
                        "land_adj_region_LOO":v["LOO_adjusted_ratio_range"]} for k,v in result.items()}},indent=2))
if __name__=="__main__":
    main()
