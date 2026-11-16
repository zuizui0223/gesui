#!/usr/bin/env python3
"""Source-only audited host moisture (EIV MOIST) already in original 303-butterfly roster.

Verify whether an ecologically defensible H7 host-moisture composition model
can even be defined: strict ACCEPTED/no used-alias conflict, N and F both numeric,
identical species roster and original 454 H7 x period geography.

CRITICAL: host_moisture_mean and host_nutrient_mean are already columns but
the roster does NOT prove they were computed using the same scored HOST species
for each butterfly. Therefore don't fit an N×F, wetland wastewater response,
or aquatic↔butterfly causal model until original per-host paired scores
and measured same site-period water-side source support are verified.
"""
from __future__ import annotations
import argparse,csv,collections,hashlib,json,math
from pathlib import Path

ORIGINAL_ROSTER=Path("results/20261115_sync/GESUI_303_HOST_NUTRIENT_SCORES.csv")
ORIGINAL_PANEL=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")
ROSTER_SHA256="00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb"
PANEL_SHA256="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
SCHEMA="GESUI_20261117_ORIGINAL_303_BUTTERFLY_LARVAL_HOST_MOISTURE_NUTRIENT_SOURCE_ONLY_V1"
FAMILIES={"Lycaenidae","Nymphalidae","Hesperiidae","Pieridae","Papilionidae"}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p,expected):
    if sha(p)!=expected:raise RuntimeError("Original source immutable SHA256 mismatch: "+str(p))
    with p.open(newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))
def numeric(row,col):
    val=str(row.get(col,"")).strip()
    if not val:return None
    try:x=float(val)
    except (TypeError,ValueError):return None
    return x if math.isfinite(x) else None
def ranked(values):
    idx=sorted(range(len(values)),key=values.__getitem__)
    rank=[0.]*len(values)
    i=0
    while i<len(values):
        j=i+1
        while j<len(values) and values[idx[j]]==values[idx[i]]:j+=1
        v=(i+1+j)/2.
        for k in idx[i:j]:rank[k]=v
        i=j
    return rank
def corr(x,y):
    if len(x)<4:raise RuntimeError("Not enough original source numeric paired butterfly scores")
    meanx=sum(x)/len(x);meany=sum(y)/len(y)
    xx=sum((v-meanx)**2 for v in x)
    yy=sum((v-meany)**2 for v in y)
    if xx<1e-12 or yy<1e-12:raise RuntimeError("Degenerate source scored traits")
    return sum((a-meanx)*(b-meany) for a,b in zip(x,y))/math.sqrt(xx*yy)

def main(args):
    args.output.mkdir(parents=True,exist_ok=True)
    a=load(ORIGINAL_ROSTER,ROSTER_SHA256)
    p=load(ORIGINAL_PANEL,PANEL_SHA256)
    if len(a)!=303 or len({q["butterflyKey"] for q in a})!=303 or {q["family"] for q in a}!=FAMILIES:
        raise RuntimeError("Original 303 speciesKey and 5 family scores source drift")
    if len(p)!=454 or len({(q["h"],q["period"]) for q in p})!=454:
        raise RuntimeError("Original same H7×period basin source 454 changed")
    if len({q["h"] for q in p})!=232 or len({q["h6"] for q in p})!=90 or len({q["country"] for q in p})!=15:
        raise RuntimeError("Original H7/H6 country source geography changed")
    if {q["period"] for q in p}!={"2017-2019","2020-2022"}:
        raise RuntimeError("Original two source GBIF periods changed")
    actual=[q for q in a if q.get("taxon_status")=="ACCEPTED" and
             q.get("actual_used_alias_conflict","").lower() not in ("true","1","yes")]
    N=[q for q in actual if numeric(q,"host_nutrient_mean") is not None]
    F=[q for q in actual if numeric(q,"host_moisture_mean") is not None]
    pairs=[q for q in actual if numeric(q,"host_nutrient_mean") is not None and numeric(q,"host_moisture_mean") is not None]
    if len(actual)!=290 or len(N)!=186 or len(F)!=184 or len(pairs)!=180:
        raise RuntimeError("Original 303 butterfly strict-source N and MOIST support drift")
    fam={f:{"original_strict_accepted_both_numeric":sum(q["family"]==f for q in pairs),
            "original_strict_accepted_moisture_numeric":sum(q["family"]==f for q in F),
            "original_exploratory_MOIST>=7_butterfly_species":sum(q["family"]==f and numeric(q,"host_moisture_mean")>=7 for q in pairs)}
         for f in sorted(FAMILIES)}
    if {k:v["original_strict_accepted_both_numeric"] for k,v in fam.items()}!={
         "Hesperiidae":23,"Lycaenidae":44,"Nymphalidae":85,"Papilionidae":10,"Pieridae":18}:
         raise RuntimeError("Original butterfly family N/MOIST strict source support drift")
    fm=[numeric(q,"host_moisture_mean") for q in pairs]
    fn=[numeric(q,"host_nutrient_mean") for q in pairs]
    lo=[q for q in pairs if numeric(q,"host_moisture_mean")>=7]
    if len(lo)!=14 or sum(q["family"]=="Lycaenidae" for q in lo)!=2:
        raise RuntimeError("Original preliminary MOIST >=7 rare subgroup source drift")
    # The N host-completeness source denominator is NOT a valid completeness
    # denominator for the host MOIST trait without original host-by-trait data.
    audit={
        "schema":SCHEMA,
        "status":"PASS_303_SPECIES_SCORE_SOURCE_MOISTURE_EXISTS_BUT_STOP_UNPAIRED_ORIGINAL_HOST_TRAITS",
        "original_Gesui_butterfly_score_source_SHA256":ROSTER_SHA256,
        "original_EU_H7_butterfly_HNS_panel_SHA256":PANEL_SHA256,
        "original_full_butterfly_speciesKey_candidate_n":len(a),
        "original_true_H7_period_n":len(p),
        "original_distinct_H7":len({q["h"] for q in p}),
        "original_distinct_parent_H6":len({q["h6"] for q in p}),
        "original_country_count":len({q["country"] for q in p}),
        "original_true_H7_periods":sorted({q["period"] for q in p}),
        "source_strict_accepted_alias_unconflicted_species":len(actual),
        "source_strict_numeric_host_N_score_species":len(N),
        "source_strict_numeric_host_MOIST_score_species":len(F),
        "source_strict_numeric_N_and_MOIST_at_butterfly_level":len(pairs),
        "source_numeric_N_MOIST_same_butterfly_species_pearson_descriptive":corr(fn,fm),
        "source_numeric_N_MOIST_same_butterfly_species_spearman_descriptive":corr(ranked(fn),ranked(fm)),
        "source_per_butterfly_family_N_MOIST_numeric":fam,
        "exploratory_host_MOIST_ge7_butterfly_species":len(lo),
        "exploratory_host_MOIST_ge7_Lycaenidae_species":sum(q["family"]=="Lycaenidae" for q in lo),
        "exploratory_ge7_binary_is_NOT_a_frozen_eco_model":True,
        "existing_trait_coverage_describes_nutrient_host_matching_not_guaranteed_joint_N_MOIST_host_coverage":True,
        "different_host_subsets_MAY_underlie_mean_N_vs_mean_F_in_same_butterfly":True,
        "source_required_next_step":"restore original per-larval-host EIV NUTRI and MOIST, score both on identical accepted host plants, then reconstruct original 2808 genuine H7×period source speciesKey butterfly records",
        "2022_actual_receiving_river_TP_benthic_EQR_same_H7_and_2022_BMS_NOT_linked":True,
        "no_NUTRI_MOIST_species_ecological_regression_estimated":True,
        "no_point_source_to_water_to_plant_to_butterfly_causal_path_fitted":True,
        "source_EIV_indicator_not_leafN_or_water_NP_delivered_or_adult_puddling":True
    }
    (args.output/"GESUI_303_ORIGINAL_HOST_MOIST_NUTRI_SOURCE_GATE.json").write_text(
        json.dumps(audit,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_ORIGINAL_303_HOST_MOISTURE_AUDIT="+json.dumps({
        "status":audit["status"],"source_all_butterflies":len(a),
        "source_strict_moisture_species":len(F),
        "source_strict_nutrient_species":len(N),
        "source_strict_both_numeric":len(pairs),
        "moisture_atleast7_species":len(lo),
        "moisture_atleast7_Lycaenidae":audit["exploratory_host_MOIST_ge7_Lycaenidae_species"],
        "butterfly_nutrient_moisture_corr":audit["source_numeric_N_MOIST_same_butterfly_species_pearson_descriptive"],
        "source_454_H7_rows":len(p),
        "same_original_scored_HOST_individuals_per_butterfly_not_proven":True,
        "no_new_ecological_or_sewage_effect_test":True
    },sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    pa=argparse.ArgumentParser()
    pa.add_argument("--output",required=True,type=Path)
    main(pa.parse_args())
