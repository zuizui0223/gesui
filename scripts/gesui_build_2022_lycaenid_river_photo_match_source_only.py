#!/usr/bin/env python3
"""Pre-photo source-only matched 2022 aquatic river receptor panel.

Source is the original NO-NUMERIC-EQR registry from gesui PR14, pinned
by Git blob SHA. Match same country + 3degree region, then hydrological
flow and geographic distance. Do not open any butterfly/photo outcomes.
A zero *registered* WWTP source is NOT a sewage-free location.
"""
import argparse,csv,hashlib,json,math
from pathlib import Path
from collections import Counter,defaultdict

EXPECTED_BLOB="16f6b531e29ed09f070f9fb86f7d8feb899108ed"
YEAR="2022"
BQE="EEA_13-01-4"
MAX_PAIRS=24
RADIUS_KM=1.0
REACH_SPACING_KM=3.0

def source_hash(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()

def distance(a,b):
    ra=math.pi/180
    dlat=(float(a["lat"])-float(b["lat"]))*ra
    dlon=(float(a["lon"])-float(b["lon"]))*ra
    q=math.sin(dlat/2)**2+math.cos(float(a["lat"])*ra)*math.cos(float(b["lat"])*ra)*math.sin(dlon/2)**2
    return 6371*2*math.atan2(math.sqrt(q),math.sqrt(max(0,1-q)))

def num(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except (TypeError,ValueError):return None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    original=a.source.read_bytes()
    if source_hash(original)!=EXPECTED_BLOB:
        raise ValueError("2022 original source-only registry Git blob drift")
    records=list(csv.DictReader(original.decode("utf-8-sig").splitlines()))
    if not records:raise ValueError("Empty source registry")
    cols=list(records[0])
    # Do not read or use numeric EQR. Reject any outcome-bearing registry.
    forbidden=["normalised_EQR","numeric_EQR","resultEQRValue","resultNormalisedEQRValue",
               "resultEcologicalStatusClassValue"]
    if any(x in cols for x in forbidden):
        raise ValueError("Outcome-containing source registry rejected")
    fields=["site_key","countryCode","lat","lon","region_3deg_block",
            "observedPropertyDeterminandBiologyEQRCode","DIS_AV_CMS","UPLAND_SKM"]
    if any(k not in cols for k in fields):
        raise ValueError("Original registry schema lacks source-only fields: "+str(cols))
    source_col=next((k for k in ("upstream_2022_unique_EEA_both_load_WWTP_count",
                               "upstream_2022_EEA_both_load_WWTP_count")
                     if k in cols),None)
    if source_col is None:raise ValueError("No actual hydrologically upstream registered WWTP source count")
    rows=[]
    for x in records:
        if x["observedPropertyDeterminandBiologyEQRCode"]!=BQE:continue
        if "phenomenonTimeReferenceYear" in x and x["phenomenonTimeReferenceYear"]!=YEAR:continue
        if any(num(x.get(k)) is None for k in ("lat","lon","DIS_AV_CMS","UPLAND_SKM",source_col)):continue
        if not(-12<=float(x["lon"])<=26 and 44<=float(x["lat"])<=64):continue
        if not x.get("site_key"):continue
        rows.append(x)
    dedup={}
    for x in rows:
        if x["site_key"] in dedup:
            raise ValueError("Duplicate same BQE source-only biological site")
        dedup[x["site_key"]]=x
    pos=[x for x in rows if num(x[source_col])>0]
    neg=[x for x in rows if num(x[source_col])==0]
    if len(pos)<10 or len(neg)<10:raise ValueError("Insufficient source-label support even before controls")

    # Source-only: hash order independent from outcome, record/date, and
    # existing wastewater-biology results. Same country×3° is mandatory.
    pos.sort(key=lambda x:hashlib.sha256((x["site_key"]+"|2026-11-08").encode()).hexdigest())
    controls=defaultdict(list)
    for x in neg:controls[(x["countryCode"],x["region_3deg_block"])].append(x)
    used_sites=[];paired=[]
    for exposed in pos:
        if any(distance(exposed,x)<REACH_SPACING_KM for x in used_sites):continue
        choice=[]
        for x in controls[(exposed["countryCode"],exposed["region_3deg_block"])]:
            if x["site_key"]==exposed["site_key"]:continue
            d=distance(exposed,x)
            if d<REACH_SPACING_KM or d>120:continue
            if any(distance(x,y)<REACH_SPACING_KM for y in used_sites):continue
            if any(distance(x,e)<REACH_SPACING_KM for e in [exposed]):continue
            lq=abs(math.log1p(float(exposed["DIS_AV_CMS"]))-math.log1p(float(x["DIS_AV_CMS"])))
            la=abs(math.log1p(float(exposed["UPLAND_SKM"]))-math.log1p(float(x["UPLAND_SKM"])))
            # tie-breaking source-only and deterministic
            score=0.05*d+lq+la
            choice.append((score,d,x["site_key"],x))
        if not choice:continue
        choice.sort(key=lambda z:z[:3])
        score,d,_,control=choice[0]
        id=len(paired)+1
        paired.append({"pair_id":id,"country":exposed["countryCode"],
                       "region":exposed["region_3deg_block"],"pair_km":d,
                       "match_score":score,"exposed":exposed,"control":control})
        used_sites.extend([exposed,control])
        if len(paired)>=MAX_PAIRS:break
    out=[]
    for p in paired:
        for label,x in (("upstream_registered_WWTP",p["exposed"]),("no_registered_qualifying_upstream_WWTP",p["control"])):
            out.append({
              "pair_id":p["pair_id"],"status":label,"site_key":x["site_key"],
              "countryCode":x["countryCode"],"region_3deg_block":x["region_3deg_block"],
              "lat":x["lat"],"lon":x["lon"],
              "DIS_AV_CMS":x["DIS_AV_CMS"],"UPLAND_SKM":x["UPLAND_SKM"],
              "upstream_wwtp_count":x[source_col],
              "matching_pair_distance_km":p["pair_km"],
              "api_search_radius_km":RADIUS_KM,
            })
    a.out.mkdir(parents=True,exist_ok=True)
    path=a.out/"2022_PREPHOTO_EEA_BQE_SOURCE_MATCHED_RIVER_SITE_PAIRS.csv"
    with path.open("w",newline="") as f:
        if out:
            writer=csv.DictWriter(f,fieldnames=list(out[0]));writer.writeheader();writer.writerows(out)
    report={
      "schema":"gesui_Lycaenidae_2022_original_EEA_source_only_matched_river_receptors_v01",
      "original_git_blob_sha1":EXPECTED_BLOB,
      "original_registry_columns":cols,
      "original_registry_total_BQE_elements":len(records),
      "benthic_selected_valid_sites":len(rows),
      "registered_source_positive":len(pos),
      "registered_qualifying_source_zero":len(neg),
      "country_3deg_matched_independent_pairs":len(paired),
      "paired_country_counts":dict(Counter(p["country"] for p in paired)),
      "all_pair_distances_km":[p["pair_km"] for p in paired],
      "paired_sample_sha256":hashlib.sha256(path.read_bytes()).hexdigest() if out else None,
      "paired_sample_CSV":path.name if out else None,
      "no_butterfly_or_other_ecological_outcomes_read":True,
      "source_zero_is_not_no_wastewater":True,
      "source_selection_might_be_biased_by_river_quality_monitoring_sites":True,
      "no_population_abundance_or_puddling_claim":True,
      "next_step":"Count official public 2022 iNaturalist butterfly photos in 1-km radius per fixed source-only matched river monitoring point; needs pair and photo effort quality gate BEFORE any effect inference",
    }
    (a.out/"2022_PREPHOTO_SOURCE_MATCH_QC.json").write_text(
       json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_LYCAENIDAE_EEA_2022_PREPHOTO_RIVER_PAIR_QC",
          json.dumps({k:report[k] for k in ("benthic_selected_valid_sites",
                    "registered_source_positive","registered_qualifying_source_zero",
                    "country_3deg_matched_independent_pairs","paired_country_counts")}),flush=True)
if __name__=="__main__":main()
