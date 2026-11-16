#!/usr/bin/env python3
"""Figure 4 pre-outcome GBIF coordinate feasibility for 8 frozen EU H7 watersheds.

The authentic 2808 GBIF speciesKey facet rows in the original H7 study have
species-level counts but NOT occurrenceKey/geographic coordinates. Reuse their
exact original, SHA-verified public search query polygons, years and filters
but retrieve at most 300 indexed occurrences per H7-period, ONLY to audit
source field presence and source drift. No EIV NUTRI, HNS, EQR or W numeric
response is ever accessed. Pilot pages do not represent the complete species
population or give river-distance source support.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import random
import time
from collections import Counter,defaultdict
from pathlib import Path
from urllib.parse import urlsplit,parse_qs,urlencode
import requests

PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
FACET_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
FACET_FILE="GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json"
PILOT_H7=8
MAX_PAGE=300
PERIODS=["2017-2019","2020-2022"]
ALLOWED_HOST="api.gbif.org"
USER_AGENT="gesui-20261116-figure4-gbif-original-geometry-coordinate-source-gate/1.0"
FIELD_NAMES=[
    "key","speciesKey","decimalLatitude","decimalLongitude",
    "coordinateUncertaintyInMeters","year","eventDate","datasetKey",
    "basisOfRecord","occurrenceStatus","issues","taxonKey","gbifID","recordedBy",
    "eventID","institutionCode","collectionCode"
]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def one(base,name):
    hits=list(base.rglob(name))
    if len(hits)!=1:raise RuntimeError("Expected one original source "+name+"; got "+str(len(hits)))
    return hits[0]
def target_panel():
    p=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")
    if sha(p)!=PANEL_SHA:raise RuntimeError("Original 454 EEA H7 panel SHA changed")
    with p.open(encoding="utf-8-sig",newline="") as f:
        r=list(csv.DictReader(f))
    if len(r)!=454 or len({(x["h"],x["period"]) for x in r})!=454 or len({x["h6"] for x in r})!=90:
        raise RuntimeError("Genuine H7-period 454/90 source missing")
    info=defaultdict(dict)
    for x in r:
        h=x["h"];per=x["period"]
        info[h][per]={"h6":str(x["h6"]),
                      "country_holdout":str(x["country_holdout"])}
    both={k:v for k,v in info.items() if set(v)==set(PERIODS)
         and v[PERIODS[0]]["h6"]==v[PERIODS[1]]["h6"]}
    if len(both)<200:raise RuntimeError("Not enough original true 2-period H7 IDs for source pilot")
    return both

def freeze_h7(both):
    candidates=defaultdict(list)
    for h,period in both.items():
        h6=period[PERIODS[0]]["h6"]
        c=period[PERIODS[0]]["country_holdout"]
        candidates[c].append((h6,h))
    # Sort actual source country groups by label; take 8 distinct groups with
    # deterministic hash-selected genuine H6 and one genuine H7 inside.
    ordered=[]
    for c,pairs in candidates.items():
        best=min(pairs,key=lambda x:hashlib.sha256(("FIG4_20261116_"+x[0]+"_"+x[1]).encode()).hexdigest())
        ordered.append((c,best[0],best[1]))
    ordered=sorted(ordered,key=lambda x:hashlib.sha256(("FIG4_COUNTRY_"+x[0]).encode()).hexdigest())
    selected=ordered[:PILOT_H7]
    if len(selected)!=PILOT_H7 or len({x[0] for x in selected})!=PILOT_H7:
        raise RuntimeError("Eight distinct original pre-outcome country groups not available")
    if len({x[1] for x in selected})!=PILOT_H7:
        raise RuntimeError("Pilot H6 parents not eight distinct original hydrological regions")
    return selected

def read_original(facets_dir,selected):
    p=one(facets_dir,FACET_FILE)
    if sha(p)!=FACET_SHA:raise RuntimeError("Original 2808 GBIF facet SHA changed")
    j=json.loads(p.read_text())
    if len(j)!=2808 or any(x.get("source_status")!="ok" for x in j):
        raise RuntimeError("Original 2808 source facets missing/status drift")
    sampled={x[2] for x in selected}
    flora={(str(x["HYBAS_ID_L07"]),str(x["period"])):x
           for x in j if x.get("group")=="Plantae" and str(x["HYBAS_ID_L07"]) in sampled}
    if len(flora)!=PILOT_H7*2:
        raise RuntimeError("Original pilot plant H7-period facet source missing")
    fields=set().union(*(set(x) for x in flora.values()))
    if "source_speciesKey_counts" not in fields:
        raise RuntimeError("Source original speciesKey counts absent")
    coord_fields=sorted(fields&set(["decimalLatitude","decimalLongitude",
          "occurrenceKey","gbifID","coordinateUncertaintyInMeters"]))
    if coord_fields:
        raise RuntimeError("Original GBIF facetted source unexpectedly has per-record geodata fields "+repr(coord_fields))
    for x in flora.values():
        if not isinstance(x.get("source_speciesKey_counts"),dict):
            raise RuntimeError("Source already lacks genuine speciesKey counts")
        if not x.get("source_request_receipts"):
            raise RuntimeError("Original GBIF query receipt absent")
    return flora

def source_query(facet):
    receipts=facet["source_request_receipts"]
    actual=[r for r in receipts if r.get("status")=="source_success" and r.get("facetOffset",999)==0]
    if not actual:raise RuntimeError("Original source has no complete 0-offset GBIF facet query")
    a=actual[0]
    url=a["url"]
    parsed=urlsplit(url)
    if parsed.scheme!="https" or parsed.hostname!=ALLOWED_HOST or parsed.path!="/v1/occurrence/search":
        raise RuntimeError("Unexpected original GBIF source URL")
    q=parse_qs(parsed.query)
    for k,expected in {
        "taxonKey":"6","hasCoordinate":"true","hasGeospatialIssue":"false",
        "occurrenceStatus":"PRESENT"}.items():
        values=q.get(k) or q.get(k.lower()) or []
        if len(values)!=1 or values[0].lower()!=expected.lower():
            raise RuntimeError("Original GBIF exact source query field mismatch "+k+" "+repr(values))
    if "geometry" not in q or len(q["geometry"])!=1 or not q["geometry"][0]:
        raise RuntimeError("Source original WKT query polygon absent")
    if len(q.get("year",[]))!=1 or q["year"][0] not in ("2017,2019","2020,2022"):
        raise RuntimeError("Source original year range absent")
    for k in ["facet","facetLimit","facetOffset","limit","offset"]:
        q.pop(k,None)
    q={k:v[0] for k,v in q.items()}
    q["limit"]=str(MAX_PAGE)
    return q

def fetch_one(session,params):
    errors=[]
    for i in range(4):
        try:
            r=session.get("https://api.gbif.org/v1/occurrence/search",params=params,
                 timeout=(15,60),headers={"User-Agent":USER_AGENT,"Accept":"application/json"})
            r.raise_for_status()
            data=r.content
            if len(data)>4000000:raise RuntimeError("GBIF source pilot per-page JSON exceeds 4MB cap")
            obj=r.json()
            if not isinstance(obj,dict) or not isinstance(obj.get("count"),int) or not isinstance(obj.get("results"),list):
                raise RuntimeError("Original GBIF source result schema invalid")
            if len(obj["results"])>MAX_PAGE:raise RuntimeError("GBIF returned >300 pilot records")
            return {"status":"ok","source_total_count":obj["count"],
                    "pilot_page_n":len(obj["results"]),"samples":obj["results"],
                    "response_SHA256":hashlib.sha256(data).hexdigest(),
                    "http":r.status_code}
        except Exception as e:
            errors.append(type(e).__name__+":"+str(e)[:135])
            if i<3:time.sleep([1,3,6][i])
    return {"status":"API_UNAVAILABLE_NOT_ECOLOGICAL_ZERO","error":errors}

def record_stats(results):
    stat=Counter()
    datasets=set()
    species=set()
    fields_union=set()
    for q in results:
        fields_union.update(q)
        for fld in FIELD_NAMES:
            if q.get(fld) not in ("",None,[],{}):stat[fld+"_nonmissing"]+=1
        lat=q.get("decimalLatitude");lon=q.get("decimalLongitude")
        good=False
        try:
            good=(-90<=float(lat)<=90 and -180<=float(lon)<=180)
        except Exception:pass
        if good:stat["valid_WGS84_coordinate"]+=1
        un=q.get("coordinateUncertaintyInMeters")
        if un is not None:
            try:
                if float(un)>=0 and float(un)<=100:stat["uncertainty_le_100m"]+=1
                if float(un)>=0 and float(un)<=250:stat["uncertainty_le_250m"]+=1
            except Exception:pass
        if q.get("datasetKey"):datasets.add(str(q["datasetKey"]))
        if q.get("speciesKey") is not None:species.add(str(q["speciesKey"]))
        if q.get("issues"):stat["any_source_issue_flag"]+=1
    return {"field_presence_counts":dict(stat),
            "all_observation_field_names":sorted(fields_union),
            "source_distinct_dataset_keys":len(datasets),
            "source_distinct_speciesKeys":len(species),
            "source_sample_records_are_not_independent_field_events":True}

def main(a):
    both=target_panel()
    selected=freeze_h7(both)
    flora=read_original(a.facets,selected)
    audits=[]
    with requests.Session() as sess:
        for country,h6,h in selected:
            for per in PERIODS:
                q=source_query(flora[(h,per)])
                if q["year"] != {"2017-2019":"2017,2019","2020-2022":"2020,2022"}[per]:
                    raise RuntimeError("Real GBIF receipt period different from source study period")
                response=fetch_one(sess,q)
                row={"original_H7":h,"original_H6":h6,"majority_country":country,
                     "period":per,"original_H7_plant_facet_source_n":
                         flora[(h,per)]["source_query_occurrence_count"],
                     "original_source_plant_index_speciesKeys":
                         flora[(h,per)]["source_speciesKey_count"],
                     "original_facets_have_no_occurrence_level_coordinates":True,
                     "original_source_exact_filters_reused":True,
                     "GBIF_retrieve_status":response["status"]}
                if response["status"]=="ok":
                    row["current_GBIF_source_query_component_count"]=response["source_total_count"]
                    row["pilot_first_page_record_n"]=response["pilot_page_n"]
                    row["source_response_SHA256"]=response["response_SHA256"]
                    row["source_field_quality"]=record_stats(response["samples"])
                else:
                    row["source_error_categories"]=response["error"]
                audits.append(row)
    if len(audits)!=16:raise RuntimeError("Original pilot requests incomplete")
    allvalid=[x for x in audits if x["GBIF_retrieve_status"]=="ok"]
    source_rows=sum(x["pilot_first_page_record_n"] for x in allvalid)
    coords=sum(x["source_field_quality"]["field_presence_counts"].get("valid_WGS84_coordinate",0) for x in allvalid)
    uncert100=sum(x["source_field_quality"]["field_presence_counts"].get("uncertainty_le_100m",0) for x in allvalid)
    out={"schema":"GESUI_20261116_GBIF_ORIGINAL_QUERY_COORDINATE_FEASIBILITY_8H7_PILOT_V1",
         "original_panel_SHA256":PANEL_SHA,
         "original_2808_speciesKey_facet_SHA256":FACET_SHA,
         "original_same_454_H7period_population_unchanged":True,
         "prespecified_source_H7":8,
         "original_source_true_H6":len({x[1] for x in selected}),
         "source_original_8_country_holdouts":len({x[0] for x in selected}),
         "source_only_study_H7":len(both),
         "pilot_source_siteperiod_queries":len(audits),
         "pilot_source_query_success":len(allvalid),
         "pilot_source_query_missing":len(audits)-len(allvalid),
         "pilot_FIRST_PAGE_records_scanned":source_rows,
         "pilot_FIRST_PAGE_valid_WGS84_coordinate":coords,
         "pilot_FIRST_PAGE_explicit_uncertainty_le_100m":uncert100,
         "pilot_FIRST_PAGE_le100m_fraction_among_first_page":uncert100/source_rows if source_rows else None,
         "original_facet_lacks_individual_occurrence_coordinate":True,
         "sample_source_is_NOT_original_full_observation_population":True,
         "current_GBIF_index_cannot_be_claimed_byte_identical_to_20261114_index":True,
         "not_enough_source_to_compute_any_river_buffer_or_EIV_HNS_effect":True,
         "no_wastewater_effect_or_EQR_NUTRI_HNS_W_values_opened":True,
         "per_site_period":audits}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"FIG4_GBIF_8H7_COORDINATE_SOURCE_PILOT_QA.json").write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    summary={k:out[k] for k in out if k!="per_site_period"}
    print("GESUI_FIG4_GBIF_GEOFEASIBILITY="+json.dumps(summary,ensure_ascii=False,sort_keys=True),flush=True)
    print("RIVER_DISTANCE_MODEL_GATE=NOT_FITTED_REQUIRED_COMPLETE_RETRIEVAL_AND_NETWORK",flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--facets",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    main(ap.parse_args())
