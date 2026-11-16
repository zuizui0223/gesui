#!/usr/bin/env python3
"""Official JRC GREEN 2014–2018 N/P modelled point source SOURCE only.

Read real public ArcGIS MapServer 3+32 country, 4+33 NUTS2, 5+34
river-basin source values (units tonnes, averaged 2014–2018). Original 2016
EEA BQE station COUNTRY/METHOD coverage is fixed BEFORE numeric 2016 EQR.
No 2016 EQR values. Not individual WWTP measured N/P, not 2016 observed mass.
All two nutrients and all three geography scales reported, no sign selection.
"""
from __future__ import annotations
import argparse,hashlib,json,math,urllib.parse,urllib.request,urllib.error,time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd

BASE="https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/Losses_to_water_MIS/MapServer"
LAYER={"country":{"N":3,"P":32,"key":["cntrCode"]},
       "NUTS2":{"N":4,"P":33,"key":["cntrCode","Nuts2"]},
       "river_basin":{"N":5,"P":34,"key":["BasinID"]}}
SOURCE=Path("results/eea2016_thirdyear_source_20261111/JRC_GREEN_2014TO2018_POINTSOURCE_NP_OFFICIAL_ARCGIS_SOURCE_LAYERS.json")
BQE=Path("results/eea2016_thirdyear_source_20261111/2016_BQE_PREOUTCOME_SOURCE_ONLY/2016_BQE_PREOUTCOME_SOURCE_QC.json")
HEAD={"User-Agent":"gesui-official-JRC-2014to18-source-real-NP-geo-v1","Accept":"application/json"}
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def query(layer,params):
    url=BASE+f"/{layer}/query?"+urllib.parse.urlencode({"f":"json",**params})
    for i in range(4):
        try:
            req=urllib.request.Request(url,headers=HEAD)
            with urllib.request.urlopen(req,timeout=85) as res:
                body=res.read(10_000_000)
                obj=json.loads(body)
                if not isinstance(obj,dict) or obj.get("error"):
                    raise ValueError("JRC original ArcGIS source error: "+str(obj.get("error"))[:250])
                return obj,hashlib.sha256(body).hexdigest(),url
        except (urllib.error.HTTPError,urllib.error.URLError,ValueError,TimeoutError) as e:
            if i==3:raise ValueError("Official GREEN MapServer SOURCE unavailable: "+str(e)[:200])
            time.sleep((i+1)*2)
    raise ValueError("Official JRC query failed")
def read(layer,fields):
    # Query number of objects independently. No "count returned = true zero".
    count,hcount,ucount=query(layer,{"where":"1=1","returnCountOnly":"true"})
    n=count.get("count")
    if type(n) is not int or n<1 or n>1_000_000:
        raise ValueError("JRC official GIS source layer impossible record count: "+repr(n))
    # Basin layers have many tens of thousands of original polygons. Do NOT
    # select a partial basin subset or claim their omitted nutrients are zero.
    # Country and NUTS2 layers are exhaustively indexed and retained.
    if n>6000:
        return None,[{"stage":"source_original_count_ONLY_no_values",
                      "SHA256":hcount,"url":ucount,"n":n}],n
    # Native OBJECTID list source allows stable batches, not trusting offset.
    ids,hids,uids=query(layer,{"where":"1=1","returnIdsOnly":"true"})
    allids=ids.get("objectIds") or []
    if len(set(allids))!=n:raise ValueError("ArcGIS JRC count and objectIds differ")
    allids=sorted(allids)
    rows=[]
    receipts=[{"stage":"count","SHA256":hcount,"url":ucount,"n":n},
              {"stage":"ids","SHA256":hids,"url":uids,"n":len(allids)}]
    for offset in range(0,n,80):
        # ArcGIS Web Application Firewall rejects overlong GET objectIds URLs.
        # Fixed small ID pages do not change the original study source population.
        chunk=allids[offset:offset+80]
        out,stamp,url=query(layer,{"objectIds":",".join(map(str,chunk)),"outFields":",".join(["OBJECTID",*fields]),
                   "returnGeometry":"false"})
        values=out.get("features",[])
        if len(values)!=len(chunk):raise ValueError("Published GREEN layer truncated objectIds page "+str(offset))
        r=[x["attributes"] for x in values]
        if len({z.get("OBJECTID") for z in r})!=len(chunk) or set(int(z["OBJECTID"]) for z in r)!=set(chunk):
            raise ValueError("Actual official ArcGIS source objectIDs not matching source inventory")
        rows.extend(r)
        receipts.append({"stage":"paged_original_attributes","sha256":stamp,
                         "first_original_objectID":chunk[0],"last_original_objectID":chunk[-1],
                         "n":len(r),
                         # API URL contains full original object ID list, no private data
                         "original_layer_query_url":url})
        time.sleep(.12)
    if len(rows)!=n:raise ValueError("Original JRC full geographic source incomplete")
    return rows,receipts,n
def normalize_key(row,fields):
    q=[str(row.get(x) or "").strip().upper() for x in fields]
    if any(not x for x in q):return None
    return "|".join(q)
def main():
    a=argparse.ArgumentParser();a.add_argument("--out",type=Path,required=True);v=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    meta=json.loads((root/SOURCE).read_text())
    bqe=json.loads((root/BQE).read_text())
    if meta["original_layer_count"]!=46 or meta["official_original_service_http"]!=200:
        raise ValueError("Prior original JRC 46-layer source schema not frozen")
    if bqe["independent_year_fixed"]!="2016" or bqe["pre_2016_Biology_EQR_values_read"] is not False:
        raise ValueError("2016 aquatic EQR source values seen or wrong year")
    indices={x["layer"]["id"]:x for x in meta["layer_schema_and_query_support"]}
    result={}
    for scale,expected in LAYER.items():
        pair={}
        receipt={}
        for nutrient in ("N","P"):
            layer=expected[nutrient]
            fields=expected["key"]+["Ps"+nutrient]
            if layer not in indices or not set(fields).issubset({z["name"] for z in indices[layer]["service_fields"]}):
                raise ValueError("Official original GREEN N/P point-source fields changed")
            rows,log,original_count=read(layer,fields)
            prepared={}
            repeated=[]
            if rows is None:
                # Original full basin source too large for bounded ALL-objectID
                # retrieval here; public source count retained but no effect data.
                pair[nutrient]={
                    "rows":{},
                    "source_data_status":"COUNT_ONLY_TOO_LARGE_NO_NUTRIENT_VALUES",
                    "all_original_features_n":original_count,
                    "distinct_nonmissing_source_geographic_ids":None,
                    "duplicate_country_region_key_count":None,
                    "source_Ps_positive_count":None,
                    "source_Ps_missing_numeric_count":None}
                receipt[nutrient]=log
                print("JRC_NP_SOURCE",scale,nutrient,"original features",original_count,
                      "original basin source COUNT ONLY, not sampled",flush=True)
                continue
            for r in rows:
                code=normalize_key(r,expected["key"])
                if code is None:continue
                if code in prepared:
                    repeated.append(code)
                    continue
                v0=r.get("Ps"+nutrient)
                q=None
                if v0 is not None:
                    try:
                        z=float(v0)
                        if math.isfinite(z) and z>=0:q=z
                    except (TypeError,ValueError):pass
                prepared[code]={"source_GIS_recorded_Ps"+nutrient+"_t":q,
                                "original_layer_OBJECTID":int(r["OBJECTID"])}
            pair[nutrient]={"rows":prepared,"source_data_status":"FULL_SOURCE_COMPLETE","duplicate_country_region_key_count":len(repeated),
                "all_original_features_n":len(rows),"distinct_nonmissing_source_geographic_ids":len(prepared),
                "source_Ps_positive_count":sum(z["source_GIS_recorded_Ps"+nutrient+"_t"] is not None and z["source_GIS_recorded_Ps"+nutrient+"_t"]>0 for z in prepared.values()),
                "source_Ps_missing_numeric_count":sum(z["source_GIS_recorded_Ps"+nutrient+"_t"] is None for z in prepared.values())}
            receipt[nutrient]=log
            print("JRC_NP_SOURCE",scale,nutrient,"n",len(rows),"unique",len(prepared),
                  "duplicates",len(repeated),flush=True)
        entire_scale=all(pair[n]["source_data_status"]=="FULL_SOURCE_COMPLETE" for n in ("N","P"))
        keys=sorted(set(pair["N"]["rows"])|set(pair["P"]["rows"])) if entire_scale else []
        joint=[]
        for key in keys:
            n=pair["N"]["rows"].get(key,{})
            p=pair["P"]["rows"].get(key,{})
            joint.append({"original_geography_ID":key,
                          "original_GREEN2014_18_pointN_modelled_t":n.get("source_GIS_recorded_PsN_t"),
                          "original_GREEN2014_18_pointP_modelled_t":p.get("source_GIS_recorded_PsP_t"),
                          "original_original_layerN_objectID":n.get("original_layer_OBJECTID"),
                          "original_original_layerP_objectID":p.get("original_layer_OBJECTID"),
                          "N_and_P_original_both_coverage":n.get("source_GIS_recorded_PsN_t") is not None and p.get("source_GIS_recorded_PsP_t") is not None})
        usable=[z for z in joint if z["N_and_P_original_both_coverage"]]
        source_country={z["original_geography_ID"].split("|")[0] for z in usable} if scale in ("country","NUTS2") else set()
        bqe_cntr=set(bqe["original_2016_BQE_by_group"]["benthic_invertebrates"]["source_reporting_by_country"])
        result[scale]={
            "official_N_and_P_layer_indices":{"N":expected["N"],"P":expected["P"]},
            "original_region_id_field_names":expected["key"],
            "original_complete_both_np_geographic_source_acquisition":entire_scale,
            "basin_source_values_not_sampled_when_original_dataset_exceeds_6000":not entire_scale,
            "original_modelled_NP_joint_geographic_records_n":len(joint) if entire_scale else None,
            "original_modelled_NP_both_nonmissing_regions_n":len(usable) if entire_scale else None,
            "N_source_quality":{k:v for k,v in pair["N"].items() if k!="rows"},
            "P_source_quality":{k:v for k,v in pair["P"].items() if k!="rows"},
            "original_2016_EEA_benthic_country_overlap_n":len(bqe_cntr&source_country) if scale!="river_basin" else None,
            "original_2016_EEA_benthic_country_overlap_codes":sorted(bqe_cntr&source_country) if scale!="river_basin" else [],
            "source_all_original_2014to2018_JRC_modelled_point_NP_t_region_aggregates":joint,
            "source_official_API_pages_original_SHA256_and_urls":receipt,
            "no_biological_2016_EQR_read":True
        }
    n=result["NUTS2"]
    returnValue={
        "schema":"gesui_OFFICIAL_JRC_GREEN_MAPSERVER_original_2014to2018_POINT_SOURCE_NP_3GEO_SCALES_SOURCE_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "official_ArcGIS_original_JRC_source":"https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/Losses_to_water_MIS/MapServer",
        "source_original_46layer_GIS_schema_SHA256":sha(root/SOURCE),
        "source_2016_EEA_BQE_site_metadata_only_SHA256":sha(root/BQE),
        "original_SOURCE_nutrients_are_modelled_GREEN_point_sources_DOMESTIC_AND_INDUSTRIAL":True,
        "original_SOURCE_time_fiveyear_average_2014_2018_NOT_single2016":True,
        "original_source_UNITS_tons_over_aggregated_area_without_2016_WWTP_specific_mass":True,
        "original_geographic_scale_is_country_NUTS2_RIVER_BASIN_NOT_individual_WWTP_reach":True,
        "original_JRC_GREEN_2014_2018_source_data":result,
        "2016_benthic_EEA_numeric_EQR_values_NOT_READ":True,
        "source_gate_to_consider_EEA2016_NUTS2_spatial_join":bool(
            n["original_modelled_NP_both_nonmissing_regions_n"]>=100 and
            n["original_2016_EEA_benthic_country_overlap_n"]>=6 and
            not n["N_source_quality"]["duplicate_country_region_key_count"] and
            not n["P_source_quality"]["duplicate_country_region_key_count"]),
        "source_do_not_equate_GREEN_aggregated_NP_with_2018_2022_reported_WWTP_actual_upstream_loads":True
    }
    v.out.parent.mkdir(parents=True,exist_ok=True)
    v.out.write_text(json.dumps(returnValue,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("JRC_ACTUAL_GREEN_MODELLED_POINT_SOURCE_2014TO18",json.dumps({
       k:{"n":x["original_modelled_NP_joint_geographic_records_n"],
          "bothNP":x["original_modelled_NP_both_nonmissing_regions_n"],
          "benthic2016_country_overlap":x["original_2016_EEA_benthic_country_overlap_n"],
          "source_duplicates":{nut:x[nut+"_source_quality"]["duplicate_country_region_key_count"] for nut in ("N","P")},
          "basin_data_FULL":x["original_complete_both_np_geographic_source_acquisition"]}
       for k,x in result.items()},indent=2),flush=True)
if __name__=="__main__":main()
