#!/usr/bin/env python3
"""Official original Austria67 monitored river HYRIV_ID -> published upstream
RiverATLAS catchment urban/crop/forest source ID coverage.

Public ArcGIS item service is third-party clipped RiverATLAS and can omit
Austria (and excludes Q<1 m3/s). This is only a source-coverage test:
does not read actual P, EQR, source wastewater loads, landcover outcomes.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,math
from pathlib import Path
from urllib.parse import urlparse
import requests

SITE_ORIG=Path("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv")
BQE_BLOB="325f199e985d9862ffee5b0edb205c61a7be3833"
EXCLUDE_SITE="ATFW54110117"
REGIONS={"15_3":9,"15_4":14,"15_5":21,"16_4":13,"16_5":10}
ITEM_ID="8ee1b502377541f0ae45f3cf4842d0d8"
ITEM_URL="https://www.arcgis.com/sharing/rest/content/items/"+ITEM_ID
SCHEMA="GESUI_20261117_AUSTRIA67_TRUE_RIVERATLAS_UPSTREAM_CATCHMENT_ARCGIS_SOURCE_ONLY_V1"
FIELDS=("HYRIV_ID","urb_pc_use","crp_pc_use","for_pc_use")
MIN_VALID=60;MIN_REGIONS=5;MIN_IN_REGION=5

def gitblob(s):return hashlib.sha1(("blob "+str(len(s))+"\0").encode()+s).hexdigest()

def source_sites():
    data=SITE_ORIG.read_bytes()
    if gitblob(data)!=BQE_BLOB:
        raise RuntimeError("Original EEA2022 benthos source Git blob changed")
    with SITE_ORIG.open(encoding="utf-8-sig",newline="") as f:rows=list(csv.DictReader(f))
    if len(rows)!=2475:
        raise RuntimeError("Original EEA2022 benthic field source coverage changed")
    at=[]
    for q in rows:
        if q["countryCode"]!="AT":continue
        parts=q["site_key"].split("|")
        if len(parts)!=3 or parts[0]!="AT" or parts[2]!="euMonitoringSiteCode":
            raise RuntimeError("Original EEA Austrian government PKINT source malformed")
        if parts[1]==EXCLUDE_SITE:continue
        at.append({"station":parts[1],"region":q["region_3deg_block"],
                   "source_reach":int(q["HYRIV_ID"]),
                   "lat":float(q["lat"]),"lon":float(q["lon"])})
    if len(at)!=67 or collections.Counter(x["region"] for x in at)!=REGIONS:
        raise RuntimeError("Authentic source Austria67 field stations or 5 regional groups changed")
    return at

def read_json(session,url,params):
    r=session.get(url,params=params,timeout=(12,70))
    if r.status_code!=200:raise RuntimeError("Public ArcGIS source HTTP "+str(r.status_code)+" "+str(r.url))
    if urlparse(r.url).hostname not in (
        "www.arcgis.com","arcgis.com","services.arcgis.com",
        "services1.arcgis.com","services2.arcgis.com","services3.arcgis.com",
        "services4.arcgis.com","services5.arcgis.com","services6.arcgis.com",
        "services7.arcgis.com","services8.arcgis.com","services9.arcgis.com"):
        raise RuntimeError("Original public ArcGIS feature source unexpectedly changed host")
    obj=r.json()
    if not isinstance(obj,dict):raise RuntimeError("Public ArcGIS did not publish JSON object source")
    if "error" in obj:raise RuntimeError("Public ArcGIS source returned error "+str(obj["error"])[:500])
    return obj,{"requested_url":r.url,"source_response_sha256":hashlib.sha256(r.content).hexdigest(),
                "source_response_bytes":len(r.content)}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    sites=source_sites()
    rid=sorted({x["source_reach"] for x in sites})
    session=requests.Session()
    session.headers.update({"Accept":"application/json",
      "User-Agent":"gesui-20261117-Austria67-original-upstream-watershed-RiverATLAS-sourceQA/1.0"})
    report={"schema":SCHEMA,"source_population":"Austria real EEA2022 benthic and BMLUK2022 river TP sites",
      "original_region_sizes":REGIONS,"original_original_2022_station_n":len(sites),
      "original_unique_HYRIV_ID_reach_n":len(rid),
      "original_source_excluded_prior_geocoordinate_mismatch_site":EXCLUDE_SITE,
      "derived_thirdparty_RiverATLAS_is_not_full_edition":True,
      "RiverATLAS_source_landcover_reference_urb2016_cropland2000_not2021":True,
      "source_only_values_no_Austria_WWTP_TP_EQR_values_read":True}
    try:
        item,receipt=read_json(session,ITEM_URL,{"f":"json"})
        url=item.get("url","")
        report["original_public_arcgis_item"]={"id":item.get("id"),
          "title":item.get("title"),"type":item.get("type"),
          "licenseInfo":str(item.get("licenseInfo"))[:1400],
          "source_public_data_URL":url,"description_excerpt":str(item.get("description"))[:1600],
          "original_source_query_receipt":receipt}
        if item.get("id")!=ITEM_ID or "FeatureServer" not in url:
            raise RuntimeError("Original public RiverATLAS ArcGIS item is not actual FeatureServer source")
        url=url.replace("http://","https://").rstrip("/")
        service,srece=read_json(session,url,{"f":"json"})
        layers=service.get("layers") or []
        if not layers:
            raise RuntimeError("Public ArcGIS RiverATLAS source contains no layer")
        report["original_public_feature_service"]={"service_url":url,
             "layers":[{"id":z.get("id"),"name":z.get("name")} for z in layers],
             "service_description_excerpt":str(service.get("description"))[:1000],
             "source_receipt":srece}
        accepted=[]
        summaries=[]
        for layer in layers:
            layer_id=layer["id"]
            layerurl=url+"/"+str(layer_id)
            m,mdrec=read_json(session,layerurl,{"f":"json"})
            names={f["name"] for f in m.get("fields",[]) if f.get("name")}
            if not set(FIELDS).issubset(names):
                summaries.append({"layer":layer_id,
                  "name":layer.get("name"),"source_required_hydroatlas_fields_missing":sorted(set(FIELDS)-names)})
                continue
            # Query exact unique native HYRIV_IDs; no aquatic outcomes;
            # a service geographically clipped outside Austria must return 0.
            where="HYRIV_ID IN ("+",".join(str(x) for x in rid)+")"
            obj,qr=read_json(session,layerurl+"/query",{
                "where":where,"outFields":",".join(FIELDS),
                "returnGeometry":"false","f":"json","outSR":"4326"})
            hits=obj.get("features",[])
            if not isinstance(hits,list):
                raise RuntimeError("Actual RiverATLAS query source features not array")
            rows=[]
            for feature in hits:
                attr=feature.get("attributes",{})
                try:original_id=int(attr["HYRIV_ID"])
                except (TypeError,ValueError,KeyError):continue
                if original_id not in rid:continue
                values={"HYRIV_ID":original_id}
                valid=True
                for field in FIELDS[1:]:
                    x=attr.get(field)
                    try:v=float(x)
                    except (TypeError,ValueError):
                        valid=False;break
                    if not math.isfinite(v) or not (0<=v<=100):
                        valid=False;break
                    values[field]=v
                if valid:rows.append(values)
            summaries.append({"layer":layer_id,
               "name":layer.get("name"),
               "source_records_returned":len(hits),
               "source_original_RiverATLAS_fields":sorted(FIELDS),
               "source_actual_valid_upstream_landuse_reach_n":len({q["HYRIV_ID"] for q in rows}),
               "source_service_field_metadata_receipt":mdrec,
               "source_original_exact_ID_query_receipt":qr,
               "source_service_extent":m.get("extent"),
               "source_max_record_count":m.get("maxRecordCount")})
            accepted.extend(rows)
        report["public_RiverATLAS_layers_tested"]=summaries
        byreach=collections.defaultdict(list)
        for r in accepted:byreach[r["HYRIV_ID"]].append(r)
        matches=[]
        for q in sites:
            candidates=byreach.get(q["source_reach"],[])
            if len(candidates)!=1:
                matches.append({"station":q["station"],"original_region":q["region"],
                  "source_actual_HYRIV_ID":q["source_reach"],
                  "match":"MISSING" if not candidates else "AMBIGUOUS_DUPLICATE_LAYER"})
                continue
            matches.append({"station":q["station"],"original_region":q["region"],
              "source_actual_HYRIV_ID":q["source_reach"],
              "match":"VALID",**{k:candidates[0][k] for k in FIELDS[1:]}})
        valid=[q for q in matches if q["match"]=="VALID"]
        reg=collections.Counter(q["original_region"] for q in valid)
        passed=len(valid)>=MIN_VALID and sum(z>=MIN_IN_REGION for z in reg.values())>=MIN_REGIONS
        report.update({"status":"PASS_AT67_ORIGINAL_RIVERATLAS_FULL_UPSTREAM_WATERSHED_COVARIATES" if passed
                       else "STOP_THIRD_PARTY_RIVERATLAS_AUSTRIA_SOURCE_NOT_COVERED",
           "original_Austrian_field_sites_matched_RiverATLAS_upstream_use":len(valid),
           "source_valid_site_region_counts":dict(reg),
           "original_source_valid_for_model":passed,
           "source_all_67_PKINT_original_reach_join":matches})
        path=a.output/"AUSTRIA67_RIVERATLAS_UPSTREAM_WATERSHED_SOURCE_COVARIATES.csv"
        with path.open("w",encoding="utf-8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(matches[0]))
            w.writeheader();w.writerows(matches)
        report["original_derived_67_station_csv_SHA256"]=hashlib.sha256(path.read_bytes()).hexdigest()
    except (requests.RequestException,RuntimeError,ValueError,KeyError) as e:
        report.update({"status":"STOP_PUBLIC_RIVERATLAS_FEATURE_SERVICE_SOURCE_ACCESS",
          "source_error_type":type(e).__name__,"source_error_message":str(e)[:550],
          "original_source_valid_for_model":False})
    (a.output/"AUSTRIA67_RIVERATLAS_REAL_UPSTREAM_CATCHMENT_ARCGIS_SOURCE_GATE.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_AUSTRIA67_REAL_RIVERATLAS_CATCHMENT_SOURCE_GATE="+json.dumps({
        "status":report["status"],
        "source_valid_sites":report.get("original_Austrian_field_sites_matched_RiverATLAS_upstream_use",0),
        "region_counts":report.get("source_valid_site_region_counts",{}),
        "original67_reaches":len(rid),
        "original_source_error":report.get("source_error_message"),
        "no_NP_EQR_numerical_model_fitted":True},sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True);main(p.parse_args())
