#!/usr/bin/env python3
"""Official Saarland WFS DescribeFeatureType schema audit, SOURCE-only.

Read the exact previously verified WFS GetCapabilities endpoint. A DescribeFeatureType
request is allowed only for the layer named in that actual official metadata and if
the operation itself is advertised. No GetFeature/field polygons or ecological outcomes.
"""
import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET

CAPABILITIES="https://geoportal.saarland.de/arcgis/services/Internet/Abfallbeseitigung/MapServer/WFSServer?request=GetCapabilities&VERSION=2.0.0&SERVICE=WFS"
EXPECTED_LAYER="Abfallbeseitigung:Klaerschlamm_Verfahren"
MAX_META_BYTES=2_500_000

def lname(tag):
    return str(tag).split("}")[-1]

def parse_capabilities(raw):
    root=ET.fromstring(raw)
    if lname(root.tag)!="WFS_Capabilities":
        raise ValueError("Source is not official WFS capabilities XML")
    operations=set()
    for node in root.iter():
        if lname(node.tag)=="Operation" and node.attrib.get("name"):
            operations.add(node.attrib["name"])
    layers=set()
    for node in root.iter():
        if lname(node.tag)!="FeatureType":
            continue
        for child in node:
            if lname(child.tag)=="Name" and child.text:
                layers.add(child.text.strip())
    return {
        "WFS_layer_names":sorted(layers),
        "WFS_operations_advertised":sorted(operations),
        "target_layer_exactly_listed":EXPECTED_LAYER in layers,
        "DescribeFeatureType_operation_advertised":"DescribeFeatureType" in operations,
    }

def parse_describe_schema(raw):
    root=ET.fromstring(raw)
    if lname(root.tag)!="schema":
        raise ValueError("Not an XSD schema")
    vals=[]
    for node in root.iter():
        if lname(node.tag) in ("element","attribute"):
            name=node.attrib.get("name")
            if name and name not in vals:
                vals.append(name[:100])
    if len(vals)>1000:raise ValueError("Unbounded schema")
    # CANDIDATES only, not verification of the semantics or any actual use
    patterns={
        "potential_time_fields":r"(jahr|year|date|datum|zeit|periode|time|beg[ui]nn|ende)",
        "potential_source_station_fields":r"(st[e]?u|klaranlage|klaeranlage|anlage|wwtp|sand[re]|quelle|herkunft|betreiber)",
        "potential_actual_use_fields":r"(ausbr[ie]ng|epand|verwert|durchgef|real|anwend|nutz|status|verfahren)",
        "potential_quantity_fields":r"(menge|quant|mass|gewicht|tonn|gewicht|kg|trocken|tms|ts|masse)",
        "potential_nutrient_fields":r"(stickstoff|phosph|n_total|p_total|naeh|naehr)",
    }
    candidates={key:[x for x in vals if re.search(expr,x,re.I)] for key,expr in patterns.items()}
    return {
        "schema_root":"schema",
        "source_field_names":vals,
        "source_field_count":len(vals),
        "source_field_candidates":candidates,
        "realized_application_year_verified":False,
        "realized_application_event_or_plan_status_verified":False,
        "WWTP_identity_verified":False,
        "date_values_checked":False,
        "geometry_or_feature_values_downloaded":False,
        "nutrient_masses_validated":False,
        "source_year_2017to22_field_exposure_GO":False,
        "biological_model_gate":"STOP"
    }

def get(url):
    if not url.startswith("https://geoportal.saarland.de/arcgis/services/Internet/Abfallbeseitigung/MapServer/WFSServer?"):
        raise ValueError("Unexpected host or WFS route")
    req=Request(url,headers={"User-Agent":"gesui-university-source-metadata-only/1.0",
                             "Accept":"application/xml,text/xml"},method="GET")
    with urlopen(req,timeout=25) as resp:
        data=resp.read(MAX_META_BYTES+1)
        if len(data)>MAX_META_BYTES:
            raise ValueError("Unbounded metadata response")
        return data,{"status":resp.status,"content_type":resp.headers.get("Content-Type",""),
                     "size_bytes":len(data),"sha256":hashlib.sha256(data).hexdigest()}

def audit(fetch=get):
    report={
        "schema":"gesui_DE_Saarland_sludge_WFS_official_DescribeFeatureType_source_only_v1",
        "utc":datetime.now(timezone.utc).isoformat(),
        "official_WFS_capabilities":CAPABILITIES,
        "target_layer_from_prior_confirmed_WFS":EXPECTED_LAYER,
        "WFS_source_metadata":None,
        "WFS_DescribeFeatureType_metadata":None,
        "source_gate":"ACCESS_STOP",
        "no_GetFeature_or_field_geometries_downloaded":True,
        "no_EQR_NUTRI_HNS_values_read":True,
        "biological_model_gate":"STOP",
    }
    try:
        c,receipt=fetch(CAPABILITIES)
        report["WFS_capabilities_receipt"]=receipt
        if receipt["status"]!=200:raise ValueError("WFS cap HTTP not 200")
        cap=parse_capabilities(c)
        report["WFS_source_metadata"]=cap
        if not cap["target_layer_exactly_listed"]:
            report["source_gate"]="NO_MATCHING_OFFICIAL_LAYER_STOP"
            return report
        if not cap["DescribeFeatureType_operation_advertised"]:
            report["source_gate"]="DESCRIBE_OPERATION_NOT_ADVERTISED_STOP"
            return report
        # Only operation/layer advertised by the official capabilities. URL
        # remains on exactly that official WFS route; no guessed endpoint.
        base=CAPABILITIES.split("?")[0]
        query=urlencode({
            "SERVICE":"WFS","VERSION":"2.0.0","REQUEST":"DescribeFeatureType",
            "TYPENAMES":EXPECTED_LAYER
        })
        url=base+"?"+query
        d,desc_receipt=fetch(url)
        report["WFS_DescribeFeatureType_url"]=url
        report["WFS_describe_receipt"]=desc_receipt
        if desc_receipt["status"]!=200:raise ValueError("Describe HTTP not 200")
        schema=parse_describe_schema(d)
        report["WFS_DescribeFeatureType_metadata"]=schema
        report["source_gate"]="WFS_SCHEMA_FIELDS_CONFIRMED_VALUES_NOT_CHECKED"
        if not schema["source_field_candidates"]["potential_time_fields"]:
            report["source_gate"]="NO_OBVIOUS_TIME_FIELD_FIELD_YEAR_STOP"
        return report
    except Exception as exc:
        report["source_error"]=type(exc).__name__+": "+str(exc)[:250]
        return report

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    r=audit()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf8")
    print("GESUI_SAARLAND_OFFICIAL_WFS_DESCRIBE_SOURCE",json.dumps(r,ensure_ascii=False,sort_keys=True))

if __name__=="__main__":main()
