#!/usr/bin/env python3
"""SOURCE-only GET probe of officially catalogued Saarland sewage sludge GIS layer.

No WFS GetFeature or biological data are fetched. Metadata, feature-layer names
and attributable date fields are audited before any real spatial exposure join.
"""
import argparse
import hashlib
import json
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree

SOURCES = {
  "official_wfs_capabilities": "https://geoportal.saarland.de/arcgis/services/Internet/Abfallbeseitigung/MapServer/WFSServer?request=GetCapabilities&VERSION=2.0.0&SERVICE=WFS",
  "official_ogc_features_collection": "https://geoportal.saarland.de/spatial-objects/321/collections/Abfallbeseitigung:Klaerschlamm_Verfahren",
}
CATALOG = "https://md.umwelt.info/dataset/metaver-sl/222f0723-3436-aaca-5f7d-778779f1d407"
LICENSE = "CC BY 4.0"
MAX_BYTES = 4_000_000


def localname(tag):
    return str(tag).split("}")[-1].split(":")[-1]


def parse_wfs_capabilities(data):
    root = ElementTree.fromstring(data)
    if "Capabilities" not in localname(root.tag):
        raise ValueError("not WFS capabilities")
    found = []
    for e in root.iter():
        if localname(e.tag) != "FeatureType":
            continue
        attrs = {}
        for node in e.iter():
            name = localname(node.tag)
            if name in ("Name", "Title", "Abstract") and node.text and name not in attrs:
                attrs[name] = node.text.strip()[:500]
        if attrs.get("Name"):
            attrs["has_klarschlamm_token"] = ("klaerschlamm" in attrs["Name"].lower()
                 or "klärschlamm" in attrs.get("Title", "").lower())
            found.append(attrs)
    return {
       "XML_root": localname(root.tag),
       "WFS_feature_type_count": len(found),
       "feature_type_names_matching_sludge": [row for row in found if row["has_klarschlamm_token"]][:20],
       "matching_layer_name_confirmed": any(row["has_klarschlamm_token"] for row in found),
       "year_field_verified": False,
       "actual_vs_plan_status_verified": False,
       "source_WWTP_ID_verified": False,
       "actual_application_mass_verified": False,
       "water_and_land_H6_biotic_model_gate": "STOP",
    }


def probe(url):
    if not url.startswith("https://geoportal.saarland.de/"):
        raise ValueError("Endpoint must be actual official geoportal")
    headers = {"User-Agent":"gesui-academic-source-metadata-QA/1.0",
               "Accept":"application/xml,application/json,text/html;q=0.8"}
    out = {"source_url":url, "request_method":"GET",
           "access_status":"SOURCE_ACCESS_STOP", "content_sha256":None,
           "raw_data_redistributed":False}
    try:
        req = urllib.request.Request(url, headers=headers, method="GET")
        with urllib.request.urlopen(req, timeout=22) as resp:
            body = resp.read(MAX_BYTES + 1)
            if len(body) > MAX_BYTES:
                out["access_note"]="source metadata exceeds bounded read size"
                return out
            out["http_status"]=resp.status
            out["content_type"]=resp.headers.get("Content-Type","")
            out["received_bytes"]=len(body)
            out["content_sha256"]=hashlib.sha256(body).hexdigest()
            if resp.status != 200:
                out["access_note"]="HTTP not 200"
                return out
        out["access_status"]="SOURCE_METADATA_HTTP200"
        if "wfs" in url.lower():
            try:
                out["wfs_capabilities"]=parse_wfs_capabilities(body)
                out["access_status"]="VALID_WFS_LAYER_METADATA" if out["wfs_capabilities"]["matching_layer_name_confirmed"] else "WFS_METADATA_NO_CONFIRMED_LAYER"
            except (ElementTree.ParseError,ValueError) as exc:
                out["access_status"]="SOURCE_SCHEMA_STOP"
                out["access_note"]="WFS XML schema invalid: " + type(exc).__name__
        else:
            # The official OGC landing page can be JSON or HTML. Do not call
            # this proven field-application data unless its schema says so.
            out["source_schema_format"]="json" if body.lstrip().startswith(b"{") else "html_or_other"
            if out["source_schema_format"]=="json":
                try:
                    v=json.loads(body)
                    out["official_catalog_response_keys"]=sorted(v)[:60]
                    out["collection_id"]=v.get("id")
                except (UnicodeDecodeError,ValueError):
                    out["access_status"]="SOURCE_SCHEMA_STOP"
                    out["access_note"]="JSON schema invalid"
        return out
    except Exception as exc:
        out["access_note"]=type(exc).__name__ + ": " + re.sub(r"https?://[^ ]+", "[url]", str(exc))[:300]
        return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    result={
      "schema":"gesui_Saarland_sewage_sludge_spreading_OGC_metadata_source_only_v1",
      "utc":datetime.now(timezone.utc).isoformat(),
      "source_catalog":CATALOG,
      "catalog_license":LICENSE,
      "scope":"Saarland Germany, prospective source only",
      "sources":{k:probe(url) for k,url in SOURCES.items()},
      "features_or_field_polygons_downloaded":False,
      "years_2017to22_for_actual_application_verified":False,
      "distinguishes_actual_application_from_permitted_areas":False,
      "sewage_station_attribution_verified":False,
      "independent_H6_for_cross_compartment_estimation":0,
      "source_GO_for_biological_outcomes":False,
      "NUTRI_HNS_EQR_values_accessed":False,
      "next_step":"Obtain WFS DescribeFeatureType only if a real sludge layer Name is returned; audit actual application years/field status and station ID BEFORE GetFeature"
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_SAARLAND_SLUDGE_OGC_METADATA_GATE",json.dumps(result,ensure_ascii=False,sort_keys=True))


if __name__=="__main__":
    main()
