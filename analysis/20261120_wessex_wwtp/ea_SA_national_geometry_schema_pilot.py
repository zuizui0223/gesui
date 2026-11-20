#!/usr/bin/env python3
"""Minimal EA WQE SA geometry structure check from verified native SA filter.

Never export sampling point notations, names or coordinates. Determine whether
source JSON-LD geometry is a BNG point and whether the documented filter stays
exclusive SA on a second page. No ChemPop ecological response accessed.
"""
from __future__ import annotations
import argparse,hashlib,json,re,urllib.parse,urllib.request
from pathlib import Path

URL="https://environment.data.gov.uk/water-quality/data/sampling-point"
HEADER={"User-Agent":"gesui-EA2026-SA-geometry-original-source-schema",
        "Accept":"application/ld+json","Content-Type":"application/json",
        "API-Version":"1"}
def normalized(value,level=0):
    if level>5:return "MAX_DEPTH"
    if isinstance(value,(int,float)):return "numeric"
    if isinstance(value,str):
        if len(value)>120:return "string>120chars"
        if any(x in value.lower() for x in ("point","coord","crs","bng")):
            return "string with geometry text"
        return "string"
    if isinstance(value,list):
        return {"type":"list","length":len(value),
                "first_value_schema":normalized(value[0],level+1) if value else None}
    if isinstance(value,dict):
        return {k:normalized(v,level+1) for k,v in value.items()}
    return str(type(value).__name__)

def main(a):
    result={"schema":"GESUI_WQE_SA_OFFICIAL_BNG_GEOMETRY_SCHEMA_ONLY_V1",
            "probes":[],"source_only":True,"no_source_station_name_or_coordinates_published":True}
    for skip in (0,250):
        query=urllib.parse.urlencode({"samplingPointType":"SA","limit":3,"skip":skip})
        req=urllib.request.Request(URL+"?"+query,data=b"",method="POST",headers=HEADER)
        try:
            with urllib.request.urlopen(req,timeout=40) as r:
                raw=r.read(45000)
            d=json.loads(raw)
            members=d["member"]
            if not members:raise ValueError("Original filtered SA registry returned empty")
            typ={str(m.get("samplingPointType",{}).get("notation")) for m in members}
            if typ!={"SA"}:raise ValueError("Native SA samplingPointType filter not exclusive")
            rec={"status":"SOURCE_METADATA_SCHEMA",
                 "skip":skip,"totalItems":d.get("totalItems"),
                 "member_n":len(members),
                 "member_field_types":normalized(members[0]),
                 "geometry_field_example_structure":normalized(members[0].get("geometry")),
                 "source_response_SHA256":hashlib.sha256(raw).hexdigest()}
        except Exception as e:
            rec={"status":"SOURCE_STOP","skip":skip,
                 "error":type(e).__name__+": "+str(e)[:320]}
        result["probes"].append(rec)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_WQE_SA_GEO_SCHEMA="+json.dumps(result,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
