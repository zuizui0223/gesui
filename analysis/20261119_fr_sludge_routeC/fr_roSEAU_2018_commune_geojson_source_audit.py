#!/usr/bin/env python3
"""Source-only 2018 French sludge municipality-string and GeoJSON-origin audit.

Only metadata and five PUBLIC municipality text examples. The GeoJSON may
locate the wastewater treatment works rather than the sludge application
field. Never infer exposure from a point without source semantics.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,requests
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime,timezone

BASES=(
 "https://www.assainissement.developpement-durable.gouv.fr",
 "https://assainissement.developpement-durable.gouv.fr",
 "https://www2.assainissement.developpement-durable.gouv.fr",
)
JS="/public/js/exportbasededonneeBoue_min.js"
PAGE="/pages/data/basededonneeboue.php"
HANDLER="/controllers/GetDataExportBoueHandler.php"
EXPECTED_JS_SHA="f4cf1e2fe84e0ccf608ece6c5d363426960e2b730e21900d48265eba4a78f0e5"
EXPECTED_CSV_SHA="bb9660f13e3dc0b7db359cfc414f14e4510db3b51d73ce0baa4bd5cca821f6c2"
MAXBYTES=30_000_000
REQ={"bassin":"ALL","region":"ALL","annee":"2018","tranche":"ALL","departement":"ALL"}

def sh(data):return hashlib.sha256(data).hexdigest()
def official_get_two_formats():
    errors=[]
    for base in BASES:
        client=requests.Session()
        client.headers.update({"User-Agent":"GESUI-FrenchROSEAU-commune-geojson-sourceQA-20261119"})
        try:
            pg=client.get(base+PAGE,timeout=(15,35))
            pg.raise_for_status()
            if "2018" not in pg.text:raise ValueError("2018 year not in official HTML")
            js=client.get(base+JS,timeout=(15,35))
            js.raise_for_status()
            if sh(js.content)!=EXPECTED_JS_SHA:raise ValueError("Official original exporter JS SHA mismatch")
            data={}
            for fmt in ("csv","geojson"):
                post={**REQ,"format":fmt}
                rr=client.post(base+HANDLER,data=post,timeout=(18,110),headers={
                    "Origin":base,"Referer":base+PAGE,"Accept":"application/json,text/csv,*/*"
                })
                rr.raise_for_status()
                if not 3500<len(rr.content)<=MAXBYTES:raise ValueError("Original public export size invalid: "+fmt)
                data[fmt]=(rr.content,rr.headers.get("Content-Type",""))
            return base,data,errors
        except Exception as ex:
            errors.append({"official_host":base,"error":str(ex)[:400]})
    raise RuntimeError("All public official hosts unavailable: "+json.dumps(errors))

def municipality_schema(raw):
    if sh(raw)!=EXPECTED_CSV_SHA:
        raise ValueError("The 2018 public CSV changed from prior captured independent SHA; source drift must be resolved")
    text=raw.decode("utf-8-sig")
    rows=csv.DictReader(io.StringIO(text),delimiter="|")
    cols=list(rows.fieldnames or [])
    commune="Commune d'épandage"
    station="Code SANDRE la station"
    mass="Quantité épandage agricole (tMS/an)"
    if not all(x in cols for x in (commune,station,mass)):
        raise ValueError("Original CSV fields changed")
    vals=Counter()
    types=defaultdict(set)
    delimtypes=Counter()
    positive=0
    for row in rows:
        s=row[commune].strip()
        if not s:continue
        vals[s]+=1
        cats=[]
        for sep,label in ((";", "semicolon"),(",","comma"),("/","slash"),("\n","newline"),(" - ","spaced_hyphen")):
            if sep in s:cats.append(label)
        if not cats:cats=["no_detected_delimiter"]
        for label in cats:delimtypes[label]+=1;types[label].add(s)
        try:
            m=float(row[mass].replace(" ","").replace(",","."))
            if m>0:positive+=1
        except (ValueError,TypeError):pass
    # Municipalities are public geographic names, not parcel IDs.
    # Keep a tiny excerpt for technical interpretation only; no station,
    # parcel, person or original full row is exported.
    return {
      "2018_original_csv_SHA_matched":True,
      "raw_distinct_commune_strings":len(vals),
      "text_delimiter_occurrences_across_nonempty_raw_rows":dict(delimtypes),
      "distinct_literals_by_delimiter":{k:len(v) for k,v in types.items()},
      "public_commune_text_examples_by_type":{
          k:[v.replace("\n"," ")[:120] for v in sorted(vals) if v in group][:4]
          for k,group in types.items()
      },
      "raw_rows_with_positive_station_level_qty_and_nonblank_commune":positive,
      "original_commune_text_NOT_proven_independent_INSEE_code":True,
    }

def geojson_schema(raw):
    j=json.loads(raw)
    if isinstance(j,list):features=j
    else:features=j.get("features",[])
    if not isinstance(features,list):
        raise ValueError("Not a valid GeoJSON FeatureCollection")
    if len(features)<1:raise ValueError("No exported GeoJSON features")
    types=Counter()
    props=set()
    sourcekey=Counter()
    missinggeometry=0
    for feature in features:
        geom=feature.get("geometry") or {}
        kinds=geom.get("type")
        if not kinds:missinggeometry+=1
        types[kinds or "missing"]+=1
        attributes=feature.get("properties") or {}
        props.update(str(k) for k in attributes.keys())
        for label in attributes:
            if "commune" in label.lower() or "epand" in label.lower():
                sourcekey[label]+=1
    return {
        "feature_count":len(features),
        "geometry_types":dict(types),
        "no_geometry_features":missinggeometry,
        "public_property_fields":sorted(props),
        "property_name_hints_commune_or_spreading":dict(sourcekey),
        "GeoJSON_geometry_does_NOT_auto_prove_parcel_or_commune_location":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    out={"schema":"gesui_France_ROSEAU_2018_source_commune_name_and_GeoJSON_geometry_grain_QA_v1",
         "UTC":datetime.now(timezone.utc).isoformat(),
         "year":2018,"no_EQR_NUTRI_HNS_outcomes_read":True,
         "official_csv_not_saved_or_republished":True,
         "official_geojson_not_saved_or_republished":True}
    try:
        host,src,failures=official_get_two_formats()
        out["official_host"]=host
        out["previous_host_attempt_errors"]=failures
        out["original_CSV"]={"bytes":len(src["csv"][0]),"sha256":sh(src["csv"][0]),
                             "metadata":municipality_schema(src["csv"][0])}
        out["original_GeoJSON"]={"bytes":len(src["geojson"][0]),"sha256":sh(src["geojson"][0]),
                                 "metadata":geojson_schema(src["geojson"][0])}
        out["status"]="SOURCE_FORMATS_INSPECTED_GEO_NOT_PARCEL_CONFIRMED"
    except Exception as e:
        out["status"]="SOURCE_FORMAT_ACCESS_OR_QA_STOP"
        out["error"]=str(e)[:700]
    (a.out/"FR_2018_ROSEAU_COMMUNE_STRINGS_AND_GEOJSON_GEOMETRY_QA.json").write_text(
        json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_FR_ROSEAU_2018_COMMUNE_GEOJSON_SOURCE",json.dumps({
        "status":out["status"],
        "commune_schema":out.get("original_CSV",{}).get("metadata"),
        "geojson_count":out.get("original_GeoJSON",{}).get("metadata",{}).get("feature_count"),
        "geojson_types":out.get("original_GeoJSON",{}).get("metadata",{}).get("geometry_types"),
        "error":out.get("error")
    },ensure_ascii=False),flush=True)
if __name__=="__main__":main()
