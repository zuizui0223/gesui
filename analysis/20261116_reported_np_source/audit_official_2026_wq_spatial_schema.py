#!/usr/bin/env python3
"""Source-only native official EEA 2026 ICM spatial-monitoring file schema.

Official public 2,797,036B WISE6_SpatialObjects_DerivedData-csv.zip SHA
locked before processing. Counts WQ spatial monitoring identifiers, coordinate
schema and nonmissingness, but does NOT rehost coordinates/IDs, download all
chemistry, or compute H7 or EQR/plant/butterfly/wastewater effects.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,subprocess,zipfile
from pathlib import Path
import requests
SOURCE_URL="https://sdi.eea.europa.eu/datashare/public.php/dav/files/ZcJ37WK6oaWonFJ/WISE6_SpatialObjects_DerivedData-csv.zip"
EXPECTED_BYTES=2797036
UA="gesui-official-2026-ICM-spatial-site-source-schema-only/20261116"
COORD_NAMES=("latitude","longitude","lat","lon","the_geom","geometry","wkt","point","easting","northing","latitudeWGS84","longitudeWGS84",
             "monitoringSiteIdentifier","monitoringSiteIdentifierScheme","monitoringSite","spatialObjectIdentifier")

def obtain(path):
    headers={"User-Agent":UA,"Accept":"application/zip"}
    with requests.get(SOURCE_URL,headers=headers,stream=True,timeout=(10,80)) as r:
        if r.status_code!=200:raise RuntimeError("EEA official spatial metadata ZIP source unavailable HTTP "+str(r.status_code))
        if r.url.split("/")[2]!="sdi.eea.europa.eu":raise RuntimeError("Official EEA spatial metadata redirected off host")
        size=0;dig=hashlib.sha256()
        with path.open("wb") as f:
            for chunk in r.iter_content(chunk_size=262144):
                if not chunk:continue
                size+=len(chunk)
                if size>EXPECTED_BYTES:raise RuntimeError("Original official spatial metadata exceeds listed bytes")
                dig.update(chunk);f.write(chunk)
        if size!=EXPECTED_BYTES or not zipfile.is_zipfile(path):
            raise RuntimeError("Original EEA official ZIP size or structure mismatch")
        return {"source_original_url":SOURCE_URL,"publisher_listed_bytes":EXPECTED_BYTES,
                "original_source_SHA256":dig.hexdigest(),
                "verified_original_zip_size":size}

def read_zip(zip_path):
    reported=[]
    z=zipfile.ZipFile(zip_path)
    entries=[f for f in z.infolist() if not f.is_dir()]
    source_meta=[{"filename":f.filename,"compressed_bytes":f.compress_size,
                  "uncompressed_bytes":f.file_size,"compression_method":f.compress_type} for f in entries]
    if not entries:raise RuntimeError("Official spatial metadata ZIP empty")
    for f in entries:
        if not f.filename.lower().endswith(".csv"):continue
        proc=None
        if f.compress_type==9:
            proc=subprocess.Popen(["7z","x","-so",str(zip_path),f.filename],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            binary=proc.stdout
        else:binary=z.open(f)
        with binary:
            stream=io.TextIOWrapper(binary,encoding="utf-8-sig",newline="",errors="strict")
            reader=csv.reader(stream)
            header=next(reader)
            if len(header)>150 or any(len(x)>250 for x in header):
                raise RuntimeError("Official spatial source header invalid")
            sample=next(reader,None)
            records=1 if sample is not None else 0
            checks={k:header.index(k) for k in header if k.lower() in set(v.lower() for v in COORD_NAMES)}
            eligible={k:int(bool(sample[ix].strip())) for k,ix in checks.items()} if sample else {k:0 for k in checks}
            for row in reader:
                records+=1
                for k,ix in checks.items():
                    if ix<len(row) and row[ix].strip():eligible[k]+=1
            if proc is not None:
                status=proc.wait(timeout=20)
                if status!=0:raise RuntimeError("Official EEA spatial CSV 7z exit nonzero")
        reported.append({"original_member_filename":f.filename,"original_source_csv_rows":records,
            "original_column_names":header,"selected_spatial_id_coordinate_field_nonempty_count":eligible,
            "real_coordinate_or_station_values_not_exported":True})
    return source_meta,reported

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    path=a.output/"TEMP_OFFICIAL_2026_SPATIAL_OBJECTS.zip"
    src=obtain(path)
    source,details=read_zip(path)
    path.unlink()
    result={"schema":"GESUI_20261116_EEA_2026_OFFICIAL_WQ_SPATIAL_OBJECTS_SOURCE_SCHEMA_V1",
       "official_2026_publisher_spatial_file":src,
       "original_ZIP_members":source,
       "original_CSV_source_structure":details,
       "original_site_identifier_and_coordinate_values_not_republished":True,
       "spatial_H7_polygon_or_flow_connection_NOT_PERFORMED":True,
       "same_h7_not_same_native_EEA_biology_station":True,
       "original_ambient_water_chemistry_NOT_JOINED":True,
       "no_wastewater_plant_butterfly_causal_claim":True}
    (a.output/"EEA2026_WQ_SPATIAL_OBJECTS_SOURCE_SCHEMA_QA.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_EEA2026_WQ_SPATIAL_SOURCE_SCHEMA="+json.dumps(result,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True);main(p.parse_args())
