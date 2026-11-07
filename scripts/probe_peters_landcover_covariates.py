#!/usr/bin/env python3
"""Verify independent catchment/population/agriculture covariate schemas in Peters v5.

Only two small published CSVs; source MD5 checked. Do not inspect concentration
values, plant outcomes, or model results.
"""
import hashlib, io, json
from pathlib import Path
import requests,pandas as pd
ARTICLE="https://api.figshare.com/v2/articles/30546371"
FILES=("TN.4_land.csv","TP.4_land.csv")
def get(url,cap):
    r=requests.get(url,stream=True,timeout=(10,90),
        headers={"User-Agent":"gesui-source-covariate-probe/0.1"})
    r.raise_for_status();b=bytearray()
    for x in r.iter_content(131072):
        b.extend(x)
        if len(b)>cap:raise ValueError("File exceeded source cap")
    r.close();return bytes(b)
meta=json.loads(get(ARTICLE,500000))
assert meta["doi"]=="10.25400/lincolnuninz.30546371.v5"
idx={x["name"]:x for x in meta["files"]}
report={"data_version":meta["doi"],"public_file_schemas":[],"used_nutrient_values":False}
for name in FILES:
    f=idx[name]
    raw=get(f["download_url"],3_000_000)
    assert hashlib.md5(raw).hexdigest()==f["supplied_md5"]
    df=pd.read_csv(io.BytesIO(raw),low_memory=False)
    cols=list(df.columns)
    record={"file":name,"file_id":f["id"],"md5":f["supplied_md5"],
            "rows":len(df),"columns":cols,
            "dtypes":{x:str(df[x].dtype) for x in cols},
            "missingness":{x:int(df[x].isna().sum()) for x in cols},
            "examples_of_identifier_fields":{x:df[x].astype(str).head(2).tolist()
                for x in cols if x.lower() in {"wq.id","id","hybas_id","basin_id","region"}},
            "numeric_covariates_candidate":[x for x in cols if any(
                q in x.lower() for q in ("pop","urban","crop","agric","land","farm","fert","soil","flow","hydro","runoff","slope","precip"))]}
    report["public_file_schemas"].append(record)
    print("GESUI_VERIFIED_LANDCOVER_COVARIATE_SOURCE",json.dumps(record,ensure_ascii=False),flush=True)
Path("/tmp/gesui_landcov_probe").mkdir(exist_ok=True)
Path("/tmp/gesui_landcov_probe/coverage.json").write_text(json.dumps(report,indent=2)+"\n")
