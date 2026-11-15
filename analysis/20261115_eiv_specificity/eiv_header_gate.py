#!/usr/bin/env python3
"""Source-only schema and six-column completeness of published European EIV.

No original plant species locations, no environmental exposures, no outcomes.
Original 8908-row Tichy et al. source is checksum-locked. The result will fix
exact printed raw column names before constructing site trait summaries.
"""
from __future__ import annotations
import hashlib, io, json, os
from pathlib import Path
import pandas as pd, requests

URL="https://www.sci.muni.cz/botany/juice/NET/sources/EIV_Europe.csv"
SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"

def main():
    response=requests.get(URL,timeout=(15,120),headers={"User-Agent":"gesui-original-source-EIV-sixtrait-probe/1.0"})
    response.raise_for_status()
    data=response.content
    digest=hashlib.sha256(data).hexdigest()
    if digest!=SHA:raise RuntimeError(f"Published original EIV SHA drift: {digest}")
    d=pd.read_csv(io.BytesIO(data),encoding="utf-8-sig")
    if len(d)!=8908 or "NUTRI" not in d.columns or "Unnamed: 0" not in d.columns:
        raise RuntimeError("Original EIV schema changed")
    candidates={}
    for column in d.columns:
        if column=="Unnamed: 0":continue
        numeric=pd.to_numeric(d[column],errors="coerce")
        candidates[column]={"numeric_rows":int(numeric.notna().sum()),"min":float(numeric.min()) if numeric.notna().any() else None,
                            "max":float(numeric.max()) if numeric.notna().any() else None}
    result={"schema":"GESUI_EU_EIV_ORIGINAL_SIXTRAIT_RAW_COLUMN_SOURCE_GATE",
            "original_EIV_SHA256":digest,
            "original_EIV_rows":len(d),
            "original_column_names":list(d.columns),
            "numeric_observation_counts":candidates,
            "no_plant_occurrence_or_W_or_butterfly_outcome_opened":True}
    p=Path(os.environ.get("EIV_SCHEMA_OUT","/tmp/eiv_six_source"))
    p.mkdir(parents=True,exist_ok=True)
    (p/"EUROPE_EIV_SIX_RAW_COLUMNS.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("EIV_SIX_SOURCE_HEADER="+json.dumps(result,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
