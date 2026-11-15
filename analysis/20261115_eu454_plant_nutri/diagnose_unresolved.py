#!/usr/bin/env python3
"""Diagnose only 181 original speciesKey source schema gaps from full EIV join.

No taxonomy imputation, no score changes, no ecology outcomes, no threshold
relaxation. Re-query exact original 181 GBIF API IDs once for metadata QA and
official field presence; test whether accepted canonicalName is genuinely
available or original GBIF speciesKey no longer resolves cleanly.
"""
from __future__ import annotations
import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import requests

def main(a):
    shards=list(a.shards.rglob("TAXONOMY_KEYS_SHARD_*.json"))
    if len(shards)!=8:raise RuntimeError("Eight source-original completed shards required")
    source=[]
    for path in shards:
        source.extend(json.loads(path.read_text()))
    if len(source)!=11511:raise RuntimeError("Original unique speciesKey count changed")
    invalid=[r for r in source if r.get("status")!="OK"]
    if len(invalid)!=181:raise RuntimeError("Unexpected original invalid key source count "+str(len(invalid)))
    fields=["key","rank","canonicalName","scientificName","taxonomicStatus","kingdom",
            "phylum","genus","species","acceptedKey","acceptedUsageKey","nubKey",
            "authorship"]
    headers={"User-Agent":"gesui-h7-source-taxonomy-QA-20261115","Accept":"application/json"}
    def source_one(row):
        k=str(row["speciesKey"])
        try:
            resp=requests.get("https://api.gbif.org/v1/species/"+k,headers=headers,timeout=(10,35))
            code=resp.status_code
            resp.raise_for_status()
            val=resp.json()
            if not isinstance(val,dict):raise ValueError("not JSON object")
            out={f:val.get(f) for f in fields}
            out["speciesKey"]=k
            out["source_http"]=code
            out["missing_fields"]=[f for f in ["rank","canonicalName","scientificName",
                      "taxonomicStatus","kingdom","phylum"] if not val.get(f)]
            return out
        except Exception as exc:
            return {"speciesKey":k,"status":"ERROR","error":type(exc).__name__+": "+str(exc)[:120]}
    with ThreadPoolExecutor(max_workers=8) as exe:
        reports=list(exe.map(source_one,invalid))
    counts={}
    for out in reports:
        if out.get("status")=="ERROR":signature="ERROR"
        else:signature=";".join(out["missing_fields"]) or "NONE"
        counts[signature]=counts.get(signature,0)+1
    out={"schema":"GESUI_EU454_181_UNRESOLVED_TAXONOMY_SOURCE_ONLY_AUDIT_V1",
         "n_original_keys":11511,"original_invalid_n":181,
         "no_trait_or_W_or_butterfly_response_used":True,
         "missing_field_signature_counts":counts,
         "source_field_examples":reports[:18],
         "can_recover_exact_canonical_rank_both":sum(bool(x.get("canonicalName")) and bool(x.get("rank")) for x in reports),
         "not_automatically_replacing_original_strict_scores":True,
         "source_taxon_raw_api_responses":reports}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"GBIF_181_UNRESOLVED_SOURCE_DIAGNOSTIC.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_181_UNRESOLVED_SOURCE_DIAG="+json.dumps({
        k:out[k] for k in ["n_original_keys","original_invalid_n",
            "missing_field_signature_counts","can_recover_exact_canonical_rank_both",
            "source_field_examples"]},ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--shards",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
