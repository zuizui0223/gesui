#!/usr/bin/env python3
"""Source-only EPA NRSA 2018-19 / 2023-24 periphyton-benthos availability.
No nutrient-biomass/diversity regressions, no WWTP causal source attribution.
No ecological response values are read into the saved results.
"""
from __future__ import annotations

import collections
import csv
import hashlib
import json
import os
import pathlib
import time
import urllib.request

BASE = pathlib.Path(os.environ.get("RUNNER_TEMP", "/tmp")) / "gesui_nrsa_biomass_source"
BASE.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "gesui-source-audit/2026-11-13 (academic research)", "Accept": "text/csv,text/plain,*/*"}
EPA_18 = "https://www.epa.gov"
SOURCES = {
    "2018_19": {
        "sites": ("https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_SiteInfo.csv", "1b5b2119c881290191edeb47f1d80a7711f10440e9851574cfa86fd57be770e3"),
        "water": ("https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_water_chemistry_chla_-_data.csv", "cb9a3349ae3708834eab533e71423d74c1778bf76aba0737f0ae1bcc8f50baa1"),
        "benthos_count": ("https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_benthic_macroinvertebrate_count_-_data.csv", "4e9b3a0798b25f67c7a96d4800b94f86544bbd623abaf119358e067b158324a7"),
        "periphyton_biomass": ("https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_pbio_0.csv", None),
        "periphyton_chla": ("https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_PeriChla.csv", None),
    },
    "2023_24": {
        "sites": ("https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_siteinfo.csv", None),
        "water": ("https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_waterchemistry.csv", None),
        "benthos_count": ("https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_benthiccount.csv", None),
        "benthos_metrics": ("https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_benthicmetrics.csv", None),
        "periphyton_biomass": ("https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_periphytonbiomass.csv", None),
        "periphyton_chla": ("https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_periphytonchlorophyll-a.csv", None),
    },
}
METADATA = {
    "2018_19_biomass": "https://www.epa.gov/system/files/other-files/2023-02/nrsa1819_pbio_metadata.txt",
    "2018_19_chla": "https://www.epa.gov/system/files/other-files/2023-01/nrsa1819_PeriChla.txt",
    "2023_24_biomass": "https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_periphytonbiomass_metadata.csv",
    "2023_24_chla": "https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_periphytonchlorophyll-a_metadata.csv",
    "2023_24_water": "https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_waterchemistry_metadata.csv",
    "2023_24_benthic": "https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_benthiccount_metadata_0.csv",
}
def retrieve(url, filename, expected=None, required=True):
    p = BASE / filename
    last = None
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=90) as resp, p.open("wb") as out:
                while True:
                    b = resp.read(1024 * 1024)
                    if not b: break
                    out.write(b)
            first = p.read_bytes()[:512].lower()
            if not p.stat().st_size or b"<html" in first or b"<!doctype" in first:
                raise RuntimeError("Not raw EPA data (empty or HTML)")
            sha = hashlib.sha256(p.read_bytes()).hexdigest()
            if expected is not None and sha != expected:
                raise RuntimeError("Pinned official 2018-19 source SHA changed: " + filename)
            return p, sha
        except Exception as exc:
            last = str(exc)
            p.unlink(missing_ok=True)
            if isinstance(exc, urllib.error.HTTPError) and exc.code in (403,404,410):
                break  # original source unavailable; do not invent zero records
            time.sleep(min(2 ** attempt, 8))
    if required: raise RuntimeError(f"EPA download FAILED ({filename}): {last}")
    return None, "ERROR: " + str(last)

def inspect_csv(p):
    # Inspect ONLY identifiers, headers, row keys, flags and availability;
    # never read or summarize ecological numerical results here.
    with p.open(encoding="utf-8-sig", newline="") as fp:
        it = csv.DictReader(fp)
        cols = it.fieldnames or []
        names = {c.strip().upper(): c for c in cols}
        uid_col = names.get("UID")
        site_col = names.get("SITE_ID")
        sample_col = names.get("SAMPLE_TYPE")
        year_col = names.get("YEAR")
        date_col = names.get("DATE_COL")
        uids = set()
        sites = set()
        sampletypes = collections.Counter()
        years = collections.Counter()
        counts = collections.Counter()
        n = 0
        nonblank_uids = 0
        for r in it:
            n += 1
            if uid_col:
                u = str(r.get(uid_col, "")).strip()
                if u:
                    nonblank_uids += 1
                    uids.add(u)
                    counts[u] += 1
            if site_col and r.get(site_col):
                sites.add(str(r[site_col]).strip())
            if sample_col: sampletypes[str(r.get(sample_col,"")).strip()]+=1
            if year_col: years[str(r.get(year_col, "")).strip()]+=1
    return ({
        "rows": n, "columns": cols, "uid_field": uid_col,
        "unique_uid": len(uids), "nonblank_uid_rows": nonblank_uids,
        "repeated_uid_count": sum(v>1 for v in counts.values()),
        "unique_site_id": len(sites), "sample_type_count": dict(sampletypes),
        "year_count": dict(years), "date_column_present": bool(date_col),
    }, uids)

def main():
    outputs = {}
    selected_sets = {}
    for period, files in SOURCES.items():
        outputs[period] = {}
        selected_sets[period] = {}
        for kind, (url, expected) in files.items():
            filename = period+"_"+kind+".csv"
            p, sha = retrieve(url,filename,expected,required=(period=="2018_19"))
            if p is None:
                outputs[period][kind]={"source_url":url,"status":"UNAVAILABLE",
                                       "download_error":sha,"columns":[]}
                selected_sets[period][kind]=set()
                print(f"SOURCE {period}/{kind}: unavailable ({sha}); NOT ZERO",flush=True)
                continue
            info, ids = inspect_csv(p)
            info.update({"source_url":url, "sha256":sha, "bytes":p.stat().st_size,
                         "pinned_from_earlier_analysis":bool(expected),"status":"OK"})
            outputs[period][kind]=info
            selected_sets[period][kind]=ids
            print(f"SOURCE {period}/{kind}: rows={info['rows']} UIDs={len(ids)} SHA={sha[:16]}", flush=True)
            p.unlink()  # raw EPA data belong in the reproducible source, not the artifact.
        core=["sites","water","benthos_count","periphyton_biomass","periphyton_chla"]
        source_complete=all(len(selected_sets[period][k])>0 for k in core)
        allkeys=set.intersection(*(selected_sets[period][k] for k in core)) if source_complete else None
        outputs[period]["shared_UID_availability"]={
            "all_five_original_tables_raw_UID_intersection":len(allkeys) if allkeys is not None else None,
            "five_table_original_source_complete":source_complete,
            "pairwise_intersections":{k:len(selected_sets[period]["sites"] & selected_sets[period][k])
                                      for k in core if k!="sites"},
            "source_only_not_numeric_eligible": True,
        }
        if not source_complete:
            print(f"SOURCE {period}: incomplete original UID availability; NO effect fitting.",flush=True)
            if period=="2018_19":
                raise RuntimeError("Previously verified 2018/19 source now missing; stop all fitting")
    meta_out = {}
    for label,url in METADATA.items():
        p, sha = retrieve(url,"metadata_"+label+".txt",required=False)
        meta_out[label]={"url":url, "sha256":sha,
                         "bytes":p.stat().st_size if p else None,
                         "first_lines":p.read_text(encoding="utf-8-sig",errors="replace").splitlines()[:12] if p else []}
        if p: p.unlink()
    report={
        "schema":"gesui_NRSA_periphyton_benthos_2018_19_2023_24_SOURCE_ONLY_v1",
        "periods":outputs,"original_metadata":meta_out,
        "source_blind_eligibility_only":True,
        "no_effect_size_calculated":True,
        "no_WWTP_exposure_assigned":True,
        "not_comparable_years_until_metadata_and_method_audit":True,
        "science":"UID overlap is an availability count, not measured N/P/AFDM/Chla/EPT complete cases or evidence of associations.",
    }
    target=BASE/"SOURCE_ONLY_NRSA_BIOMASS_BENTHOS_2018_2024.json"
    target.write_text(json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False))
    print("SOURCE_QA_JSON_START")
    print(json.dumps({"overlap":{t:outputs[t]["shared_UID_availability"] for t in outputs},
                      "columns":{t:{k:v["columns"] for k,v in d.items() if k!="shared_UID_availability"} for t,d in outputs.items()},
                      "metadata":meta_out},ensure_ascii=False))
    print("SOURCE_QA_JSON_END")

if __name__=="__main__":main()
