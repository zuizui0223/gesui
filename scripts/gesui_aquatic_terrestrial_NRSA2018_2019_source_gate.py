#!/usr/bin/env python3
"""2026-11-11 SOURCE-ONLY EPA NRSA2018-19 river benthos/TN/TP joint support.
Reads no ecological outcome relationships, computes no slopes, and commits no raw records.
Only official EPA original CSV URLs defined before any source response is inspected.
"""
import collections
import csv
import hashlib
import json
import math
import os
import pathlib
import sys
import time
import urllib.request

URLS = {
    "benthos": "https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_benthic_macroinvertebrate_count_-_data.csv",
    "sites": "https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_SiteInfo.csv",
    "water": "https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_water_chemistry_chla_-_data.csv",
}
OUT = pathlib.Path(os.environ.get("RUNNER_TEMP", "/tmp")) / "gesui_nrsa1819_source_gate"
OUT.mkdir(parents=True,exist_ok=True)
FILL = {"", "NA", "N/A", "NULL", "NAN", "NONE", "-9999", "MISSING"}
def valid(s):
    return str(s or "").strip().upper() not in FILL

def get(k, url):
    p=OUT/(k+".csv")
    if p.exists() and p.stat().st_size>300:
        return p
    exc=None
    for i in range(5):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"gesui-source-quality-research/1.0 (public academic reproducibility)","Accept":"text/csv,*/*"})
            with urllib.request.urlopen(req,timeout=75) as resp, p.open("wb") as fh:
                while True:
                    chunk=resp.read(1024*1024)
                    if not chunk:break
                    fh.write(chunk)
            with p.open("rb") as f:
                header=f.readline(1024).decode("utf8-sig","replace")
            if p.stat().st_size<1000 or "SITE_ID" not in header.upper() or "," not in header:
                raise ValueError("Official CSV content/schema failed: "+k+", header "+repr(header[:100]))
            return p
        except Exception as e:
            exc=e
            p.unlink(missing_ok=True)
            time.sleep(min(12,2**i))
    raise RuntimeError("Official source inaccessible after five attempts: "+k+", "+str(exc))

def rows(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        h=[s.strip().upper() for s in reader.fieldnames]
        for line in reader:
            yield {str(k).strip().upper():v for k,v in line.items() if k is not None}

def header(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        return [s.strip().upper() for s in next(csv.reader(f))]

def key(row):
    sid=str(row.get("SITE_ID") or "").strip().upper()
    visit=str(row.get("VISIT_NO") or "").strip().upper()
    if not sid: return None
    # The canonical first-visit key is a site-level unit. UID/UNIQUE_ID audited as alternatives.
    if not visit: visit="1"
    return (sid,visit)

def numeric_gt0(x):
    if not valid(x):return False
    try:
        y=float(x)
        return math.isfinite(y) and y>0
    except (ValueError,TypeError):return False

out={"phase":"SOURCE_ONLY_PRE_NUMERIC_BIOLOGICAL_ASSOCIATION","not_population_benefit":True,"sources":{}}
paths={}
for k,url in URLS.items():
    print("Fetching official EPA "+k+" original...",flush=True)
    paths[k]=get(k,url)
    raw=paths[k].read_bytes()
    out["sources"][k]={"url":url,"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest(),"columns":header(paths[k])}

# Observed site samples: valid river site/visit metadata. Country, state, ecoregion and
# lat/long column availability are determined independently of nutrient outcomes.
site_counts=collections.Counter()
site_rows={}
for row in rows(paths["sites"]):
    k=key(row)
    if not k:continue
    site_counts[k]+=1
    if k not in site_rows:site_rows[k]=row

# EPA chemistry visit has no zero/imputation rule; require positive finite N/P numbers
# solely for SOURCE SUPPORT. Their values are not stored or plotted.
chem_counts=collections.Counter()
chem_ok=collections.Counter()
chem_units=collections.Counter()
chem_n_ok, chem_p_ok=set(),set()
for row in rows(paths["water"]):
    k=key(row)
    if not k:continue
    chem_counts[k]+=1
    if numeric_gt0(row.get("NTL_RESULT")):chem_n_ok.add(k)
    if numeric_gt0(row.get("PTL_RESULT")):chem_p_ok.add(k)
    if numeric_gt0(row.get("NTL_RESULT")) and numeric_gt0(row.get("PTL_RESULT")):
        chem_ok[k]+=1
        chem_units[(str(row.get("NTL_UNITS","")).strip(),str(row.get("PTL_UNITS","")).strip())]+=1

# Taxonomic SOURCE metadata only: sample type, positive numeric subsample,
# existence of Order / Family and functional feeding group.
type_rows=collections.Counter()
type_sites=collections.defaultdict(set)
type_countpositive=collections.Counter()
type_ffg_countpositive=collections.Counter()
type_orderfamily_countpositive=collections.Counter()
protocol_rows=collections.defaultdict(list)
for row in rows(paths["benthos"]):
    s=str(row.get("SAMPLE_TYPE") or "").strip().upper()
    k=key(row)
    if not k:continue
    type_rows[s]+=1;type_sites[s].add(k)
    if numeric_gt0(row.get("TOTAL300")):
        type_countpositive[s]+=1
        if valid(row.get("FFG")):type_ffg_countpositive[s]+=1
        if valid(row.get("ORDER")) or valid(row.get("FAMILY")):type_orderfamily_countpositive[s]+=1
        if s=="BERW":protocol_rows[k].append(row)
# The agreement assessment uses row support first, never biological sign.
primary="BERW"
berw_sites=set(protocol_rows)
common=set(site_rows)&set(chem_ok)&berw_sites
# Keep case identities internal to pipeline artifact (no site IDs in stdout or repo).
unique_states=set()
unique_ecos=set()
for k in common:
    s=site_rows[k]
    one=protocol_rows[k][0]
    state=str(s.get("STATE") or one.get("STATE") or "").strip()
    eco=str(s.get("AG_ECO9") or one.get("AG_ECO9") or "").strip()
    if valid(state):unique_states.add(state)
    if valid(eco):unique_ecos.add(eco)
total=type_countpositive[primary]
ffg_frac=type_ffg_countpositive[primary]/total if total else 0
taxon_frac=type_orderfamily_countpositive[primary]/total if total else 0

# SOURCE support gate: could be >300 if measured chemistry and one subtype share keys.
# Keeping potential many-to-one joins from being interpreted as independent samples.
out["source_support"]={
 "official_source_gate_only":True,
 "n_sites_metadata_unique":len(site_rows),
 "n_duplicate_site_metadata_visit_keys":sum(v>1 for v in site_counts.values()),
 "n_water_unique_visit_keys":len(chem_counts),
 "n_water_visit_keys_with_both_N_and_P_positive":len(chem_ok),
 "n_water_duplicate_visit_keys":sum(v>1 for v in chem_counts.values()),
 "n_primary_berw_unique_sites":len(berw_sites),
 "n_primary_berw_site_metadata_water_joint":len(common),
 "joint_support_distinct_states":len(unique_states),
 "joint_support_distinct_ecoregions":len(unique_ecos),
 "source_types":{k:{"rows":v,"sites":len(type_sites[k]),"positive_total300_taxon_rows":type_countpositive[k]} for k,v in type_rows.items()},
 "berw_positive_taxon_rows_FFG_covered_fraction":round(ffg_frac,6),
 "berw_positive_taxon_rows_order_or_family_covered_fraction":round(taxon_frac,6),
 "observed_NTL_PTL_positive_field_units": [{"TN":n,"TP":p,"rows":m} for (n,p),m in sorted(chem_units.items())],
 "site_geography_candidate_column_names":[c for c in out["sources"]["sites"]["columns"] if any(x in c for x in ("LAT","LON","ECO","STATE","UID","COORD","LOC","REG"))],
}
out["source_support"]["precommitted_gate_pass"]=bool(
 len(common)>=300 and len(unique_states)>=25 and len(unique_ecos)>=6
 and ffg_frac>=.60 and taxon_frac>=.80
 and all(v==1 for k,v in site_counts.items() if k in common)
)
out["decision"]="PASS_TO_FROZEN_SOURCE_MANIFEST_ONLY" if out["source_support"]["precommitted_gate_pass"] else "STOP_NO_RESPONSE_ASSOCIATION"
(out["source_support"]["source_note"])="An NRSA positive result is chemistry source availability, not a wastewater discharge or delivered nutrient dose"
with (OUT/"NRSA_1819_SOURCE_ONLY_QA.json").open("w") as f:
    json.dump(out,f,indent=2,ensure_ascii=False)
print("NRSA_1819_SOURCE_ONLY_QA_JSON_START",flush=True)
print(json.dumps(out,ensure_ascii=False,sort_keys=True),flush=True)
print("NRSA_1819_SOURCE_ONLY_QA_JSON_END",flush=True)
print("Output artifact: "+str(OUT/"NRSA_1819_SOURCE_ONLY_QA.json"),flush=True)
