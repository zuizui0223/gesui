#!/usr/bin/env python3
"""Step 1 of the WWTP-upgrade source gate: acquire and inventory UWWTD plant tables.

Source only. Never reads butterfly, plant or benthic outcomes.
For each reporting cycle: list the public WebDAV folder, download plant-level
(UWWTP) tables, SHA-256 lock them, and record whether the fields the frozen
protocol needs are present. Field names are NOT assumed equal across cycles.
"""
import csv, hashlib, io, json, re, sys, time, urllib.request, zipfile
from pathlib import Path
from urllib.parse import quote, unquote
from xml.etree import ElementTree as ET

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
ROOT = "https://sdi.eea.europa.eu/webdav/datastore/public/"
CYCLES = {  # data call -> (catalogue record, datastore id, temporal coverage)
    "2014": ("d6dbe6a3-49d1-4987-b55e-94dd680d0254", "eea_t_uwwtd-waterbase_p_2013-2014_v05_r00", "2013-2014"),
    "2016": ("36115534-1c21-4c5e-9ec8-2ad22e5ffbb7", "eea_t_uwwtd-waterbase_p_2015-2016_v06_r00", "2015-2016"),
    "2019": ("dbbccc77-ec41-4b19-8ec5-093bba1ba7a8", "eea_t_uwwtd-waterbase_p_2017-2018_v07_r01", "2017-2018"),
    "2021": ("a55c14cd-aa37-4045-892d-2c3ea36abaff", "eea_t_uwwtd-final_p_2019-2020_v08_r00", "2019-2020"),
    "2023": ("52b2e779-a146-414f-bf00-d63dfbf9c4f1", "eea_t_uwwtd-final_p_2023_v10_r00", "2021-2022/2022-2023"),
}
NEED = {
    "plant_code": ["uwwCode"], "reporter": ["repCode", "rptMStateKey"],
    "N_removal": ["uwwNRemoval"], "P_removal": ["uwwPRemoval"],
    "load_entering_pe": ["uwwLoadEnteringUWWTP"], "capacity": ["uwwCapacity"],
    "secondary": ["uwwSecondaryTreatment"], "other_treatment": ["uwwOtherTreatment"],
    "lat": ["uwwLatitude"], "lon": ["uwwLongitude"], "state": ["uwwState"],
    "closing": ["uwwDateClosing"], "history": ["uwwHistorie", "uwwHistory"],
    "inspire_id": ["uwwInspireIDFacility"], "begin_life": ["uwwBeginLife"],
}
UA = {"User-Agent": "gesui-wwtp-upgrade-source-gate"}

def open_url(url, method="GET", headers=None, timeout=120, tries=3):
    last = None
    for k in range(tries):
        try:
            h = dict(UA); h.update(headers or {})
            return urllib.request.urlopen(urllib.request.Request(url, headers=h, method=method), timeout=timeout)
        except Exception as ex:
            last = ex; time.sleep(10 * (k + 1))
    raise last

def listing(dsid):
    with open_url(ROOT + dsid + "/", "PROPFIND", {"Depth": "infinity"}) as r:
        body = r.read()
    ns = {"d": "DAV:"}; out = []
    for e in ET.fromstring(body).findall("d:response", ns):
        href = unquote(e.find("d:href", ns).text)
        size = e.find(".//d:getcontentlength", ns)
        if size is not None and size.text:
            out.append({"href": href, "bytes": int(size.text)})
    return out

def is_plant_table(name):
    n = name.lower()
    return bool(re.search(r"(uwwtps?|treatmentplant|t_uwwtps)", n)) and not re.search(r"(agglo|disch|dcp|uwwtps_agg|emission)", n)

def header_check(cols):
    low = {c.strip().lower(): c.strip() for c in cols}
    return {k: next((low[v.lower()] for v in vs if v.lower() in low), None) for k, vs in NEED.items()}

def inspect_csv(raw, name):
    txt = raw.decode("utf-8-sig", errors="replace")
    dialect = csv.Sniffer().sniff(txt[:20000], delimiters=",;\t|")
    rows = list(csv.DictReader(io.StringIO(txt), dialect=dialect))
    cols = list(rows[0].keys()) if rows else []
    hc = header_check(cols); res = {"file": name, "rows": len(rows), "columns": cols, "fields": hc}
    for k in ("N_removal", "P_removal"):
        if hc[k]:
            vals = [str(r.get(hc[k], "")).strip() for r in rows]
            res[k + "_values"] = dict(sorted(__import__("collections").Counter(vals).most_common(12)))
    if hc["reporter"]:
        res["reporter_suffix"] = dict(__import__("collections").Counter(
            (str(r.get(hc["reporter"], ""))[-4:]) for r in rows).most_common(10))
    return res, rows

summary = {"protocol": "PROTOCOL.md v0.1", "outcomes_read": False, "cycles": {}}
for call, (rec, dsid, cov) in CYCLES.items():
    c = {"catalogue_record": rec, "datastore_id": dsid, "coverage": cov}
    try:
        files = listing(dsid); c["listing"] = files
        targets = [f for f in files if is_plant_table(f["href"].rsplit("/", 1)[-1]) or f["href"].lower().endswith(".zip")]
        c["inspected"] = []
        for f in targets[:12]:
            name = f["href"].rsplit("/", 1)[-1]
            url = "https://sdi.eea.europa.eu" + quote(f["href"]) if f["href"].startswith("/") else f["href"]
            with open_url(url, timeout=600) as r: raw = r.read()
            sha = hashlib.sha256(raw).hexdigest()
            (OUT / f"{call}_{name}").write_bytes(raw)
            entry = {"file": name, "bytes": len(raw), "sha256": sha}
            if name.lower().endswith(".zip"):
                with zipfile.ZipFile(io.BytesIO(raw)) as z:
                    entry["members"] = z.namelist()
                    for m in z.namelist():
                        if is_plant_table(m.rsplit("/", 1)[-1]) and m.lower().endswith((".csv", ".txt")):
                            res, _ = inspect_csv(z.read(m), m); entry.setdefault("plant_tables", []).append(res)
            elif name.lower().endswith((".csv", ".txt")):
                res, _ = inspect_csv(raw, name); entry["plant_tables"] = [res]
            c["inspected"].append(entry)
    except Exception as ex:
        c["error"] = repr(ex)[:400]
    summary["cycles"][call] = c
    print(call, dsid, "ERROR " + c["error"] if "error" in c else f"files={len(c['listing'])} inspected={len(c['inspected'])}", flush=True)
    for e in c.get("inspected", []):
        for t in e.get("plant_tables", []):
            print("  ", t["file"], t["rows"], {k: v for k, v in t["fields"].items()}, t.get("N_removal_values"), flush=True)

(OUT / "UWWTD_CYCLE_INVENTORY.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
