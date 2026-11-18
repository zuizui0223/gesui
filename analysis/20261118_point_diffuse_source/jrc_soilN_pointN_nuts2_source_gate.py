#!/usr/bin/env python3
"""Audit official JRC modelled water-loss sources versus soil N inputs at NUTS2.

This is a source-only spatial coverage / collinearity census, not a biotic model
or measured land/water nutrient delivery. No EQR, EIV or butterfly data read.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

WATER = "https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/Losses_to_water_MIS/MapServer"
SOIL = "https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/MIS_Inputs_to_soil/MapServer"
COUNTRY = "giscoEsriNuts2016_NUTS2_3035_CNTR_CODE"
NUTS = "giscoEsriNuts2016_NUTS2_3035_NUTS_ID"
SOURCES = {
    "water_point_N_t": (WATER, 4, ("cntrCode", "Nuts2"), "PsN"),
    "water_point_P_t": (WATER, 33, ("cntrCode", "Nuts2"), "PsP"),
    "land_atmos_N_mg_m2": (SOIL, 4, (COUNTRY, NUTS), "MEAN_mgsqmy"),
    "land_mineral_N_kg_ha": (SOIL, 9, (COUNTRY, NUTS), "MEAN_kgha"),
    "land_organic_N_kg_ha": (SOIL, 14, (COUNTRY, NUTS), "MEAN_kgha"),
    "land_fixation_N_kg_ha": (SOIL, 19, (COUNTRY, NUTS), "MEAN_kgha"),
}
REQHEAD = {"User-Agent": "gesui-nuts2-soil-input-source-census-20261118", "Accept": "application/json"}


def retrieve(url, params):
    full = url + "?" + urllib.parse.urlencode({"f": "json", **params})
    for i in range(4):
        try:
            req = urllib.request.Request(full, headers=REQHEAD)
            with urllib.request.urlopen(req, timeout=65) as res:
                raw = res.read(2_500_000)
            data = json.loads(raw)
            if not isinstance(data, dict) or "error" in data:
                raise ValueError(f"Source error for {url}: {data}")
            return data, hashlib.sha256(raw).hexdigest()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, ValueError) as exc:
            if i == 3:
                raise RuntimeError(f"Official source unavailable at {url}: {exc}") from exc
            time.sleep(i + 1)
    raise RuntimeError("unreachable")


def numeric(v):
    try:
        x = float(v)
    except (ValueError, TypeError):
        return None
    return x if math.isfinite(x) and x >= 0 else None


def source_key(r, country_field, nuts_field):
    country = r.get(country_field)
    nuts = r.get(nuts_field)
    if not isinstance(country, str) or not isinstance(nuts, str):
        return None
    country, nuts = country.strip().upper(), nuts.strip().upper()
    return f"{country}|{nuts}" if country and nuts else None


def parse_rows(rows, country_field, nuts_field, value_field):
    byregion = {}
    missing_ids = 0
    for row in rows:
        key = source_key(row, country_field, nuts_field)
        if key is None:
            missing_ids += 1
            continue
        if key in byregion:
            raise ValueError(f"Duplicated native NUTS2 source key: {key}")
        byregion[key] = numeric(row.get(value_field))
    return byregion, {"rows": len(rows), "regions": len(byregion), "missing_key_rows": missing_ids,
                       "numeric_regions": sum(x is not None for x in byregion.values()),
                       "true_zero_regions": sum(x == 0 for x in byregion.values())}


def source(layer, cfields, value):
    meta, hashmeta = retrieve(layer, {})
    needed = ["OBJECTID", *cfields, value]
    present = {x.get("name") for x in meta.get("fields", [])}
    if not set(needed).issubset(present) or meta.get("type") != "Feature Layer":
        raise ValueError(f"JRC layer schema changed: {layer}, missing {set(needed)-present}")
    idsobj, hashids = retrieve(layer + "/query", {"where": "1=1", "returnIdsOnly": "true"})
    countobj, hashcount = retrieve(layer + "/query", {"where": "1=1", "returnCountOnly": "true"})
    ids = idsobj.get("objectIds")
    n = countobj.get("count")
    if not isinstance(ids, list) or type(n) is not int or n < 1 or n > 5000:
        raise ValueError(f"Invalid JRC NUTS2 full source count: {layer}, {n}")
    if len(ids) != n or len(set(ids)) != n:
        raise ValueError(f"JRC objectID discrepancy: {layer}")
    rows, pages = [], []
    for pos in range(0, n, 80):
        batch = sorted(ids)[pos:pos+80]
        obj, h = retrieve(layer + "/query", {"objectIds": ",".join(map(str,batch)),
                                            "outFields": ",".join(needed),
                                            "returnGeometry": "false"})
        chunk = [x.get("attributes", {}) for x in obj.get("features", [])]
        if len(chunk) != len(batch) or {x.get("OBJECTID") for x in chunk} != set(batch):
            raise ValueError(f"JRC original ObjectID page incomplete: {layer}, offset {pos}")
        pages.append(h)
        rows.extend(chunk)
    out, qc = parse_rows(rows, *cfields, value)
    return out, {**qc, "layer": layer, "source_meta_sha256": hashmeta,
                 "source_ids_sha256": hashids, "source_count_sha256": hashcount,
                 "page_response_sha256": hashlib.sha256("\n".join(pages).encode()).hexdigest()}


def correlations(x, y):
    pts = [(math.log1p(a), math.log1p(b)) for a,b in zip(x,y) if a is not None and b is not None]
    if len(pts) <= 2:
        return None
    mx = sum(a for a,_ in pts)/len(pts)
    my = sum(b for _,b in pts)/len(pts)
    xx = sum((a-mx)**2 for a,_ in pts)
    yy = sum((b-my)**2 for _,b in pts)
    return sum((a-mx)*(b-my) for a,b in pts)/math.sqrt(xx*yy) if xx>0 and yy>0 else None


def join(sources):
    # Do not include soil-only units and silently fill unavailable water sources with zero.
    ids = sorted(set().union(*(set(v) for v in sources.values())))
    rows = []
    for key in ids:
        row = {"NUTS2_country_id": key}
        row.update({name: sources[name].get(key) for name in SOURCES})
        row["complete_all_six_modelled_sources"] = all(row[name] is not None for name in SOURCES)
        rows.append(row)
    complete = [r for r in rows if r["complete_all_six_modelled_sources"]]
    countries = {r["NUTS2_country_id"].split("|")[0] for r in complete}
    cvariability = {}
    for name in SOURCES:
        per_country = defaultdict(set)
        for r in complete:
            per_country[r["NUTS2_country_id"].split("|")[0]].add(r[name])
        cvariability[name] = sum(len(vals)>1 for vals in per_country.values())
    correlations_by_land_source = {
        name: correlations([r["water_point_N_t"] for r in complete], [r[name] for r in complete])
        for name in SOURCES if name.startswith("land_")
    }
    gate = {
        "all_unique_unioned_nuts2":len(ids),
        "exact_source_supported_6_variables_nuts2":len(complete),
        "supported_country_n":len(countries),
        "supported_countries":sorted(countries),
        "countries_with_within_country_variation":cvariability,
        "source_only_log1p_pointN_vs_land_input_correlations":correlations_by_land_source,
        "source_gate_thresholds": {"six_source_NUTS2_at_least_100":len(complete)>=100,
                                  "supported_countries_at_least_6":len(countries)>=6,
                                  "each_land_source_varies_within_at_least_5_countries":all(
                                      cvariability[k]>=5 for k in SOURCES if k.startswith("land_"))},
        "crosswalk_source_only_pass": len(complete)>=100 and len(countries)>=6 and all(
            cvariability[k]>=5 for k in SOURCES if k.startswith("land_")),
        "biotic_source_overlap_and_ecological_effect_NOT_ASSESSED": True,
    }
    return rows, gate


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    values,receipts={},{}
    for name,(base,lid,fields,val) in SOURCES.items():
        values[name],receipts[name]=source(f"{base}/{lid}",fields,val)
        print(name,receipts[name]["regions"],receipts[name]["numeric_regions"],flush=True)
    rows,gate=join(values)
    csv_path=args.output_dir/"jrc_nuts2_soilN_vs_pointN_source_only.csv"
    with csv_path.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "schema": "gesui_JRC_NUTS2_original_water_point_soilN_separate_modelled_source_v1",
        "UTC": datetime.now(timezone.utc).isoformat(),
        "source_only_GEO_unit": "Native 2016 GISCO NUTS2 country-code|NUTS_ID; not HydroBASINS H7",
        "water_point_NP": "GREEN modelled domestic PLUS industrial losses to water, 2014-2018 average, tonnes over each NUTS2",
        "land_atmos_N": "EMEP modelled 2015 atmospheric deposition averaged mg/m2/year all land; not received leaf N",
        "land_agricultural_N": "DAYCENT modelled 2010-2019 mineral/organic fertilization, biological fixation averaged kg/ha/year on agricultural land",
        "cannot_sum_different_spatial_denominators_or_call_landed_loads_direct_NP_delivery": True,
        "source_qc":receipts,
        "source_gate":gate,
        "csv_sha256":hashlib.sha256(csv_path.read_bytes()).hexdigest(),
        "NO_EEA_EQR_EIV_GBIF_BIOLOGICAL_VALUES_READ":True,
    }
    (args.output_dir/"jrc_nuts2_soilN_pointN_source_only_report.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("JRC_NUTS2_SOURCE_GATE",json.dumps(gate,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
