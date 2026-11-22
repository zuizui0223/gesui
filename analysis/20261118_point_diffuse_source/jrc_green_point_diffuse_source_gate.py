#!/usr/bin/env python3
"""Source-only census of JRC GREEN modelled point and diffuse N/P.

No biological outcome, WWTP-specific dose, deposition-specific exposure, or
causal effect is estimated. Download the entire native attribute layers before
constructing any matched source panel. Hard-fail on truncated pages/duplicates.
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
from datetime import datetime, timezone
from pathlib import Path

BASE = "https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/Losses_to_water_MIS/MapServer"
LAYERS = {
    "point_N": (5, "PsN", "point"),
    "point_P": (34, "PsP", "point"),
    "diffuse_N": (16, "De", "diffuse"),
    "diffuse_P": (45, "De", "diffuse"),
}
HEADERS = {"User-Agent": "gesui-JRC-GREEN-source-first-point-diffuse-20261118", "Accept": "application/json"}
PAGE = 1000


def request(url: str, params: dict) -> tuple[dict, str]:
    full = url + "?" + urllib.parse.urlencode({"f": "json", **params})
    for attempt in range(4):
        try:
            req = urllib.request.Request(full, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=70) as response:
                raw = response.read(12_000_000)
            obj = json.loads(raw)
            if not isinstance(obj, dict) or obj.get("error"):
                raise ValueError(f"ArcGIS source error: {obj.get('error') if isinstance(obj, dict) else type(obj)}")
            return obj, hashlib.sha256(raw).hexdigest()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, ValueError) as exc:
            if attempt == 3:
                raise RuntimeError(f"Source retrieval failed ({url}): {exc}") from exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("Unreachable")


def fetch_layer(key: str) -> tuple[list[dict], dict]:
    layer_id, field, grain = LAYERS[key]
    url = f"{BASE}/{layer_id}"
    meta, meta_hash = request(url, {})
    fields = {x.get("name") for x in meta.get("fields", [])}
    required = {"OBJECTID", "BasinID", field}
    if grain == "point":
        required.add("AreaSqKm")
    else:
        required.add("HydroID")
    if not required.issubset(fields) or meta.get("type") != "Feature Layer":
        raise ValueError(f"Original JRC source schema changed for {key}: missing {required - fields}")
    count_response, count_hash = request(url + "/query", {"where": "1=1", "returnCountOnly": "true"})
    ids_response, ids_hash = request(url + "/query", {"where": "1=1", "returnIdsOnly": "true"})
    count = count_response.get("count")
    ids = ids_response.get("objectIds")
    if type(count) is not int or count <= 0 or count > 2_000_000 or not isinstance(ids, list):
        raise ValueError(f"Unreliable JRC source count for {key}")
    if len(ids) != count or len(set(ids)) != count:
        raise ValueError(f"Original source count/ID inventory mismatch for {key}: {count} vs {len(ids)}")
    sorted_ids = sorted(ids)
    rows = []
    page_hashes = []
    for offset in range(0, count, PAGE):
        q = {
            "where": "1=1",
            "outFields": ",".join(sorted(required)),
            "returnGeometry": "false",
            "orderByFields": "OBJECTID ASC",
            "resultOffset": offset,
            "resultRecordCount": min(PAGE, count - offset),
        }
        data, stamp = request(url + "/query", q)
        page = [x.get("attributes", {}) for x in data.get("features", [])]
        expected_ids = sorted_ids[offset: offset + len(page)]
        if len(page) != min(PAGE, count - offset) or [x.get("OBJECTID") for x in page] != expected_ids:
            raise ValueError(f"Truncated/reordered source page for {key} at offset {offset}")
        rows.extend(page)
        page_hashes.append(stamp)
        print(f"{key}: {len(rows)}/{count}", flush=True)
        time.sleep(0.08)
    if len(rows) != count:
        raise ValueError(f"Incomplete JRC layer {key}")
    return rows, {
        "layer_id": layer_id,
        "original_count": count,
        "meta_sha256": meta_hash,
        "count_sha256": count_hash,
        "object_ids_sha256": ids_hash,
        "attribute_pages": len(page_hashes),
        "all_attribute_page_hashes_sha256": hashlib.sha256("\n".join(page_hashes).encode()).hexdigest(),
    }


def valid_mass(x):
    try:
        v = float(x)
    except (ValueError, TypeError):
        return None
    return v if math.isfinite(v) and v >= 0 else None


def normalize(rows: list[dict], kind: str, nutrient: str) -> tuple[dict, dict, dict]:
    """Return {BasinID: mass or None}, {BasinID: AreaSqKm}, diagnostics."""
    mass_field = f"Ps{nutrient}" if kind == "point" else "De"
    vals = {}
    areas = {}
    seen_hydro = set()
    missing_mass = set()
    source_basins = set()
    for r in rows:
        basin_id = r.get("BasinID")
        if type(basin_id) is not int:
            raise ValueError("BasinID not a valid original integer")
        source_basins.add(basin_id)
        if kind == "point":
            if basin_id in vals:
                raise ValueError(f"Duplicate native point basin ID {basin_id} (no implicit summing)")
            vals[basin_id] = valid_mass(r.get(mass_field))
            area = valid_mass(r.get("AreaSqKm"))
            areas[basin_id] = area if area is not None and area > 0 else None
        else:
            hydro_id = r.get("HydroID")
            if type(hydro_id) is not int or hydro_id in seen_hydro:
                raise ValueError("Missing or duplicated diffuse HydroID")
            seen_hydro.add(hydro_id)
            v = valid_mass(r.get(mass_field))
            if v is None:
                missing_mass.add(basin_id)
            else:
                vals[basin_id] = vals.get(basin_id, 0.0) + v
    if kind == "diffuse":
        for bid in missing_mass:
            vals[bid] = None
        for bid in source_basins:
            vals.setdefault(bid, None)
    return vals, areas, {
        "native_rows": len(rows),
        "unique_BasinID": len(source_basins),
        "missing_mass_basins": sum(vals[b] is None for b in source_basins),
        "true_reported_zero_basins": sum(vals[b] == 0 for b in source_basins),
    }


def pearson_log(x, y):
    pairs = [(math.log1p(a), math.log1p(b)) for a, b in zip(x, y) if a is not None and b is not None]
    if len(pairs) < 3:
        return None
    n = len(pairs)
    ax = sum(a for a, _ in pairs) / n
    by = sum(b for _, b in pairs) / n
    xx = sum((a - ax) ** 2 for a, _ in pairs)
    yy = sum((b - by) ** 2 for _, b in pairs)
    return sum((a - ax) * (b - by) for a, b in pairs) / math.sqrt(xx * yy) if xx > 0 and yy > 0 else None


def make_panel(layers: dict[str, tuple[dict, dict, dict]]) -> tuple[list[dict], dict]:
    all_basin_ids = sorted(set().union(*(set(x[0]) for x in layers.values())))
    panel = []
    for bid in all_basin_ids:
        row = {"BasinID": bid}
        for key in LAYERS:
            row[key + "_modelled_tonnes_2014_2018"] = layers[key][0].get(bid)
        narea = layers["point_N"][1].get(bid)
        parea = layers["point_P"][1].get(bid)
        if narea is not None and parea is not None and abs(narea - parea) > max(0.01, 0.005 * narea):
            raise ValueError(f"Point N/P basin area mismatch {bid}")
        area = narea if narea is not None else parea
        row["native_JRC_basin_area_km2"] = area
        row["full_four_source_support"] = all(row[key + "_modelled_tonnes_2014_2018"] is not None for key in LAYERS) and area is not None
        panel.append(row)
    complete = [r for r in panel if r["full_four_source_support"]]
    four_both_positive = sum(all(r[k + "_modelled_tonnes_2014_2018"] > 0 for k in LAYERS) for r in complete)
    positive_pairs = {
        n: sum(r[f"point_{n}_modelled_tonnes_2014_2018"] > 0 and r[f"diffuse_{n}_modelled_tonnes_2014_2018"] > 0 for r in complete)
        for n in ("N", "P")
    }
    def xs(key):
        return [r[key + "_modelled_tonnes_2014_2018"] for r in complete]
    corrs = {
        "point_vs_diffuse_N": pearson_log(xs("point_N"), xs("diffuse_N")),
        "point_vs_diffuse_P": pearson_log(xs("point_P"), xs("diffuse_P")),
        "point_N_vs_point_P": pearson_log(xs("point_N"), xs("point_P")),
        "diffuse_N_vs_diffuse_P": pearson_log(xs("diffuse_N"), xs("diffuse_P")),
    }
    criteria = {
        "at_least_100_complete_distinct_native_basins": len(complete) >= 100,
        "at_least_30_basins_with_both_sources_positive_for_N": positive_pairs["N"] >= 30,
        "at_least_30_basins_with_both_sources_positive_for_P": positive_pairs["P"] >= 30,
    }
    return panel, {
        "all_union_basins": len(panel),
        "native_basins_with_all_4_modelled_source_fields_and_area": len(complete),
        "source_only_N_and_P_both_positive_basins": positive_pairs,
        "source_only_all_4_positive_basins": four_both_positive,
        "log1p_mass_pearson_correlations_descriptive": corrs,
        "source_only_gate_criteria": criteria,
        "source_only_native_BasinID_crosswalk_pass": all(criteria.values()),
        "EQR_butterfly_plant_source_match_pass": False,
        "no_ecological_effect_test_permitted": True,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    originals = {}
    receipts = {}
    for key, (_, _, grain) in LAYERS.items():
        rows, receipt = fetch_layer(key)
        originals[key] = normalize(rows, grain, key[-1])
        receipts[key] = {**receipt, **originals[key][2]}
    panel, support = make_panel(originals)
    out = args.output_dir / "jrc_GREEN_2014_18_native_basin_source_only.csv"
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(panel[0]), extrasaction="raise")
        writer.writeheader()
        writer.writerows(panel)
    report = {
        "schema": "gesui_jrc_GREEN_2014to18_native_BasinID_point_vs_diffuse_source_first_v1",
        "run_UTC": datetime.now(timezone.utc).isoformat(),
        "official_JRC_MapServer": BASE,
        "model_period": "2014-2018 average; not annual measured WWTP discharge",
        "point_definition": "GREEN modelled domestic plus industrial; NOT WWTP-only",
        "diffuse_definition": "GREEN modelled composite diffuse losses; NOT atmospheric-only or agriculture-only",
        "spatial_definition": "JRC native BasinID; NOT HydroBASINS H7 or local receiving reach",
        "units": "modelled tonnes averaged over 2014-2018, not concentrations or delivered doses",
        "layers": receipts,
        "source_only_support": support,
        "panel_sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
        "NO_outcome_or_species_or_EEA_EQR_values_read": True,
        "NEXT_GATE": "Independent source-only spatial alignment of exact JRC BasinID polygons to native EEA 2016 and GBIF H7 cohorts, with counts by year/country/independent basin BEFORE any outcome fit",
    }
    dst = args.output_dir / "jrc_GREEN_point_diffuse_source_only_gate.json"
    dst.write_text(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(support, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
