"""Aggregate source support and TP availability; never inspect TP values.

The candidate mapping is historical UWWTD outlets to current EA SA points.
Its hash is a reproducibility fingerprint, not a verified facility-ID join.
Every support count is descriptive: no eligibility threshold, weight, model,
outcome summary, or missing-value imputation is introduced here.
"""
from __future__ import annotations

import hashlib
import json
import math
import statistics


CAPACITY_BANDS = (
    "<2000", "2000_to_9999", "10000_to_99999", "100000plus",
    "missing_or_nonpositive",
)
NUTS1_REGIONS = frozenset("UK" + suffix for suffix in "CDEFGHIJKLMN")
SITE_FATES = (
    "no_returned_TP_record", "no_eligible_numeric_TP", "eligible_numeric_TP",
)
FLAGS = ("0", "1")


def _finite_number(value):
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return number if math.isfinite(number) else None


def _capacity_band(capacity):
    if capacity is None:
        return "missing_or_nonpositive"
    if capacity < 2000:
        return "<2000"
    if capacity < 10000:
        return "2000_to_9999"
    if capacity < 100000:
        return "10000_to_99999"
    return "100000plus"


def _counts(candidate_n, observed_n):
    return {
        "candidate_n": candidate_n,
        "eligible_numeric_TP_n": observed_n,
        "not_eligible_numeric_TP_n": candidate_n - observed_n,
        "eligible_numeric_TP_rate": observed_n / candidate_n if candidate_n else None,
    }


def _table(rows, flag):
    return {
        value: _counts(
            sum(row[flag] == value for row in rows),
            sum(row[flag] == value and row["observed"] for row in rows),
        )
        for value in FLAGS
    }


def _summary(values):
    valid = [value for value in values if value is not None]
    return {
        "count": len(valid),
        "missing_n": len(values) - len(valid),
        "min": min(valid) if valid else None,
        "median": float(statistics.median(valid)) if valid else None,
        "max": max(valid) if valid else None,
    }


def _summaries(rows, flag, field):
    return {
        population: {
            value: _summary([
                row[field] for row in rows
                if row[flag] == value and (population == "candidates" or row["observed"])
            ])
            for value in FLAGS
        }
        for population in ("candidates", "eligible_numeric_TP")
    }


def _support(cells):
    result = {}
    for population, count_field in (
        ("candidates", "candidate_n"),
        ("eligible_numeric_TP", "eligible_numeric_TP_n"),
    ):
        covered = {value: 0 for value in FLAGS}
        lacking = {value: 0 for value in FLAGS}
        any_site, both_flags = 0, 0
        for cell in cells:
            counts = {value: cell[value][count_field] for value in FLAGS}
            any_site += any(counts.values())
            has_both = all(counts.values())
            both_flags += has_both
            destination = covered if has_both else lacking
            for value in FLAGS:
                destination[value] += counts[value]
        result[population] = {
            "cells_with_any_site": any_site,
            "cells_with_both_flags": both_flags,
            "covered_n": sum(covered.values()),
            "lacking_support_n": sum(lacking.values()),
            "covered_by_flag": covered,
            "lacking_support_by_flag": lacking,
        }
    return result


def build_source_support(matches, plants, tps, site_fates=None) -> dict:
    """Return only counts, availability rates, source summaries and a hash.

    ``tps`` is used only for its key set and membership. A site with no
    eligible numeric TP remains in every candidate denominator; its TP is
    neither read nor imputed. ``site_fates``, when supplied, must cover all
    candidates and agree with TP membership. Per-site record counts in that
    optional mapping are not interpreted by this helper.

    Both treatment flags are audited independently. Region and capacity
    support means at least one source with each flag in the same displayed
    stratum, assessed separately among all candidates and eligible TP sites.
    Missing region/capacity categories remain explicit descriptive strata.
    """
    candidate_keys = set(matches)
    observed_keys = set(tps)
    if not observed_keys <= candidate_keys:
        raise ValueError("eligible TP keys outside candidate mapping")
    if site_fates is not None and set(site_fates) != candidate_keys:
        raise ValueError("site fate keys must equal candidate mapping keys")

    seen_plants, seen_outlets = set(), set()
    mapping_rows, rows = [], []
    for sa, match in matches.items():
        native = match.get("outlet_id")
        if not isinstance(sa, str) or not sa or not isinstance(native, tuple) or len(native) < 2:
            raise ValueError("invalid native candidate identity")
        if not all(isinstance(part, str) and part for part in native):
            raise ValueError("invalid native outlet identity")
        if native in seen_outlets:
            raise ValueError("duplicate native outlet among candidates")
        plant_id = native[1]
        if plant_id in seen_plants:
            raise ValueError("duplicate native plant among candidates")
        if plant_id not in plants or plants[plant_id] is None:
            raise ValueError("missing native plant for candidate outlet")
        seen_outlets.add(native)
        seen_plants.add(plant_id)
        plant = plants[plant_id]
        for field in ("PRemoval", "NRemoval"):
            if plant.get(field) not in (0, 1):
                raise ValueError("unrecognized treatment flag")

        coords = match.get("outlet_coord")
        if not isinstance(coords, (tuple, list)) or len(coords) != 2:
            raise ValueError("invalid outlet coordinates")
        coords = [_finite_number(value) for value in coords]
        if any(value is None for value in coords):
            raise ValueError("nonfinite outlet coordinates")
        distance = _finite_number(match.get("original2018_FW_outlet_to_SA_distance_m"))
        if distance is None or distance < 0:
            raise ValueError("invalid match distance")
        capacity = _finite_number(plant.get("capacity"))
        if capacity is not None and capacity <= 0:
            capacity = None
        region = plant.get("original_region")
        if region is not None and region not in NUTS1_REGIONS:
            raise ValueError("unrecognized original NUTS1 region")
        observed = sa in observed_keys
        fate = None
        if site_fates is not None:
            fate = site_fates[sa].get("fate")
            if fate not in SITE_FATES:
                raise ValueError("unrecognized site fate")
            if (fate == "eligible_numeric_TP") != observed:
                raise ValueError("site fate inconsistent with eligible TP membership")
        rows.append({
            "PRemoval": str(int(plant["PRemoval"])),
            "NRemoval": str(int(plant["NRemoval"])),
            "region": region if region is not None else "missing",
            "capacity_band": _capacity_band(capacity),
            "capacity": capacity, "distance": distance,
            "observed": observed, "fate": fate,
        })
        mapping_rows.append([sa, list(native), coords])

    canonical_mapping = json.dumps(
        sorted(mapping_rows, key=lambda row: row[0]),
        ensure_ascii=False, separators=(",", ":"), allow_nan=False,
    ).encode("utf-8")
    result = {
        "schema": "GESUI_PREMOVAL_SOURCE_SUPPORT_V1",
        "source_key_counts": {
            "matches": len(candidate_keys), "plants": len(plants),
            "tps": len(observed_keys),
            "site_fates": len(site_fates) if site_fates is not None else 0,
            "candidate_native_plants": len(seen_plants),
            "candidate_native_outlets": len(seen_outlets),
        },
        "totals": _counts(len(candidate_keys), len(observed_keys)),
        "pair_mapping_sha256": hashlib.sha256(canonical_mapping).hexdigest(),
        "pair_mapping_sha256_canonicalization": (
            "UTF-8 compact JSON, ensure_ascii=False, sorted by native SA string; "
            "each row [SA, native outlet tuple as list, [BNG easting float, northing float]]"
        ),
        "support_definition": (
            "Descriptive cells containing at least one source with each flag, "
            "separately among candidates and eligible numeric TP sites; "
            "missing categories retained; no new source gate or eligibility threshold"
        ),
        "by_treatment_flag": {},
    }
    regions = sorted({row["region"] for row in rows})
    for flag in ("PRemoval", "NRemoval"):
        by_region = {
            region: _table([row for row in rows if row["region"] == region], flag)
            for region in regions
        }
        by_capacity = {
            band: _table([row for row in rows if row["capacity_band"] == band], flag)
            for band in CAPACITY_BANDS
        }
        by_joint = {
            region: {
                band: _table([
                    row for row in rows
                    if row["region"] == region and row["capacity_band"] == band
                ], flag)
                for band in CAPACITY_BANDS
            }
            for region in regions
        }
        group = {
            "by_flag": _table(rows, flag),
            "by_flag_and_NUTS1": by_region,
            "by_flag_and_capacity_band": by_capacity,
            "by_flag_and_NUTS1_and_capacity_band": by_joint,
            "both_flag_support": {
                "NUTS1": _support(list(by_region.values())),
                "capacity_band": _support(list(by_capacity.values())),
                "NUTS1_x_capacity_band": _support([
                    cell for bands in by_joint.values() for cell in bands.values()
                ]),
            },
            "capacity_PE_summary": _summaries(rows, flag, "capacity"),
            "outlet_to_SA_distance_m_summary": _summaries(rows, flag, "distance"),
        }
        if site_fates is not None:
            group["site_fates_by_flag"] = {
                value: {
                    fate: sum(row[flag] == value and row["fate"] == fate for row in rows)
                    for fate in SITE_FATES
                }
                for value in FLAGS
            }
        result["by_treatment_flag"][flag] = group
    return result
