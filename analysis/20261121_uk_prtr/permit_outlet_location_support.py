"""Version-linked outlet grid cells; public outputs contain aggregates only.

The grid width describes notation resolution, never surveyed accuracy or
historical coordinate validity. Spatial agreement never proves native identity.
"""
from __future__ import annotations

import collections
import math
import re

from permit_sample_identity_support import permit_key, text

FIELDS = ("OUTLET_GRID_REF", "DISCHARGE_NGR", "EFFLUENT_GRID_REF")
RADII = (100, 250, 500)


def grid_cell(value):
    """Decode an OS National Grid reference into its southwest corner and width."""
    raw = re.sub(r"\s+", "", text(value)).upper()
    if not re.fullmatch(r"[A-HJ-Z]{2}(?:[0-9]{2}){1,5}", raw):
        raise ValueError("Missing or invalid National Grid reference")
    alphabet = "ABCDEFGHJKLMNOPQRSTUVWXYZ"
    first, second = (alphabet.index(c) for c in raw[:2])
    east_tile = ((first - 2) % 5) * 5 + second % 5
    north_tile = 19 - (first // 5) * 5 - second // 5
    if not (0 <= east_tile < 7 and 0 <= north_tile < 13):
        raise ValueError("National Grid tile outside declared GB domain")
    digits = raw[2:]
    half = len(digits) // 2
    width = 10 ** (5 - half)
    east = east_tile * 100000 + int(digits[:half]) * width
    north = north_tile * 100000 + int(digits[half:]) * width
    return east, north, width


def distance_bounds(point, cell):
    try:
        x, y = (float(v) for v in point)
        east, north, width = (float(v) for v in cell)
    except (TypeError, ValueError):
        raise ValueError("Invalid geometric input") from None
    if not all(math.isfinite(v) for v in (x, y, east, north, width)) or width <= 0:
        raise ValueError("Nonfinite coordinate or invalid cell width")
    near_x = max(east - x, x - east - width, 0.0)
    near_y = max(north - y, y - north - width, 0.0)
    far_x = max(abs(x - east), abs(x - east - width))
    far_y = max(abs(y - north), abs(y - north - width))
    return math.hypot(near_x, near_y), math.hypot(far_x, far_y)


def distance_class(point, cell, radius):
    if not math.isfinite(radius) or radius < 0:
        raise ValueError("Invalid comparison radius")
    minimum, maximum = distance_bounds(point, cell)
    if maximum <= radius:
        return "within"
    if minimum > radius:
        return "outside"
    return "boundary"


def _cell_token(value):
    raw = re.sub(r"\s+", "", text(value)).upper()
    if not raw:
        return "missing", None
    try:
        return "valid", grid_cell(raw)
    except ValueError:
        # Retain distinct invalid strings privately so conflicts are not hidden.
        return "invalid", raw


def _state(tokens):
    if not tokens:
        return {"status": "absent_key"}
    if len(tokens) != 1:
        return {"status": "conflicting"}
    status, value = next(iter(tokens))
    return {"status": status, **({"cell": value} if status == "valid" else {})}


def location_index(rows, eligible_keys):
    eligible_keys = set(eligible_keys)
    outlets = {key[:4] for key in eligible_keys}
    key_values = collections.defaultdict(lambda: {field: set() for field in FIELDS})
    outlet_values = collections.defaultdict(set)
    source_rows = 0
    for row in rows:
        source_rows += 1
        key = permit_key(row, "consents")
        if key is None or key[:4] not in outlets:
            continue
        # Every effluent type on the same versioned outlet is checked, even if
        # it is not one of the final-treated eligible five-keys.
        outlet_values[key[:4]].add(_cell_token(row.get("OUTLET_GRID_REF")))
        if key in eligible_keys:
            for field in FIELDS:
                key_values[key][field].add(_cell_token(row.get(field)))
    outlet_states = {key: _state(outlet_values[key]) for key in outlets}
    index = {}
    for key in eligible_keys:
        fields = {field: _state(key_values[key][field]) for field in FIELDS}
        if (fields["OUTLET_GRID_REF"]["status"] == "valid"
                and outlet_states[key[:4]]["status"] != "valid"):
            fields["OUTLET_GRID_REF"] = {"status": "conflicting_versioned_outlet"}
        index[key] = fields
    qa = {
        "source_rows": source_rows,
        "eligible_distinct_keys": len(eligible_keys),
        "distinct_versioned_outlet_keys": len(outlets),
        "versioned_outlet_status_counts": dict(sorted(collections.Counter(
            value["status"] for value in outlet_states.values()).items())),
        "field_status_counts": {field: dict(sorted(collections.Counter(
            index[key][field]["status"] for key in eligible_keys).items())) for field in FIELDS},
        "all_effluent_types_checked_on_same_versioned_outlet": True,
        "coordinate_field_has_no_separate_survey_or_effective_date": True,
    }
    return index, qa


def summarize_locations(pairs, index, matching, registry):
    pairs = set(pairs)
    points = {point for point, _ in pairs}
    if len(points) != len(pairs):
        raise ValueError("A sampling point has multiple permit identities")
    if not points <= set(matching) or not points <= set(registry):
        raise ValueError("Original matched point geometry missing")
    if any(key not in index for _, key in pairs):
        raise ValueError("Required permit key missing from location index")
    statuses, widths = collections.Counter(), collections.Counter()
    comparisons = {label: collections.Counter() for label in ("outlet_vs_site_grid", "outlet_vs_effluent_grid")}
    distance = {label: {str(radius) + "m": {k: 0 for k in ("within", "outside", "boundary")}
                       for radius in RADII}
                for label in ("original_EEA_outlet", "current_WQE_samplepoint")}
    valid_outlets = set()
    points_by_outlet = collections.defaultdict(set)
    both_end_zero, both_end_double_zero = 0, 0
    for point, key in pairs:
        fields = index[key]
        outlet = fields["OUTLET_GRID_REF"]
        statuses[outlet["status"]] += 1
        points_by_outlet[key[:4]].add(point)
        if outlet["status"] != "valid":
            continue
        cell = outlet["cell"]
        valid_outlets.add(key[:4])
        widths[str(cell[2])] += 1
        if cell[2] == 1:
            both_end_zero += int(cell[0] % 10 == 0 and cell[1] % 10 == 0)
            both_end_double_zero += int(cell[0] % 100 == 0 and cell[1] % 100 == 0)
        for field, label in (("DISCHARGE_NGR", "outlet_vs_site_grid"),
                             ("EFFLUENT_GRID_REF", "outlet_vs_effluent_grid")):
            other = fields[field]
            fate = ("same_cell" if other.get("cell") == cell else "different_cells") if other["status"] == "valid" else other["status"]
            comparisons[label][fate] += 1
        for label, target in (("original_EEA_outlet", matching[point]["outlet_coord"]),
                              ("current_WQE_samplepoint", registry[point])):
            for radius in RADII:
                distance[label][str(radius) + "m"][distance_class(target, cell, radius)] += 1
    return {
        "candidate_points": len(points), "candidate_five_key_pairs": len(pairs),
        "distinct_versioned_outlets": len(points_by_outlet),
        "versioned_outlets_with_multiple_candidate_points": sum(len(v) > 1 for v in points_by_outlet.values()),
        "outlet_status_counts": dict(sorted(statuses.items())),
        "valid_location_points": statuses["valid"], "valid_distinct_versioned_outlets": len(valid_outlets),
        "recorded_grid_width_m_counts": dict(sorted(widths.items())),
        "one_metre_notation_with_both_axes_trailing_zero": both_end_zero,
        "one_metre_notation_with_both_axes_trailing_double_zero": both_end_double_zero,
        **{label: dict(sorted(value.items())) for label, value in comparisons.items()},
        "distances": distance,
        "distance_counts_are_points_not_independent_physical_outlets": True,
        "whole_cell_bounds_reflect_notation_not_survey_accuracy": True,
        "current_publication_attaches_geometry_to_historical_permit_version": True,
        "does_not_establish_2018_coordinate_survey_or_sample_date_validity": True,
        "EEA_native_outlet_to_permit_identity_still_unverified": True,
    }
