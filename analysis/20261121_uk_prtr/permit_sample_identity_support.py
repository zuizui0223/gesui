"""Exact historical identity support; all returned public summaries are aggregate.

Ofwat supplies complete sample notations in the same rows as permit keys.
The documented Midlands spelling difference is the only region alias used.
No source names, approximate IDs, coordinates, chemistry or flow are accepted.
"""
from __future__ import annotations

import collections
import datetime as dt
import re

EA_TO_WQE = {"AN": "AN", "MI": "MD", "NE": "NE", "NW": "NW",
             "SO": "SO", "SW": "SW", "TH": "TH"}
KEY_FIELDS = {
    "ofwat": ("EA_REGION", "permit_number", "version", "outlet_number", "effluent_number"),
    "consents": ("EA_REGION", "PERMIT_NUMBER", "PERMIT_VERSION", "OUTLET_NUMBER", "EFFLUENT_NUMBER"),
    "effluents": ("EA_REGION", "PERMIT_REF", "VERSION", "OUTLET_NUMBER", "EFFLUENT_NUMBER"),
}


def text(value):
    return "" if value is None else str(value).strip()


def permit_key(row, source, *, normalize_numbers=True):
    """Preserve opaque IDs; canonicalize only components that are integers.

    Original EA outlet identifiers include letters despite the column name
    OUTLET_NUMBER. A non-integer identifier is not a missing identifier.
    """
    values = [text(row.get(k)) for k in KEY_FIELDS[source]]
    if values[0] not in EA_TO_WQE or not all(values):
        return None
    for i in (2, 3, 4):
        if normalize_numbers and re.fullmatch(r"[0-9]+(?:\.0+)?", values[i]):
            values[i] = str(int(values[i].split(".")[0]))
    return tuple(values)


def ofwat_pairs(rows, year=2018):
    """Deduplicate a reported annual identity across all source sheets."""
    pairs = set()
    counts = collections.Counter()
    excluded = collections.Counter()
    prefix_pairs = collections.Counter()
    for row in rows:
        counts["input_rows"] += 1
        if text(row.get("year")) not in {str(year), str(year) + ".0"}:
            continue
        counts["year_rows"] += 1
        region = text(row.get("EA_REGION"))
        point = text(row.get("sample_samplingpoint_notation"))
        if region not in EA_TO_WQE:
            excluded["missing_or_unknown_English_region"] += 1
            continue
        if not point:
            excluded["missing_point_notation"] += 1
            continue
        prefix, separator, local = point.partition("-")
        if not separator or not local:
            excluded["not_a_complete_point_notation"] += 1
            continue
        if prefix != EA_TO_WQE[region]:
            excluded["point_region_mismatch"] += 1
            continue
        prefix_pairs[region + "->" + prefix] += 1
        key = permit_key(row, "ofwat")
        if key is None:
            excluded["missing_or_invalid_full_permit_key"] += 1
            continue
        counts["complete_identity_rows"] += 1
        pairs.add((point, key))
    by_point = collections.defaultdict(set)
    physical = collections.defaultdict(set)
    for point, key in pairs:
        by_point[point].add(key)
        physical[point].add((key[0], key[1], key[3], key[4]))
    qa = {k: int(counts[k]) for k in ("input_rows", "year_rows", "complete_identity_rows")}
    qa.update(distinct_pairs=len(pairs), distinct_points=len(by_point),
        duplicate_identity_rows=counts["complete_identity_rows"] - len(pairs),
        points_with_multiple_keys=sum(len(v) > 1 for v in by_point.values()),
        points_with_multiple_keys_ignoring_version=sum(len(v) > 1 for v in physical.values()),
        excluded=dict(sorted(excluded.items())),
        full_notation_region_correspondence=dict(sorted(prefix_pairs.items())))
    return pairs, qa


def _date(value):
    if isinstance(value, dt.datetime):
        return value.date()
    if isinstance(value, dt.date):
        return value
    return dt.date.fromisoformat(text(value))


def date_overlap_status(row, year=2018):
    """Retrospective date screen, never proof of complete historical activity."""
    if not text(row.get("EFFECTIVE_DATE")):
        return "missing_effective_date"
    try:
        effective = _date(row["EFFECTIVE_DATE"])
    except (ValueError, TypeError):
        return "invalid_effective_date"
    revocation = None
    if text(row.get("REVOCATION_DATE")):
        try:
            revocation = _date(row["REVOCATION_DATE"])
        except (ValueError, TypeError):
            return "invalid_revocation_date"
    if revocation is not None and revocation <= effective:
        return "reversed_or_empty_interval"
    if effective >= dt.date(year + 1, 1, 1):
        return "effective_after_year"
    if revocation is not None and revocation <= dt.date(year, 1, 1):
        return "revoked_before_year"
    return "overlaps_by_recorded_dates"


def summarize_candidates(candidates, pairs, corroborated, active):
    if not active <= corroborated <= pairs:
        raise ValueError("Corroboration must be a subset of source-reported identities")
    by_point = collections.defaultdict(set)
    for point, key in pairs:
        by_point[point].add(key)
    covered = set(by_point) & candidates
    confirmed = {p for p, _ in corroborated} & candidates
    current = {p for p, _ in active} & candidates
    unique = {p for p in covered if len(by_point[p]) == 1}
    return {
        "candidate_points": len(candidates),
        "with_historical_ofwat_pair": len(covered),
        "without_ofwat_pair": len(candidates - covered),
        "with_exactly_one_ofwat_key": len(unique),
        "with_multiple_ofwat_keys": len(covered - unique),
        "ofwat_distinct_pairs": sum(len(by_point[p]) for p in covered),
        "with_ea_all_corroborated_pair": len(confirmed),
        "ea_all_corroborated_distinct_pairs": sum(p in candidates for p, _ in corroborated),
        "with_ea_active_corroborated_pair": len(current),
        "ea_active_corroborated_distinct_pairs": sum(p in candidates for p, _ in active),
        "ea_all_corroborated_distinct_pairs_absent_from_active": sum(
            p in candidates for p, _ in corroborated - active),
        "with_corroboration_only_in_all_not_active": len(confirmed - current),
        "with_unique_ofwat_key_and_ea_corroboration": len(unique & confirmed),
        "source_reported_annual_association_not_sample_date_version_assignment": True,
        "does_not_verify_UWWTD_dcpCode_to_permit_outlet": True,
    }
