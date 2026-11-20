"""Audit source keys, coverage and design; never estimate ecological effects.

Downloads are kept outside the repository. Outputs contain provenance/schema
and source-level candidate records, not redistributable raw ecological rows.
"""
import argparse
import collections
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import urllib.request
import zipfile
import xml.etree.ElementTree as ET

PIN = "d6520b530f9310391596cd0095ecc0ad8ca07f0e"
BASE = f"https://raw.githubusercontent.com/juangallegozamorano/NImpact_Invertebrates_MetaAnalysis/{PIN}/"
FILES = ["README.md", "Databases/AbundanceDataset_csv.csv",
         "Databases/RichnessDataset_csv.csv", "Databases/References_AllPapersAdded.xlsx"]
SEARCH_FIELDS = ["Source", "Site_Name", "Fertilizer_Type", "Vegetation_Type", "Land_Use_Described", "Notes"]
PATTERNS = {
    "wastewater_terms": r"waste.?water|sewage|effluent|biosolid|sludge",
    "water_margin_terms": r"riparian|river|stream|floodplain|wetland|marsh|swamp",
}


def fetch(url, cache, name):
    request = urllib.request.Request(url, headers={"User-Agent": "gesui-source-opportunity-audit/1"})
    with urllib.request.urlopen(request, timeout=60) as response:
        content = response.read()
    (cache / name).write_bytes(content)
    return content


def decode_csv(content):
    try:
        text, encoding = content.decode("utf-8-sig"), "utf-8-sig"
    except UnicodeDecodeError:
        text, encoding = content.decode("cp1252"), "cp1252"
    return list(csv.DictReader(io.StringIO(text))), encoding


def xlsx_rows(content):
    ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(io.BytesIO(content)) as archive:
        shared = ["".join(node.itertext()) for node in ET.fromstring(
            archive.read("xl/sharedStrings.xml")).findall("s:si", ns)]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    records = []
    for row in sheet.findall(".//s:row", ns):
        record = {}
        for cell in row.findall("s:c", ns):
            column = re.sub(r"\d", "", cell.attrib["r"])
            value = cell.find("s:v", ns)
            if value is not None:
                record[column] = shared[int(value.text)] if cell.get("t") == "s" else value.text
            elif cell.get("t") == "inlineStr":
                record[column] = "".join(cell.find("s:is", ns).itertext())
        if record:
            records.append(record)
    header = records[0]
    return [{name: row.get(col, "") for col, name in header.items()} for row in records[1:]]


def csv_audit(content):
    rows, encoding = decode_csv(content)
    records = [row for row in rows if row.get("Source", "").strip()]
    result = {
        "encoding": encoding, "raw_rows": len(rows), "source_populated_rows": len(records),
        "rows_without_source": len(rows) - len(records), "columns": list(rows[0]),
        "unique_source_labels": len({row["Source"] for row in records}),
        "source_study_keys_NOT_independent_sources": len({(r["Source"], r["Study"]) for r in records}),
        "exact_herbivore_label_rows": sum(r["Feeding_Guild"] == "Herbivore" for r in records),
        "outcomes_or_effect_sizes_recomputed": False, "candidate_search": {},
    }
    for label, pattern in PATTERNS.items():
        matches = [r for r in records if re.search(pattern, " ".join(r.get(c, "") for c in SEARCH_FIELDS), re.I)]
        groups = collections.defaultdict(list)
        for row in matches:
            groups[(row["Source"], row["Study"])].append(row)
        result["candidate_search"][label] = {
            "matched_rows": len(matches), "unique_source_labels": sorted({r["Source"] for r in matches}),
            "source_study_candidates": [{
                "source": key[0], "study": key[1], "rows": len(values),
                "fertilizer": sorted({r["Fertilizer_Type"] for r in values}),
                "vegetation": sorted({r["Vegetation_Type"] for r in values}),
                "guilds": sorted({r["Feeding_Guild"] for r in values}),
            } for key, values in sorted(groups.items())],
        }
    return result


def run(cache, out):
    cache.mkdir(parents=True, exist_ok=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    provenance, audits = [], {}
    for path in FILES:
        url = BASE + path
        content = fetch(url, cache, "nimpact_" + Path(path).name)
        provenance.append({"url": url, "bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()})
        if path.endswith(".csv"):
            audits[path] = csv_audit(content)
        elif path.endswith(".xlsx"):
            refs = xlsx_rows(content)
            audits["references"] = {"rows": len(refs), "candidate_search": {}}
            for label, pattern in PATTERNS.items():
                selected = [r for r in refs if re.search(pattern, " ".join([r.get("Title", ""), r.get("Abstract", "")]), re.I)]
                audits["references"]["candidate_search"][label] = [
                    {k: r.get(k, "") for k in ("Source", "Year", "Title", "DOI")} for r in selected]
    meta_url = "https://zenodo.org/api/records/14510967"
    metadata = json.loads(fetch(meta_url, cache, "hammond_metadata.json"))
    hammond = {"record": meta_url, "files": [], "verified_schemas": {}}
    for item in metadata["files"]:
        name = item["key"]
        entry = {"name": name, "bytes": item["size"], "checksum": item["checksum"]}
        if name in ("traits.csv", "data column descriptions.csv"):
            content = fetch(item["links"]["self"], cache, "hammond_" + name)
            actual = "md5:" + hashlib.md5(content).hexdigest()
            if actual != item["checksum"]:
                raise ValueError("Hammond source checksum mismatch: " + name)
            entry.update({"checksum_verified": True, "sha256": hashlib.sha256(content).hexdigest()})
            rows, encoding = decode_csv(content)
            info = {"rows": len(rows), "columns": list(rows[0]), "encoding": encoding}
            if name == "traits.csv":
                if not {"site", "tree", "N", "P"}.issubset(rows[0]):
                    raise ValueError("Unexpected Hammond schema")
                info["site_tree_counts"] = dict(collections.Counter(r["site"] for r in rows))
                info["tree_ids"] = len({r["tree"] for r in rows})
                info["N_P_complete_rows"] = sum(bool(r["N"]) and bool(r["P"]) and r["N"] != "NA" and r["P"] != "NA" for r in rows)
                info["independent_wwtp_sources_from_original_design"] = 1
                info["aquatic_nutrients_row_table_present_in_this_deposit"] = False
                info["herbivore_table_present_in_this_deposit"] = False
                info["full_chain_identified"] = False
            hammond["verified_schemas"][name] = info
        hammond["files"].append(entry)
    result = {
        "purpose": "Source-opportunity audit, not ecological effect or meta-analysis",
        "nimpact_git_commit": PIN, "provenance": provenance, "nimpact": audits, "hammond": hammond,
        "limitations": [
            "Keyword candidates require original-paper verification; no match does not prove absence.",
            "Source/Study labels can share sites, dates and controls; these are not independent replicate counts.",
            "Raw extraction sheets differ from the 126-publication/4365-observation published analysis cohort.",
            "N addition, sludge, animal effluent and municipal river effluent cannot be pooled as one exposure.",
        ],
    }
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"audit_file": str(out), "files_pinned": len(provenance),
                      "hammond_schemas": hammond["verified_schemas"],
                      "new_ecological_effects": 0}, ensure_ascii=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.cache, args.out)
