#!/usr/bin/env python3
"""Verify byte-exact synced result artifacts and common hydrological grain."""
import csv, hashlib, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/"results/20261115_sync/SOURCE_BYTE_EXACT_63_FILES.json"
def source_sha(b):
    return hashlib.sha1(f"blob {len(b)}".encode()+bytes([0])+b).hexdigest()
def read_rows(name):
    with (ROOT/name).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))
def check():
    manifest=json.loads(MANIFEST.read_text(encoding="utf-8"))
    rows=manifest["entries"]
    assert len(rows)==63,len(rows)
    keys=set()
    bytes_verified=0
    for r in rows:
        name=r["path"]
        assert name not in keys,name
        keys.add(name)
        assert name.startswith(("docs/RESULTS_","analysis/","results/20261115_sync/"))
        b=(ROOT/name).read_bytes()
        assert len(b)==r["bytes"],name
        assert hashlib.sha256(b).hexdigest()==r["sha256"],"SHA256 "+name
        assert source_sha(b)==r["git_blob_sha1"],"Git blob "+name
        bytes_verified+=len(b)
    assert bytes_verified==4820502,bytes_verified
    panel=read_rows("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")
    assert len(panel)==454
    rowkeys={(r["h"],r["period"]) for r in panel}
    assert len(rowkeys)==454
    assert len({r["h"] for r in panel})==232
    assert len({r["h6"] for r in panel})==90
    assert {r["period"] for r in panel}=={"2017-2019","2020-2022"}
    pair=read_rows("results/20261115_sync/hns_macro/paired_222_H7_HNS_changes.csv")
    assert len(pair)==222 and len({r["h"] for r in pair})==222
    assert all(r["period"]=="2017-2019_to_2020-2022" for r in pair)
    for x in ("network/H7_candidate_host_network_panel.csv",
              "robustness/H7_454_geographic_null_network_panel.csv"):
        other=read_rows("results/20261115_sync/"+x)
        assert len(other)==454
        assert {(r["h"],r["period"]) for r in other}==rowkeys
    assert len(read_rows("results/20261115_sync/GESUI_303_HOST_NUTRIENT_SCORES.csv"))==303
    assert len(read_rows("results/20261115_sync/CHOCHO_230_HOST_NUTRIENT_SCORES.csv"))==230
    audit=json.loads((ROOT/"results/20261115_sync/VERIFICATION_20261115.json").read_text())
    assert audit["sha_mismatches"]==0 and audit["unit_tests"]["total_pass"]==42
    print(json.dumps({"status":"PASS","source_sha256_and_git_blob_verified":len(rows),"bytes":bytes_verified,
    "source_unit_h7_period":454,"H7":232,"H6":90,"repeated_H7":222,
    "gesui_species_keys":303,"chocho_catalogue_species_names":230,
    "does_not_refit_field_effects":True},sort_keys=True))
if __name__=="__main__":check()
