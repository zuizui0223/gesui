"""Offline response-blind tests for fixed 16-outfall source geometry."""
import importlib.util
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("directional",ROOT/"scripts/gesui_directed_river_plantae.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def test_real_geometry_selected_before_gbif():
    root=ROOT/"results/hydrorivers_eu_source_gate_20261107"
    selected=m.strict_source_manifest(json.loads((root/"candidate_site_design.json").read_text()),
                                      json.loads((root/"source_qc.json").read_text()))
    sites=m.site_registry(selected)
    assert len(selected)==16
    assert len(sites)==48
    assert len({s["site_id"] for s in sites})==48
    assert all(s["observed_N_mass_kg_year"] is None for s in sites)
    assert {s["river_role"] for s in sites}=={"u2","u1","d1"}

def test_exact_rarefaction():
    assert m.expected_rarefaction([30],30)==1
    assert m.expected_rarefaction([15,15],30)==2
    assert m.expected_rarefaction([2,3],30) is None
    assert m.expected_rarefaction([],10) is None

def test_jaccard_species_sets():
    assert abs(m.jaccard({1:4,2:1},{2:2,3:1})-2/3)<1e-12
    assert m.jaccard({1:1},{1:12})==0
    assert m.jaccard({}, {}) is None

def test_original_source_has_unequal_control_distances():
    root=ROOT/"results/hydrorivers_eu_source_gate_20261107"
    records=json.loads((root/"candidate_site_design.json").read_text())["selected"]
    assert all(x["river_hydraulic_chain"]=="U2>U1>OUTFALL_REACH>D1" for x in records)
    assert any(x["distance_u1_d1_km"]>2*x["distance_u2_u1_km"] for x in records)
