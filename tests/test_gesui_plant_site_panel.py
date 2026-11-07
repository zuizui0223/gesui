"""Synthetic quality-control tests, not new ecological measurements."""
import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "gesui_validate", ROOT / "scripts" / "validate_plant_site_panel.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def valid():
    row = {key: "" for key in mod.REQUIRED}
    row.update(site_id="site1", period="2011-2024", taxon_id="taxon1",
        catchment_id="basin1", habitat_type="riparian", origin_status="native",
        occurrence_records="3", site_total_records="100",
        downstream_status="downstream", outfall_id="wwtp4",
        network_distance_km="0.5", ww_volume_m3_yr="1000",
        effluent_N_kg_yr="20", effluent_N_source="reported_facility_mass_load",
        trait_N="6", trait_source="Ellenberg")
    return row

def problems(rows):
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "panel.csv"
        with path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=mod.REQUIRED)
            w.writeheader()
            w.writerows(rows)
        return mod.check(path)

class PanelContractTests(unittest.TestCase):
    def test_valid(self):
        self.assertEqual(problems([valid()]), [])

    def test_duplicate_key(self):
        self.assertTrue(any("duplicate" in x for x in problems([valid(), valid()])))

    def test_ww_volume_is_not_n(self):
        row = valid()
        row["effluent_N_source"] = "derived_from_ww_volume_only"
        self.assertTrue(any("N mass" in x for x in problems([row])))

    def test_downstream_requires_connectivity(self):
        row = valid()
        row["outfall_id"] = ""
        self.assertTrue(any("downstream" in x for x in problems([row])))

    def test_no_true_absence_from_zero_records(self):
        row = valid()
        row["site_total_records"] = "0"
        self.assertTrue(any("denominator" in x for x in problems([row])))

if __name__ == "__main__":
    unittest.main()
