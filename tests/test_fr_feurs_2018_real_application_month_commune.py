import csv
import importlib.util
import unittest
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "analysis/20261119_fr_sludge_routeC/fr_feurs_2018_actual_application_month_commune_audit.py"
SOURCE = ROOT / "analysis/20261119_fr_sludge_routeC/FR_FEURS_2018_ANONYMIZED_ACTUAL_SPREADING_43_ROWS.csv"
spec = importlib.util.spec_from_file_location("actual2018", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def rows():
    with SOURCE.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


class FeursRealizedApplicationTest(unittest.TestCase):
    def test_independently_reconcile_public_table6_and_table7(self):
        result = m.audit_csv(SOURCE)
        self.assertEqual(result["source_table7_actual_application_rows"], 43)
        self.assertEqual(result["municipalities_with_actual_application"], 8)
        self.assertEqual(result["original_annual_applied_area_ha"], 161.53)
        self.assertEqual(result["original_annual_applied_wet_t"], 1918.0)
        self.assertEqual({k:int(v["actual_wet_sludge_t"]) for k,v in result["month"].items()},
                         {"3":761,"4":50,"7":613,"9":214,"10":208,"11":72})
        self.assertEqual(result["march_july_wet_t"], 1374.0)
        self.assertAlmostEqual(result["march_july_share_of_2018_wet_t"], 1374 / 1918, places=5)
        self.assertEqual(result["independent_realized_stations_years_in_current_case"], 1)
        self.assertEqual(result["biological_analysis_gate"], "STOP")
        self.assertTrue(result["no_NUTRI_HNS_EQR_values_read"])

    def test_commune_reconciliation_fixed(self):
        result = m.audit_csv(SOURCE)
        expect = {
            "42066":(14,50.11,623),"42174":(8,18.83,218),
            "42251":(5,19.69,231),"42155":(3,20.67,244),
            "42319":(4,20.09,234),"42094":(4,17.13,192),
            "42143":(3,6.64,79),"42135":(2,8.37,97),
        }
        for code, (n, ha, tonnes) in expect.items():
            c=result["commune"][code]
            self.assertEqual(c["application_log_rows"], n)
            self.assertAlmostEqual(c["actual_applied_area_ha"], ha)
            self.assertEqual(c["actual_wet_sludge_t"], tonnes)

    def test_month_relabel_even_with_annual_total_same_fails(self):
        data=rows()
        data[0]["month"]="4"
        with self.assertRaises(ValueError):
            m.summary(data)

    def test_2018_only_other_year_fails(self):
        data=rows()
        data[0]["year"]="2021"
        with self.assertRaises(ValueError):
            m.summary(data)

    def test_actual_annual_mass_tamper_fails(self):
        data=rows()
        data[0]["wet_sludge_t"]="37"
        with self.assertRaises(ValueError):
            m.summary(data)

    def test_actual_annual_area_tamper_fails(self):
        data=rows()
        data[0]["applied_area_ha"]="3.01"
        with self.assertRaises(ValueError):
            m.summary(data)

    def test_private_farmer_field_fails_closed(self):
        data=rows()
        data[0]["farmer_name"]="Do not export"
        with self.assertRaises(ValueError):
            m.summary(data)

    def test_identical_valid_application_rows_not_deduplicated(self):
        # Two 2018 Nervieux events really share the same month, ha and wet tonnes.
        data=rows()
        same=[x for x in data if x["INSEE_code"]=="42155" and x["month"]=="9"]
        self.assertEqual(len(same),2)
        self.assertEqual(same[0],same[1])
        self.assertEqual(m.summary(data)["source_table7_actual_application_rows"],43)

    def test_dropping_one_row_fails(self):
        with self.assertRaises(ValueError):
            m.summary(rows()[:-1])


if __name__ == "__main__":
    unittest.main()
