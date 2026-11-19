import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

P = Path(__file__).resolve().parents[1] / 'analysis/20261119_fr_sludge_routeC/fr_external_performed_spreading_coverage.py'
spec = importlib.util.spec_from_file_location('performed_coverage', P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class ExternalPositiveCoverageTest(unittest.TestCase):
    def write(self, name, rows, fields):
        path = Path(self.tmp.name) / name
        with path.open('w', newline='', encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        return path

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.candidate_fields = ['INSEE_code', 'candidate_source_years', 'source_interpretation']
        self.ext_fields = ['INSEE_code','commune_name','application_year','source_station_SANDRE','evidence_class','primary_report_url','source_locator']
        self.candidate = {'INSEE_code':'42192','candidate_source_years':'2018;2020','source_interpretation':m.SOURCE_TAG}
        self.external = {'INSEE_code':'42066','commune_name':'Cleppé','application_year':'2018','source_station_SANDRE':'0442094S0002', 'evidence_class':m.EVIDENCE_TAG,'primary_report_url':'https://government.example/document.pdf','source_locator':'Table 7'}

    def tearDown(self):
        self.tmp.cleanup()

    def test_confirmed_positive_missing_never_becomes_negative_control(self):
        a = self.write('a.csv', [self.candidate], self.candidate_fields)
        b = self.write('b.csv', [self.external], self.ext_fields)
        result = m.audit(a,b)
        self.assertEqual(result['external_cases_not_in_2018_ROSEAU_candidate_list'], 1)
        self.assertEqual(result['external_cases_not_in_2017to22_union'], 1)
        self.assertTrue(result['absence_from_list_is_not_absence_of_field_application'])
        self.assertEqual(result['biological_model_source_gate'], 'STOP')

    def test_same_year_verified_candidate_is_found(self):
        a = self.write('a.csv', [dict(self.candidate, INSEE_code='42066')], self.candidate_fields)
        b = self.write('b.csv', [self.external], self.ext_fields)
        result = m.audit(a,b)
        self.assertEqual(result['external_cases_not_in_2018_ROSEAU_candidate_list'], 0)
        self.assertTrue(result['cases'][0]['in_ROSEAU_positive_station_commune_candidate_same_year'])

    def test_other_year_cannot_be_claimed_same_year(self):
        a = self.write('a.csv', [dict(self.candidate,INSEE_code='42066',candidate_source_years='2020')], self.candidate_fields)
        b = self.write('b.csv', [self.external], self.ext_fields)
        result = m.audit(a,b)
        self.assertTrue(result['cases'][0]['in_ROSEAU_candidate_union_2017to22'])
        self.assertFalse(result['cases'][0]['in_ROSEAU_positive_station_commune_candidate_same_year'])

    def test_planned_spreading_not_accepted(self):
        a = self.write('a.csv', [self.candidate], self.candidate_fields)
        b = self.write('b.csv', [dict(self.external, evidence_class='PLANNED')], self.ext_fields)
        with self.assertRaises(ValueError):
            m.audit(a,b)

    def test_duplicates_fail_closed(self):
        a = self.write('a.csv', [self.candidate,self.candidate], self.candidate_fields)
        b = self.write('b.csv', [self.external], self.ext_fields)
        with self.assertRaises(ValueError):
            m.audit(a,b)


if __name__ == '__main__':
    unittest.main()
