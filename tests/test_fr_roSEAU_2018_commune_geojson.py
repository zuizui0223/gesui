import importlib.util,unittest
from pathlib import Path

P=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_sludge_routeC/fr_roSEAU_2018_commune_geojson_source_audit.py"
spec=importlib.util.spec_from_file_location("geo_roSEAU",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
class FranceSludgeGeoJSONTests(unittest.TestCase):
    def test_geojson_points_not_farms(self):
        blob=b'{"type":"FeatureCollection","features":[{"type":"Feature","geometry":{"type":"Point","coordinates":[2,48]},"properties":{"Code SANDRE la station":"123","Commune d\\u0027epandage":"Example"}},{"type":"Feature","geometry":{"type":"Point","coordinates":[3,47]},"properties":{}}]}'
        q=m.geojson_schema(blob)
        self.assertEqual(q["feature_count"],2)
        self.assertEqual(q["geometry_types"],{"Point":2})
        self.assertTrue(q["GeoJSON_geometry_does_NOT_auto_prove_parcel_or_commune_location"])
    def test_missing_geojson_is_stop(self):
        with self.assertRaises(ValueError):m.geojson_schema(b'{"type":"FeatureCollection","features":[]}')
    def test_original_year_unmodified(self):
        self.assertEqual(m.REQ["annee"],"2018")
        self.assertEqual(m.REQ["tranche"],"ALL")
if __name__=="__main__": unittest.main()
