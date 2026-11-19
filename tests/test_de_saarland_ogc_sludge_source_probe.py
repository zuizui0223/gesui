import importlib.util
import unittest
from pathlib import Path

P=Path(__file__).resolve().parents[1]/"analysis/20261119_de_sludge_saarland/sl_saarland_official_ogc_SOURCE_probe.py"
spec=importlib.util.spec_from_file_location("saarland",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class SaarlandOfficialSourceMetadataTest(unittest.TestCase):
    def test_official_wfs_layer_identification_only(self):
        xml='''<?xml version="1.0" encoding="UTF-8"?>
        <wfs:WFS_Capabilities xmlns:wfs="http://www.opengis.net/wfs/2.0">
         <wfs:FeatureTypeList>
           <wfs:FeatureType><wfs:Name>ABC:Forest</wfs:Name><wfs:Title>Forest</wfs:Title></wfs:FeatureType>
           <wfs:FeatureType><wfs:Name>Abfallbeseitigung:Klaerschlamm_Verfahren</wfs:Name>
            <wfs:Title>Ausbringungsflächen von Klärschlamm</wfs:Title></wfs:FeatureType>
         </wfs:FeatureTypeList>
        </wfs:WFS_Capabilities>'''.encode('utf-8')
        r=m.parse_wfs_capabilities(xml)
        self.assertEqual(r["WFS_feature_type_count"],2)
        self.assertTrue(r["matching_layer_name_confirmed"])
        self.assertFalse(r["year_field_verified"])
        self.assertFalse(r["actual_vs_plan_status_verified"])
        self.assertFalse(r["source_WWTP_ID_verified"])
        self.assertEqual(r["water_and_land_H6_biotic_model_gate"],"STOP")

    def test_empty_layer_does_not_become_positive(self):
        r=m.parse_wfs_capabilities(b'<WFS_Capabilities><FeatureTypeList><FeatureType><Name>Waste</Name></FeatureType></FeatureTypeList></WFS_Capabilities>')
        self.assertFalse(r["matching_layer_name_confirmed"])

    def test_invalid_document_rejected(self):
        with self.assertRaises(ValueError):
            m.parse_wfs_capabilities(b'<html><body>Not a WFS document</body></html>')

    def test_untrusted_domain_rejected_before_network(self):
        with self.assertRaises(ValueError):
            m.probe("https://example.org/data")
        with self.assertRaises(ValueError):
            m.probe("http://geoportal.saarland.de/")

if __name__=="__main__":
    unittest.main()
