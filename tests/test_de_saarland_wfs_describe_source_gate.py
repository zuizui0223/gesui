import importlib.util
import unittest
from pathlib import Path

P=Path(__file__).resolve().parents[1]/"analysis/20261119_de_sludge_saarland/sl_saarland_official_wfs_describe_feature_SOURCE_gate.py"
spec=importlib.util.spec_from_file_location("saarland_describe",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

CAP=b'''<wfs:WFS_Capabilities xmlns:wfs="http://www.opengis.net/wfs/2.0"
  xmlns:ows="http://www.opengis.net/ows/1.1">
  <ows:OperationsMetadata><ows:Operation name="GetCapabilities"/>
     <ows:Operation name="DescribeFeatureType"/></ows:OperationsMetadata>
  <wfs:FeatureTypeList>
    <wfs:FeatureType><wfs:Name>Abfallbeseitigung:Klaerschlamm_Verfahren</wfs:Name>
    </wfs:FeatureType></wfs:FeatureTypeList>
  </wfs:WFS_Capabilities>'''
XSD=b'''<xsd:schema xmlns:xsd="http://www.w3.org/2001/XMLSchema">
  <xsd:complexType name="SewageSludgeType">
    <xsd:sequence>
      <xsd:element name="ausbringungsjahr" type="xsd:int"/>
      <xsd:element name="klaeranlage_id" type="xsd:string"/>
      <xsd:element name="masse" type="xsd:double"/>
    </xsd:sequence>
  </xsd:complexType></xsd:schema>'''

def response(blob):
    return blob,{"status":200,"content_type":"text/xml","size_bytes":len(blob),"sha256":"fixture"}

class WFSDescribeSchemaTest(unittest.TestCase):
    def test_official_capability_advertises_exact_layer_and_schema_operation(self):
        c=m.parse_capabilities(CAP)
        self.assertTrue(c["target_layer_exactly_listed"])
        self.assertTrue(c["DescribeFeatureType_operation_advertised"])

    def test_schema_field_names_do_not_elevate_to_confirmed_application(self):
        r=m.parse_describe_schema(XSD)
        self.assertIn("ausbringungsjahr",r["source_field_candidates"]["potential_time_fields"])
        self.assertIn("klaeranlage_id",r["source_field_candidates"]["potential_source_station_fields"])
        self.assertFalse(r["realized_application_year_verified"])
        self.assertFalse(r["source_year_2017to22_field_exposure_GO"])
        self.assertEqual(r["biological_model_gate"],"STOP")

    def test_official_only_operation_sequence(self):
        calls=[]
        def fake(url):
            calls.append(url)
            return response(CAP if len(calls)==1 else XSD)
        r=m.audit(fetch=fake)
        self.assertEqual(len(calls),2)
        self.assertEqual(calls[0],m.CAPABILITIES)
        self.assertIn("DescribeFeatureType",calls[1])
        self.assertIn("Abfallbeseitigung%3AKlaerschlamm_Verfahren",calls[1])
        self.assertEqual(r["source_gate"],"WFS_SCHEMA_FIELDS_CONFIRMED_VALUES_NOT_CHECKED")
        self.assertTrue(r["no_GetFeature_or_field_geometries_downloaded"])

    def test_missing_operation_does_not_fetch_unknown_endpoint(self):
        calls=[]
        def fake(url):
            calls.append(url)
            return response(CAP.replace(b'<ows:Operation name="DescribeFeatureType"/>',b''))
        r=m.audit(fetch=fake)
        self.assertEqual(len(calls),1)
        self.assertEqual(r["source_gate"],"DESCRIBE_OPERATION_NOT_ADVERTISED_STOP")

    def test_invalid_schema_preserves_stop_and_error(self):
        calls=[]
        def fake(url):
            calls.append(url)
            return response(CAP if len(calls)==1 else b'<html/>')
        r=m.audit(fetch=fake)
        self.assertIn("source_error",r)
        self.assertEqual(r["source_gate"],"ACCESS_STOP")
        self.assertEqual(r["biological_model_gate"],"STOP")

    def test_wrong_domain_not_accepted(self):
        with self.assertRaises(ValueError):
            m.get("https://example.com/arcgis/services/Internet/Abfallbeseitigung/MapServer/WFSServer?REQUEST=DescribeFeatureType")

if __name__=="__main__":unittest.main()
