import importlib.util
from pathlib import Path
import unittest

SRC=Path(__file__).resolve().parents[1]/"analysis/20261119_jrc_native_basins/jrc_soilN_waterpoint_native_inventory.py"
spec=importlib.util.spec_from_file_location("jrcnative",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class NativeJRCSourceTest(unittest.TestCase):
    def test_four_soil_layers_share_hydroid_schema_not_point_key(self):
        soil=[r for k,r in m.CONFIG.items() if k.startswith("soil_")]
        self.assertEqual(len(soil),4)
        self.assertTrue(all(v[2]=="Basin_3035_HydroID" for v in soil))
        self.assertEqual(m.CONFIG["water_pointN_2014_18"][2],"BasinID")
    def test_native_layer_schema_rejects_missing_key(self):
        m0={"type":"Feature Layer","geometryType":"esriGeometryPolygon",
            "fields":[{"name":"OBJECTID"},{"name":"Basin_3035_HydroID"},
                      {"name":"Basin_3035_AreaSqKm"},{"name":"MEAN_mgsqmy"}]}
        fields=m.validate_meta(m0,"deposition","Basin_3035_HydroID","MEAN_mgsqmy","Basin_3035_AreaSqKm")
        self.assertIn("MEAN_mgsqmy",fields)
        with self.assertRaises(ValueError):
            m.validate_meta(m0,"point","BasinID","PsN","AreaSqKm")
    def test_polygon_only(self):
        m0={"type":"Feature Layer","geometryType":"esriGeometryPoint",
            "fields":[{"name":n} for n in ["OBJECTID","Basin_3035_HydroID","Basin_3035_AreaSqKm","MEAN_kgha"]]}
        with self.assertRaises(ValueError):
            m.validate_meta(m0,"mineral","Basin_3035_HydroID","MEAN_kgha","Basin_3035_AreaSqKm")

if __name__=="__main__":
    unittest.main()
