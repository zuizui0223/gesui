"""Pure source-join invariants; no external queries or biological outcomes."""
from __future__ import annotations
import importlib.util
from pathlib import Path
import unittest
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon

P=Path(__file__).resolve().parents[1]/"analysis/20261118_point_diffuse_source/biotic_2x2_2018_geographic_source_gate.py"
spec=importlib.util.spec_from_file_location("gisgate",P)
g=importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)


class GeographicGateTest(unittest.TestCase):
    def test_h7_nuts2_area_threshold_is_strict(self):
        nuts=gpd.GeoDataFrame([
            {"geo":"AA|AA11","CNTR_CODE":"AA","geometry":Polygon([(0,0),(10,0),(10,10),(0,10)])},
            {"geo":"AA|AA12","CNTR_CODE":"AA","geometry":Polygon([(10,0),(20,0),(20,10),(10,10)])}
        ],geometry="geometry",crs="EPSG:3035")
        bas=gpd.GeoDataFrame([
            {"h":"h1","geometry":Polygon([(1,1),(9,1),(9,9),(1,9)])},
            {"h":"h2","geometry":Polygon([(5,1),(15,1),(15,9),(5,9)])}
        ],geometry="geometry",crs="EPSG:3035")
        land=pd.DataFrame({"h":["h1","h2"],"h6":["p1","p2"],"period":["2017-2019"]*2,"country":["AA"]*2})
        mapped, qc=g.map_land(land,bas,nuts)
        self.assertEqual(mapped.loc[mapped.h=="h1","assigned_NUTS2"].iloc[0],"AA|AA11")
        self.assertIsNone(mapped.loc[mapped.h=="h2","assigned_NUTS2"].iloc[0])
        self.assertEqual(qc["H7_strict_assigned_over90pct"],1)

    def test_nuts2_h6_geographic_components_protect_against_leakage(self):
        frame=pd.DataFrame([
            {"h":"A","h6":"P1","assigned_NUTS2":"AA|A"},
            {"h":"B","h6":"P1","assigned_NUTS2":"AA|B"},
            {"h":"C","h6":"P2","assigned_NUTS2":"AA|B"},
            {"h":"D","h6":"P3","assigned_NUTS2":"AA|C"}
        ])
        n, comp, load=g.components(frame)
        self.assertEqual(n,2)
        self.assertEqual(sum(comp),2)
        self.assertEqual(sum(load),4)

    def test_no_false_positive_on_tiny_shared_cohort(self):
        source=pd.DataFrame({
            "water_point_N_t":[2,3],
            "land_atmos_N_mg_m2":[10,11]
        },index=["AA|AA11","AA|AA12"])
        water=pd.DataFrame({"site_key":["site1","site2"],"geo":["AA|AA11","AA|AA12"]})
        land=pd.DataFrame({"h":["h1","h2"],"h6":["p1","p2"],"assigned_NUTS2":["AA|AA11","AA|AA12"]})
        result=g.gate(source,water,land)
        self.assertEqual(result["common_water_and_land_strict_NUTS2_source_upper_bound"],2)
        self.assertFalse(result["geographic_upper_bound_thresholds_pass"])
        self.assertFalse(result["final_biotic_source_gate_pass"])

if __name__=="__main__":
    unittest.main()
