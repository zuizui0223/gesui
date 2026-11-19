"""Geo-only test: 2018 ROSEAU candidate communes are NOT sludge field exposure."""
from __future__ import annotations
import importlib.util
from pathlib import Path
import unittest

import geopandas as gpd
from shapely.geometry import Polygon
P=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_sludge_routeC/fr_roSEAU_2018_communes_to_H7_opportunity_geo.py"
spec=importlib.util.spec_from_file_location("fr2018geo",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class ScreeningOnlyTest(unittest.TestCase):
    def test_intersection_is_geographic_not_true_application(self):
        # EPSG3035 metres: one municipality 1000m^2 by 1000m; two H7 basins
        communes=gpd.GeoDataFrame([{"insee_code":"12345","geometry":Polygon([(0,0),(1000,0),(1000,1000),(0,1000)])}],
                                  crs="EPSG:3035")
        h7=gpd.GeoDataFrame([
            {"h":"hA","h6":"h6A","geometry":Polygon([(0,0),(500,0),(500,1000),(0,1000)])},
            {"h":"hB","h6":"h6B","geometry":Polygon([(500,0),(1500,0),(1500,1000),(500,1000)])}
        ],crs="EPSG:3035")
        pairs,counts=m.overlap(communes,h7)
        self.assertEqual(len(pairs),2)
        self.assertEqual(counts["positive_area_intersection"]["unique_H7_with_candidate_commune_overlap"],2)
        self.assertTrue(all(x["minimum_10ha_opportunity_screen"] for x in pairs))
        self.assertTrue(all("NOT confirmed field application" in x["source_interpretation"] for x in pairs))
        self.assertAlmostEqual(sum(x["intersection_ha"] for x in pairs),100.)
    def test_no_contact_only_false_positive(self):
        communes=gpd.GeoDataFrame([{"insee_code":"12345","geometry":Polygon([(0,0),(1000,0),(1000,1000),(0,1000)])}],crs="EPSG:3035")
        h7=gpd.GeoDataFrame([{"h":"hC","h6":"h6C","geometry":Polygon([(1000,0),(1500,0),(1500,1000),(1000,1000)])}],crs="EPSG:3035")
        pairs,counts=m.overlap(communes,h7)
        self.assertEqual(len(pairs),0)
        self.assertEqual(counts["positive_area_intersection"]["unique_H7_with_candidate_commune_overlap"],0)
    def test_predefined_area_floor_10ha(self):
        self.assertEqual(m.MIN_MEANINGFUL_HA,10.)
        self.assertEqual(m.ORIGINAL_2018_CANDIDATE_CODES,72)
        self.assertEqual(m.ORIGINAL_FIRST_PERIOD_H7,226)

if __name__=="__main__":unittest.main()
