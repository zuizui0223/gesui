"""No-result tests for 226-H7 JRC native SOURCE-only eligibility rules."""
import importlib.util
from pathlib import Path
import unittest

P=Path(__file__).resolve().parents[1]/"analysis/20261119_JRC_native_basin_N_source/jrc_original_226h7_source_completeness_go_stop.py"
spec=importlib.util.spec_from_file_location("nativefull",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def make(h,coverage=100.,single=100.,same=True,repair=False):
    soil_id="11" if same else None
    rows=[]
    for kind in m.SOURCE_TYPES:
        id_="11" if kind!="pointN" else "22"
        if kind=="organicN" and not same:id_="33"
        rows.append({
            "h7":h,"source":kind,
            "status":"COMPLETE_OFFICIAL_JRC_SOURCE_ONLY",
            "true_h6":"p1","original_source_country":"FR",
            "numeric_geom_coverage_pct":coverage,
            "max_numeric_single_native_coverage_pct":single,
            "max_single_native_id":id_,
            "soil_model_basin_weighted_mean":12. if kind!="pointN" else None,
            "model_pointN_t_per_native_area_km2_proxy":4. if kind=="pointN" else None,
            "material_original_jrc_ring_repair":repair
        })
    return rows

class GateTests(unittest.TestCase):
    def test_union_90_is_not_unique_90(self):
        rows,qa=m.source_support(make("h1",single=60))
        self.assertEqual(qa["H7_with_all4_model_union_source_coverage_at_least90pct"],1)
        self.assertEqual(qa["H7_with_all4_single_native_polygon_90pct_and_same_soil_id"],0)
        self.assertEqual(qa["H7_conservative_geometry_qualified"],0)
        self.assertFalse(qa["new_biotic_model_permitted_by_this_source_gate"])

    def test_single_native_requires_consistent_soil_id(self):
        rows,qa=m.source_support(make("h1",same=False))
        self.assertEqual(qa["H7_with_all4_single_native_polygon_90pct_and_same_soil_id"],0)
        self.assertFalse(qa["new_biotic_model_permitted_by_this_source_gate"])

    def test_material_original_ring_repair_requires_independent_qa(self):
        rows,qa=m.source_support(make("h1",repair=True))
        self.assertEqual(qa["H7_with_all4_single_native_polygon_90pct_and_same_soil_id"],1)
        self.assertEqual(qa["H7_conservative_geometry_qualified"],0)
        self.assertEqual(qa["H7_material_native_Esri_ring_repair_unverified"],1)

    def test_shared_point_basin_and_H6_dont_count_independent_exposures(self):
        r=make("h1")+make("h2")
        _,qa=m.source_support(r)
        self.assertEqual(qa["H7_conservative_geometry_qualified"],2)
        self.assertEqual(qa["number_native_water_point_BasinID_independent"],1)
        self.assertEqual(qa["independent_H6_x_soil_native_x_point_native_connected_components"],1)

if __name__=="__main__":unittest.main()
