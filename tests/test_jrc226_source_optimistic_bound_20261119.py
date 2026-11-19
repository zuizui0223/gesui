"""Independent optimistic upperbound tests; zero biological outputs."""
import importlib.util
from pathlib import Path
import unittest
P=Path(__file__).resolve().parents[1]/"analysis/20261119_JRC_native_basin_N_source/jrc226_source_optimistic_repair_upperbound.py"
spec=importlib.util.spec_from_file_location("upper",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class SourceBoundTests(unittest.TestCase):
    def records(self):
        out=[]
        for i in range(226):
            strict=i<120
            good=i<8
            out.append({
                "original_H7":str(i),"true_H6":str(i//3),
                "original_source_country":"FR",
                "soil_native_source_ID":str(i//10) if strict else "",
                "water_native_point_BasinID":str(i//10) if strict else "",
                "point_model_N_density_proxy_t_km2":str(1+i//10),
                "atmos_model_N_mean":str(1+i//10),
                "union_90_all4_core":"True",
                "single_model_source_strict_geom_eligible":str(strict),
                "conservative_geometry_quality_pass":str(good),
            })
        return out
    def test_optimistic_superset_fails_unique_basin_even_if_120_h7(self):
        data=m.summarize(self.records())
        self.assertEqual(data["optimistic_IF_ALL_ESRI_RING_REPAIRS_ACCEPTED"]["H7"],120)
        self.assertEqual(data["optimistic_IF_ALL_ESRI_RING_REPAIRS_ACCEPTED"]["unique_soil_native_HydroID"],12)
        self.assertTrue(data["independent_native_source_basins_FAIL_even_if_all_repairs_pass"])
    def test_strict_cannot_be_outside_optimistic(self):
        data=self.records()
        data[150]["conservative_geometry_quality_pass"]="True"
        with self.assertRaises(ValueError):
            m.summarize(data)
if __name__=="__main__":unittest.main()
