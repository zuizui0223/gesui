import importlib.util
from pathlib import Path
import unittest

P=Path(__file__).resolve().parents[1]/"analysis/20261119_jrc_native_basins/jrc_native_point_soil_polygon30_gate.py"
spec=importlib.util.spec_from_file_location("geom30",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class NativeGeometryTests(unittest.TestCase):
    def test_esri_exterior_clockwise_rectangle(self):
        ring=[[0,0],[0,10],[10,10],[10,0],[0,0]]
        self.assertLess(m.signed_area(ring),0)
        shape=m.esri_polygon({"rings":[ring]})
        self.assertAlmostEqual(shape.area,100.0)
    def test_hole_and_intersection(self):
        outer=[[0,0],[0,10],[10,10],[10,0],[0,0]]
        hole=[[2,2],[8,2],[8,8],[2,8],[2,2]]
        p=m.esri_polygon({"rings":[outer,hole]})
        self.assertAlmostEqual(p.area,64.0)
    def test_no_exterior_must_fail(self):
        wrong=[[0,0],[10,0],[10,10],[0,10],[0,0]]
        with self.assertRaises(ValueError):
            m.esri_polygon({"rings":[wrong]})

if __name__=="__main__":
    unittest.main()
