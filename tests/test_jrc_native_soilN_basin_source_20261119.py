"""Offline invariants for JRC original soil-versus-water native basin ID audit."""
from pathlib import Path
import importlib.util
import unittest
SRC=Path(__file__).resolve().parents[1]/"analysis/20261119_JRC_native_basin_N_source/jrc_native_soilN_layers_source_gate.py"
spec=importlib.util.spec_from_file_location("jrcNative",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
class NativeBasinSourceTests(unittest.TestCase):
    def test_missing_not_equal_true_zero(self):
        self.assertIsNone(m.nonnegative(None))
        self.assertIsNone(m.nonnegative(float("nan")))
        self.assertIsNone(m.nonnegative(-1))
        self.assertEqual(m.nonnegative(0),0)
    def test_objects_must_have_integer_oid(self):
        with self.assertRaises(ValueError):m.validate_id("12")
        with self.assertRaises(ValueError):m.validate_id(True)
        self.assertEqual(m.validate_id(12),12)
    def test_full_soil_source_gate(self):
        names=list(m.SOURCES)
        rows={k:{str(i):float(i) for i in range(1,102)} for k in names}
        x=m.source_gate(rows)
        self.assertTrue(x["soil_source_first_native_basin_pass"])
        self.assertTrue(x["soil_source_hydroID_not_equated_to_point_source_BasinID"])
    def test_missing_original_not_imputed(self):
        names=list(m.SOURCES)
        rows={k:{str(i):float(i) for i in range(1,102)} for k in names}
        for k in names[:3]:
            rows[k]["1"]=None
        self.assertEqual(m.source_gate(rows)["native_soil_basins_3layers_complete_nonmissing"],100)
        for k in names[:3]:
            rows[k]["2"]=None
        self.assertFalse(m.source_gate(rows)["soil_source_first_native_basin_pass"])
if __name__=="__main__":unittest.main()
