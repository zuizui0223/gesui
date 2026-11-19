import importlib.util
import json
import tempfile
from pathlib import Path
import unittest

P=Path(__file__).resolve().parents[1]/"analysis/20261119_point_land_N_source/eea2016_water_only_pointN_landN_source_gate.py"
spec=importlib.util.spec_from_file_location("watergate",P)
g=importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)

class FrozenBenthicSourceTest(unittest.TestCase):
    def fixture(self,folder):
        codes=["AT","BE","DK","FR","HR","HU","IT","LT","LU","LV","NL","NO","SE"]
        rows=[]
        folds={}
        for i in range(79):
            country=codes[i%13]
            nuts=country+str(i).zfill(3)
            rows.append({"country":country,"original_GISCO_NUTS2":nuts,
                         "source_original_num_2016_EEA_sites":31 if i<78 else 51,
                         "source_NUTS2_GREEN_pointN_t":float(i+1),
                         "source_NUTS2_area_km2":10.})
            folds[nuts]=i%5+1
        src={}
        for r in rows:
            if r["country"]=="NO":continue
            key=r["country"]+"|"+r["original_GISCO_NUTS2"]
            src[key]={"water_point_N_t":r["source_NUTS2_GREEN_pointN_t"],
                      "land_atmos_N_mg_m2":float(len(key)+int(r["original_GISCO_NUTS2"][2:]))}
        j={"original_2016_EEA_numeric_EQR_outcome_VALUES_not_read":True,
           "source_biological_classes":{"benthic_invertebrates":{
               "original_NUTS2_source_inventory":rows,
               "fixed_region_to_fivefold_source_only":folds,
               "n_original_independent_NUTS2":79,
               "n_sites_source_fixed_before_EQR":2469}}}
        f=Path(folder)/"source.json"
        f.write_text(json.dumps(j))
        return f,src
    def test_original_year_and_country_stop(self):
        with tempfile.TemporaryDirectory() as t:
            f,src=self.fixture(t)
            records,gate=g.gate(f,src)
            self.assertEqual(gate["source_complete_NUTS2_shared_n"],73)
            self.assertEqual(gate["source_shared_country_n"],12)
            self.assertTrue(gate["never_executed_biotic_numeric_model"])
    def test_eqr_access_not_permitted(self):
        with tempfile.TemporaryDirectory() as t:
            f,src=self.fixture(t)
            j=json.loads(f.read_text())
            j["original_2016_EEA_numeric_EQR_outcome_VALUES_not_read"]=False
            f.write_text(json.dumps(j))
            with self.assertRaises(ValueError):
                g.gate(f,src)

if __name__=="__main__":unittest.main()
