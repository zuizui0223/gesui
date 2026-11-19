"""Offline source-only tests: native Esri rings and single-versus-union H7 coverage."""
from pathlib import Path
from unittest import TestCase,main,mock
from importlib.util import spec_from_file_location,module_from_spec
from shapely.geometry import Polygon
from shapely.ops import transform

FILE=Path(__file__).resolve().parents[1]/"analysis/20261119_JRC_native_basin_N_source/jrc_original_226h7_sharded_source_census.py"
spec=spec_from_file_location("source_full",FILE)
m=module_from_spec(spec)
spec.loader.exec_module(m)

def rect(x0,x1):
    # Esri exterior ring convention: clockwise.
    return [[x0,0],[x0,1],[x1,1],[x1,0],[x0,0]]

class ExactOriginalSourceGates(TestCase):
    def test_ring_original_codec(self):
        p,changed=m.esri_shape([rect(0,1)])
        self.assertFalse(p.is_empty)
        self.assertTrue(p.is_valid)
        self.assertEqual(changed,0)
        self.assertAlmostEqual(p.area,1.)

    def test_union_cover_is_not_single_native_basin(self):
        def fake(endpoint,params):
            if params.get("returnCountOnly")=="true":
                return {"count":2},"count_sha"
            data={"features":[
                {"attributes":{"OBJECTID":1,"Basin_3035_HydroID":11,"MEAN_mgsqmy":3000},
                 "geometry":{"rings":[rect(0,.6)]}},
                {"attributes":{"OBJECTID":2,"Basin_3035_HydroID":12,"MEAN_mgsqmy":2500},
                 "geometry":{"rings":[rect(.6,1)]}}
            ]}
            return data,"page_sha"
        geom=Polygon([(0,0),(1,0),(1,1),(0,1)])
        with mock.patch.object(m,"req",side_effect=fake):
            out=m.h7_layer("h",geom,transform(m.P,geom),"atmos2015",m.BASE,5,"MEAN_mgsqmy","Basin_3035_HydroID")
        self.assertEqual(out["status"],"COMPLETE_OFFICIAL_JRC_SOURCE_ONLY")
        self.assertGreater(out["numeric_geom_coverage_pct"],99.9)
        self.assertLess(out["max_numeric_single_native_coverage_pct"],90)
        self.assertEqual(out["number_geom_features_intersecting_H7"],2)
        self.assertAlmostEqual(out["soil_model_basin_weighted_mean"],2800,delta=10)

    def test_zero_source_not_imputed_missing(self):
        self.assertEqual(m.numeric(0),0)
        self.assertIsNone(m.numeric(None))
        self.assertIsNone(m.numeric(float("nan")))
        self.assertIsNone(m.numeric(-1))

    def test_original_native_point_total_is_converted_to_density(self):
        geom=Polygon([(0,0),(1,0),(1,1),(0,1)])
        def fake(endpoint,params):
            if params.get("returnCountOnly")=="true":return {"count":1},"countsha"
            return {"features":[
                {"attributes":{"OBJECTID":1,"BasinID":5,"PsN":50,"AreaSqKm":10},
                 "geometry":{"rings":[rect(0,1)]}}
            ]},"sha"
        with mock.patch.object(m,"req",side_effect=fake):
            out=m.h7_layer("h",geom,transform(m.P,geom),"pointN",m.WATER,5,"PsN","BasinID")
        self.assertAlmostEqual(out["model_pointN_t_per_native_area_km2_proxy"],5)
        self.assertGreater(out["max_numeric_single_native_coverage_pct"],99)

if __name__=="__main__": main()
