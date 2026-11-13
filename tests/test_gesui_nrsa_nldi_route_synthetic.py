"""Synthetic-only route classification: NLDI exact permit+outfall and same-reach.
Never tests real biological effects or uses official source observations."""
import importlib.util
import pathlib
import unittest
from unittest.mock import patch

SCRIPT=pathlib.Path(__file__).resolve().parents[1]/"scripts"/"gesui_nrsa_nldi_upstream_registered_POTW_source_shards.py"
spec=importlib.util.spec_from_file_location("gesui_source_nldi",SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class RouteContract(unittest.TestCase):
 def test_icis_outfall_three_digit_padding(self):
  self.assertEqual(mod.npdes_key("IA0033723","1"),"IA0033723001")
  self.assertEqual(mod.npdes_key("IA0033723","12"),"IA0033723012")
  self.assertEqual(mod.npdes_key("IA0033723","001"),"IA0033723001")
  self.assertEqual(mod.npdes_key("IA0033723","A01"),"IA0033723A01")
 def test_true_up_down_and_same_reach_separate(self):
  site={"uid":"S1","site_id":"NRS_SAMPLE","comid":"100","state":"AA","huc8":"01234567",
        "candidate_outfalls_25km":[
          {"permit":"AA0000001","feature":"001","geodesic_km":3.0,"outfall_lat":39.0,"outfall_lon":-90.0},
          {"permit":"AA0000002","feature":"002","geodesic_km":7.0,"outfall_lat":39.0,"outfall_lon":-90.0},
          {"permit":"AA0000003","feature":"003","geodesic_km":7.0,"outfall_lat":39.0,"outfall_lon":-90.0},
          {"permit":"AA0000004","feature":"004","geodesic_km":7.0,"outfall_lat":39.0,"outfall_lon":-90.0}]}
  def mocked(url):
   if "/UT/" in url:
    items=[{"id":"AA0000001001","comid":"99"},{"id":"AA0000003003","comid":"100"}]
   elif "/DM/" in url:
    items=[{"id":"AA0000002002","comid":"101"},{"id":"AA0000003003","comid":"100"}]
   else:raise AssertionError("Unexpected NLDI mode")
   return {"status":200,"records":items,"n_source_features":len(items),"source_sha256":"source_test"}
  with patch.object(mod,"call",side_effect=mocked):
   result=mod.network_one(site)
  self.assertEqual(result["status"],"NETWORK_SOURCE_API_OK")
  vals={x["feature"]:x["nldi_route_status"] for x in result["matched_NLDI_exact_registered_pairs"]}
  self.assertEqual(vals["001"],"UPSTREAM_DIFFERENT_REACH")
  self.assertEqual(vals["002"],"DOWNSTREAM_DIFFERENT_REACH")
  self.assertEqual(vals["003"],"AMBIGUOUS_IN_BOTH_UT_AND_DM")
  self.assertNotIn("004",vals)
  self.assertEqual(result["n_verified_different_reach_upstream_candidate"],1)
 def test_http_failure_is_not_no_outfall(self):
  site={"uid":"X","site_id":"X","comid":"100","state":"AA","huc8":"01234567",
        "candidate_outfalls_25km":[{"permit":"AA0000001","feature":"001",
        "geodesic_km":1.0,"outfall_lat":39.0,"outfall_lon":-90.0}]}
  with patch.object(mod,"call",return_value={"status":503,"error":"failed","records":None}):
   out=mod.network_one(site)
  self.assertEqual(out["status"],"NETWORK_API_UNAVAILABLE")
  self.assertNotIn("n_verified_different_reach_upstream_candidate",out)

if __name__=="__main__":unittest.main()
