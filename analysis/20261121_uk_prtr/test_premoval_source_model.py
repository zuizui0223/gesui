"""Offline regressions for the frozen UK2018 source model and record audit."""
import argparse
import contextlib
import io
import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

import eea2018_PRemoval_to_EA_TP_effluent_source_test as model


class SourceModelTests(unittest.TestCase):
    def frame(self):
        rng = np.random.default_rng(177)
        n = 180
        flag = np.tile([0, 1], n // 2)
        east = 200000 + 150000 * flag + rng.uniform(0, 100000, n)
        north = rng.uniform(100000, 600000, n)
        capacity = rng.uniform(1000, 100000, n)
        nr = rng.integers(0, 2, n)
        region = np.tile(['UKC', 'UKD', 'UKE'], n // 3)
        y = (1.5 + .08 * flag + .13 * east / 100000 - .09 * north / 100000
             + .02 * np.log1p(capacity) + .04 * nr
             + .1 * (region == 'UKD') - .05 * (region == 'UKE'))
        return pd.DataFrame(dict(P_removal=flag, N_removal=nr,
            capacity=capacity, region=region, annual_TP=np.expm1(y),
            outlet_easting=east, outlet_northing=north,
            n_days=12, n_months=12, nearest_outlet_SA_distance=20))

    def test_known_effect_recovered_with_required_historical_coordinates(self):
        # Removing either required spatial column must bias this known effect.
        result = model.modelling(self.frame(), {'source_gate': True})
        adjusted = result['adjusted_model']
        self.assertAlmostEqual(adjusted['conditional_log1p_TP_PRemoval_beta'], .08, places=8)
        self.assertIn('historical_outlet_easting_100km', adjusted['original_source_controls'])
        self.assertIn('historical_outlet_northing_100km', adjusted['original_source_controls'])
        self.assertIn('original_uwwNUTS1_prefix', adjusted['original_source_controls'])

    def test_descriptive_population_survives_missing_model_covariate(self):
        d = self.frame()
        d.loc[0, 'capacity'] = np.nan
        result = model.modelling(d, {'source_gate': True})
        self.assertEqual(result['n_source_PRemoval_0'], 90)
        self.assertEqual(result['n_source_PRemoval_1'], 90)
        self.assertEqual(result['adjusted_model']['n_geographically_matched_unique_2018_SA'], 179)

    def test_failed_adjustment_gate_keeps_eligible_descriptive_population(self):
        d = self.frame()
        d.loc[d.index[d.P_removal == 0][:35], 'capacity'] = np.nan
        result = model.modelling(d, {'capacity_ge70pct_in_each_PRemoval_source_group': False})
        self.assertEqual(result.get('n_source_PRemoval_0'), 90)
        self.assertEqual(result.get('n_source_PRemoval_1'), 90)
        self.assertIn('STOP', result['adjusted_model']['status'])
        self.assertNotIn('conditional_log1p_TP_PRemoval_beta', result['adjusted_model'])

    def test_make_frame_uses_original_outlet_coordinates(self):
        matches = {'synthetic_sa': {'outlet_id': ('outlet', 'synthetic_plant'),
            'outlet_coord': (400000, 300000),
            'original2018_FW_outlet_to_SA_distance_m': 12}}
        plants = {'synthetic_plant': {'PRemoval': 1, 'NRemoval': 0,
            'capacity': 10000, 'original_region': 'UKC'}}
        tps = {'synthetic_sa': {'annual_median_TP_mgL': 0,
            'n_source_dates': 2, 'n_source_months': 2}}
        frame, _, _ = model.make_frame(matches, plants, tps)
        self.assertIn('outlet_easting', frame)
        self.assertEqual(frame.loc[0, 'outlet_easting'], 400000)
        self.assertEqual(frame.loc[0, 'outlet_northing'], 300000)
        self.assertEqual(frame.loc[0, 'annual_TP'], 0)

    def test_source_exception_writes_stop_and_returns_failure(self):
        with tempfile.TemporaryDirectory() as folder:
            args = argparse.Namespace(out=Path(folder) / 'result.json')
            with patch.object(model.eea_geo, 'exact_source', side_effect=RuntimeError('PRIVATE_SITE_123 synthetic failure')):
                with contextlib.redirect_stdout(io.StringIO()):
                    status = model.main(args)
            self.assertEqual(status, 1)
            self.assertEqual(json.loads(args.out.read_text())['status'], 'SOURCE_OR_MODEL_STOP')
            self.assertNotIn('PRIVATE_SITE_123', args.out.read_text())

    def test_replaced_station_at_same_candidate_count_stops_before_TP_fetch(self):
        matches = {'synthetic_' + str(i): {} for i in range(930)}
        with tempfile.TemporaryDirectory() as folder, contextlib.ExitStack() as stack:
            args = argparse.Namespace(out=Path(folder) / 'result.json')
            stack.enter_context(patch.object(model.eea_geo, 'exact_source', return_value=({}, {}, {})))
            stack.enter_context(patch.object(model, 'plant_records', return_value=({i: {} for i in range(1799)}, {})))
            stack.enter_context(patch.object(model.ea, 'collect', return_value=({}, {})))
            stack.enter_context(patch.object(model.geo, 'select_current_SA_2018_UWWTD_historic_outlets', return_value=(matches, {})))
            fetch = stack.enter_context(patch.object(model, 'original_TP_dates', side_effect=AssertionError('should not open changed cohort')))
            with contextlib.redirect_stdout(io.StringIO()):
                status = model.main(args)
            fetch.assert_not_called()
            self.assertEqual(status, 1)

    def test_nonfinite_output_still_writes_safe_stop_artifact(self):
        import premoval_source_support
        matches = {'synthetic_' + str(i): {} for i in range(930)}
        expected = hashlib.sha256('\n'.join(sorted(matches)).encode()).hexdigest()
        with tempfile.TemporaryDirectory() as folder, contextlib.ExitStack() as stack:
            args = argparse.Namespace(out=Path(folder) / 'result.json')
            stack.enter_context(patch.object(model, 'EXPECTED_SA_COHORT_SHA256', expected, create=True))
            stack.enter_context(patch.object(model.eea_geo, 'exact_source', return_value=({}, {}, {})))
            stack.enter_context(patch.object(model, 'plant_records', return_value=({i: {} for i in range(1799)}, {})))
            stack.enter_context(patch.object(model.ea, 'collect', return_value=({}, {})))
            stack.enter_context(patch.object(model.geo, 'select_current_SA_2018_UWWTD_historic_outlets', return_value=(matches, {})))
            stack.enter_context(patch.object(model, 'original_TP_dates', return_value=({}, {}, {})))
            stack.enter_context(patch.object(model, 'original_sampling_selection_QA', return_value={}))
            stack.enter_context(patch.object(premoval_source_support, 'build_source_support', return_value={}))
            stack.enter_context(patch.object(model, 'make_frame', return_value=(None, {}, {})))
            stack.enter_context(patch.object(model, 'modelling', return_value={'status': 'EXPLORATORY_SOURCE_TP_DIFFERENCE_CALCULATED', 'CI': [float('nan'), 1]}))
            with contextlib.redirect_stdout(io.StringIO()):
                status = model.main(args)
            self.assertEqual(status, 1)
            result = json.loads(args.out.read_text())
            self.assertEqual(result['status'], 'SOURCE_OR_MODEL_STOP')
            self.assertEqual(result['source_error_stage'], 'serialize_aggregate')

    def test_record_fates_receipts_duplicates_and_daily_weighting(self):
        # Bound strings must not become exact values; export duplicates must
        # not change the daily median or give one sample day extra weight.
        def row(site, day, value, observation):
            return {'samplingPoint.notation': site,
                'samplingPoint.samplingPointType': model.original_wqe.SA,
                'sampleMaterialType': 'FINAL SEWAGE EFFLUENT',
                'determinand.notation': '0348', 'phenomenonTime': day,
                'result': value, 'unit': 'MILLIGRAM PER LITRE',
                'samplingPurpose': 'synthetic', 'observation_id': observation}
        a = row('a', '2018-01-01', '0', 'a1')
        rows = [a, dict(a), row('a', '2018-01-01', '4', 'a2'),
            row('a', '2018-02-01', '10', 'a3'),
            row('b', '2018-03-01', '<0.05', 'b1'),
            row('b', '2018-03-02', 'bad-value', 'b2')]
        receipt = dict(year=2018, n_requested_native_SA_points=3,
            skip=0, status=200, n_observations=6, sha256='a' * 64)
        with patch.object(model.wq, 'page', return_value=(rows, receipt)):
            output = model.original_TP_dates(['a', 'b', 'c'])
        self.assertEqual(len(output), 3)
        tps, qa, fates = output
        self.assertEqual(tps['a']['annual_median_TP_mgL'], 6)
        self.assertNotIn('b', tps)
        self.assertEqual(qa['exact_export_duplicate_rows_excluded'], 1)
        self.assertEqual(qa['bounded_censored_results'], 1)
        self.assertEqual(qa['malformed_or_nonnumeric_results'], 1)
        self.assertEqual(qa['numeric_zero_values'], 1)
        self.assertEqual(qa['source_page_receipts'][0]['sha256'], 'a' * 64)
        self.assertEqual(fates['a']['fate'], 'eligible_numeric_TP')
        self.assertEqual(fates['b']['fate'], 'no_eligible_numeric_TP')
        self.assertEqual(fates['c']['fate'], 'no_returned_TP_record')


if __name__ == '__main__':
    unittest.main()
