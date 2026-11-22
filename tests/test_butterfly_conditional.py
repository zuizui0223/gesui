"""Independent numerical checks of the pair-stratified conditional likelihood."""
import importlib
import math
import unittest

import numpy as np


MODULE = "analysis.20261122_butterfly_within_genus.conditional"


class ConditionalLikelihoodTests(unittest.TestCase):
    def setUp(self):
        try:
            self.evaluate = importlib.import_module(MODULE).conditional_loglike
        except ModuleNotFoundError as exc:
            if MODULE.startswith(exc.name):
                self.fail("the conditional likelihood kernel is missing")
            raise
        self.x = np.array([
            [[0., 1.], [1., -2.], [2., 3.], [3., 1.]],
            [[-1., 2.], [0.5, -1.], [2., 1.], [np.nan, np.nan]],
        ])
        self.y = np.array([[1, 0, 1, 0], [0, 1, 1, 0]])
        self.lengths = np.array([4, 3])
        self.theta = np.array([0.5, -0.3])

    def test_two_of_four_hand_enumeration(self):
        # Odds at the four sites are 1, 2, 4, 8. The six possible pairs
        # have masses 2, 4, 8, 8, 16, 32 and respective x sums 1..5.
        value, gradient = self.evaluate(
            [math.log(2)], self.x[:1, :, :1], self.y[:1], [4])
        self.assertAlmostEqual(value, math.log(70 / 4))
        self.assertAlmostEqual(gradient[0], (2 + 8 + 24 + 24 + 64 + 160) / 70 - 2)

    def test_gradient_matches_centered_finite_difference(self):
        weights = [0.75, 1.25]
        def objective(theta):
            return self.evaluate(theta, self.x, self.y, self.lengths, weights)
        _, gradient = objective(self.theta)
        steps = np.eye(2) * 1e-5
        differences = [(objective(self.theta + step)[0] - objective(self.theta - step)[0])
                       / 2e-5 for step in steps]
        np.testing.assert_allclose(gradient, differences, rtol=1e-8, atol=1e-9)

    def test_weighted_scores_and_loglikelihood(self):
        weights = [0.25, 3.]
        # At theta=0 every subset of the observed size is equally likely.
        value, gradient, scores = self.evaluate(
            [0., 0.], self.x, self.y, self.lengths, weights, return_scores=True)
        expected_scores = np.array([[-1., 2.5], [1.5, -4. / 3.]]) * np.array(weights)[:, None]
        self.assertAlmostEqual(value, 0.25 * math.log(6) + 3 * math.log(3))
        np.testing.assert_allclose(scores, expected_scores)
        np.testing.assert_allclose(gradient, -expected_scores.sum(axis=0))

    def test_stratum_specific_constants_and_padded_cells_drop_out(self):
        shifted = self.x + np.array([[[128., -256.]], [[-64., 256.]]])
        shifted[1, 3] = [np.inf, -np.inf]
        labels = self.y.copy()
        labels[1, 3] = -999
        reference = self.evaluate(self.theta, self.x, self.y, self.lengths, return_scores=True)
        changed = self.evaluate(self.theta, shifted, labels, self.lengths, return_scores=True)
        for actual, expected in zip(changed, reference):
            np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)

    def test_reversing_pair_orientation(self):
        reference = self.evaluate(self.theta, self.x, self.y, self.lengths, return_scores=True)
        reversed_pair = self.evaluate(
            self.theta, -self.x, 1 - self.y, self.lengths, return_scores=True)
        for actual, expected in zip(reversed_pair, reference):
            np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)

    def test_degenerate_and_zero_weight_strata_contribute_zero(self):
        x = np.array([[[0.], [5.]], [[7.], [4.]], [[3.], [-8.]], [[np.nan], [np.nan]]])
        y = [[0, 0], [1, 1], [0, 1], [0, 0]]
        value, gradient, scores = self.evaluate(
            [2.], x, y, [2, 2, 2, 0], [1., 1., 0., 1.], return_scores=True)
        self.assertEqual(value, 0.)
        np.testing.assert_array_equal(gradient, [0.])
        np.testing.assert_array_equal(scores, np.zeros((4, 1)))

    def test_extreme_predictors_and_unlikely_observed_subset(self):
        x = np.array([[[1000.], [900.], [-800.], [-900.]]])
        y = [[0, 0, 1, 1]]
        with np.errstate(over="raise", invalid="raise"):
            value, gradient = self.evaluate([1.], x, y, [4])
        self.assertAlmostEqual(value, 3600.)
        self.assertAlmostEqual(gradient[0], 3600.)


if __name__ == "__main__":
    unittest.main()
