"""Numerical checks independent of the normalizer's analytical integral."""

import unittest

import numpy as np
from scipy.integrate import quad

from copa2026.betting import (
    _histogram_bet, _matlab_hist_centers, combination,
    evaluate_knn, prepare_knn,
)


def independently_integrated_mass(history, k, z):
    x = np.sort(history)
    m = min(x.size, k)
    knots = np.unique(np.r_[0, 1, (x[: -m] + x[m:]) / 2,
                                (x[: x.size - m + 1] + x[m - 1 :]) / 2])
    left, right = knots[:-1], knots[1:]
    midpoint = left + (right - left) / 2
    order = np.argsort(np.abs(x[:, None] - midpoint[None, :]), axis=0, kind="stable")
    anchors = x[order[m - 1]]
    direction = np.sign(midpoint - anchors)
    r0, r1 = np.abs(left - anchors), np.abs(right - anchors)
    low_log, high_log = np.log(np.minimum(r0, r1)), np.log(np.maximum(r0, r1))
    width = high_log - low_log

    def integrand(t):
        radius = np.exp(low_log + width * t)
        u = np.clip(anchors + direction * radius, 0, 1)
        independent_distances = np.sort(np.abs(x[:, None] - u[None, :]), axis=0)
        dist = 2 * independent_distances[m - 1]
        dist[dist == 0] = 1
        raw = min(x.size - 1, k - 1) / (x.size * dist)
        return np.sum(raw * radius * width) / z

    return quad(integrand, 0, 1, epsabs=2e-10, epsrel=2e-9, limit=1000)[0]


class NormalizationTests(unittest.TestCase):
    def test_multiple_histories(self):
        rng = np.random.default_rng(74129)
        histories = [np.array([0., 1.]), np.array([0.2, 0.7]), np.linspace(.1, .9, 5),
                     np.linspace(0, 1, 100), rng.random(100),
                     np.r_[np.linspace(0, .01, 50), np.linspace(.99, 1, 50)],
                     np.sort(np.r_[.15+.002*rng.random(45), .75+.005*rng.random(55)]),
                     np.repeat(np.linspace(0, 1, 25), 4), .5+np.arange(40)*1e-9,
                     np.array([0., np.nextafter(0., 1.)]),
                     np.array([0., 4*np.nextafter(0., 1.)]), np.full(30, .4),
                     np.r_[np.zeros(15), np.linspace(.1, 1, 15)], np.array([.4]), np.empty(0)]
        for history in histories:
            for k in (5, 10, 15):
                with self.subTest(n=history.size, k=k):
                    model = prepare_knn(history, k)
                    query = np.unique(np.r_[0, 1, np.linspace(0, 1, 101), history])
                    bet, raw = evaluate_knn(model, query)
                    self.assertTrue(np.all(np.isfinite(bet) & (bet >= 0)))
                    if model.fallback:
                        np.testing.assert_array_equal(bet, np.ones(query.size))
                    else:
                        mass = independently_integrated_mass(history, k, model.Z)
                        self.assertLess(abs(mass - 1), 2e-8)
                        sorted_distances = np.sort(np.abs(history[:, None] - query[None, :]), axis=0)
                        dist = 2 * sorted_distances[min(history.size, k)-1]
                        dist[dist == 0] = 1
                        expected = min(history.size-1, k-1) / (history.size * dist)
                        np.testing.assert_array_equal(raw, expected)

    def test_closed_form_and_invalid(self):
        for k in (5, 10, 15):
            self.assertAlmostEqual(prepare_knn([0, 1], k).Z, np.log(2)/2, places=15)
        for bad in ([.2, np.nan], [.2, np.inf], [-1e-15, .4], [.2, 1+1e-15]):
            with self.assertRaises(ValueError):
                prepare_knn(bad, 5)
        for bad_k in (0, 1, 2.2, np.inf):
            with self.assertRaises(ValueError):
                prepare_knn([0, 1], bad_k)


class CautiousBettingTests(unittest.TestCase):
    def test_matlab_histogram_midpoint_ties_and_colon(self):
        # R2025b oracle: HIST([.25 .5 .75],[.125 .375 .625 .875])
        # returns [1 1 1 0], so exact midpoint ties belong to the left bin.
        h = np.array([.125, .375, .625, .875, .25, .5, .75])
        expected = [8/7, 8/7, 8/7, 4/7, 6/7, 8/7, 4/7]
        for q, bet in zip([.125, .375, .625, .875, .75, 0, 1], expected):
            self.assertAlmostEqual(_histogram_bet(h, q, 4), bet, places=14)
        # Oracle values exercise both symmetric colon filling and rounding.
        np.testing.assert_array_equal(_matlab_hist_centers(3), [1/6, .5, .8333333333333333])
        np.testing.assert_array_equal(_matlab_hist_centers(5), [.1, .30000000000000004, .5, .7, .9])
        np.testing.assert_array_equal(_matlab_hist_centers(10),
            [.05, .15000000000000002, .25, .35000000000000003, .45,
             .55, .65, .75, .8500000000000001, .9500000000000001])
        for bins in range(1, 16):
            centers = _matlab_hist_centers(bins)
            edges = (centers[:-1] + centers[1:])/2
            history = np.r_[centers, edges]
            # Every bin is occupied, so no decrement occurs. Every midpoint
            # adds to its left bin: counts=[2,...,2,1] for bins>1.
            expected_density = np.r_[np.full(bins-1, 2.), 1.] * bins/(2*bins-1)
            for q, value in zip(centers, expected_density):
                self.assertAlmostEqual(_histogram_bet(history, q, bins), value, places=13)

    def test_histograms_candidate_histories_and_predictability(self):
        rng = np.random.default_rng(781)
        p = np.r_[rng.random(160), rng.random(160)**3]
        wealth, log_wealth, d = combination(p, return_diagnostics=True)
        for i in range(2, len(p)):
            history = p[max(0, i-200):i]
            expected = [_histogram_bet(history, p[i], bins) for bins in (5, 10, 15)]
            np.testing.assert_allclose(d['candidateValues'][i, 4:], expected, rtol=0, atol=1e-14)
        np.testing.assert_array_equal(d['candidateValues'][:, 0], np.ones(len(p)))
        np.testing.assert_allclose(np.cumsum(np.log(d['candidateValues']), axis=0), d['candidateLog'], atol=1e-12, rtol=0)
        changed = p[:130].copy()
        changed[-1] = 1-changed[-1]
        _, _, dc = combination(changed, return_diagnostics=True)
        np.testing.assert_array_equal(d['Z'][:130], dc['Z'])
        np.testing.assert_array_equal(d['selected'][:130], dc['selected'])
        np.testing.assert_allclose(np.log(wealth), log_wealth, atol=1e-12, rtol=0)
        self.verify_selection(d)

    def verify_selection(self, d):
        for i in range(2, d['candidateLog'].shape[0]):
            gains = d['candidateLog'][i-1] - np.min(d['candidateLog'][max(0, i-5000):i], axis=0)
            best = np.argmax(gains)
            expected = best+1 if gains[best] > np.log(100) else 0
            self.assertEqual(d['selected'][i], expected)

    def test_window_rollover_and_degenerates(self):
        p = np.random.default_rng(912).random(5010)**2
        _, log_wealth, d = combination(p, return_diagnostics=True)
        self.verify_selection(d)
        self.assertTrue(np.all(np.isfinite(log_wealth)))
        wealth, logs, d = combination(np.full(32, .5), return_diagnostics=True)
        self.assertTrue(np.all(d['fallback'][2:]))
        np.testing.assert_array_equal(wealth, np.ones(32))
        np.testing.assert_array_equal(logs, np.zeros(32))
        for n in range(3):
            wealth, logs, d = combination(np.full(n, .5), return_diagnostics=True)
            self.assertEqual(d['Z'].shape, (n, 3))
            np.testing.assert_array_equal(wealth, np.ones(n))


if __name__ == '__main__':
    unittest.main()
