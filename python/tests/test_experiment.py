import unittest
import numpy as np

from copa2026.experiment import (conformal_pvalues, ellipsoid_density, summarize,
                                 run_experiment, run_knn_cm_experiment, pvalues)


class ExperimentTests(unittest.TestCase):
    def test_sequential_tie_randomization(self):
        actual = conformal_pvalues([2, 1, 2, 1], [.2, .4, .6, .8])
        np.testing.assert_allclose(actual, [.2, .7, .4, .9], rtol=0, atol=1e-15)

    def test_density_one_dimensional_manual(self):
        # Xi={-1,1}: population variance=1; regularized variance=1.25.
        # Radius to cover query0 is 1/sqrt(1.25), giving volume2 exactly.
        actual = ellipsoid_density(np.array([[-1.], [1.]]), [[0.]], [[0, 1]], .25)
        np.testing.assert_allclose(actual, [0.5], rtol=1e-14, atol=0)

    def test_detected_only_sample_standard_deviation(self):
        summary = summarize([0, 2, 4], [True, False, False])
        self.assertEqual(summary['DetectionRate'], 2/3)
        self.assertEqual(summary['FalseAlarmRate'], 1/3)
        self.assertEqual(summary['MeanDelay'], 3)
        self.assertEqual(summary['MedianDelay'], 3)
        self.assertAlmostEqual(summary['StdDelay'], np.sqrt(2))
        self.assertEqual(summarize([0, 2], [False, False])['StdDelay'], 0)

    def test_k1_legacy_scalar_covariance(self):
        # MATLAB cov([1,3],1)=1: scalar addition gives S=[[1.5,1],[1,1.5]],
        # rather than a zero local covariance. D*S^-1*D'=7.2, det(S)=1.25.
        actual = ellipsoid_density([[1., 3.], [4., 8.]], [[0., 0.]], [[0]], .5)
        expected = 1/(2*np.pi*7.2*np.sqrt(1.25))
        np.testing.assert_allclose(actual, [expected], rtol=1e-14, atol=0)

    def test_native_repeat_and_no_global_rng_change(self):
        x = np.random.RandomState(7).normal(size=(80, 2))
        global_state = np.random.get_state()
        one = run_experiment(x, reps=2, TR_size=30, nu=8, K=5)
        two = run_experiment(x, reps=2, TR_size=30, nu=8, K=5)
        np.testing.assert_array_equal(one['pvalues'], two['pvalues'])
        np.testing.assert_array_equal(one['point'], two['point'])
        after = np.random.get_state()
        self.assertEqual(global_state[0], after[0])
        np.testing.assert_array_equal(global_state[1], after[1])
        self.assertEqual(global_state[2:], after[2:])

    def test_shared_rng_inherits_reset_and_final_state(self):
        x = np.random.RandomState(7).normal(size=(80, 2))
        shared = np.random.RandomState(16)
        first = run_experiment(x, reps=2, TR_size=30, nu=8, K=5, rng=shared)
        after_first = shared.get_state()
        shared.standard_normal(17)
        second = run_experiment(x, reps=2, TR_size=30, nu=8, K=5, rng=shared)
        np.testing.assert_array_equal(first['pvalues'], second['pvalues'])
        after_second = shared.get_state()
        np.testing.assert_array_equal(after_first[1], after_second[1])
        self.assertEqual(after_first[2:], after_second[2:])
        expected = np.random.RandomState(42)
        for _ in range(2):
            expected.permutation(80)
            expected.permutation(50)
            expected.standard_normal(2)
            expected.random_sample(16)
        np.testing.assert_array_equal(expected.random_sample(7), shared.random_sample(7))

    def test_pvalue_only_replay_and_strict_threshold(self):
        from copa2026.betting import combination
        streams = np.random.RandomState(6).random_sample((2, 32))**2
        state = np.random.get_state()
        out = run_knn_cm_experiment([], reps=2, nu=16, pvalue_streams=streams)
        after = np.random.get_state()
        np.testing.assert_array_equal(state[1], after[1])
        self.assertEqual(state[2:], after[2:])
        for rep in range(2):
            wealth, logs, diag = combination(streams[rep], return_diagnostics=True)
            np.testing.assert_array_equal(out['wealth'][rep], wealth)
            np.testing.assert_array_equal(out['logwealth'][rep], logs)
            np.testing.assert_array_equal(out['fallback'][rep], diag['fallback'])
            np.testing.assert_array_equal(out['rawKnn'][rep], diag['rawKnn'])
        out = run_knn_cm_experiment([], reps=1, nu=2, significance_l=1,
                                    pvalue_streams=np.full((1, 4), 0.5))
        self.assertFalse(out['fa'][0])
        self.assertEqual(out['points'][0], 0)
        self.assertTrue(np.isnan(out['summary']['MeanDelay']))

    def test_pvalues_numeric_column_order_and_nonfinite_behavior(self):
        scores = np.array([[2, 1], [3, 2]])
        expected = conformal_pvalues([2, 3, 1, 2], [.2, .4, .6, .8])
        np.testing.assert_array_equal(pvalues(scores, [.2, .4, .6, .8]), expected)
        np.testing.assert_array_equal(pvalues([np.nan, 1, np.nan], [.2, .4, .6]), [.2, .2, 0])
        scores = np.array([2**63+3, 2**63+1, 2**63+2], dtype=np.uint64)
        np.testing.assert_allclose(pvalues(scores, [.2, .4, .6]), [.2, .7, 1.6/3])

    def test_pvalues_cell_struct_filtering_and_draw_count(self):
        items = [None, {}, {'ncm': []}, {'ncm': [np.nan, np.inf]},
                 {'ncm': [2]}, {'ncm': [1, 2, 3]}, {'ncm': [3, np.nan]}]
        values = pvalues(items, [.2, .4, .6])
        np.testing.assert_array_equal(values[:4], np.full(4, np.nan))
        np.testing.assert_allclose(values[4:], [.2, .4/3, .6])
        structured = np.empty(2, dtype=[('ncm', object)])
        structured['ncm'] = [np.array([2]), np.array([2, 1])]
        np.testing.assert_allclose(pvalues(structured, [.2, .4]), [.2, .7])
        global_state = np.random.get_state()
        try:
            np.random.seed(72)
            expected_rng = np.random.RandomState(72)
            pvalues(items)
            expected_rng.random_sample(3)
            np.testing.assert_array_equal(np.random.random_sample(5), expected_rng.random_sample(5))
        finally:
            np.random.set_state(global_state)


if __name__ == '__main__':
    unittest.main()
