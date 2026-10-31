"""Independent checks for the batched candidate path and reference fallback."""

import numpy as np
import pytest

from copa2026.betting import _combination, combination


def compare_reference(pvalues, epsilon=100.0):
    expected_m, expected_log, expected = _combination(
        pvalues, epsilon, 100.0, True
    )
    actual_m, actual_log, actual = combination(
        pvalues, epsilon=epsilon, return_diagnostics=True
    )
    np.testing.assert_allclose(actual_m, expected_m, rtol=1e-8, atol=1e-10)
    np.testing.assert_allclose(actual_log, expected_log, rtol=0, atol=1e-8)
    for key in ("Z", "candidateValues", "candidateLog"):
        np.testing.assert_allclose(
            actual[key], expected[key], rtol=1e-9, atol=1e-10, equal_nan=True
        )
    for key in ("selected", "fallback", "statusCode", "rawKnn"):
        np.testing.assert_array_equal(actual[key], expected[key])
    return actual


@pytest.mark.parametrize("length", [2, 3, 14, 15, 16, 99, 100, 101, 512, 5000, 5001, 5010])
def test_batch_boundaries_and_reference_fallback(length):
    pvalues = np.random.RandomState(197 + length).random_sample(length) ** 2
    compare_reference(pvalues)


def test_degenerate_histories_and_midpoint_ties():
    compare_reference(np.full(180, 0.5))
    compare_reference(np.tile([0.0, 0.25, 0.5, 0.75, 1.0, 0.125, 0.375], 30))


def test_tightly_clustered_variable_length_histories():
    for length in (101, 512):
        compare_reference(0.5 + np.arange(length) * 1e-9)


def test_history_only_normalizers_and_selection():
    pvalues = np.random.RandomState(761).random_sample(512)
    _, _, original = combination(pvalues, return_diagnostics=True)
    changed = pvalues.copy()
    changed[220:] = 1 - changed[220:]
    _, _, after = combination(changed, return_diagnostics=True)
    # At position 220, the current p-value affects its factor but neither its
    # normalization nor the candidate chosen from earlier candidate histories.
    np.testing.assert_array_equal(original["Z"][:221], after["Z"][:221])
    np.testing.assert_array_equal(original["selected"][:221], after["selected"][:221])


def test_strict_cautious_comparison():
    # All candidates are exactly uniform for a repeated atom. Gain==log(1)
    # must not bet when epsilon==1; epsilon<1 chooses the first tied candidate.
    pvalues = np.full(120, 0.5)
    at_boundary = compare_reference(pvalues, epsilon=1.0)
    np.testing.assert_array_equal(at_boundary["selected"], np.zeros(120))
    below_boundary = compare_reference(pvalues, epsilon=0.5)
    np.testing.assert_array_equal(below_boundary["selected"][2:], np.ones(118))


def test_matlab_tiny_distance_underflow_is_retained():
    # KNNSEARCH's Euclidean squaring underflows at this scale. The current
    # MATLAB implementation patches the resulting zero distance to one;
    # preserve that known source behavior rather than silently fixing it.
    pvalues = np.array([0., 1e-200, 5e-201])
    _, _, diagnostics = combination(pvalues, return_diagnostics=True)
    np.testing.assert_array_equal(diagnostics['rawKnn'][2], np.full(3, 0.5))
    assert not diagnostics['fallback'][2].any()
    np.testing.assert_allclose(diagnostics['Z'][2], 115.47582823998226, rtol=1e-14)
    compare_reference(np.arange(180)*1e-200)
