"""MATLAB-compatible extended-real ranks at the Credit K=1 corner."""
import numpy as np
import pytest

from copa2026.experiment import conformal_pvalues, ncms_from_arrays


def test_infinite_ties_match_matlab_numeric_pvalues():
    # Verified by a fresh R2025b pvalues.m call with rng(42,'twister').
    alpha = [-np.inf, 0, np.inf, np.inf, -np.inf, 0]
    uniforms = np.random.RandomState(42).random_sample(6)
    expected = [uniforms[0], uniforms[1]/2, uniforms[2]/3,
                uniforms[3]/2, (3+2*uniforms[4])/5,
                (2+2*uniforms[5])/6]
    np.testing.assert_array_equal(conformal_pvalues(alpha, uniforms), expected)
    with pytest.raises(ValueError, match="NaN"):
        conformal_pvalues([0, np.nan], [.5, .5])


def test_zero_volume_high_dimension_remains_a_rankable_score():
    # radius=0 is floored at 1e-12 by the original NCM; its 29th power
    # underflows, producing density=Inf and score=-Inf exactly as in MATLAB.
    reference = np.vstack([np.arange(29.), np.arange(29.)+1])
    scores = ncms_from_arrays(reference, reference, K=1)
    assert np.isneginf(scores).all()
    np.testing.assert_array_equal(conformal_pvalues(scores, [.25, .75]), [.25, .75])
