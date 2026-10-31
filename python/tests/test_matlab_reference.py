"""Independent fixtures exported from the current MATLAB source package."""
from pathlib import Path

import numpy as np
import pytest

from copa2026.betting import combination
from copa2026.experiment import conformal_pvalues, ellipsoid_density, ncms_from_arrays

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name):
    with np.load(FIXTURES / (name + ".npz"), allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


@pytest.mark.parametrize("name", ["boundaries", "ties", "alternating", "rolling_window", "tiny_distances"])
def test_complete_betting_history_matches_matlab(name):
    expected = fixture("betting_" + name)
    wealth, logs, diagnostics = combination(expected["pv"].ravel(), return_diagnostics=True)
    np.testing.assert_allclose(wealth, expected["wealth"].ravel(), rtol=2e-10, atol=2e-11)
    np.testing.assert_allclose(logs, expected["logwealth"].ravel(), rtol=2e-10, atol=2e-11)
    for field in ("Z", "candidateValues", "candidateLog", "rawKnn"):
        np.testing.assert_allclose(diagnostics[field], expected[field], rtol=2e-10, atol=2e-11)
    for field in ("selected", "fallback", "statusCode"):
        np.testing.assert_array_equal(np.asarray(diagnostics[field]).squeeze(), expected[field].squeeze())


@pytest.mark.parametrize("name", ["protein", "particle", "sensorless", "credit"])
def test_shared_benchmark_inputs_match_matlab(name):
    expected = fixture("experiment_" + name)
    scores, indices = ncms_from_arrays(expected["reference"], expected["queries"],
                                       K=50, cov_reg=1e-3, return_indices=True)
    np.testing.assert_array_equal(indices, expected["indices"])
    np.testing.assert_allclose(scores, expected["scores"].ravel(), rtol=2e-10, atol=2e-10)
    pv = conformal_pvalues(scores, expected["uniforms"].ravel())
    np.testing.assert_array_equal(pv, expected["pv"].ravel())
    wealth, logs, diag = combination(pv, return_diagnostics=True)
    np.testing.assert_allclose(wealth, expected["wealth"].ravel(), rtol=2e-10, atol=2e-11)
    np.testing.assert_allclose(logs, expected["logwealth"].ravel(), rtol=2e-10, atol=2e-11)
    np.testing.assert_array_equal(diag["selected"], expected["selected"].ravel())
    distance = ncms_from_arrays(expected["reference"], expected["queries"],
                                K=50, ncm_type="knn_distance")
    np.testing.assert_allclose(distance, expected["distance_scores"].ravel(), rtol=2e-13, atol=2e-13)
    np.testing.assert_array_equal(conformal_pvalues(distance, expected["uniforms"].ravel()),
                                  expected["distance_pvalues"].ravel())


def test_scalar_covariance_convention_matches_matlab():
    expected = fixture("ellipsoid_k1")
    actual = ellipsoid_density(expected["reference"], expected["queries"],
                               expected["idx"].astype(int) - 1, cov_reg=1e-3)
    np.testing.assert_allclose(actual, expected["density"].ravel(), rtol=2e-10, atol=2e-11)


@pytest.mark.parametrize("k", [1, 2, 3, 4])
def test_neighbor_tie_order_matches_matlab(k):
    expected = fixture("neighbor_ties")
    scores, indices = ncms_from_arrays(expected["reference"], expected["queries"],
                                       K=k, ncm_type="knn_distance", return_indices=True)
    np.testing.assert_array_equal(indices, expected["indices" + str(k)])
    np.testing.assert_allclose(scores, expected["distances" + str(k)][:, -1])
