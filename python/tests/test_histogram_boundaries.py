"""Legacy HIST boundary regressions against independent MATLAB outputs."""
from pathlib import Path

import numpy as np
import pytest

from copa2026.betting import _HIST_CENTERS, _combination, _histogram_bet, combination

FIXTURE = Path(__file__).parent / "fixtures" / "histogram_boundaries.npz"


@pytest.mark.parametrize("bins", range(1, 16))
def test_histogram_boundary_memberships_counts_and_factors_match_matlab(bins):
    with np.load(FIXTURE, allow_pickle=False) as expected:
        centers = expected[f"centers_{bins:02d}"]
        probes = expected[f"probes_{bins:02d}"]
        memberships = expected[f"memberships_{bins:02d}"]
        counts = expected[f"counts_{bins:02d}"]
        factors = expected[f"factors_{bins:02d}"]
    np.testing.assert_array_equal(_HIST_CENTERS[bins - 1], centers)

    # Every center is present, so the histogram keeps the requested bin count.
    # At centers its interpolated factors recover the counts after rounding.
    # This checks membership through the runtime function, without duplicating
    # its edge arithmetic in the test.
    for probe, matlab_bin in zip(probes, memberships):
        history = np.r_[centers, probe]
        actual = np.array([_histogram_bet(history, query, bins) for query in centers])
        actual_counts = np.rint(actual * len(history) / bins).astype(np.int64)
        expected_counts = np.ones(bins, dtype=np.int64)
        expected_counts[matlab_bin - 1] += 1
        np.testing.assert_array_equal(actual_counts, expected_counts,
                                      err_msg=f"bins={bins}, p={float(probe).hex()}")

    actual_factors = np.array([_histogram_bet(probes, query, bins) for query in centers])
    np.testing.assert_allclose(actual_factors, factors, rtol=2e-13, atol=2e-13)
    np.testing.assert_array_equal(np.rint(actual_factors * len(probes) / bins).astype(np.int64), counts)


@pytest.mark.parametrize("implementation", [combination, _combination], ids=["batched", "rolling"])
@pytest.mark.parametrize("case", [1032, 1033])
def test_boundary_stream_candidates_and_cautious_selection_match_matlab(implementation, case):
    with np.load(FIXTURE, allow_pickle=False) as expected:
        prefix = f"case{case}_"
        pvalues = expected[prefix + "pvalues"]
        reference = {name: expected[prefix + name] for name in
                     ("wealth", "logwealth", "candidateValues", "candidateLog", "selected")}
    wealth, logs, diagnostics = implementation(pvalues, 100, 100, True)
    # Case1033's final position exercises a real cautious-choice discrepancy
    # from the former bin arithmetic, beyond a small numerical factor error.
    np.testing.assert_array_equal(diagnostics["selected"], reference["selected"])
    np.testing.assert_array_equal(wealth > 100, reference["wealth"] > 100)
    for name in ("candidateValues", "candidateLog"):
        np.testing.assert_allclose(diagnostics[name], reference[name], rtol=2e-10, atol=2e-10)
    np.testing.assert_allclose(wealth, reference["wealth"], rtol=2e-10, atol=2e-10)
    np.testing.assert_allclose(logs, reference["logwealth"], rtol=2e-10, atol=2e-10)
