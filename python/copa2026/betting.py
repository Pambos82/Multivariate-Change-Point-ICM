"""History-only kNN normalization and the MATLAB cautious betting strategy.

The operational integrals use analytical logarithmic antiderivatives evaluated
in float64. They are not an exact-arithmetic or interval certificate. Histories
whose original density has no usable finite integral use a complete uniform
betting function. The histogram formulas follow the legacy MATLAB strategy.
"""

from dataclasses import dataclass

import numpy as np


CANDIDATE_NAMES = ("uniform", "knn15", "knn10", "knn5", "hist5", "hist10", "hist15")
STATUS_NAMES = ("startup", "analytic", "uniform_small_sample", "uniform_divergent", "uniform_numerical_failure")


def _probability_array(values, *, vector=False):
    raw = np.asarray(values)
    if np.iscomplexobj(raw):
        raise ValueError("p-values must be real")
    values = np.asarray(raw, dtype=np.float64)
    if vector and (values.ndim > 2 or (values.ndim == 2 and min(values.shape, default=0) > 1)):
        raise ValueError("p-values must be a vector")
    if np.any(~np.isfinite(values)) or np.any((values < 0) | (values > 1)):
        raise ValueError("p-values must be finite and lie in [0, 1]")
    return values


@dataclass(frozen=True)
class KnnBettingModel:
    history_sorted: np.ndarray
    N: int
    k: int
    m: int
    coefficient: float
    Z: float
    fallback: bool
    status: str
    transitions: np.ndarray


def prepare_knn(history, k):
    """Prepare the integral of the original 1D kNN estimate from history only.

    The original m-th nearest distance is used without boundary reflection or
    radius truncation. Repeated atoms of multiplicity >= m cause a divergent
    integral even though the legacy estimator patches its value at radius zero.
    """
    if not np.isscalar(k) or not np.isfinite(k) or k < 2 or int(k) != k:
        raise ValueError("k must be an integer >= 2")
    k = int(k)
    x = np.sort(_probability_array(history, vector=True).reshape(-1))
    n = x.size
    m = min(n, k)
    coefficient = min(n - 1, k - 1) / (2 * n) if n else np.nan
    transitions = np.empty(0)

    def model(z, status, fallback=True):
        return KnnBettingModel(x, n, k, m, coefficient, float(z), fallback, status, transitions)

    if n < 2:
        return model(np.nan if n == 0 else 0.0, "uniform_small_sample")
    transitions = x[: n - m] + (x[m:] - x[: n - m]) / 2
    first = x[: n - m + 1]
    last = x[m - 1 :]
    if np.any(first == last):
        return model(np.inf, "uniform_divergent")
    with np.errstate(over="ignore", divide="ignore", invalid="ignore", under="ignore"):
        minimum_radius = np.min(last - first) / 2
        maximum_raw = coefficient / minimum_radius
        if minimum_radius <= 0 or not np.isfinite(maximum_raw):
            return model(np.nan, "uniform_numerical_failure")
        left = np.r_[0.0, transitions]
        right = np.r_[transitions, 1.0]
        middle = first + (last - first) / 2
        low_end = np.minimum(right, middle)
        low_width = np.maximum(0, low_end - left)
        low_radius = last - low_end
        high_start = np.maximum(left, middle)
        high_width = np.maximum(0, right - high_start)
        high_radius = high_start - first
        low_active = low_width > 0
        high_active = high_width > 0
        terms = np.r_[
            np.log1p(low_width[low_active] / low_radius[low_active]),
            np.log1p(high_width[high_active] / high_radius[high_active]),
        ]
        z = coefficient * np.sum(terms)
        if not np.isfinite(z) or z <= 0 or np.any(~np.isfinite(terms)) or not np.isfinite(maximum_raw / z):
            return model(z, "uniform_numerical_failure")
    return model(z, "analytic", False)


def evaluate_knn(model, p):
    """Return (normalized bet, original raw value) at the supplied queries.

    Normalizers and fallback decisions were already determined by prepare_knn.
    The dist=1 convention at zero radius is retained literally in the raw value.
    """
    query = _probability_array(p)
    shape = query.shape
    flat = query.reshape(-1)
    if model.N == 0:
        raw = np.full(flat.shape, np.nan)
    elif model.N == 1:
        raw = np.zeros(flat.shape)
    else:
        # MATLAB KNNSEARCH computes Euclidean distances through squares.
        # Retain its tiny-distance underflow behavior for source parity.
        with np.errstate(under="ignore"):
            distances = np.sqrt(np.square(model.history_sorted[:, None] - flat[None, :]))
        radius = np.partition(distances, model.m - 1, axis=0)[model.m - 1]
        dist = 2 * radius
        dist[dist == 0] = 1
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            raw = min(model.N - 1, model.k - 1) / (model.N * dist)
    if model.fallback:
        value = np.ones(flat.shape)
    else:
        value = raw / model.Z
        if np.any(~np.isfinite(value) | (value < 0)):
            raise FloatingPointError("normalized kNN value cannot be represented as a finite nonnegative float64")
    return value.reshape(shape), raw.reshape(shape)


def _matlab_hist_centers(bins):
    """Reproduce ``1/(2*bins):1/bins:1`` in MATLAB float64 arithmetic.

    MATLAB's colon fills symmetrically from both endpoints. Direct division of
    (index + 0.5) by bins changes several centers by one ulp and can change
    exact-midpoint bin membership. There are exactly ``bins`` centers here.
    """
    first = 1 / (bins * 2)
    step = 1 / bins
    last = first + (bins - 1) * step
    centers = np.empty(bins)
    for j in range(bins // 2):
        centers[j] = first + j * step
        centers[-j - 1] = last - j * step
    if bins % 2:
        centers[bins // 2] = first + (bins // 2) * step
    return centers


_HIST_CENTERS = tuple(_matlab_hist_centers(bins) for bins in range(1, 16))


def _histogram_bet(history, p, maximum_bins):
    """Legacy HIST centers, decrement-until-nonempty rule and interpolation."""
    for bins in range(maximum_bins, 0, -1):
        centers = _HIST_CENTERS[bins - 1]
        # Reproduce HIST's edge arithmetic and eps shift before its left-
        # inclusive HISTC count. An algebraic midpoint changes one-ulp ties.
        edges = centers[:-1] + np.diff(centers) / 2
        edges = edges + np.spacing(edges)
        counts = np.bincount(np.searchsorted(edges, history, side="right"), minlength=bins)
        if np.all(counts > 0):
            break
    density = counts / np.sum(counts) * bins
    p = max(float(p), 0.00001)
    distances = np.abs(centers - p)
    active = distances < 1 / bins
    weights = 1 / bins - distances[active]
    return float((weights / np.sum(weights)) @ density[active])


def _combination(pvalues, epsilon, threshold, return_diagnostics):
    pvalues = _probability_array(pvalues, vector=True).reshape(-1)
    for name, value in (("epsilon", epsilon), ("threshold", threshold)):
        if not np.isscalar(value) or not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    n = pvalues.size
    wealth = np.ones(n)
    log_wealth = np.zeros(n)
    candidate_log = np.zeros((n, 7))
    diagnostics = None
    if return_diagnostics:
        diagnostics = {
            "k": np.array([5, 10, 15]),
            "Z": np.full((n, 3), np.nan),
            "fallback": np.zeros((n, 3), dtype=bool),
            "statusCode": np.zeros((n, 3), dtype=np.uint8),
            "statusNames": STATUS_NAMES,
            "selected": np.zeros(n, dtype=np.int64),
            "rawKnn": np.full((n, 3), np.nan),
            "candidateValues": np.ones((n, 7)),
            "candidateLog": candidate_log,
            "candidateNames": CANDIDATE_NAMES,
            "epsilon": epsilon,
            "W": 5000,
            "threshold": threshold,
        }
    log_epsilon = np.log(epsilon)
    for i in range(2, n):
        history = pvalues[max(0, i - 100) : i]
        models = [prepare_knn(history, k) for k in (5, 10, 15)]
        values = np.ones(7)
        raw_knn = np.empty(3)
        for j, model in enumerate(models):
            value, raw = evaluate_knn(model, pvalues[i])
            raw_knn[j] = float(raw)
            values[3 - j] = float(value)
        hist_history = pvalues[max(0, i - 200) : i]
        for j, bins in enumerate((5, 10, 15)):
            values[4 + j] = _histogram_bet(hist_history, pvalues[i], bins)
        gains = candidate_log[i - 1] - np.min(candidate_log[max(0, i - 5000) : i], axis=0)
        best = int(np.argmax(gains))
        selected = best + 1 if gains[best] > log_epsilon else 0
        bet = values[best] if selected else 1.0
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            candidate_log[i] = candidate_log[i - 1] + np.log(values)
            wealth[i] = wealth[i - 1] * bet
            log_wealth[i] = log_wealth[i - 1] + np.log(bet)
        if diagnostics is not None:
            diagnostics["Z"][i] = [model.Z for model in models]
            diagnostics["fallback"][i] = [model.fallback for model in models]
            diagnostics["statusCode"][i] = [STATUS_NAMES.index(model.status) for model in models]
            diagnostics["selected"][i] = selected
            diagnostics["rawKnn"][i] = raw_knn
            diagnostics["candidateValues"][i] = values
    return wealth, log_wealth, diagnostics


def combination(pvalues, epsilon=100.0, threshold=100.0, return_diagnostics=False):
    """Return (wealth, log_wealth, diagnostics_or_None) for normalized betting.

    Candidate indexes in diagnostics follow MATLAB: 1..7, with 0 meaning no
    cautious bet. The first two positions are startup wealth=1. Selection uses
    only past candidate log histories; the current p-value affects factors only.
    """
    # Import lazily so the batch kernel can reuse the validated reference
    # helpers without a module-initialization cycle. Long streams retain the
    # reference rolling-window path; all normalizers remain history-only.
    from ._fast_betting import fast_combination
    return fast_combination(pvalues, epsilon, threshold, return_diagnostics)
