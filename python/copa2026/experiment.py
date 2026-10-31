"""Feature-based conformal change detection, ported from the MATLAB driver.

The replay API accepts exported random inputs.  This avoids conflating numerical
agreement with the different random-permutation/normal-draw algorithms used by
MATLAB and NumPy, even when both are initialized with seed 42.
"""
from __future__ import annotations

from collections.abc import Mapping
import math
import numpy as np
from scipy.spatial import cKDTree

from .betting import combination


def _uniform_draws(count, uniforms=None, rng=None):
    if uniforms is None:
        rng = np.random if rng is None else rng
        draw = rng.random_sample if hasattr(rng, "random_sample") else rng.random
        uniforms = draw(count)
    values = np.asarray(uniforms, dtype=float).reshape(-1, order="F")
    if values.size != count or np.any(~np.isfinite(values)) or np.any((values < 0) | (values > 1)):
        raise ValueError("One finite uniform draw in [0,1] is required per computed p-value")
    return values


def pvalues(x, uniforms=None, rng=None):
    """Python equivalent of the MATLAB numeric/cell/struct pvalues helper.

    Numeric arrays are flattened in MATLAB column order and retain the
    original extended-real and NaN comparison behavior. A mapping with an
    ``ncm`` field, or a list/object array of such mappings, follows MATLAB's
    cell/struct branch: nonfinite scores are removed before choosing the last
    retained score. Missing/empty entries produce NaN without consuming a
    random draw. With no supplied ``rng``, NumPy's global stream is used.

    Explicit ``uniforms`` must supply one draw per computed p-value. The strict
    numeric experiment interface ``conformal_pvalues`` rejects NaN scores.
    """
    structural = isinstance(x, Mapping)
    if not structural:
        try:
            array = np.asarray(x)
        except ValueError:
            array = np.asarray(x, dtype=object)
        structural = array.dtype == object and any(isinstance(v, Mapping) or v is None for v in array.ravel())
        if array.dtype.names and "ncm" in array.dtype.names:
            structural = True
    if structural:
        if isinstance(x, Mapping):
            entries = [x]
        elif array.dtype.names:
            entries = array.ravel(order="F")
        else:
            entries = np.asarray(x, dtype=object).ravel(order="F")
        retained = []
        for item in entries:
            if isinstance(item, Mapping):
                scores = item.get("ncm")
            elif isinstance(item, np.void) and item.dtype.names and "ncm" in item.dtype.names:
                scores = item["ncm"]
            else:
                scores = None
            if scores is None:
                retained.append(None)
                continue
            scores = np.asarray(scores, dtype=float).reshape(-1, order="F")
            scores = scores[np.isfinite(scores)]
            retained.append(scores if scores.size else None)
        draws = iter(_uniform_draws(sum(scores is not None for scores in retained), uniforms, rng))
        result = np.full(len(retained), np.nan)
        for index, scores in enumerate(retained):
            if scores is not None:
                newest = scores[-1]
                result[index] = (np.count_nonzero(scores[:-1] > newest) + next(draws) * np.count_nonzero(scores == newest)) / scores.size
        return result
    if np.iscomplexobj(array) or array.dtype.kind not in "iuf":
        raise TypeError("Input must be a real numeric array or entries with an ncm field")
    scores = array.ravel(order="F")
    draws = _uniform_draws(scores.size, uniforms, rng)
    result = np.empty(scores.size)
    if scores.size:
        result[0] = draws[0]
    for index in range(1, scores.size):
        score = scores[index]
        result[index] = (np.count_nonzero(scores[:index] > score) + draws[index] * np.count_nonzero(scores[:index+1] == score)) / (index+1)
    return result


def conformal_pvalues(ncms, uniforms=None, rng=None):
    """Sequential smoothed p-values, including MATLAB's first-value rule."""
    alpha = np.asarray(ncms, dtype=float).reshape(-1, order="F")
    # The MATLAB numeric pvalues path compares extended-real scores directly.
    # A zero floating-point ellipsoid volume can legitimately produce -Inf;
    # its strict-order and equality comparisons still define valid tie ranks.
    if np.any(np.isnan(alpha)):
        raise ValueError("Nonconformity scores must not contain NaN")
    if uniforms is None:
        rng = np.random if rng is None else rng
        uniforms = rng.random_sample(alpha.size) if hasattr(rng, "random_sample") else rng.random(alpha.size)
    uniforms = np.asarray(uniforms, dtype=float).reshape(-1)
    if uniforms.size != alpha.size or np.any(~np.isfinite(uniforms)) or np.any((uniforms < 0) | (uniforms > 1)):
        raise ValueError("One finite uniform draw in [0,1] is required per score")
    result = np.empty(alpha.size)
    for i, score in enumerate(alpha):
        result[i] = (np.count_nonzero(alpha[:i] > score) + uniforms[i] * np.count_nonzero(alpha[:i+1] == score)) / (i+1)
    return result


def ellipsoid_density(reference, queries, indices, cov_reg=1e-3):
    """Preserve cov(Xi,1), upper Cholesky and MATLAB's D/R convention.

    This deliberately retains the source formula, including its orientation.
    Indices are zero-based and have shape (number of queries, K).
    """
    tr = np.asarray(reference, dtype=float)
    q = np.asarray(queries, dtype=float)
    idx = np.asarray(indices, dtype=np.int64)
    if tr.ndim != 2 or q.ndim != 2 or tr.shape[1] != q.shape[1] or idx.ndim != 2 or idx.shape[0] != q.shape[0]:
        raise ValueError("Incompatible reference/query/neighbor shapes")
    if idx.shape[1] == 0 or np.any(idx < 0) or np.any(idx >= tr.shape[0]) or not np.isfinite(cov_reg) or cov_reg <= 0:
        raise ValueError("Invalid neighbor indices or regularization")
    n_ref, d = tr.shape
    k = idx.shape[1]
    neighbors = tr[idx]
    # MATLAB cov(vector,1) is scalar for K=1; its scalar is expanded across
    # the matrix when lambda*eye(d) is added. Preserve that legacy convention.
    if k == 1:
        scalar = np.var(neighbors[:, 0, :], axis=1, ddof=0)
        covariance = np.broadcast_to(scalar[:, None, None], (len(q), d, d)).copy()
    else:
        centered = neighbors - neighbors.mean(axis=1, keepdims=True)
        covariance = np.matmul(centered.transpose(0, 2, 1), centered) / k
    covariance[:, np.arange(d), np.arange(d)] += cov_reg
    lower = np.linalg.cholesky(covariance)
    delta = neighbors - q[:, None, :]
    # R = lower.T, D/R = solve(lower,D.T).T -- do not transpose this solve.
    transformed = np.linalg.solve(lower, delta.transpose(0, 2, 1)).transpose(0, 2, 1)
    radius = np.sqrt(np.max(np.sum(transformed**2, axis=2), axis=1))
    unit_volume = math.pi**(d/2) / math.gamma(d/2 + 1)
    volume = unit_volume * np.maximum(radius, 1e-12)**d * np.prod(np.diagonal(lower, axis1=1, axis2=2), axis=1)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return k / (n_ref * volume)


def _matlab_neighbors(reference, queries, k):
    """KD-tree neighbors, resolving exact distance ties by source row index.

    MATLAB KNNSEARCH selects the lowest reference index at a tied boundary.
    SciPy does not promise that tie order. Requesting one extra neighbor
    detects a tie at K; only those rows need a radius query for all tied
    candidates. Interior ties are sorted by distance, then reference index.
    """
    tree = cKDTree(reference)
    search_k = min(k+1, len(reference))
    distances, indices = tree.query(queries, k=search_k, workers=1)
    if search_k == 1:
        distances, indices = distances[:, None], indices[:, None]
    order = np.lexsort((indices, distances), axis=1)
    distances = np.take_along_axis(distances, order, axis=1)
    indices = np.take_along_axis(indices, order, axis=1)
    if search_k > k:
        tied_rows = np.flatnonzero(distances[:, k-1] == distances[:, k])
        for row in tied_rows:
            radius = np.nextafter(distances[row, k-1], np.inf)
            candidates = np.unique(np.r_[tree.query_ball_point(queries[row], radius), indices[row]])
            delta = reference[candidates] - queries[row]
            candidate_distances = np.sqrt(np.sum(delta*delta, axis=1))
            keep = np.lexsort((candidates, candidate_distances))[:k]
            indices[row, :k] = candidates[keep]
            distances[row, :k] = candidate_distances[keep]
    return distances[:, :k], indices[:, :k]


def ncms_from_arrays(TR, Y, cov_reg=1e-3, K=50, ncm_type="ellipsoid_pdf", return_indices=False):
    """Compute independent SciPy nearest-neighbor queries and source NCMs."""
    tr, queries = np.asarray(TR, dtype=float), np.asarray(Y, dtype=float)
    if tr.ndim != 2 or queries.ndim != 2 or tr.shape[1] != queries.shape[1] or not np.all(np.isfinite(tr)) or not np.all(np.isfinite(queries)):
        raise ValueError("Reference and query arrays must be finite matrices with equal feature counts")
    if int(K) != K or not 1 <= K <= len(tr):
        raise ValueError("K must be an integer between one and reference size")
    distances, indices = _matlab_neighbors(tr, queries, int(K))
    ncm_type = str(ncm_type).lower()
    if ncm_type == "ellipsoid_pdf":
        density = ellipsoid_density(tr, queries, indices, cov_reg)
        values = -np.log(np.maximum(density, 1e-300))
    elif ncm_type == "knn_distance":
        values = distances[:, -1]
    else:
        raise ValueError(f"Unknown NCM: {ncm_type}")
    return (values, indices) if return_indices else values


def summarize(points, fa):
    """Rates and detected-only delay summaries; MATLAB uses sample SD."""
    points, fa = np.asarray(points), np.asarray(fa, dtype=bool)
    if points.ndim != 1 or len(points) == 0 or fa.shape != points.shape:
        raise ValueError("At least one paired per-stream point/false-alarm flag is required")
    delays = points[points > 0]
    return {
        "DetectionRate": float(np.mean(points > 0)),
        "FalseAlarmRate": float(np.mean(fa)),
        "MissRate": float(np.mean(points == 0)),
        "MedianDelay": float(np.median(delays)) if len(delays) else float("nan"),
        "MeanDelay": float(np.mean(delays)) if len(delays) else float("nan"),
        "StdDelay": float(np.std(delays, ddof=1)) if len(delays) > 1 else (0.0 if len(delays) else float("nan")),
    }


def _matrix(mapping, name, dtype=None):
    result = np.asarray(mapping[name], dtype=dtype)
    if result.ndim == 1:
        result = result[None, :]
    if result.ndim != 2:
        raise ValueError(f"{name} must be a matrix")
    return result


def run_replay(Xprepared, plan: Mapping, epsilon=100, threshold=100,
               ncm_type="ellipsoid_pdf", reference=None, K=50, cov_reg=1e-3,
               standardize=True, progress=None):
    """Run all stages independently with supplied random experimental inputs.

    `plan` keys: reference_indices, query_indices, shifts, uniforms. Indices are
    zero-based original-dataset row numbers. Query columns contain nu pre-change
    observations followed by nu post-change observations. Shifts are in
    standardized coordinates. No MATLAB score or betting output is used to
    compute Python scores; optional `reference` is used only after computation.
    """
    x = np.asarray(Xprepared, dtype=float)
    if x.ndim != 2 or not np.all(np.isfinite(x)):
        raise ValueError("Dataset must be a finite numeric matrix")
    tr_idx = _matrix(plan, "reference_indices", np.int64)
    q_idx = _matrix(plan, "query_indices", np.int64)
    shifts = _matrix(plan, "shifts", float)
    uniforms = _matrix(plan, "uniforms", float)
    reps, steps = q_idx.shape
    if steps % 2 or reps == 0 or tr_idx.shape[0] != reps or shifts.shape != (reps, x.shape[1]) or uniforms.shape != (reps, steps):
        raise ValueError("Inconsistent replay plan dimensions")
    if np.any(tr_idx < 0) or np.any(tr_idx >= len(x)) or np.any(q_idx < 0) or np.any(q_idx >= len(x)):
        raise ValueError("Replay row index out of range")
    nu = steps // 2
    points, false_alarms = np.zeros(reps, dtype=np.int64), np.zeros(reps, dtype=bool)
    scores = np.empty((reps, steps))
    ps, wealth, logs = np.empty_like(scores), np.empty_like(scores), np.empty_like(scores)
    zs = np.full((reps, steps, 3), np.nan)
    selected = np.zeros((reps, steps), dtype=np.int64)
    candidate_values, candidate_logs = np.ones((reps, steps, 7)), np.zeros((reps, steps, 7))
    fallbacks = np.zeros((reps, steps, 3), dtype=bool)
    statuses = np.zeros((reps, steps, 3), dtype=np.uint8)
    raw_knn = np.full((reps, steps, 3), np.nan)
    for rep in range(reps):
        tr, queries = x[tr_idx[rep]].copy(), x[q_idx[rep]].copy()
        if standardize:
            mu = tr.mean(axis=0)
            sigma = tr.std(axis=0, ddof=1) if len(tr) > 1 else np.zeros(tr.shape[1])
            sigma[sigma < 1e-12] = 1.0
            tr, queries = (tr-mu)/sigma, (queries-mu)/sigma
        queries[nu:] += shifts[rep]
        scores[rep] = ncms_from_arrays(tr, queries, cov_reg, K, ncm_type)
        ps[rep] = conformal_pvalues(scores[rep], uniforms[rep])
        current, log_current, diag = combination(ps[rep], epsilon, threshold, return_diagnostics=True)
        wealth[rep], logs[rep] = current, log_current
        zs[rep], selected[rep] = diag["Z"], diag["selected"]
        fallbacks[rep], statuses[rep], raw_knn[rep] = diag["fallback"], diag["statusCode"], diag["rawKnn"]
        candidate_values[rep], candidate_logs[rep] = diag["candidateValues"], diag["candidateLog"]
        hits = np.flatnonzero(current[nu:] > threshold)
        points[rep] = int(hits[0]+1) if len(hits) else 0
        false_alarms[rep] = bool(np.any(current[:nu] > threshold))
        if progress is not None:
            progress(rep+1, reps)
    result = {
        "point": points, "points": points, "fa": false_alarms,
        "summary": summarize(points, false_alarms), "ncms": scores,
        "pvalues": ps, "wealth": wealth, "logwealth": logs, "Z": zs,
        "selected": selected, "candidateValues": candidate_values,
        "candidateLog": candidate_logs, "fallback": fallbacks,
        "statusCode": statuses, "rawKnn": raw_knn,
        "params": {"reps": reps, "TR_size": tr_idx.shape[1], "nu": nu,
                   "K": K, "epsilon": epsilon, "threshold": threshold,
                   "cov_reg": cov_reg, "ncm_type": ncm_type,
                   "standardize": standardize, "betting_mode": "normalized"},
    }
    errors = {}
    if reference is not None:
        for name in ("ncms", "pvalues", "wealth", "logwealth", "Z", "candidateValues", "candidateLog"):
            if name in reference:
                a, b = np.asarray(result[name]), np.asarray(reference[name])
                if a.shape != b.shape:
                    raise ValueError(f"Reference shape differs for {name}: {a.shape} != {b.shape}")
                finite = np.isfinite(a) & np.isfinite(b)
                errors[name+"_max_abs"] = float(np.max(np.abs(a[finite]-b[finite]))) if np.any(finite) else 0.0
                errors[name+"_max_scaled"] = float(np.max(np.abs(a[finite]-b[finite]) / np.maximum(1.0, np.abs(b[finite])))) if np.any(finite) else 0.0
                errors[name+"_nonfinite_mismatch"] = int(np.count_nonzero(np.isnan(a) != np.isnan(b)) + np.count_nonzero(np.isposinf(a) != np.isposinf(b)) + np.count_nonzero(np.isneginf(a) != np.isneginf(b)))
        for name in ("selected", "point", "fa"):
            if name in reference:
                a, b = np.asarray(result[name]), np.asarray(reference[name])
                if name in ("point", "fa"):
                    b = b.reshape(-1)
                errors[name+"_mismatches"] = int(np.count_nonzero(a != b))
    result["comparison_errors"] = errors
    return result


def run_experiment(X, *, reps=100, seed=42, TR_size=4096, nu=256, K=50,
                   sigma_delta=0.1, epsilon=100, threshold=100, cov_reg=1e-3,
                   ncm_type="ellipsoid_pdf", standardize=True,
                   change_type="mean_shift", progress=None, rng=None):
    """Native NumPy experiment; its seed is not a MATLAB stream replay."""
    x = np.asarray(X, dtype=float)
    if x.ndim != 2 or len(x) < TR_size+2*nu:
        raise ValueError("Need at least TR_size + 2*nu dataset rows")
    if int(reps) != reps or reps < 1 or int(TR_size) != TR_size or TR_size < 1 or int(nu) != nu or nu < 1:
        raise ValueError("Repetitions, reference size and nu must be positive integers")
    reps, TR_size, nu = int(reps), int(TR_size), int(nu)
    change_type = str(change_type).lower()
    if change_type not in ("mean_shift", "none"):
        raise ValueError(f"Unknown change type: {change_type}")
    if not np.isscalar(seed) or not np.isfinite(seed) or int(seed) != seed:
        raise ValueError("seed must be an integer")
    seed = int(seed)
    if rng is None:
        rng = np.random.RandomState(seed)
    elif isinstance(rng, np.random.RandomState):
        # MATLAB rng(params.seed) resets the inherited generator, whose final
        # state is then used by the benchmark's next dataset preparation.
        rng.seed(seed)
    else:
        raise TypeError("rng must be a numpy.random.RandomState")
    plan = {"reference_indices": np.empty((reps, TR_size), dtype=np.int64),
            "query_indices": np.empty((reps, 2*nu), dtype=np.int64),
            "shifts": np.zeros((reps, x.shape[1])), "uniforms": np.empty((reps, 2*nu))}
    for rep in range(reps):
        permutation = rng.permutation(len(x))
        plan["reference_indices"][rep] = permutation[:TR_size]
        pool = permutation[TR_size:]
        plan["query_indices"][rep] = pool[rng.permutation(len(pool))[:2*nu]]
        if change_type == "mean_shift":
            plan["shifts"][rep] = rng.standard_normal(x.shape[1])*sigma_delta
        plan["uniforms"][rep] = rng.random_sample(2*nu)
    result = run_replay(x, plan, epsilon, threshold, ncm_type, K=K, cov_reg=cov_reg,
                        standardize=standardize, progress=progress)
    result["params"].update(seed=seed, sigma_delta=sigma_delta, change_type=change_type,
                            random_input_mode="numpy.RandomState(MT19937)")
    return result


def _run_pvalue_streams(streams, parameters):
    """Replay the same p-value matrix accepted by the MATLAB engine."""
    values = np.asarray(streams)
    if np.iscomplexobj(values):
        raise ValueError("pvalue_streams must be real")
    values = np.asarray(values, dtype=float)
    reps, nu = parameters["reps"], parameters["nu"]
    if values.shape != (reps, 2*nu) or np.any(~np.isfinite(values)) or np.any((values < 0) | (values > 1)):
        raise ValueError("pvalue_streams must have shape (reps,2*nu) and lie in [0,1]")
    shape = values.shape
    wealth, logs = np.ones(shape), np.zeros(shape)
    zs, raw_knn = np.full(shape+(3,), np.nan), np.full(shape+(3,), np.nan)
    fallbacks, statuses = np.zeros(shape+(3,), bool), np.zeros(shape+(3,), np.uint8)
    selected = np.zeros(shape, np.int64)
    candidate_values, candidate_logs = np.ones(shape+(7,)), np.zeros(shape+(7,))
    points, false_alarms = np.zeros(reps, np.int64), np.zeros(reps, bool)
    for rep, stream in enumerate(values):
        wealth[rep], logs[rep], diag = combination(stream, parameters["epsilon"], parameters["threshold"], True)
        zs[rep], raw_knn[rep] = diag["Z"], diag["rawKnn"]
        fallbacks[rep], statuses[rep] = diag["fallback"], diag["statusCode"]
        selected[rep] = diag["selected"]
        candidate_values[rep], candidate_logs[rep] = diag["candidateValues"], diag["candidateLog"]
        hits = np.flatnonzero(wealth[rep, nu:] > parameters["threshold"])
        points[rep] = int(hits[0]+1) if hits.size else 0
        false_alarms[rep] = bool(np.any(wealth[rep, :nu] > parameters["threshold"]))
    return {"point": points, "points": points, "fa": false_alarms,
            "summary": summarize(points, false_alarms), "pvalues": values.copy(),
            "wealth": wealth, "logwealth": logs, "Z": zs, "rawKnn": raw_knn,
            "fallback": fallbacks, "statusCode": statuses, "selected": selected,
            "candidateValues": candidate_values, "candidateLog": candidate_logs,
            "params": dict(parameters), "replayed": True}


def run_knn_cm_experiment(X, params=None, **kwargs):
    """Run or replay an experiment with the MATLAB engine's parameter names.

    Returns the same flat result dictionary as ``run_experiment``. ``params``
    is an optional mapping. ``significance_l`` is the capital threshold;
    ``pvalue_streams`` replays a repetitions-by-time matrix without sampling
    data or consuming an RNG. Only the supplied betting implementation is
    supported. ``rng`` may be a shared NumPy RandomState for dataset drivers.
    """
    options = {} if params is None else dict(params)
    duplicated = options.keys() & kwargs.keys()
    if duplicated:
        raise TypeError("Repeated parameters: "+", ".join(sorted(duplicated)))
    options.update(kwargs)
    if "significance_l" in options:
        if "threshold" in options:
            raise TypeError("Supply significance_l or threshold, not both")
        options["threshold"] = options.pop("significance_l")
    implementation = options.pop("implementation", "normalized")
    if implementation.lower() != "normalized":
        raise ValueError("Only the supplied betting implementation is supported")
    streams = options.pop("pvalue_streams", None)
    defaults = {"reps": 100, "seed": 42, "TR_size": 4096, "nu": 256,
                "K": 50, "sigma_delta": 0.1, "epsilon": 100,
                "threshold": 100, "standardize": True,
                "change_type": "mean_shift", "cov_reg": 1e-3,
                "ncm_type": "ellipsoid_pdf"}
    unknown = options.keys() - defaults.keys() - {"rng", "progress"}
    if unknown:
        raise TypeError("Unknown parameters: "+", ".join(sorted(unknown)))
    parameters = {**defaults, **{k: v for k, v in options.items() if k in defaults}}
    for name in ("reps", "nu"):
        value = parameters[name]
        if not np.isscalar(value) or not np.isfinite(value) or int(value) != value or value < 1:
            raise ValueError(name+" must be a positive integer")
        parameters[name] = int(value)
    for name in ("epsilon", "threshold"):
        value = parameters[name]
        if not np.isscalar(value) or not np.isfinite(value) or value <= 0:
            raise ValueError(name+" must be finite and positive")
    if streams is not None and np.size(streams):
        result = _run_pvalue_streams(streams, parameters)
    else:
        result = run_experiment(X, **parameters,
                                **{k: options[k] for k in ("rng", "progress") if k in options})
        result["replayed"] = False
    result["params"]["implementation"] = "normalized"
    return result
