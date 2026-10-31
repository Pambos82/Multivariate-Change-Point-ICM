"""Four-dataset benchmark, p-value replay and MATLAB-compatible summaries."""
from __future__ import annotations

from datetime import datetime, timezone
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy

from . import experiment
from .data import DATASETS, ROOT, dataset_directory, load_dataset, verify_datasets
from .io import load_benchmark_replay, output_directory, save_benchmark_archive, save_csv, save_json

DISPLAY_COLUMNS = ("Dataset", "Repetitions", "DetectionRate", "FalseAlarmRate",
                   "MissRate", "MedianDelay", "MeanDelay", "StdDelay")
PAPER = {
    "protein": (.84, .01, 83.9), "particle": (1., .01, 15.1),
    "sensorless": (.99, 0., 60.7), "credit": (.99, .01, 83.5),
}
CANDIDATE_NAMES = ("cautious/no bet", "uniform", "knn15", "knn10", "knn5",
                   "histogram5", "histogram10", "histogram15")


def _sample_sd(x):
    return float(np.std(x, ddof=1)) if len(x) > 1 else (0.0 if len(x) else float("nan"))


def _results_mapping(outputs):
    if isinstance(outputs, dict):
        return outputs
    return {item.get("Dataset", item.get("dataset", "")): item for item in outputs}


def summarize_results(outputs):
    """Return the same five detailed tables as the MATLAB summary function.

    Each table is a list of dictionaries, suitable for CSV or DataFrame use.
    Startup positions one and two are omitted from betting diagnostics.
    """
    tables = {name: [] for name in ("Performance", "NormalizationConstants",
                                    "SelectionFrequencies", "PerRepetition", "NumericalIssues")}
    for dataset, result in _results_mapping(outputs).items():
        points = np.asarray(result["point"], dtype=int).reshape(-1)
        fa = np.asarray(result["fa"], dtype=bool).reshape(-1)
        reps = len(points)
        summary = experiment.summarize(points, fa)
        paper_detection, paper_fa, paper_mean = PAPER.get(dataset, (np.nan, np.nan, np.nan))
        row = {"Dataset": dataset, "Implementation": "normalized", "Repetitions": reps, **summary,
               "PaperDetectionRate": paper_detection, "PaperFalseAlarmRate": paper_fa,
               "PaperMeanDelay": paper_mean,
               "DetectionMinusPaper": summary["DetectionRate"] - paper_detection,
               "FalseAlarmMinusPaper": summary["FalseAlarmRate"] - paper_fa,
               "MeanDelayMinusPaper": summary["MeanDelay"] - paper_mean}
        tables["Performance"].append(row)
        for rep in range(reps):
            tables["PerRepetition"].append({
                "Dataset": dataset, "Implementation": "normalized", "Repetition": rep + 1,
                "Detected": bool(points[rep] > 0), "FalseAlarm": bool(fa[rep]),
                "DelayZeroIfMiss": int(points[rep]),
            })
        wealth = np.asarray(result["wealth"])
        logwealth = np.asarray(result["logwealth"])
        tables["NumericalIssues"].append({
            "Dataset": dataset, "Implementation": "normalized",
            "InfiniteLinearValues": int(np.count_nonzero(np.isinf(wealth))),
            "NaNLinearValues": int(np.count_nonzero(np.isnan(wealth))),
            "ZeroLinearValues": int(np.count_nonzero(wealth == 0)),
            "NonfiniteLogValues": int(np.count_nonzero(~np.isfinite(logwealth))),
        })
        selected = np.asarray(result["selected"], dtype=int)
        steps = wealth.shape[1]
        nu = int(result["params"]["nu"])
        if selected.shape != (reps, steps) or np.any((selected < 0) | (selected > 7)):
            raise ValueError("Unexpected candidate selection codes or dimensions")
        for phase, indices in (("all", slice(2, steps)), ("pre_change", slice(2, nu)),
                               ("post_change", slice(nu, steps))):
            decisions = selected[:, indices]
            denominator = decisions.size
            active = int(np.count_nonzero(decisions > 0))
            for code, label in enumerate(CANDIDATE_NAMES):
                count = int(np.count_nonzero(decisions == code))
                tables["SelectionFrequencies"].append({
                    "Dataset": dataset, "Implementation": "normalized", "Phase": phase,
                    "CandidateCode": code, "Candidate": label, "Count": count,
                    "DecisionCount": denominator, "ActiveBetCount": active,
                    "FrequencyAllDecisions": count / denominator if denominator else float("nan"),
                    "FrequencyAmongActiveBets": count / active if code > 0 and active else float("nan"),
                })
        z = np.asarray(result["Z"], dtype=float)[:, 2:, :].reshape(-1, 3)
        fallback = np.asarray(result["fallback"], dtype=bool)[:, 2:, :].reshape(-1, 3)
        for col, k in enumerate((5, 10, 15)):
            current = z[:, col]
            finite = current[np.isfinite(current)]
            # MATLAB's default quantiles use the empirical midpoint/Hazen rule.
            q = np.quantile(finite, [0, .01, .05, .25, .5, .75, .95, .99, 1],
                            method="hazen") if len(finite) else np.full(9, np.nan)
            row = {"Dataset": dataset, "K": k, "DecisionCount": len(current),
                   "FiniteCount": len(finite), "NaNCount": int(np.count_nonzero(np.isnan(current))),
                   "InfiniteCount": int(np.count_nonzero(np.isinf(current))),
                   "FallbackCount": int(np.count_nonzero(fallback[:, col]))}
            row.update(zip(("Min", "Q01", "Q05", "Q25", "Median", "Q75", "Q95", "Q99", "Max"),
                           (float(value) for value in q)))
            row.update(Mean=float(np.mean(finite)) if len(finite) else float("nan"),
                       Std=_sample_sd(finite))
            tables["NormalizationConstants"].append(row)
    return tables


def _print_table(rows):
    values = [[str(row[column]) if column == "Dataset" else
               (str(int(row[column])) if column == "Repetitions" else f"{row[column]:.8g}")
               for column in DISPLAY_COLUMNS] for row in rows]
    widths = [max(len(column), *(len(row[i]) for row in values))
              for i, column in enumerate(DISPLAY_COLUMNS)]
    print("  ".join(column.ljust(widths[i]) for i, column in enumerate(DISPLAY_COLUMNS)))
    for row in values:
        print("  ".join(value.ljust(widths[i]) for i, value in enumerate(row)))


def _write_report(path, tables, metadata):
    lines = ["# kNN betting benchmark", "", f"Generated UTC: {metadata['CreatedUTC']}. Executed in Python.", "",
             f"Cautious betting; epsilon={metadata['Options']['Epsilon']:g}, "
             f"alarm threshold={metadata['Options']['Threshold']:g}, W=5000.",
             "Reference size=4096, two segments of 256, NCM K=50, covariance diagonal adjustment=1e-3, seed=42.", "",
             "| Dataset | Repetitions | Detection | False alarm | Mean delay | Median delay | SD delay |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for row in tables["Performance"]:
        lines.append(f"| {row['Dataset']} | {row['Repetitions']} | {row['DetectionRate']:.4f} | "
                     f"{row['FalseAlarmRate']:.4f} | {row['MeanDelay']:.8g} | {row['MedianDelay']:.8g} | "
                     f"{row['StdDelay']:.8g} |")
    lines += ["", "Delay statistics are conditional on post-change detection. False alarms are strict "
              "threshold crossings in the first segment. The first two betting positions are startup.",
              "Fewer than 100 repetitions is a smoke or exploratory run.", "",
              "## Output files", "",
              "The benchmark.npz archive preserves exact p-values, candidate histories, wealth, diagnostics "
              "and settings. Dataset archives contain the corresponding individual results.",
              "CSV files report performance and paper-reference differences, betting diagnostics, candidate "
              "selection frequencies, per-repetition decisions, and numerical issues.", "",
              "Diagnostics omit startup positions. Uniform fallback is selected from history when needed.",
              "Random-generator and library metadata are retained in metadata.json. Native Python and MATLAB "
              "random streams do not become identical merely by using the same numeric seed.", ""]
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def _replay_result(streams, parameters, epsilon, threshold):
    streams = np.asarray(streams, dtype=float)
    if streams.ndim != 2 or streams.shape[1] % 2:
        raise ValueError("Saved p-value streams must be repetitions by twice nu")
    options = dict(reps=streams.shape[0], seed=int(parameters.get("seed", 42)),
                   TR_size=int(parameters.get("TR_size", 4096)), nu=streams.shape[1] // 2,
                   K=int(parameters.get("K", 50)), sigma_delta=float(parameters.get("sigma_delta", .1)),
                   epsilon=epsilon, significance_l=threshold,
                   standardize=bool(parameters.get("standardize", True)),
                   change_type=str(parameters.get("change_type", "mean_shift")),
                   cov_reg=float(parameters.get("cov_reg", 1e-3)),
                   ncm_type=str(parameters.get("ncm_type", "ellipsoid_pdf")), pvalue_streams=streams)
    return experiment.run_knn_cm_experiment(np.empty((0, 0)), **options)


def run_benchmarks(*, reps=100, epsilon=100.0, threshold=100.0, data_dir=None,
                   output_root=None, replay_file=None):
    """Execute four benchmarks or replay saved p-values; print only the final table."""
    if isinstance(reps, bool) or int(reps) != reps or reps < 1:
        raise ValueError("Reps must be a positive integer")
    for label, value in (("epsilon", epsilon), ("threshold", threshold)):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{label} must be finite and positive")
    directory = output_directory("benchmark", output_root or ROOT / "results" / "benchmarks")
    started = time.perf_counter()
    metadata = {
        "CreatedUTC": datetime.now(timezone.utc).isoformat(),
        "PythonVersion": platform.python_version(), "NumPyVersion": np.__version__,
        "SciPyVersion": scipy.__version__, "Generator": "numpy.RandomState(MT19937)",
        "Options": {"Reps": int(reps), "Epsilon": float(epsilon), "Threshold": float(threshold),
                    "DataDir": str(dataset_directory(data_dir)), "ReplayFile": str(replay_file or "")},
        "LookBackWindow": 5000, "KnnHistoryLength": 100, "HistogramHistoryLength": 200,
        "StartupExcludedFromDiagnostics": 2, "Status": "running",
    }
    outputs = {}
    save_json(directory / "metadata.json", metadata)
    try:
        if replay_file:
            for name, parameters, streams in load_benchmark_replay(replay_file):
                if name in outputs:
                    raise ValueError("Replay archive contains duplicate dataset identifiers")
                result = _replay_result(streams, parameters, epsilon, threshold)
                outputs[name] = result
                save_benchmark_archive(directory / f"{name}.npz", {name: result}, metadata)
        else:
            metadata["DatasetsVerified"] = verify_datasets(data_dir)
            # MATLAB also interleaves preparation with each experiment. Each
            # engine reseeds this same generator, leaving its final state for
            # the following dataset's preparation.
            rng = np.random.RandomState(42)
            for name in DATASETS:
                x = load_dataset(name, data_dir=data_dir)
                if name in ("particle", "sensorless"):
                    x += 1e-3 * rng.standard_normal(x.shape)
                result = experiment.run_experiment(x, reps=int(reps), seed=42, TR_size=4096,
                                                  nu=256, K=50, sigma_delta=.1,
                                                  epsilon=epsilon, threshold=threshold,
                                                  cov_reg=1e-3, ncm_type="ellipsoid_pdf", rng=rng)
                outputs[name] = result
                save_benchmark_archive(directory / f"{name}.npz", {name: result}, metadata)
        tables = summarize_results(outputs)
        metadata.update(Status="complete", ElapsedSeconds=time.perf_counter() - started)
        save_benchmark_archive(directory / "benchmark.npz", outputs, metadata)
        for name, rows in tables.items():
            save_csv(directory / f"{name}.csv", rows)
        save_json(directory / "metadata.json", metadata)
        _write_report(directory / "benchmark.md", tables, metadata)
    except Exception as error:
        metadata.update(Status="failed", Error=repr(error), ElapsedSeconds=time.perf_counter() - started)
        save_json(directory / "metadata.json", metadata)
        raise
    _print_table(tables["Performance"])
    return outputs, tables, directory
