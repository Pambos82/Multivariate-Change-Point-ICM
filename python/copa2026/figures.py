"""The MATLAB package's sensitivity-sweep, figure and Air Quality workflows."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import date, datetime, time as daytime, timezone
from pathlib import Path
from uuid import uuid4

for _thread_variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_thread_variable] = "1"

import numpy as np
import scipy

from .betting import combination
from .data import FILES, load_dataset
from .experiment import conformal_pvalues, ncms_from_arrays, run_experiment, summarize
from .io import save_csv, save_json

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DATASET_NAMES = ("protein", "particle", "sensorless", "credit")
PANEL_NAMES = ("protein", "particle", "credit", "sensorless")
DEFAULT_LAMBDAS = (1e-5, 1e-4, 1e-3, 1e-2, 1e-1)
_WORKER_CACHE = {}
_SOURCES = ("figures.py", "betting.py", "_fast_betting.py", "experiment.py", "data.py", "benchmarks.py", "io.py")


def _hash(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=True,
                      default=lambda item: item.tolist() if isinstance(item, np.ndarray) else item.item())


def _utc():
    return datetime.now(timezone.utc).isoformat()


def _output_path(path):
    return Path(path).expanduser().resolve()


def new_output_dir(prefix, parent=None):
    parent = _output_path(PACKAGE_ROOT / "results" if parent is None else parent)
    parent.mkdir(parents=True, exist_ok=True)
    path = parent / f"{prefix}_{datetime.now(timezone.utc):%Y%m%d_%H%M%S_%f}_{uuid4().hex[:8]}"
    path.mkdir()
    return path


def _atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + f".{uuid4().hex}.tmp")
    save_json(temporary, value)
    # Windows readers/antivirus can briefly deny replacement of an open file.
    # Write once; retry only this atomic replacement. No random delays.
    for attempt in range(50):
        try:
            temporary.replace(path)
            return
        except PermissionError as error:
            if getattr(error, "winerror", None) not in (5, 32, 33) or attempt == 49:
                raise
            time.sleep(0.05)


def _source_hashes():
    directory = Path(__file__).resolve().parent
    return {name: _hash(directory / name) for name in _SOURCES}


def _deduplicate(values):
    return tuple(dict.fromkeys(values))


def _positive_integer(value, name):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value) or not np.isfinite(value) or value < 1 or int(value) != value:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _jobs(datasets, ks, lambdas, include_distance, epsilon):
    jobs = []
    for dataset in datasets:
        for k in ks:
            settings = [(value, "ellipsoid_pdf") for value in lambdas]
            if include_distance:
                settings.append((1e-3, "knn_distance"))
            for value, ncm in settings:
                jobs.append(dict(Number=len(jobs)+1, Dataset=dataset, K=k,
                                 Lambda=value, NCM=ncm, Epsilon=epsilon))
    return jobs


def _checkpoint_path(output, job):
    return Path(output) / "jobs" / f"job_{job['Number']:04d}.npz"


def _load_checkpoint(path, job, metadata):
    with np.load(path, allow_pickle=False) as stored:
        expected_keys = {"PValues", "points", "fa", "metadata_json"}
        if not expected_keys.issubset(stored.files):
            raise ValueError(f"Incomplete checkpoint: {path}")
        contents = json.loads(str(stored["metadata_json"].item()))
        pvalues = np.array(stored["PValues"], copy=True)
        points = np.array(stored["points"], copy=True)
        fa = np.array(stored["fa"], copy=True)
    reps = metadata["Scientific"]["Reps"]
    if contents.get("SchemaVersion") != 1 or contents.get("Signature") != metadata["Signature"] or contents.get("Job") != job:
        raise ValueError(f"Checkpoint signature/job differs: {path}")
    if pvalues.shape != (reps, 512) or not np.all(np.isfinite(pvalues)) or np.any((pvalues < 0) | (pvalues > 1)):
        raise ValueError(f"Invalid checkpoint p-values: {path}")
    if points.shape != (reps,) or fa.shape != (reps,) or fa.dtype != np.bool_ or np.any(points != np.floor(points)) or np.any((points < 0) | (points > 256)):
        raise ValueError(f"Invalid checkpoint decisions: {path}")
    summary = summarize(points, fa)
    for key, value in summary.items():
        recorded = contents["Performance"][key]
        if not (value == recorded or (math.isnan(value) and math.isnan(recorded))):
            raise ValueError(f"Checkpoint summaries disagree with decisions: {path}")
    tables = contents["Tables"]
    if (len(tables["NormalizationConstants"]) != 3 or len(tables["SelectionFrequencies"]) != 24
            or len(tables["PerRepetition"]) != reps or len(tables["NumericalIssues"]) != 1):
        raise ValueError(f"Incomplete checkpoint diagnostic tables: {path}")
    if any(row["DecisionCount"] != reps * 510 for row in tables["NormalizationConstants"]):
        raise ValueError(f"Checkpoint diagnostic counts disagree: {path}")
    parameters = contents["Parameters"]
    expected_parameters = dict(reps=reps, K=job["K"], cov_reg=job["Lambda"], epsilon=job["Epsilon"],
                               threshold=100, TR_size=4096, nu=256, sigma_delta=.1, seed=42,
                               ncm_type=job["NCM"], standardize=True, change_type="mean_shift")
    if any(parameters.get(name) != value for name, value in expected_parameters.items()):
        raise ValueError(f"Checkpoint engine parameters disagree: {path}")
    return contents


def _execute_job(output, job, metadata):
    """Top-level process target; prepared matrices are cached per process."""
    from .benchmarks import summarize_results
    started = time.perf_counter()
    if _source_hashes() != metadata["Sources"]:
        raise RuntimeError("Scientific source changed during the sweep")
    output = Path(output)
    destination = _checkpoint_path(output, job)
    if destination.exists():
        return _load_checkpoint(destination, job, metadata)
    cache_key = (str(output), job["Dataset"], metadata["Signature"])
    if cache_key not in _WORKER_CACHE:
        prepared_path = output / "prepared" / f"{job['Dataset']}.npy"
        if _hash(prepared_path) != metadata["PreparedFiles"][prepared_path.name]:
            raise ValueError("Prepared data changed during the sweep")
        _WORKER_CACHE[cache_key] = np.load(prepared_path, mmap_mode="r", allow_pickle=False)
    result = run_experiment(_WORKER_CACHE[cache_key], reps=metadata["Scientific"]["Reps"],
                            seed=42, TR_size=4096, nu=256, K=job["K"], sigma_delta=0.1,
                            epsilon=job["Epsilon"], threshold=100.0,
                            cov_reg=job["Lambda"], ncm_type=job["NCM"])
    tables = summarize_results({job["Dataset"]: result})
    performance = dict(Dataset=job["Dataset"], Implementation="normalized", K=job["K"],
                       Lambda=job["Lambda"], NCM=job["NCM"], Epsilon=job["Epsilon"], **result["summary"])
    tables.pop("Performance", None)
    contents = dict(SchemaVersion=1, Signature=metadata["Signature"], Job=job,
                    CompletedUTC=_utc(), ElapsedSeconds=time.perf_counter()-started,
                    Performance=performance, Tables=tables, Parameters=result["params"])
    temporary = destination.with_name(destination.stem + f"_{uuid4().hex}.inprogress.npz")
    with temporary.open("xb") as stream:
        np.savez_compressed(stream, PValues=result["pvalues"], points=result["points"], fa=result["fa"],
                            metadata_json=np.asarray(_json(contents)))
    _load_checkpoint(temporary, job, metadata)
    if destination.exists():
        raise FileExistsError(f"Completed checkpoint appeared during computation: {destination}")
    temporary.rename(destination)
    return contents


def _aggregate(output, jobs, completed):
    tables = {name: [] for name in ("Performance", "SelectionFrequencies", "NormalizationConstants", "NumericalIssues", "JobIndex")}
    for job in jobs:
        if job["Number"] not in completed:
            continue
        item = completed[job["Number"]]
        tables["Performance"].append(item["Performance"])
        for name in ("SelectionFrequencies", "NormalizationConstants", "NumericalIssues"):
            tables[name].extend(dict(row, Job=job["Number"], NCMK=job["K"], Lambda=job["Lambda"],
                                     NCM=job["NCM"], Epsilon=job["Epsilon"]) for row in item["Tables"][name])
        tables["JobIndex"].append(dict(Job=job["Number"], Dataset=job["Dataset"], K=job["K"],
                                       Lambda=job["Lambda"], NCM=job["NCM"], Epsilon=job["Epsilon"],
                                       ElapsedSeconds=item["ElapsedSeconds"], File=f"jobs/job_{job['Number']:04d}.npz"))
    for name, rows in tables.items():
        if rows:
            temporary = Path(output) / f"{name}.csv.tmp"
            save_csv(temporary, rows)
            temporary.replace(Path(output) / f"{name}.csv")
    return tables


def _progress(output, jobs, completed, failed, status, workers):
    _atomic_json(Path(output) / "progress.json", dict(Status=status, CompletedJobs=len(completed),
                 TotalJobs=len(jobs), FailedJobs=list(failed), ConfiguredWorkers=workers, UpdatedUTC=_utc()))


def run_full_figure_sweeps(*, data_dir=None, dataset_names=DATASET_NAMES, k_list=range(1, 101),
                           lambda_list=DEFAULT_LAMBDAS, include_distance=True, reps=100,
                           epsilon=100.0, workers=16, output_root=None, max_jobs=math.inf):
    """Compute/resume the joint 2,400-setting Figure 1/2 grid; plotting is separate.

    Workers and max_jobs can change on resume. Scientific arguments, source,
    input bytes, and prepared cache bytes must match the immutable metadata.
    Native NumPy random streams are recorded rather than labeled MATLAB-identical.
    """
    reps, workers = _positive_integer(reps, "reps"), _positive_integer(workers, "workers")
    selected = tuple(dataset_names)
    if not selected or any(name not in DATASET_NAMES for name in selected):
        raise ValueError("dataset_names must select protein, particle, sensorless or credit")
    datasets = tuple(name for name in DATASET_NAMES if name in selected)
    ks = _deduplicate(_positive_integer(k, "k_list") for k in k_list)
    lambdas = _deduplicate(float(value) for value in lambda_list)
    if not ks or any(k > 4096 for k in ks) or not lambdas or any(not np.isfinite(v) or v <= 0 for v in lambdas):
        raise ValueError("Require K in 1:4096 and positive finite lambda values")
    if not np.isfinite(epsilon) or epsilon <= 0:
        raise ValueError("epsilon must be positive and finite")
    if not (max_jobs == math.inf or (np.isfinite(max_jobs) and max_jobs >= 0 and int(max_jobs) == max_jobs)):
        raise ValueError("max_jobs must be a nonnegative integer or positive infinity")
    data_directory = Path(data_dir or os.environ.get("COPA_DATA_DIR") or PACKAGE_ROOT / "datasets").resolve()
    scientific = dict(DataDir=str(data_directory), DatasetNames=list(datasets), KList=list(ks),
                      LambdaList=list(lambdas), IncludeDistance=bool(include_distance), Reps=reps,
                      Epsilon=float(epsilon), Seed=42, Threshold=100, TRSize=4096, Nu=256,
                      SigmaDelta=0.1, LookBackWindow=5000, ExperimentGenerator="NumPy RandomState/MT19937",
                      PreprocessingGenerator="NumPy RandomState/MT19937", PreprocessingOrder=list(DATASET_NAMES),
                      WorkerMaxNumCompThreads=1, K1CovarianceBehavior="Preserve scalar covariance convention",
                      PythonVersion=platform.python_version(), NumPyVersion=np.__version__, SciPyVersion=scipy.__version__)
    inputs = {FILES[name]: _hash(data_directory / FILES[name]) for name in DATASET_NAMES}
    sources = _source_hashes()
    signature = hashlib.sha256(_json(dict(Scientific=scientific, Sources=sources, Inputs=inputs)).encode()).hexdigest()
    jobs = _jobs(datasets, ks, lambdas, bool(include_distance), float(epsilon))
    output = (new_output_dir("full", PACKAGE_ROOT / "results" / "figure_sweeps") if output_root is None
              else _output_path(output_root))
    metadata_path = output / "metadata.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("Signature") != signature or metadata.get("Jobs") != jobs:
            raise ValueError("Resume rejected: scientific options, runtime version, sources or inputs changed")
        actual_prepared = {name: _hash(output / "prepared" / name) for name in metadata["PreparedFiles"]}
        if actual_prepared != metadata["PreparedFiles"]:
            raise ValueError("Resume rejected: prepared data changed")
    else:
        if output.exists() and any(output.iterdir()):
            raise ValueError("OutputRoot is nonempty without metadata; choose a fresh directory")
        for child in ("jobs", "prepared", "errors"):
            (output / child).mkdir(parents=True, exist_ok=True)
        prep_rng = np.random.RandomState(42)
        k1_rows = []
        for name in DATASET_NAMES:
            x = load_dataset(name, data_dir=data_directory)
            if name in ("particle", "sensorless"):
                x = x + 1e-3 * prep_rng.standard_normal(x.shape)
            np.save(output / "prepared" / f"{name}.npy", x, allow_pickle=False)
            k1_rows.append(dict(Dataset=name, Dimension=x.shape[1], CovarianceIsScalar=True,
                                CovarianceValue=float(np.var(x[0], ddof=0))))
        metadata = dict(SchemaVersion=1, CreatedUTC=_utc(), Scientific=scientific, Signature=signature,
                        Sources=sources, Inputs=inputs, Jobs=jobs, TotalJobs=len(jobs),
                        MatlabStreamEquivalence=False,
                        PreparedFiles={f"{name}.npy": _hash(output / "prepared" / f"{name}.npy") for name in DATASET_NAMES})
        _atomic_json(metadata_path, metadata)
        save_csv(output / "PlannedJobs.csv", jobs)
        save_csv(output / "K1CovarianceDiagnostic.csv", k1_rows)
    completed, failed = {}, []
    for job in jobs:
        path = _checkpoint_path(output, job)
        if path.exists():
            completed[job["Number"]] = _load_checkpoint(path, job, metadata)
    pending = [job for job in jobs if job["Number"] not in completed]
    scheduled = pending if max_jobs == math.inf else pending[:int(max_jobs)]
    _progress(output, jobs, completed, failed, "running", workers)
    try:
        if workers == 1:
            for job in scheduled:
                try:
                    completed[job["Number"]] = _execute_job(str(output), job, metadata)
                except BaseException as error:
                    failed.append(job["Number"])
                    _atomic_json(output / "errors" / f"job_{job['Number']:04d}_{uuid4().hex[:8]}.json",
                                 dict(Job=job, Error=repr(error), CreatedUTC=_utc()))
                    raise
                _progress(output, jobs, completed, failed, "running", workers)
        elif scheduled:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(_execute_job, str(output), job, metadata): job for job in scheduled}
                for future in as_completed(futures):
                    job = futures[future]
                    try:
                        completed[job["Number"]] = future.result()
                    except BaseException as error:
                        failed.append(job["Number"])
                        _atomic_json(output / "errors" / f"job_{job['Number']:04d}_{uuid4().hex[:8]}.json",
                                     dict(Job=job, Error=repr(error), CreatedUTC=_utc()))
                        for pending_future in futures:
                            pending_future.cancel()
                        if hasattr(pool, "terminate_workers"):
                            pool.terminate_workers()
                        raise
                    _progress(output, jobs, completed, failed, "running", workers)
    except BaseException:
        _aggregate(output, jobs, completed)
        _progress(output, jobs, completed, failed, "failed", workers)
        raise
    tables = _aggregate(output, jobs, completed)
    status = "complete" if len(completed) == len(jobs) else "partial"
    _progress(output, jobs, completed, failed, status, workers)
    _atomic_json(output / "summary.json", dict(Signature=signature, Status=status, Tables=tables))
    return tables["Performance"], output


def _validated_plot_rows(performance_input):
    source = Path(performance_input).resolve()
    if source.is_dir():
        source = source / "Performance.csv"
    with source.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        required = {"Dataset", "Implementation", "K", "Lambda", "NCM", "Epsilon", "DetectionRate", "FalseAlarmRate"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Performance.csv is missing required columns")
        rows = [dict(row) for row in reader if row["Implementation"].lower() == "normalized"]
    if len(rows) != 2400:
        raise ValueError(f"Expected exactly 2400 rows; found {len(rows)}. Incomplete grids cannot be plotted.")
    curves = {}
    for row in rows:
        row["Dataset"], row["NCM"] = row["Dataset"].lower(), row["NCM"].lower()
        for field in ("K", "Lambda", "Epsilon", "DetectionRate", "FalseAlarmRate"):
            row[field] = float(row[field])
        if row["Epsilon"] != 100 or row["Dataset"] not in PANEL_NAMES or row["NCM"] not in ("ellipsoid_pdf", "knn_distance"):
            raise ValueError("Unexpected dataset, NCM or epsilon in plotting grid")
        if not np.isfinite(row["K"]) or row["K"] != int(row["K"]) or not 1 <= row["K"] <= 100:
            raise ValueError("K must be an integer in 1:100")
        row["K"] = int(row["K"])
        if any(not np.isfinite(row[field]) or not 0 <= row[field] <= 1 for field in ("DetectionRate", "FalseAlarmRate")):
            raise ValueError("Plot rates must be finite probabilities")
        if row["NCM"] == "ellipsoid_pdf":
            matches = [value for value in DEFAULT_LAMBDAS if abs(row["Lambda"]-value) <= 8*np.spacing(value)]
            if len(matches) != 1:
                raise ValueError("Unexpected ellipsoid lambda")
            value = matches[0]
        else:
            value = None
        key = (row["Dataset"], row["NCM"], value, row["K"])
        if key in curves:
            raise ValueError("Duplicate grid setting")
        curves[key] = (row["DetectionRate"], row["FalseAlarmRate"])
    for dataset in PANEL_NAMES:
        for k in range(1, 101):
            for value in DEFAULT_LAMBDAS:
                if (dataset, "ellipsoid_pdf", value, k) not in curves:
                    raise ValueError("Missing ellipsoid grid setting")
            if (dataset, "knn_distance", None, k) not in curves:
                raise ValueError("Missing distance grid setting")
    return source, rows, curves


def reconstruct_copa_figures(performance_input, *, export_pdf=False):
    """Render complete Figures 1 and 2 from supplied rows, without recomputation."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    source, rows, curves = _validated_plot_rows(performance_input)
    output = new_output_dir("figures", source.parent)
    save_csv(output / "PerformanceUsed.csv", rows)
    colors = [(0, .447, .741), (.85, .325, .098), (.929, .694, .125), (.494, .184, .556), (.466, .674, .188)]
    markers = ["o", "s", "d", "^", "v"]
    maximum = max(row["FalseAlarmRate"] for row in rows)
    fa_max = .05 if maximum <= .05 else math.ceil(100*(maximum+.005))/100
    ks = np.arange(1, 101)
    records = []
    for number in (1, 2):
        fig = plt.figure(figsize=(15, 16.8), facecolor="white")
        for d, dataset in enumerate(PANEL_NAMES):
            for metric in range(2):
                left, bottom = .070+metric*.490, .800-d*.235
                ax = fig.add_axes([left, bottom, .405, .158])
                settings = [("ellipsoid_pdf", .001), ("knn_distance", None)] if number == 1 else [("ellipsoid_pdf", value) for value in DEFAULT_LAMBDAS]
                labels = ["Ellipsoid NCM", "kNN-Euclidean NCM"] if number == 1 else [rf"$\lambda = 10^{{{round(math.log10(value))}}}$" for value in DEFAULT_LAMBDAS]
                for index, ((ncm, value), label) in enumerate(zip(settings, labels)):
                    values = [curves[(dataset, ncm, value, int(k))][metric] for k in ks]
                    ax.plot(ks, values, "-", color=colors[index], marker=markers[index], markersize=3,
                            markerfacecolor="none", linewidth=.9, label=label)
                metric_name = "Detection rate" if metric == 0 else "False alarm rate"
                ax.set(xlim=(0, 100), ylim=(0, 1.05 if metric == 0 else fa_max),
                       xlabel="Number of nearest neighbours K", ylabel=metric_name,
                       title=dataset.capitalize() if number == 1 else f"{metric_name} vs K: {dataset}")
                ax.set_xticks(np.arange(0, 101, 10))
                if metric == 0:
                    ax.set_yticks(np.arange(0, 1.01, .1))
                elif fa_max == .05:
                    ax.set_yticks(np.arange(0, .0501, .005))
                ax.tick_params(labelsize=9, direction="in")
                ax.grid(alpha=.15)
                ax.legend(fontsize=8, loc="best", frameon=True)
                letter = chr(ord("a")+d*2+metric)
                fig.text(left+.2025, bottom-.046, f"({letter}) {dataset.capitalize()}: {metric_name.lower()}",
                         ha="center", va="center", fontsize=14, family="serif")
                records.append(dict(Figure=number, Dataset=dataset, Metric=metric_name,
                                    CurveCount=len(settings), PointsPerCurve=100, XLim=[0, 100],
                                    YLim=[0, 1.05 if metric == 0 else fa_max]))
        caption = ("Figure 1. Detection and false alarm rates for kNN betting (cautious epsilon = 100).\n"
                   "Ellipsoid lambda = 0.001; K = 1 retains the original covariance convention."
                   if number == 1 else
                   "Figure 2. Ellipsoid NCM sensitivity with kNN betting (cautious epsilon = 100).\n"
                   "lambda = 0.00001, 0.0001, 0.001, 0.01, 0.1; K = 1 retains the original covariance convention.")
        fig.text(.07, .016, caption, va="center", fontsize=11, family="serif")
        fig.savefig(output / f"figure{number}.png", dpi=220)
        if export_pdf:
            fig.savefig(output / f"figure{number}.pdf")
        plt.close(fig)
    artifacts = ["figure1.png", "figure2.png", "PerformanceUsed.csv"]
    if export_pdf:
        artifacts += ["figure1.pdf", "figure2.pdf"]
    save_json(output / "reconstruction_metadata.json", dict(CreatedUTC=_utc(), SourcePerformanceCSV=str(source),
              SourceSHA256=_hash(source), FullGridValidated=True, RowsUsed=2400,
              ExperimentsRecomputedByRenderer=False, K=list(range(1,101)), EllipsoidLambda=list(DEFAULT_LAMBDAS),
              Figure1EllipsoidLambda=.001, CautiousEpsilon=100, FalseAlarmAxisUpperLimit=fa_max,
              DatasetsInPanelOrder=list(PANEL_NAMES), PlotRecords=records, ArtifactNames=artifacts,
              EditableData="PerformanceUsed.csv plus this reconstruction function"))
    return output


def prepare_airquality_input(xlsx_path, output_dir):
    """Prepare the original XLSX locally, retaining date-only week indexing."""
    from openpyxl import load_workbook

    source, destination = Path(xlsx_path).resolve(), Path(output_dir).resolve() / "airquality.npz"
    if source.suffix.lower() != ".xlsx" or not source.is_file():
        raise ValueError("Use the original AirQualityUCI.xlsx, not the rounded CSV")
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite {destination}")
    workbook = load_workbook(source, read_only=True, data_only=True)
    try:
        rows = list(workbook.active.iter_rows(max_col=15, values_only=True))
    finally:
        workbook.close()
    while rows and all(value is None for value in rows[-1]):
        rows.pop()
    expected = ["Date", "Time", "CO(GT)", "PT08.S1(CO)", "NMHC(GT)", "C6H6(GT)", "PT08.S2(NMHC)",
                "NOx(GT)", "PT08.S3(NOx)", "NO2(GT)", "PT08.S4(NO2)", "PT08.S5(O3)", "T", "RH", "AH"]
    if len(rows) != 9358 or list(rows[0]) != expected:
        raise ValueError("Expected the original 9357-row, 15-column Air Quality spreadsheet")
    dates, values = [], []
    for record in rows[1:]:
        row_date = record[0]
        if isinstance(row_date, datetime):
            if row_date.time() != daytime():
                raise ValueError("Expected date-only spreadsheet row dates")
            row_date = row_date.date()
        if not isinstance(row_date, date):
            raise ValueError("Expected original spreadsheet date cells")
        row_time = record[1]
        if isinstance(row_time, daytime):
            row_time = (row_time.hour*3600+row_time.minute*60+row_time.second+row_time.microsecond/1e6)/86400
        if not isinstance(row_time, (int, float)) or not np.isfinite(row_time) or not 0 <= row_time < 1:
            raise ValueError("Unexpected fractional-day Time value")
        dates.append(np.datetime64(row_date, "D"))
        values.append(record[2:15])
    numeric = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(numeric)):
        raise ValueError("Unexpected nonfinite raw spreadsheet values")
    x = np.delete(numeric, 2, axis=1)
    x[x == -200] = np.nan
    keep = np.all(np.isfinite(x), axis=1)
    x, dates = x[keep], np.asarray(dates)[keep]
    if x.shape != (6941, 12) or np.any(np.diff(dates).astype("timedelta64[D]").astype(int) < 0):
        raise ValueError("Unexpected cleaned Air Quality shape or date order")
    week = (dates-dates[0]).astype("timedelta64[D]").astype(np.int64)//7+1
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("xb") as stream:
        np.savez_compressed(stream, X=x, week=week, dates=dates.astype("U10"),
                            feature_names=np.asarray(expected[2:4]+expected[5:]), InputSHA256=np.asarray(_hash(source)))
    return destination


def run_real_stream(x, week, *, TR_size=3000, K=50, epsilon=100.0, threshold=100.0,
                    cov_reg=1e-3, seed=42, uniforms=None):
    """Compute the fixed-reference Air stream; optional draws support paired checks."""
    x = np.asarray(x, dtype=np.float64)
    week = np.asarray(week, dtype=np.float64).reshape(-1)
    if (x.ndim != 2 or len(x) <= TR_size or len(week) != len(x)
            or not np.all(np.isfinite(x)) or not np.all(np.isfinite(week))):
        raise ValueError("Require finite data and aligned weeks beyond the reference set")
    reference, queries = x[:TR_size].copy(), x[TR_size:].copy()
    mean, std = reference.mean(axis=0), reference.std(axis=0, ddof=1)
    std[std < 1e-12] = 1
    reference, queries = (reference-mean)/std, (queries-mean)/std
    scores = ncms_from_arrays(reference, queries, cov_reg=cov_reg, K=K)
    if not np.all(np.isfinite(scores)):
        raise ValueError("Air Quality nonconformity scores must be finite")
    pvalues = conformal_pvalues(scores, uniforms=uniforms, rng=np.random.RandomState(seed))
    wealth, logs, diagnostics = combination(pvalues, epsilon, threshold, return_diagnostics=True)
    if not np.all(np.isfinite(logs)) or np.any(np.isnan(wealth)):
        raise FloatingPointError("Unexpected nonfinite Air Quality log capital or NaN capital")
    hits = np.flatnonzero(logs > np.log(threshold))
    alarm = int(hits[0]+1) if hits.size else None
    monitored_week = week[TR_size:]
    return dict(M=wealth, logM=logs, pval=pvalues, ncms=scores, week_test=monitored_week,
                diagnostics=diagnostics, alarm_index=alarm, alarm_time=TR_size+alarm if alarm else None,
                alarm_week=float(monitored_week[alarm-1]) if alarm else None)


def _air_tables(stream):
    diagnostics = stream["diagnostics"]
    selection, constants = [], []
    candidate_names = ["cautious/no bet", *diagnostics["candidateNames"]]
    decisions = len(diagnostics["selected"])-2
    for code in range(8):
        count = int(np.count_nonzero(diagnostics["selected"][2:] == code))
        selection.append(dict(CandidateCode=code, Candidate=candidate_names[code], Count=count,
                              DecisionCount=decisions, Frequency=count/decisions if decisions else np.nan))
    for j, k in enumerate((5,10,15)):
        z = diagnostics["Z"][2:, j]
        finite = z[np.isfinite(z)]
        q = np.quantile(finite, [0,.01,.05,.5,.95,.99,1], method="hazen") if len(finite) else np.full(7,np.nan)
        constants.append(dict(K=k, DecisionCount=len(z), FiniteCount=len(finite),
                              FallbackCount=int(np.count_nonzero(diagnostics["fallback"][2:,j])),
                              InfiniteCount=int(np.count_nonzero(np.isinf(z))), NaNCount=int(np.count_nonzero(np.isnan(z))),
                              Mean=float(np.mean(finite)) if len(finite) else np.nan,
                              Std=float(np.std(finite,ddof=1)) if len(finite)>1 else (0. if len(finite) else np.nan),
                              **dict(zip(("Min","Q01","Q05","Median","Q95","Q99","Max"),map(float,q)))))
    return selection, constants


def run_real(*, data_dir=None, output_root=None):
    """Run the prepared Air Quality input and return (stream, summary, output_dir)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    directory = Path(data_dir or os.environ.get("COPA_DATA_DIR") or PACKAGE_ROOT / "datasets").resolve()
    source = directory / "airquality.npz"
    if not source.is_file():
        raise FileNotFoundError("Air Quality is excluded. Prepare the original XLSX and pass its local folder as data_dir.")
    with np.load(source, allow_pickle=False) as stored:
        x, week = np.array(stored["X"], copy=True), np.array(stored["week"], copy=True)
        dates = np.asarray(stored["dates"], dtype="datetime64[D]")
    if x.ndim != 2 or x.shape[1] != 12 or dates.shape != (len(x),) or np.any(np.isnat(dates)) or np.any(np.diff(dates).astype(int) < 0):
        raise ValueError("Expected 12 measurements and ordered date-only row dates")
    expected_week = (dates-dates[0]).astype("timedelta64[D]").astype(np.int64)//7+1
    if not np.array_equal(week,expected_week):
        raise ValueError("Week indexing must use cleaned date-only row dates")
    output = new_output_dir("run", output_root or PACKAGE_ROOT / "results" / "airquality")
    started = time.perf_counter()
    stream = run_real_stream(x, week)
    stream["rowTimes"] = dates[3000:].astype("U10")
    summary = dict(Implementation="normalized", Epsilon=100, MonitoredSamples=len(stream["pval"]),
                   FirstAlarmMonitoredIndex=stream["alarm_index"], FirstAlarmCleanDataIndex=stream["alarm_time"],
                   FirstAlarmWeek=stream["alarm_week"], FinalLogM=float(stream["logM"][-1]),
                   MaxLogM=float(np.max(stream["logM"])), LinearOverflowCount=int(np.isinf(stream["M"]).sum()),
                   LinearNaNCount=int(np.isnan(stream["M"]).sum()), FiniteLogM=bool(np.isfinite(stream["logM"]).all()),
                   MonitoringStartWeek=float(np.min(stream["week_test"])), MonitoringEndWeek=float(np.max(stream["week_test"])))
    selection, constants = _air_tables(stream)
    metadata = dict(CreatedUTC=_utc(), PythonVersion=platform.python_version(), NumPyVersion=np.__version__,
                    SciPyVersion=scipy.__version__, Seed=42,
                    Generator="NumPy RandomState/MT19937", MatlabStreamEquivalence=False, ReferenceSize=3000,
                    K=50, CovarianceDiagonalAdjustment=.001, Epsilon=100, W=5000, Threshold=100,
                    KnnHistoryLength=100, HistogramHistoryLength=200, CleanRows=len(x), Features=12,
                    InputSHA256=_hash(source), Sources=_source_hashes(), ElapsedSeconds=time.perf_counter()-started,
                    FirstCleanTimestamp=str(dates[0]), LastCleanTimestamp=str(dates[-1]), GroundTruthChangePoint=None)
    arrays = {name:value for name,value in stream.items() if isinstance(value,np.ndarray)}
    arrays.update({"diagnostics_"+key:value for key,value in stream["diagnostics"].items() if isinstance(value,np.ndarray)})
    np.savez_compressed(output / "airquality_results.npz", **arrays, metadata_json=np.asarray(_json(metadata)))
    save_csv(output / "Summary.csv", [summary])
    save_csv(output / "SelectionFrequencies.csv", selection)
    save_csv(output / "NormalizationConstants.csv", constants)
    save_json(output / "provenance.json", metadata)
    fig, ax = plt.subplots(figsize=(11,6), facecolor="white")
    ax.plot(stream["week_test"], stream["logM"], linewidth=1.5, color=(.12,.35,.65), label="Cautious martingale")
    ax.axhline(np.log(100), color="black", linestyle="--", label="Alarm threshold: log(100)")
    ax.axvline(33, color=".5", linestyle=":", label="Week 33 reference")
    if stream["alarm_index"] is not None:
        i=stream["alarm_index"]-1
        ax.scatter(stream["week_test"][i],stream["logM"][i],s=40,color=(.78,.24,.12),label="First threshold crossing")
    ax.set(xlim=(min(stream["week_test"]),max(stream["week_test"])), xlabel="Week since first cleaned observation",
           ylabel="Natural log martingale", title="Air Quality: cautious betting, epsilon = 100")
    ax.grid(True)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(output / "figure3.png",dpi=200)
    fig.savefig(output / "figure3.pdf")
    plt.close(fig)
    report=("# Air Quality experiment\n\nPython. Seed 42, reference 3000, K=50, covariance diagonal adjustment 0.001, "
            "epsilon=100, W=5000, alarm threshold=100.\n\n"
            f"The {len(x)} complete observations have 12 measurements; {len(stream['pval'])} are monitored. "
            "Figure 3 shows the complete monitoring interval using accumulated log capital. "
            f"Linear overflow count: {summary['LinearOverflowCount']}. First alarm: monitored index "
            f"{stream['alarm_index']}, clean-data index {stream['alarm_time']}, week {stream['alarm_week']}.\n\n"
            "The alarm uses the strict rule log(M)>log(100). Week 33 is a visual reference rather than an "
            "assumed ground-truth change point. No detection delay or false-alarm rate is inferred from this unlabeled stream.\n")
    (output / "airquality_results.md").write_text(report,encoding="utf8")
    return stream, summary, output
