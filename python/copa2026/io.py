"""Unique result directories and portable, pickle-free experiment archives."""
from __future__ import annotations

import csv
import json
import uuid
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from .data import ROOT


def output_directory(name: str, output_root=None) -> Path:
    if not name.replace("_", "").replace("-", "").isalnum():
        raise ValueError("Output name must contain only letters, digits, '_' and '-'")
    parent = Path(output_root) if output_root is not None else ROOT / "results"
    directory = parent / f"{name}_{datetime.now(timezone.utc):%Y%m%d_%H%M%S_%f}_{uuid.uuid4().hex[:8]}"
    directory.mkdir(parents=True, exist_ok=False)
    return directory


def json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def save_json(file: Path, value: object) -> None:
    Path(file).write_text(json.dumps(value, indent=2, default=json_value) + "\n", encoding="utf-8")


def save_csv(file: Path, rows: list[dict], fieldnames=None) -> None:
    if fieldnames is None:
        if not rows:
            raise ValueError("Cannot infer columns from an empty table")
        fieldnames = list(rows[0])
    with Path(file).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


_ARCHIVE_ARRAYS = ("point", "fa", "ncms", "pvalues", "wealth", "logwealth", "Z",
                   "selected", "candidateValues", "candidateLog", "fallback", "statusCode", "rawKnn")


def save_benchmark_archive(path, outputs, metadata=None):
    """Preserve exact p-values, candidate histories and settings without pickle."""
    archive = {"ArchiveVersion": np.array(1, dtype=np.int64),
               "Datasets": np.asarray(list(outputs), dtype="U"),
               "MetadataJSON": np.asarray(json.dumps(metadata or {}, default=json_value))}
    for name, result in outputs.items():
        if not name or not name.replace("_", "").replace("-", "").isalnum():
            raise ValueError("Invalid dataset identifier in archive")
        archive[name + "__ParametersJSON"] = np.asarray(json.dumps(result["params"], default=json_value))
        for key in _ARCHIVE_ARRAYS:
            if key in result:
                value = np.asarray(result[key])
                if value.dtype.hasobject:
                    raise TypeError("Object arrays are not accepted in experiment archives")
                archive[name + "__" + key] = value
    with Path(path).open("xb") as stream:
        np.savez_compressed(stream, **archive)


def load_benchmark_replay(path):
    """Read p-value streams and parameters from Python or MATLAB benchmarks."""
    path = Path(path)
    if path.suffix.lower() == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            if int(archive["ArchiveVersion"]) != 1:
                raise ValueError("Unsupported experiment archive version")
            items = []
            for value in archive["Datasets"]:
                name = str(value)
                params = json.loads(str(archive[name + "__ParametersJSON"]))
                items.append((name, params, np.array(archive[name + "__pvalues"], dtype=float, copy=True)))
            return items
    if path.suffix.lower() != ".mat":
        raise ValueError("ReplayFile must be a benchmark .npz or MATLAB .mat archive")
    try:
        from scipy.io import loadmat
        saved = loadmat(path, simplify_cells=True)
    except NotImplementedError:
        return _read_matlab_hdf5(path)
    if "All" not in saved:
        raise ValueError("MATLAB replay archive must contain All")
    values = saved["All"]
    items = [values] if isinstance(values, dict) else np.asarray(values, dtype=object).reshape(-1)
    return [(str(item["Dataset"]), item["Parameters"], _stream_matrix(item["PValues"])) for item in items]


def _stream_matrix(values):
    values = np.asarray(values, dtype=float)
    return values[None, :] if values.ndim == 1 else values


def _read_matlab_hdf5(path):
    try:
        import h5py
    except ImportError as exc:
        raise ImportError("Reading MATLAB v7.3 replay files requires h5py") from exc

    with h5py.File(path, "r") as store:
        if "All" not in store:
            raise ValueError("MATLAB replay archive must contain All")

        def read(node):
            if isinstance(node, h5py.Group):
                return {key: read(value) for key, value in node.items()}
            values = np.asarray(node)
            if h5py.check_dtype(ref=node.dtype) is not None:
                decoded = [read(store[ref]) for ref in values.reshape(-1, order="F")]
                return decoded[0] if len(decoded) == 1 else decoded
            cls = node.attrs.get("MATLAB_class", b"")
            if isinstance(cls, bytes):
                cls = cls.decode("ascii")
            if cls == "char":
                return "".join(chr(int(c)) for c in values.reshape(-1, order="F") if c)
            if values.ndim > 1:
                values = values.transpose(tuple(reversed(range(values.ndim))))
            return values.item() if values.size == 1 else values

        group = store["All"]
        fields = {key: read(group[key]) for key in ("Dataset", "Parameters", "PValues")}
        dataset_values = fields["Dataset"]
        count = len(dataset_values) if isinstance(dataset_values, list) else 1
        result = []
        for index in range(count):
            item = {key: (values[index] if isinstance(values, list) else values) for key, values in fields.items()}
            result.append((str(item["Dataset"]), item["Parameters"], _stream_matrix(item["PValues"])))
        return result
