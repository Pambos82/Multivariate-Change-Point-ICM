"""Load the four numeric benchmark inputs without MATLAB or pickle."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATASETS = ("protein", "particle", "sensorless", "credit")
FILES = {
    "protein": "CASP.npz", "particle": "MiniBooNEPID.npz",
    "sensorless": "Sensorlessdrivediagnosis.npz", "credit": "creditcard_X.npz",
}


def dataset_directory(data_dir=None) -> Path:
    """Explicit directory, then COPA_DATA_DIR, then the bundled inputs."""
    value = data_dir or os.environ.get("COPA_DATA_DIR")
    return Path(value).expanduser().resolve() if value else ROOT / "datasets"


def file_sha256(path) -> str:
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def matrix_sha256(x) -> str:
    """Digest little-endian, MATLAB column-major float64 values independently of ZIP bytes."""
    return hashlib.sha256(np.asarray(x, dtype="<f8", order="F").tobytes(order="F")).hexdigest()


def load_dataset(name: str, prepared: bool = False, data_dir=None) -> np.ndarray:
    """Load a raw feature matrix; experiment-time jitter is applied by runners."""
    name = str(name).lower()
    if name not in FILES:
        raise ValueError(f"Unknown dataset {name!r}; choose from {DATASETS}")
    if prepared:
        raise ValueError("Frozen jittered inputs are not distributed; experiment runners prepare jitter")
    file = dataset_directory(data_dir) / FILES[name]
    with np.load(file, allow_pickle=False) as stored:
        if "X" not in stored:
            raise ValueError(f"Dataset archive must contain X: {file}")
        x = np.array(stored["X"], dtype=np.float64, order="C", copy=True)
    if x.ndim != 2 or min(x.shape) == 0 or not np.isfinite(x).all():
        raise ValueError(f"Dataset must be a nonempty finite numeric matrix: {file}")
    return x


def verify_datasets(data_dir=None) -> list[dict]:
    """Validate hashes, schemas and optional value digests against the manifest."""
    directory = dataset_directory(data_dir)
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8-sig"))
    specs = manifest["files"]
    if len(specs) != 4 or {s["dataset"] for s in specs} != set(DATASETS):
        raise ValueError("Dataset manifest must describe exactly the four benchmark inputs")
    checked = []
    for spec in specs:
        file = directory / spec["file"]
        if spec["file"] != FILES[spec["dataset"]]:
            raise ValueError("Dataset filename does not match its declared identifier")
        digest = file_sha256(file)
        if digest != spec["sha256"]:
            raise ValueError(f"Dataset hash mismatch: {file.name}")
        x = load_dataset(spec["dataset"], data_dir=directory)
        if list(x.shape) != spec["shape"] or str(x.dtype) != spec["dtype"]:
            raise ValueError(f"Dataset schema mismatch: {file.name}")
        values_digest = matrix_sha256(x)
        if "matrix_sha256" in spec and values_digest != spec["matrix_sha256"]:
            raise ValueError(f"Dataset value digest mismatch: {file.name}")
        checked.append({"dataset": spec["dataset"], "sha256": digest,
                        "matrix_sha256": values_digest, "shape": list(x.shape)})
    return checked
