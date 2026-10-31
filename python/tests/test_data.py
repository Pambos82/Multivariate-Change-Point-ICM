from pathlib import Path
import json

import numpy as np
import pytest

from copa2026.data import DATASETS, FILES, ROOT, load_dataset, matrix_sha256, verify_datasets


def test_bundled_datasets_have_verified_values_and_matlab_provenance():
    checked = verify_datasets()
    assert [row["dataset"] for row in checked] == list(DATASETS)
    specs = json.loads((ROOT / "datasets" / "manifest.json").read_text())["files"]
    assert [row["shape"] for row in checked] == [[45730, 9], [130064, 50], [58509, 48], [284807, 29]]
    for spec, actual in zip(specs, checked):
        assert actual["matrix_sha256"] == spec["matrix_sha256"]
        assert len(spec["source_mat_sha256"]) == 64


def test_matrix_digest_matches_matlab_column_major_byte_order():
    import hashlib
    x = np.array([[1., 2.], [3., 4.]])
    expected = hashlib.sha256(np.array([1., 3., 2., 4.], dtype="<f8").tobytes()).hexdigest()
    assert matrix_sha256(x) == expected


def test_data_directory_override_and_safe_loading(tmp_path, monkeypatch):
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    np.savez(first / FILES["protein"], X=np.arange(12.).reshape(4, 3))
    np.savez(second / FILES["protein"], X=np.ones((5, 2)))
    monkeypatch.setenv("COPA_DATA_DIR", str(first))
    assert load_dataset("protein").shape == (4, 3)
    assert load_dataset("protein", data_dir=second).shape == (5, 2)
    changed = load_dataset("protein", data_dir=second)
    changed[0, 0] = 9
    assert load_dataset("protein", data_dir=second)[0, 0] == 1


def test_nonfinite_data_and_frozen_input_mode_are_rejected(tmp_path):
    np.savez(tmp_path / FILES["protein"], X=np.array([[np.nan]]))
    with pytest.raises(ValueError, match="finite"):
        load_dataset("protein", data_dir=tmp_path)
    with pytest.raises(ValueError, match="Frozen"):
        load_dataset("protein", prepared=True)
