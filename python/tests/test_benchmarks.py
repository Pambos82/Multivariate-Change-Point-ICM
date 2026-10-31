import json

import numpy as np
import pytest
from scipy.io import savemat

from copa2026.benchmarks import DISPLAY_COLUMNS, run_benchmarks, summarize_results
from copa2026.experiment import run_knn_cm_experiment
from copa2026.io import load_benchmark_replay, save_benchmark_archive


def small_result():
    streams = np.random.RandomState(12).uniform(size=(2, 32))
    return run_knn_cm_experiment(np.empty((0, 0)), reps=2, nu=16, pvalue_streams=streams)


def test_summary_definitions_and_matlab_quantile_convention():
    result = {
        "point": np.array([0, 2]), "fa": np.array([True, False]),
        "wealth": np.ones((2, 6)), "logwealth": np.zeros((2, 6)),
        "selected": np.array([[7, 7, 1, 2, 0, 7], [7, 7, 4, 0, 4, 4]]),
        "Z": np.full((2, 6, 3), np.nan), "fallback": np.zeros((2, 6, 3), dtype=bool),
        "params": {"nu": 3},
    }
    result["Z"][:, 2:, 0] = np.arange(1., 9.).reshape(2, 4)
    result["Z"][:, 2:, 1:] = 2.
    result["Z"][0, 3, 1] = np.inf
    result["fallback"][0, 3, 1] = True
    tables = summarize_results({"protein": result})
    p = tables["Performance"][0]
    assert p["DetectionRate"] == p["FalseAlarmRate"] == p["MissRate"] == .5
    assert p["MedianDelay"] == p["MeanDelay"] == 2
    assert p["StdDelay"] == 0
    k5 = tables["NormalizationConstants"][0]
    assert k5["Q25"] == 2.5
    assert k5["Median"] == 4.5
    k10 = tables["NormalizationConstants"][1]
    assert k10["InfiniteCount"] == k10["FallbackCount"] == 1
    choices = {r["CandidateCode"]: r for r in tables["SelectionFrequencies"] if r["Phase"] == "all"}
    assert choices[7]["Count"] == 1  # startup histogram choices are excluded
    assert choices[4]["FrequencyAllDecisions"] == 3 / 8
    assert choices[4]["FrequencyAmongActiveBets"] == .5


def test_replay_archive_roundtrip_and_quiet_console(tmp_path, capsys):
    result = small_result()
    source = tmp_path / "source.npz"
    save_benchmark_archive(source, {"protein": result})
    with np.load(source, allow_pickle=False) as a:
        assert a["protein__candidateValues"].shape == (2, 32, 7)
        assert a["protein__fallback"].shape == (2, 32, 3)
    outputs, tables, directory = run_benchmarks(replay_file=source, output_root=tmp_path / "runs")
    lines = capsys.readouterr().out.strip().splitlines()
    assert len(lines) == 2
    assert lines[0].split() == list(DISPLAY_COLUMNS)
    for key in ("point", "fa", "pvalues", "wealth", "logwealth", "selected", "candidateValues", "candidateLog"):
        np.testing.assert_equal(outputs["protein"][key], result[key])
    assert set(tables) == {"Performance", "NormalizationConstants", "SelectionFrequencies", "PerRepetition", "NumericalIssues"}
    assert (directory / "benchmark.npz").is_file()
    assert (directory / "benchmark.md").is_file()
    for name in tables:
        assert (directory / f"{name}.csv").is_file()
    _, _, next_directory = run_benchmarks(replay_file=directory / "benchmark.npz", output_root=tmp_path / "runs")
    assert next_directory != directory


def test_matlab_v7_pvalue_archive_can_be_replayed(tmp_path):
    result = small_result()
    path = tmp_path / "matlab.mat"
    savemat(path, {"All": {"Dataset": "protein", "Parameters": {"nu": 16, "seed": 42},
                           "PValues": result["pvalues"]}})
    items = load_benchmark_replay(path)
    assert len(items) == 1 and items[0][0] == "protein"
    np.testing.assert_equal(items[0][2], result["pvalues"])


def test_matlab_hdf5_decoder_reads_only_replay_fields(tmp_path):
    h5py = pytest.importorskip("h5py")
    from copa2026.io import _read_matlab_hdf5
    path = tmp_path / "sample_v73.mat"
    with h5py.File(path, "w") as f:
        refs = f.create_group("#refs#")
        name = refs.create_dataset("name", data=np.asarray([ord(c) for c in "protein"], dtype=np.uint16)[:, None])
        name.attrs["MATLAB_class"] = np.bytes_("char")
        params = refs.create_group("params")
        params.create_dataset("nu", data=np.array([[16.]]))
        streams = np.arange(64., dtype=float).reshape(2, 32) / 64
        p = refs.create_dataset("streams", data=streams.T)
        group = f.create_group("All")
        for label, value in (("Dataset", name), ("Parameters", params), ("PValues", p)):
            dataset = group.create_dataset(label, (1, 1), dtype=h5py.ref_dtype)
            dataset[0, 0] = value.ref
        # Unneeded MATLAB result objects may have arbitrary/cyclic references.
        unused = group.create_dataset("Normalized", (1, 1), dtype=h5py.ref_dtype)
        unused[0, 0] = group.ref
    loaded = _read_matlab_hdf5(path)
    assert loaded[0][0] == "protein"
    assert loaded[0][1]["nu"] == 16
    np.testing.assert_equal(loaded[0][2], streams)


def test_invalid_benchmark_settings_are_rejected(tmp_path):
    for kwargs in ({"reps": 0}, {"epsilon": 0}, {"threshold": np.inf}):
        with pytest.raises(ValueError):
            run_benchmarks(output_root=tmp_path, **kwargs)


@pytest.mark.parametrize("version", ["v7", "v73"])
def test_actual_matlab_archives_preserve_decisions(tmp_path, version, capsys):
    from pathlib import Path
    fixtures = Path(__file__).parent / "fixtures"
    outputs, _, _ = run_benchmarks(replay_file=fixtures / f"benchmark_matlab_{version}.mat",
                                  output_root=tmp_path / "replay")
    assert len(capsys.readouterr().out.strip().splitlines()) == 5
    assert len(outputs) == 4
    for name, result in outputs.items():
        with np.load(fixtures / f"experiment_{name}.npz", allow_pickle=False) as reference:
            np.testing.assert_array_equal(result["pvalues"], reference["pv"].reshape(1, -1))
            np.testing.assert_array_equal(result["selected"], reference["selected"].reshape(1, -1))
            wealth = reference["wealth"].reshape(-1)
            hits = np.flatnonzero(wealth[256:] > 100)
            np.testing.assert_array_equal(result["point"], [hits[0] + 1 if len(hits) else 0])
            np.testing.assert_array_equal(result["fa"], [bool(np.any(wealth[:256] > 100))])
