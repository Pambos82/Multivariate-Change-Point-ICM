# COPA 2026: Python implementation

Python benchmark and figure-generation workflows associated with *Multivariate Change-Point Detection Using Feature-Based Inductive Conformal Martingales*, by Charalambos Eliades and Harris Papadopoulos.

This package implements the same experiment workflows and scientific settings as the sibling MATLAB package: four benchmarks, saved-p-value replay, the joint Figures 1 and 2 sensitivity grid, figure reconstruction, and the Air Quality preparation and Figure 3 workflow. Python execution does not require MATLAB.

## Install and run

Use Python 3.11 or later. The recorded verification used Python 3.14.4. From this `python` directory, on Windows:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe run_benchmarks.py --reps 1
```

On Linux or macOS, use `python3` to create the environment and `.venv/bin/python` for subsequent commands. Run from the source checkout with its included `datasets` directory. The pinned dependencies and optional project metadata are supplied in `requirements.txt` and `pyproject.toml`.

The one-repetition command is a smoke test. The full benchmark is:

```powershell
.\.venv\Scripts\python.exe run_benchmarks.py
```

Its default console output is only the final table: `Dataset`, `Repetitions`, `DetectionRate`, `FalseAlarmRate`, `MissRate`, `MedianDelay`, `MeanDelay`, `StdDelay`. Detailed arrays, parameters and diagnostic tables are saved in a new `results/benchmarks/benchmark_*` directory.

## Workflows

| MATLAB entry point | Python counterpart |
|---|---|
| `combination` | `copa2026.betting.combination` |
| `pvalues` | `copa2026.experiment.pvalues` |
| `knn_betting_prepare` | `copa2026.betting.prepare_knn` |
| `knn_ellipsoid_pdf` | `copa2026.experiment.ellipsoid_density` |
| `run_knn_cm_experiment` | `copa2026.experiment.run_knn_cm_experiment` |
| `summarize_results` | `copa2026.benchmarks.summarize_results` |
| `run_benchmarks` | `run_benchmarks.py` |
| `run_full_figure_sweeps` | `run_full_figure_sweeps.py` |
| `reconstruct_copa_figures` | `reconstruct_copa_figures.py` |
| `prepare_airquality_input` | `prepare_airquality_input.py` |
| `run_real` | `run_real.py` |

The Python functions use zero-based array indices. Saved detection delays and candidate-selection codes retain the MATLAB conventions. Python saves pickle-free NPZ arrays with JSON/CSV metadata and PNG/PDF figures. The MATLAB-specific MAT/FIG output formats are replaced by these Python formats; the plotted rows are retained for editing and reconstruction.

## Figures 1 and 2

```powershell
.\.venv\Scripts\python.exe run_full_figure_sweeps.py --workers 16
.\.venv\Scripts\python.exe reconstruct_copa_figures.py results/figure_sweeps/full_YOUR_RUN --export-pdf
```

The full grid contains 2,400 settings with 100 repetitions each. The sweep creates checkpoints and `Performance.csv`; the separate reconstruction command creates Figures 1 and 2 only after validating the complete grid. Resume instructions and output names are in [EXPERIMENTS.md](EXPERIMENTS.md).

## Air Quality and Figure 3

Obtain the original XLSX independently using the [source and preparation instructions](datasets/PREPROCESSING.md#air-quality-local-preparation-only). Keep downloaded and prepared Air Quality input outside this publication folder.

```powershell
.\.venv\Scripts\python.exe prepare_airquality_input.py YOUR_DOWNLOAD/AirQualityUCI.xlsx YOUR_LOCAL_AIR_DIRECTORY
.\.venv\Scripts\python.exe run_real.py --data-dir YOUR_LOCAL_AIR_DIRECTORY
```

The input is excluded from the repository because the original source page contains conflicting licensing statements. The workflow writes a new Figure 3 run under `results/airquality/`.

## Reproducibility and tests

The defaults retain cautious epsilon = 100, look-back W = 5000, alarm threshold = 100, the same seven candidates, and the original interpolation, p-value, covariance and detection rules. The four supplied numeric dataset matrices have been checked against the MATLAB inputs value for value.

MATLAB and NumPy use different random-draw and permutation algorithms. Using seed 42 in both languages does not guarantee identical samples. Both native Python runs and their random-generator settings are saved; shared-input fixtures separately check the calculations against MATLAB. Replay accepts Python benchmark archives and MATLAB v7/v7.3 benchmark files.

The complete native sensitivity grid (2,400 settings, 100 repetitions each) and Figures 1–3 were run locally. Across 500 experiments using identical inputs, MATLAB and Python produced identical p-values, candidate selections, detection delays, and false-alarm flags. Independently sampled native runs have the differences documented in the verification summary.

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

See [VERIFICATION.md](VERIFICATION.md) for executed checks, numerical tolerances and coverage. Generated results, logs, figures, checkpoints and local environments are excluded by `.gitignore`.

## Data and licenses

The four prepared benchmark inputs are in `datasets`. Their [sources and licenses](datasets/SOURCES_AND_LICENSES.md), [preprocessing](datasets/PREPROCESSING.md), [citations](datasets/CITATIONS.bib), [checksums](datasets/SHA256SUMS.txt), and [manifest](datasets/manifest.json) are included.

The [MIT code license](LICENSE) remains separate from dataset licenses. Keep the [license notices](LICENSE-NOTICE.md) and [data attribution](DATA_ATTRIBUTION.md), including the Credit database's ODbL and Database Contents requirements.
