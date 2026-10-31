# COPA 2026: Multivariate Change-Point Detection

MATLAB benchmark and figure-generation workflows associated with *Multivariate Change-Point Detection Using Feature-Based Inductive Conformal Martingales*, by **Charalambos Eliades** and **Harris Papadopoulos**.

The commands below generate outputs using the included implementation and recorded parameters.

## Experiment workflows

| Article output | MATLAB entry point | Scope |
|---|---|---|
| Table 1 — proposed-method rows | `run_benchmarks()` | Four datasets; detection, false-alarm and miss rates; detection-delay statistics |
| Figures 1 and 2 | `run_full_figure_sweeps()` and `reconstruct_copa_figures()` | Sensitivity analysis for the ellipsoid and distance nonconformity measures |
| Figure 3 — Air Quality | `prepare_airquality_input()` and `run_real()` | Martingale using separately obtained Air Quality data |
| Table 1 — QuantTree/QT-EWMA baselines; Figure 4 | **Not included** | The complete published baseline experiment and calibration pipelines are not available in this package |

## Requirements

- MATLAB with **Statistics and Machine Learning Toolbox**.
- **Parallel Computing Toolbox** for the full Figures 1–2 sweep; the documented run uses **16 process workers**.
- MATLAB's Java runtime for provenance and checksum operations.

Each run records its MATLAB version.

## Quick start

Set the repository root as MATLAB's current folder, then run a small smoke test:

```matlab
restoredefaultpath;
addpath(pwd);
setenv('COPA_DATA_DIR','');
[All, Tables, smokeDir] = run_benchmarks('Reps',1);
```

This runs **one repetition per dataset**, without Air Quality input. It checks the benchmark workflow but **does not** reproduce the full Table 1 results. Outputs are written under `results/benchmarks/`.

### Run the Table 1 benchmark setup (proposed method)

```matlab
[All, Tables, benchmarkDir] = run_benchmarks();
```

The default benchmark uses **100 repetitions per dataset**, seed `42`, reference size `4096`, two monitored segments of `256` observations, mean-shift scale `0.1`, NCM `K = 50`, covariance adjustment `lambda = 1e-3`, cautious `epsilon = 100`, look-back window `W = 5000`, and alarm threshold `100`.

`Tables.Performance` reports detection, false-alarm and miss rates, plus the mean, median and standard deviation of detection delay **among detected repetitions**. CSV and MAT outputs, including per-repetition diagnostics, are saved in a new benchmark run directory.

### Generate Figures 1 and 2

```matlab
[Performance, sweepDir] = run_full_figure_sweeps('Workers',16);
figureDir = reconstruct_copa_figures(sweepDir,'ExportPDF',true);
```

The complete sweep covers **2,400 settings × 100 repetitions**: four datasets, `K = 1:100`, five ellipsoid covariance adjustments (`1e-5`, `1e-4`, `1e-3`, `1e-2`, `1e-1`), and one distance-NCM setting per `K`. This is a computationally intensive experiment. **Only reconstruct figures after the sweep finishes successfully.** For resuming interrupted runs and checking completion, see [Experiment details](EXPERIMENTS.md).

### Generate Figure 3 (Air Quality)

The Air Quality dataset is **not distributed** with this repository because its source page contains conflicting licensing statements. Obtain the original XLSX for your own authorized use and follow [local preparation instructions](datasets/PREPROCESSING.md#air-quality-local-preparation-only). Keep both the download and prepared input outside the repository.

```matlab
localAirDir = 'YOUR_LOCAL_AIR_DIRECTORY';
prepare_airquality_input('YOUR_DOWNLOAD/AirQualityUCI.xlsx',localAirDir);
[Stream, Summary, airDir] = run_real('DataDir',localAirDir);
```

This experiment produces the martingale and Figure 3. Because the real stream is **unlabeled**, it does not establish a ground-truth detection delay or false-alarm rate.

## Betting functions

The cautious betting procedure in `combination.m` uses a uniform candidate, three interpolated histogram candidates, and three kNN candidates with `k = 5, 10, 15`. Candidate selection uses past candidate histories. `knn_betting_prepare.m` prepares the history-dependent quantities before the current p-value is evaluated. Unsuitable histories or numerical failures use a uniform-betting fallback. The implementation details and interval-boundary conventions are documented in the MATLAB source.

## Repository contents

| Files | Purpose |
|---|---|
| `combination.m`, `knn_betting_prepare.m` | Cautious betting and kNN preparation |
| `pvalues.m`, `knn_ellipsoid_pdf.m`, `run_knn_cm_experiment.m` | Conformal p-values, NCMs and experiment engine |
| `run_benchmarks.m`, `summarize_results.m` | Table 1 benchmarks and summaries |
| `run_full_figure_sweeps.m`, `reconstruct_copa_figures.m` | Figures 1–2 sweeps and plots |
| `prepare_airquality_input.m`, `run_real.m` | Air Quality preprocessing and Figure 3 |
| `datasets/` | Four prepared benchmark inputs and their source/licensing metadata |

The 11 MATLAB runtime files cover the workflows above.

## Data, licensing and citations

The supplied benchmark inputs are `CASP.mat`, `MiniBooNEPID.mat`, `Sensorlessdrivediagnosis.mat` and `creditcard_X.mat`. Consult [dataset sources and licenses](datasets/SOURCES_AND_LICENSES.md), [preprocessing](datasets/PREPROCESSING.md), [checksums](datasets/SHA256SUMS.txt), and [dataset BibTeX entries](datasets/CITATIONS.bib).

The code's [MIT license](LICENSE) **does not replace dataset licenses**. The three UCI datasets require CC BY 4.0 attribution; the Credit database is subject to ODbL 1.0, with individual contents under DbCL 1.0. Retain [license notices](LICENSE-NOTICE.md) and [data attribution](DATA_ATTRIBUTION.md) when redistributing data or publishing Credit-derived tables and figures.

When using this code, cite the associated article and the relevant datasets. See [dataset citations](datasets/CITATIONS.bib).

## Results and further documentation

Generated results, figures, logs, and checkpoints belong in `results/` and are ignored by Git; `results/README.md` is retained. **Do not commit local Air Quality files, temporary outputs, or private research results.**

- [Experiment details, outputs, and resuming runs](EXPERIMENTS.md)
- [Publication verification history](PUBLICATION_VERIFICATION.md)
- [Data preprocessing](datasets/PREPROCESSING.md)
