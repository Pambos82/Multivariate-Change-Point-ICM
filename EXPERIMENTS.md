# Experiment details and reproducibility notes

This document supplements the [main README](README.md) with full-run parameters, output artifacts, checkpoint behavior, and reproducibility notes. It describes the included MATLAB workflows; it is not a record of newly executed tests.

## Smoke test

From the MATLAB repository root:

```matlab
restoredefaultpath;
addpath(pwd);
setenv('COPA_DATA_DIR','');
[All, Tables, smokeDir] = run_benchmarks('Reps',1);
```

Runs one repetition on each benchmark dataset. This is **not** a full benchmark or figure sweep, and does not require Air Quality. Dataset hashes and earlier verification findings are recorded in the supplied metadata and [publication verification history](PUBLICATION_VERIFICATION.md).

## Table 1: proposed-method benchmarks

```matlab
[All, Tables, benchmarkDir] = run_benchmarks();
```

Defaults:

| Parameter | Value |
|---|---|
| Benchmark datasets | Protein (CASP), Particle (MiniBooNEPID), Sensorless, Credit |
| Repetitions | 100 per dataset |
| Seed | 42 |
| Reference size | 4,096 |
| Monitored segments | Two segments of 256 |
| Mean-shift scale | 0.1 |
| NCM neighbours, `K` | 50 |
| Covariance diagonal adjustment, `lambda` | `1e-3` |
| Cautious `epsilon` | 100 |
| Look-back `W` | 5,000 |
| Alarm threshold | 100 |
| Original Gaussian tie-scrambling noise | `1e-3` for Particle and Sensorless |

Standardization, nonconformity measures, conformal p-values and detection rules retain the source implementation's definitions.

Default console output is a final table with these columns, in order: `Dataset`, `Repetitions`, `DetectionRate`, `FalseAlarmRate`, `MissRate`, `MedianDelay`, `MeanDelay`, `StdDelay`. Progress/debug output is suppressed by default.

Each run creates a new `results/benchmarks/benchmark_*` directory containing `benchmark.mat`, per-dataset MAT files and CSV summary tables. `Tables.Performance` contains detection, false-alarm and miss rates; delay statistics are calculated **among detected repetitions**. Additional tables preserve per-repetition outcomes, betting diagnostics, candidate-selection frequencies, and numerical issues. The performance summary includes historical paper Table 1 reference values and their differences from the new run.

To rerun saved p-values with the cautious betting procedure:

```matlab
[All, Tables, replayDir] = run_benchmarks( ...
    'ReplayFile',fullfile(benchmarkDir,'benchmark.mat'));
```

## Figures 1 and 2: full sensitivity sweep

```matlab
[Performance, sweepDir] = run_full_figure_sweeps('Workers',16);
figureDir = reconstruct_copa_figures(sweepDir,'ExportPDF',true);
```

The full sweep covers **2,400 distinct settings with 100 repetitions each**:

- Four benchmark datasets.
- Ellipsoid NCM for `K = 1:100` at each `lambda` in `{1e-5, 1e-4, 1e-3, 1e-2, 1e-1}`.
- Distance NCM for `K = 1:100`.

Figure 1 uses the ellipsoid subset at `lambda = 1e-3` and the distance NCM. Figure 2 compares all five ellipsoid `lambda` settings. Both use the same cautious betting procedure.

Jobs, input preparation and checkpoints are saved beneath `results/figure_sweeps/full_*`. To resume **the same** run:

```matlab
[Performance, sweepDir] = run_full_figure_sweeps( ...
    'Workers',16,'OutputRoot',sweepDir);
```

Use only one controller per run directory. The resume mechanism validates input/source hashes and saved checkpoints. If a MATLAB parallel pool already exists, its worker count must match the requested `Workers` value.

**Before plotting**, verify that `progress.json` reports `complete` and all **2,400 jobs** are present. The figure reconstructor checks for a complete, nonduplicate CSV grid. It writes `figures_<timestamp>/` containing `figure1.png`, `figure1.fig`, `figure2.png`, `figure2.fig`, optionally `figure1.pdf` and `figure2.pdf`, `PerformanceUsed.csv` and provenance metadata.

### Reproduction details

For the full sweep, the four datasets are prepared in Protein, Particle, Sensorless, Credit order with MATLAB Twister/Ziggurat at seed `42`. Each experiment then uses Threefry/Inversion at seed `42`. The benchmark runner retains its original interleaved preparation/experiment random stream, so **the benchmark and sweep random streams are not identical**. The source implementation's scalar-covariance behavior at `K = 1` is retained.

Revised output names and labels affect new runs only; the underlying calculations, parameters and internal method fields are unchanged. Renaming the betting function or revising source lists changes provenance hashes. **Start a fresh sweep after such changes**, rather than editing older manifests to bypass source-hash checks. Existing results/checkpoints should not be rewritten.

## Figure 3: Air Quality

Air Quality input is intentionally excluded because of conflicting license statements on the original dataset page. Obtain `AirQualityUCI.xlsx` independently for authorized use and keep it (and all prepared input) outside the repository. See [data preparation](datasets/PREPROCESSING.md#air-quality-local-preparation-only).

```matlab
localAirDir = 'YOUR_LOCAL_AIR_DIRECTORY';
prepare_airquality_input('YOUR_DOWNLOAD/AirQualityUCI.xlsx',localAirDir);
[Stream, Summary, airDir] = run_real('DataDir',localAirDir);
```

The XLSX preserves measurement precision that is lost in the CSV. Preparation preserves the historical date-only timetable row times.

Default real-stream parameters are 3,000 reference observations, NCM `K = 50`, covariance adjustment `lambda = 1e-3`, seed `42`, cautious `epsilon = 100`, and alarm threshold `100`.

Each run creates `results/airquality/run_<timestamp>/` containing `airquality_results.mat`, `airquality_results.md`, summary and diagnostic CSV files, provenance and `figure3.png`, `figure3.fig`, `figure3.pdf`. Log martingale values are used for stable plotting when linear martingale values overflow. **Unlabeled data cannot establish a false-alarm rate or a ground-truth detection delay.**

## Betting implementation

`combination.m` uses a uniform candidate, three interpolated histogram candidates, and three kNN candidates with `k = 5, 10, 15`. Cautious selection uses past candidate histories. The helper `knn_betting_prepare.m` prepares history-dependent quantities before evaluating the next p-value. Insufficient or unsuitable histories and numerical failures use a uniform-betting fallback.

Empirical tests alone do not establish an anytime-valid false-alarm guarantee; the theoretical assumptions on the p-values and predictable betting procedure must also hold.

## Publication packaging

Generated results, figures, checkpoints, logs and verification reports belong under separate run folders in `results/`. Git ignores runtime output but retains `results/README.md`. Research data and saved historical runs remain separately maintained.

`PACKAGE_MANIFEST.json` records the publication files and their checksums, excluding generated outputs and the manifest itself. Only `results/README.md` is tracked from the results directory.

Keep the dataset license notices and attribution files with redistributed data and public research outputs. Generated results and local Air Quality inputs remain excluded from the publication files.
