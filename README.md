# COPA 2026 experiments with normalized betting

MATLAB code for the proposed method in *Multivariate Change-Point Detection Using Feature-Based Inductive Conformal Martingales*, by Charalambos Eliades and Harris Papadopoulos. This package runs the article's experiment setups using **normalized kNN betting and cautious epsilon = 100 throughout**.

The complete proposed-method workflow is included: four benchmark datasets, Table 1 metrics, Figures 1 and 2 sensitivity experiments, and the Figure 3 Air Quality workflow. The generated values describe the normalized variant; they are not expected to equal the historical unnormalized results. The original source code and saved research results remain in the separate research folder.

## What this package produces

| Article output | Command | Coverage |
|---|---|---|
| Table 1, proposed-method rows | `run_benchmarks()` | Four datasets, 100 repetitions each; detection/false-alarm rates and mean/median/SD detection delay |
| Figures 1 and 2 | `run_full_figure_sweeps`, then `reconstruct_copa_figures` | Full K = 1:100 grid, five choices of lambda for the ellipsoid NCM and the distance NCM, normalized betting in every experiment |
| Figure 3 | `prepare_airquality_input`, then `run_real` | Normalized martingale for independently obtained Air Quality input |
| Baseline Table 1 rows and Figure 4 | Not regenerated here | The supplied research material lacks the complete published QuantTree/QT-EWMA experiment and calibration pipelines |

Both NCM choices in Figure 1 feed the same normalized betting function. Adding lambda*I to the ellipsoid NCM covariance and normalizing the betting density are separate operations.

## MATLAB requirements

Checked using MATLAB with Statistics and Machine Learning Toolbox. The full sensitivity sweep additionally requires Parallel Computing Toolbox; the documented configuration uses **16 process workers**. MATLAB's JVM is used for checksums and provenance. Each run records the MATLAB version in its metadata.

Open this folder as MATLAB's current folder:

```matlab
restoredefaultpath;
addpath(pwd);
setenv('COPA_DATA_DIR','');
[All, Tables, smokeDir] = run_benchmarks('Reps',1);
```

This optional **smoke test** runs one normalized repetition per dataset and writes a new `results/benchmarks` folder. It is not a full benchmark or full figure sweep, and does not require Air Quality. Dataset hashes and earlier verification findings are recorded in the supplied metadata and [publication verification history](PUBLICATION_VERIFICATION.md).

## Table 1 experiments

```matlab
[All, Tables, benchmarkDir] = run_benchmarks();
```

Defaults are 100 repetitions, seed 42, reference size 4,096, two monitored segments of 256, mean-shift scale 0.1, NCM K = 50, covariance diagonal adjustment lambda = 1e-3, cautious epsilon = 100, look-back W = 5,000 and alarm threshold = 100. Particle and Sensorless receive the original 1e-3 Gaussian tie-scrambling noise. Data standardization, NCMs, conformal p-values and detection rules retain the source implementation's definitions.

The default console output is only the final table with these columns, in order: `Dataset`, `Repetitions`, `DetectionRate`, `FalseAlarmRate`, `MissRate`, `MedianDelay`, `MeanDelay`, `StdDelay`. Progress and debug text are suppressed. The runner saves a new `results/benchmarks/benchmark_*` directory containing `benchmark.mat`, per-dataset MAT files and CSV tables. `Tables.Performance` includes detection, false alarms, miss rate, and mean/median/standard deviation of delay **among detected repetitions**. The other tables preserve per-repetition outcomes, normalization constants, candidate selection frequencies and numerical issues. Paper Table 1 reference values and differences are included in the performance summary.

For a short runner check, use `run_benchmarks('Reps',1)` and label it a smoke test. To rerun saved p-values with the same normalized betting:

```matlab
[All, Tables, replayDir] = run_benchmarks( ...
    'ReplayFile',fullfile(benchmarkDir,'benchmark.mat'));
```

## Full Figures 1 and 2

```matlab
[Performance, sweepDir] = run_full_figure_sweeps('Workers',16);
figureDir = reconstruct_copa_figures(sweepDir,'ExportPDF',true);
```

The complete sweep contains **2,400 distinct settings with 100 repetitions each**: four datasets, K = 1:100, five choices of lambda in the covariance S + lambda*I (1e-5, 1e-4, 1e-3, 1e-2, 1e-1), and one distance-NCM setting per K. Figure 1 uses the lambda = 1e-3 ellipsoid subset together with the distance NCM. Figure 2 uses all five lambda settings. This is a long experiment.

Jobs and input preparation are saved under a new `results/figure_sweeps/full_*` directory. Resume the same run with:

```matlab
[Performance, sweepDir] = run_full_figure_sweeps( ...
    'Workers',16,'OutputRoot',sweepDir);
```

Use one controller per run directory. Resume validates source/input hashes and saved checkpoints. An existing MATLAB pool must have the requested worker count. Confirm that `progress.json` reports `complete` with all 2,400 jobs before calling the plotter. The plotter validates the complete CSV grid, rejects incomplete or duplicate grids, and saves a new `figures_<timestamp>` subfolder. It contains `figure1.png`, `figure1.fig`, `figure2.png`, `figure2.fig`, optional `figure1.pdf` and `figure2.pdf`, the exact plotted rows in `PerformanceUsed.csv`, and provenance.

The revised output names and figure labels apply only to future runs. The calculations, parameters and internal method fields are unchanged. Existing results and checkpoint metadata are preserved. The renamed betting function and revised source lists change provenance hashes, so start a fresh sweep for this revision. Do not rewrite older manifests to force their checkpoints to resume with different source hashes.

The sweep prepares all four datasets using Twister/Ziggurat at seed 42 in Protein, Particle, Sensorless, Credit order; each experiment uses Threefry/Inversion at seed 42. The benchmark retains its original interleaved preparation/experiment stream, so benchmark and sweep random streams differ. The source implementation's scalar covariance behavior at K = 1 is preserved.

## Figure 3: Air Quality

Air Quality input is excluded because its original source page contains conflicting licensing statements. Obtain the original XLSX independently for your authorized use, using the [original download and preparation instructions](datasets/PREPROCESSING.md#air-quality-local-preparation-only). Keep the downloaded and prepared input outside the repository.

```matlab
localAirDir = 'YOUR_LOCAL_AIR_DIRECTORY';
prepare_airquality_input('YOUR_DOWNLOAD/AirQualityUCI.xlsx',localAirDir);
[Stream, Summary, airDir] = run_real('DataDir',localAirDir);
```

The XLSX preserves precision lost in the CSV. Preparation preserves the historical date-only timetable row times. The real-stream experiment uses 3,000 reference observations, K = 50, covariance diagonal adjustment lambda = 1e-3, seed 42, normalized betting, epsilon = 100 and alarm threshold = 100. New outputs are written to a `results/airquality/run_<timestamp>` folder: `airquality_results.mat`, `airquality_results.md`, summary/diagnostic CSV files, provenance, and `figure3.png`, `figure3.fig` and `figure3.pdf`. Stable log martingale values support plotting when linear values overflow. This unlabeled real stream does not establish a false-alarm rate or ground-truth detection delay.

## Betting implementation and interpretation

`combination.m` retains the original uniform candidate, three interpolated histograms and cautious selection. The three kNN betting densities, k = 5, 10 and 15, are divided by their integrals over [0,1], computed from historical p-values before evaluating the current value. The corrected values update the candidate histories used by cautious selection. The histogram interpolation formulas are unchanged.

`knn_betting_prepare.m` uses piecewise analytical logarithmic antiderivatives in double precision, preserving the original nearest-neighbour and boundary definitions. Nonintegrable repeated histories, insufficient history and numerical failure use a uniform betting fallback. Such fallback functions are not claimed to be normalized versions of a divergent density.

Analytical formulas evaluated in floating-point arithmetic do not certify a rigorous anytime-valid guarantee. Such a guarantee additionally requires the appropriate conditional null law of the p-values, predictable selection using past information, nonnegative bets, and the required normalization/integrability condition for every selectable candidate, including the unchanged histograms. Numerical tests alone do not establish those assumptions.

The original Table 1 values matched the historical unnormalized implementation at epsilon = 10, whereas the manuscript states 100. This package consistently uses the requested normalized variant at 100. Historical Figure 1 jitter was not seeded, so differences from archived curves cannot be attributed solely to betting normalization. See [reproduction scope](REPRODUCTION_SCOPE.md) for the remaining article-level limitations.

## Included files

| Files | Role |
|---|---|
| `combination.m`, `knn_betting_prepare.m` | Normalized cautious betting |
| `pvalues.m`, `knn_ellipsoid_pdf.m`, `run_knn_cm_experiment.m` | P-values, NCMs and benchmark experiment engine |
| `run_benchmarks.m`, `summarize_results.m` | Four-dataset benchmark and reported metrics |
| `run_full_figure_sweeps.m`, `reconstruct_copa_figures.m` | Full Figures 1 and 2 |
| `prepare_airquality_input.m`, `run_real.m` | Figure 3 preparation and normalized real-stream experiment |

These 11 MATLAB runtime files support the benchmark and figure workflows. The publication package omits development verification scripts and the unused standalone betting evaluator. No unnormalized betting implementation or original-versus-normalized comparison runner is included.

## Data, licenses and GitHub

The four prepared inputs in `datasets` are `CASP.mat`, `MiniBooNEPID.mat`, `Sensorlessdrivediagnosis.mat` and `creditcard_X.mat`. No conversion is required to run them. Their [sources, citations and licenses](datasets/SOURCES_AND_LICENSES.md), [preprocessing](datasets/PREPROCESSING.md), [BibTeX](datasets/CITATIONS.bib) and [checksums](datasets/SHA256SUMS.txt) are included.

The [MIT software license](LICENSE) is separate from dataset licenses. The three UCI inputs retain CC BY 4.0 attribution. Credit's database retains ODbL 1.0, and its individual contents retain DbCL 1.0. Keep [license notices](LICENSE-NOTICE.md) and [data attribution](DATA_ATTRIBUTION.md), including with public Credit tables and figures.

All generated results, figures, logs and verification reports belong under `results` in separate run folders. Git ignores these runtime outputs and retains `results/README.md`. The source research results are preserved elsewhere. `PACKAGE_MANIFEST.json` records the current publication files and hashes, excluding generated outputs and the manifest itself. Run `Prepare_GitHub_Package.bat` in this folder to copy the manifest-listed publication files into a fresh folder under `results/exports`; generated results, logs, caches and local Air Quality input are excluded. Its optional destination argument must also stay under `results/exports`. The exporter always reads this publication folder, refuses existing destinations, and performs no upload.

In GitHub Desktop, add this folder as a local repository and choose **create a repository here** if needed. Review that the MATLAB files, documentation, notices, four datasets and dataset metadata are included; generated results and local Air Quality inputs should be ignored. Commit locally, then use **Publish repository** when ready. No upload is performed by preparation or verification.
