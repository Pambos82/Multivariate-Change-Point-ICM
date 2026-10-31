# Experiment commands and saved outputs

Run these commands from the `python` package directory with the dependencies installed. Examples below use `python` to mean the interpreter from that environment.

## Four benchmarks

```sh
python run_benchmarks.py
python run_benchmarks.py --reps 1
```

The full benchmark runs 100 repetitions on Protein, Particle, Sensorless and Credit. The scientific defaults are seed 42, 4,096 reference observations, two monitored segments of 256, mean-shift scale 0.1, NCM K = 50, covariance diagonal adjustment lambda = 1e-3, cautious epsilon = 100, W = 5000 and capital threshold = 100. Particle and Sensorless receive Gaussian jitter of scale 1e-3 before the experiment.

Reference-set standardization uses sample standard deviation. Local covariance uses divisor K, including the source implementation's scalar-covariance convention at K = 1. Detection uses the first strict capital crossing after the change; a false alarm is any strict crossing in the pre-change segment. Delay summaries include detected repetitions only, and the standard deviation uses the sample convention. The two startup positions are excluded from candidate diagnostics.

The benchmark preserves the preparation/experiment ordering: experiment initialization resets the shared random stream, and the following dataset's jitter consumes its ending state. Python's random stream is native NumPy MT19937; its realized samples differ from MATLAB's.

Each run creates `results/benchmarks/benchmark_*` with `benchmark.npz`, per-dataset archives, `metadata.json`, `benchmark.md` and five CSV tables covering performance, per-repetition decisions, candidate choices, betting constants and numerical issues. The archive contains the p-values, linear/log capital, full candidate histories and diagnostic flags. Load NPZ files with `allow_pickle=False`.

Replay the exact saved p-values:

```sh
python run_benchmarks.py --replay-file results/benchmarks/benchmark_YOUR_RUN/benchmark.npz
python run_benchmarks.py --replay-file PATH_TO_MATLAB_BENCHMARK/benchmark.mat
```

Both MATLAB v7 and v7.3 benchmark archives are supported. Replay does not resample data or generate new p-values. Optional `--epsilon`, `--threshold`, `--data-dir` and `--output-root` arguments correspond to the MATLAB runner options. `COPA_DATA_DIR` overrides the default dataset location when no explicit directory is supplied.

## Joint Figures 1 and 2 grid

```sh
python run_full_figure_sweeps.py --workers 16
```

The default grid is four datasets × K = 1:100 × six settings: the ellipsoid NCM at lambda = 1e-5, 1e-4, 1e-3, 1e-2 and 1e-1, plus the distance NCM. Each setting runs 100 repetitions. Figure 1 uses the shared lambda = 1e-3 ellipsoid subset and the distance NCM. Figure 2 uses all five ellipsoid settings.

The sweep prepares the four datasets in a fixed order and starts each setting from seed 42. Its preparation stream therefore differs from the interleaved benchmark stream. Each worker uses one numerical-library thread. Use a single controller per run directory.

The default output is `results/figure_sweeps/full_*`. Completed jobs are checkpointed individually. Metadata records parameters, source hashes, dataset hashes and prepared-input hashes. Resume the same run with:

```sh
python run_full_figure_sweeps.py --workers 16 --output-root results/figure_sweeps/full_YOUR_RUN
```

Repeat any nondefault scientific arguments when resuming. Worker count and `--max-jobs` are execution settings; changing scientific settings, source or input contents is rejected. Never edit old checkpoint metadata to bypass validation. `--max-jobs` permits a bounded check; it does not make a partial grid suitable for full figure reconstruction.

After `progress.json` reports completion of all 2,400 jobs:

```sh
python reconstruct_copa_figures.py results/figure_sweeps/full_YOUR_RUN --export-pdf
```

The renderer rejects missing, duplicate or unexpected grid settings and requires epsilon = 100. It creates a fresh `figures_*` directory beside the input CSV with `figure1.png`, `figure2.png`, optional `figure1.pdf` and `figure2.pdf`, `PerformanceUsed.csv`, and reconstruction metadata. It does not rerun experiments or fill missing points. The saved rows and renderer retain the editable plot data.

## Air Quality and Figure 3

```sh
python prepare_airquality_input.py YOUR_DOWNLOAD/AirQualityUCI.xlsx YOUR_LOCAL_AIR_DIRECTORY
python run_real.py --data-dir YOUR_LOCAL_AIR_DIRECTORY
```

Use the original XLSX and the [documented preparation steps](datasets/PREPROCESSING.md#air-quality-local-preparation-only). The helper preserves date-only row times, removes the specified measurement, replaces the -200 sentinel, and removes incomplete rows, leaving 6,941 × 12 measurements. It refuses to replace an existing prepared archive.

The first 3,000 clean observations form the reference set and the remaining 3,941 are monitored. Defaults are K = 50, lambda = 1e-3, seed 42, epsilon = 100 and threshold = 100. The strict alarm comparison uses log capital, matching the MATLAB real-stream runner when linear capital overflows.

Outputs are written to a fresh `results/airquality/run_*` directory: Figure 3 PNG/PDF, `airquality_results.npz`, a report, diagnostic CSV files and provenance. The unlabeled real stream does not provide a ground-truth detection delay or an empirical false-alarm rate.

## Scope

Only the current betting procedure and its benchmark/figure workflows are included. The repository does not add baseline experiment or calibration pipelines absent from the MATLAB package. Native file formats differ between languages, while scientific parameters and decision definitions are retained.
