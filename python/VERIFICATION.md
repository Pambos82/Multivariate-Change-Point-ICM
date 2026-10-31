# Python verification

Verified locally on 2026-10-31. The MATLAB package and existing research results were preserved. This summary distinguishes full native runs, shared-input comparisons, rendering checks, and smoke tests.

## Environment and commands

The clean environment used Python 3.14.4, NumPy 2.4.4, SciPy 1.17.1, Matplotlib 3.11.2, openpyxl 3.1.5, h5py 3.16.0 and pytest 9.2.0 on Windows. MATLAB was used to export independent numerical fixtures; it is not required to execute the Python package or its bundled tests.

Executed commands included:

```sh
python -m pip install -r requirements.txt
python -m pip check
python -m pytest -q
python run_benchmarks.py
python run_benchmarks.py --reps 1
python run_benchmarks.py --replay-file tests/fixtures/benchmark_matlab_v7.mat
python run_benchmarks.py --replay-file tests/fixtures/benchmark_matlab_v73.mat
python run_full_figure_sweeps.py --workers 16 --datasets protein --k-list 50 --lambda-list 0.001 --reps 1 --max-jobs 1
```

The bounded sweep was resumed with its recorded output directory. Full native sensitivity and Air Quality workflows were also executed, using the command-line runners and equivalent Python APIs described below. Detailed logs, checkpoints, figures, and machine-readable reports are retained in separate local `results` directories and are excluded from publication.

## Test results

**78 tests and 45 subtests passed** in the clean environment. Coverage includes:

- All four dataset checksums, shapes and exact numeric values against newly computed MATLAB matrix digests.
- Betting boundaries, repeated values, candidate histories, 100/200/5000-position windows, startup behavior, uniform fallbacks and very small distances.
- Numeric, cell-like and structured p-value inputs, tie handling and random-draw consumption.
- Both NCMs, the K = 1 covariance convention and deterministic nearest-neighbour tie ordering.
- Native RNG repeatability and the benchmark's interleaved preparation/experiment sequence.
- Python NPZ replay and actual MATLAB v7/v7.3 benchmark replay.
- Complete-grid and duplicate-setting checks, checkpoint integrity, changed-input/source/version rejection and resume behavior.
- Original XLSX schema, precision, cleaning, date-only week indexing and Figure 3 output.

The implementation retains the MATLAB source's extreme-distance floating-point behavior rather than changing its scientific procedure. Numerical tests establish agreement on their tested inputs; they are not a proof of all theoretical martingale assumptions.

## Fresh MATLAB fixtures

The current MATLAB functions independently generated five complete betting-history cases, four shared-input dataset cases, a K = 1 case and a nearest-neighbour tie case. All 14 parameterized reference tests passed. Neighbour indices, p-values and candidate selections matched; continuous quantities passed the declared floating-point tolerances. The small portable inputs and expected arrays are in `tests/fixtures`.

## 400 shared-input MATLAB repetitions

The current Python code independently recomputed all NCMs, p-values and betting histories for 100 archived MATLAB random-input plans per dataset. The archived inputs and reference outputs were reused; the Python calculations were executed anew.

Across all **400 repetitions**, p-values, candidate selections, detection delays and false-alarm flags matched exactly. Every summary below agrees between the paired MATLAB references and Python:

| Dataset | Repetitions | Detection | False alarm | Median delay | Mean delay | SD delay |
|---|---:|---:|---:|---:|---:|---:|
| Protein | 100 | 0.84 | 0.00 | 61.5 | 70.845238 | 33.837490 |
| Particle | 100 | 1.00 | 0.00 | 15.0 | 14.480000 | 3.056637 |
| Sensorless | 100 | 1.00 | 0.00 | 45.0 | 52.700000 | 27.345710 |
| Credit | 100 | 1.00 | 0.00 | 67.0 | 73.880000 | 26.441962 |

The largest score discrepancy was 1.65e-9 in absolute units, or 9.90e-11 after scaling by `max(1,abs(reference))`. The largest scaled linear-capital discrepancy was 1.64e-14; the largest absolute log-capital discrepancy was 8.53e-14. The cross-language assertion tolerance was 2e-10 relative plus 2e-10 absolute. Large absolute capital values are therefore assessed with a relative tolerance; raw floating-point arrays are not claimed to be bitwise identical.

## Full native Python benchmark

A separate **full benchmark** used native NumPy random inputs: 100 repetitions on each of the four datasets, seed 42 and the default scientific settings. It completed successfully:

| Dataset | Repetitions | Detection | False alarm | Median delay | Mean delay | SD delay |
|---|---:|---:|---:|---:|---:|---:|
| Protein | 100 | 0.82 | 0.00 | 67.0 | 74.865854 | 47.044702 |
| Particle | 100 | 1.00 | 0.00 | 14.0 | 14.560000 | 3.309979 |
| Sensorless | 100 | 0.99 | 0.00 | 51.0 | 56.757576 | 32.199978 |
| Credit | 100 | 1.00 | 0.01 | 68.5 | 72.590000 | 26.330109 |

Native Python and MATLAB random generators do not produce identical samples from equal seed numbers. This native run is distinct from the paired shared-input verification above. Delay statistics include detected repetitions only. A second native one-repetition smoke run passed in the fresh environment with only the final eight-column table on the console.

## Additional 100-repetition shared-input check

The largest detection-rate difference in the native sensitivity grids occurred for Sensorless at K = 81 and covariance diagonal adjustment = 0.01: MATLAB detected 69 of 100 changes; Python detected 51 of 100. A targeted check regenerated the actual Python random-input plan and exactly reproduced its saved checkpoint. MATLAB then independently recomputed NCMs and betting histories using the same prepared matrix, sampled reference/query rows, shifts, and tie-breaking uniforms.

MATLAB using those shared inputs also detected **51 of 100 changes**, with median delay 114, mean delay 125.392157, and delay SD 48.825435. All 51,200 p-values and candidate selections, and all 100 detection delays and false-alarm flags, matched exactly. The maximum absolute score discrepancy was 8.64e-12 and the maximum absolute log-capital discrepancy was 1.95e-14. Fallback flags and status codes matched as well.

Together with the earlier four-dataset comparison, this gives **500 shared-input experiments with identical p-values, candidate selections, detection delays, and false-alarm flags**. The controlled check identifies different realized native inputs as the explanation for that particular 69% versus 51% difference; it is not a claim of bitwise identity for every possible input.

## Full native sensitivity experiment and Figures 1–2

The complete native Python experiment finished with **2,400 settings × 100 repetitions = 240,000 experiments**, using 16 configured workers. It covers all four datasets, K = 1–100, five ellipsoid covariance diagonal adjustments (1e-5, 1e-4, 1e-3, 1e-2, 1e-1), and the distance NCM. Figures 1 and 2 were reconstructed to PNG/PDF from that complete native Python grid and visually checked.

Both saved MATLAB and Python grids contain exactly the same 2,400 unique setting keys, 100 repetitions per setting, and scientific parameters. Descriptive differences between their independently sampled native runs are:

| Figure | Settings | Mean absolute detection-rate difference | Maximum detection-rate difference | Mean absolute false-alarm-rate difference | Maximum false-alarm-rate difference |
|---|---:|---:|---:|---:|---:|
| Figure 1 | 800 | 0.910 percentage points | 12 percentage points | 0.166 percentage points | 2 percentage points |
| Figure 2 | 2,000 | 1.214 percentage points | 18 percentage points | 0.170 percentage points | 2 percentage points |

The figures share 400 ellipsoid settings, so their setting counts overlap. Across all 2,400 unique settings, the mean absolute detection-rate difference is 1.128 percentage points. These are descriptive comparisons, not paired statistical significance tests. Delay summaries exclude misses; at settings with very few detections, a conditional mean or median can vary substantially. Undefined delay summaries are retained.

The MATLAB sweep uses Threefry with the Inversion normal transform for experiments and Twister/Ziggurat for preprocessing. Python uses NumPy RandomState/MT19937. Equal seed numbers therefore do not specify identical permutations, Gaussian shifts, or jitter. Prepared Protein and Credit matrices match exactly; Particle and Sensorless use independently generated 0.001-scale jitter. The benchmark and sweep also have different generator/preparation setups within MATLAB, so their K = 50 summaries need not coincide.

Rendering was checked separately from experiment generation. A fresh read-only MATLAB session extracted actual saved FIG curve data. Rendering the same MATLAB CSV with the production Python renderer reproduced **all 5,600 plotted values across 56 curves exactly**, with identical axis limits. Differences in native curves are present in their source experiment tables; fonts, legend placement, and text rendering can differ.

## Full native Air Quality experiment and Figure 3

A fresh preparation of the original XLSX and a full native Python Air Quality run completed. The cleaned 6,941 × 12 matrix and cleaned/monitored row dates match the MATLAB inputs exactly by SHA-256. There are 3,000 reference observations and 3,941 monitored observations; all 11 checked numeric settings agree.

All 3,941 native p-values and candidate selections match the saved MATLAB run exactly. The maximum absolute NCM difference is 2.03e-13 and log-capital difference is 7.99e-15. Both first cross the threshold at monitored observation **85**, cleaned observation 3,085, in **week 27**. The overlaid curves coincide visually. Fallback and status flags also agree.

The only random inputs in this workflow are the uniform tie-breaking draws. A fresh check verified that MATLAB `rng(42,'twister'); rand(3941,1)` and independent NumPy `RandomState(42).random_sample(3941)` produce identical float64 bytes here. No MATLAB draws or expected outputs were supplied to the native Python run. This specific equality does not imply equality of the permutation and Gaussian-draw mappings used in the simulated experiments.

Air Quality input files remain outside the publication package. Figure 3 uses natural log capital. Both runs have 1,020 identical positions where linear capital overflows to infinity; every saved log-capital value is finite. The full native Python sensitivity run emitted intermediate high-dimensional volume overflow warnings, while its saved capital and log-capital diagnostics contain no nonfinite values. The original scientific formulas were retained.

## Smoke and checkpoint checks

A four-job K = 1/50 check covered both NCMs using 16 configured workers. Two-job K = 50 checks passed partial execution, resume, and completed-run reuse with 16 configured workers. Following a transient Windows progress-file replacement error in the full run, all saved scientific checkpoints were recovered without recomputing the experiment jobs. The file writer's bounded retry behavior passed seven additional tests; a fresh partial/resume check matched the corresponding full-run inputs and decisions exactly. These smoke checks are separate from the completed full sensitivity run.

## Coverage limits

Full native workflows were run; the 500 shared-input experiments cover the settings described above, not every repetition of the complete sensitivity grids. Numerical agreement on these inputs does not establish every theoretical martingale assumption. Other operating systems and Python versions were not execution-tested. MATLAB-specific MAT/FIG output formats are represented by Python NPZ/JSON/CSV arrays and PNG/PDF figures, with the plotted rows retained for reconstruction. The detailed local run and comparison reports are excluded from publication; this summary and the small portable test fixtures are included.

Credit results contain information from [Credit Card Fraud Detection](https://www.kaggle.com/datasets/mlg-ulb/creditcardfraud), available under the [Open Database License](https://opendatacommons.org/licenses/odbl/1-0/). Keep the [data attribution](DATA_ATTRIBUTION.md) and [license notices](LICENSE-NOTICE.md).
