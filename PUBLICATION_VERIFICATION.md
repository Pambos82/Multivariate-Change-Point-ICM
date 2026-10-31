# Publication verification: normalized betting

Checked on 2026-10-30 using MATLAB, Statistics and Machine Learning Toolbox, Parallel Computing Toolbox and the MATLAB JVM.

Only normalized kNN betting is used. Defaults are epsilon = 100, W = 5000, alarm threshold = 100 and betting k = 5, 10, 15. The original histogram formulas, NCM definitions, p-value calculations and experiment parameters are retained.

## Current publication revision

The runtime package now contains 11 MATLAB files. Its betting entry point is `combination.m`; normalized betting and all scientific defaults are retained. The benchmark console is restricted to the requested eight-column final table, while detailed file outputs remain available. The development verification scripts and unused standalone evaluator are not distributed.

A short MATLAB check of this revision passed: `run_benchmarks('Reps',1)` completed one repetition on each of the four datasets. Captured console output contained only the requested final eight-column table. Bare saved-stream replay produced the same results and no `ans` structure dump. Detailed MAT/CSV outputs remained available, all 11 runtime dependencies resolved within this package and MATLAB/toolboxes, and the renamed betting function, uniform fallback and epsilon = 100 / W = 5000 / threshold = 100 defaults passed checks.

The local check report is stored at `results/publication_revision_20261030_204630/smoke_check.json`; this generated artifact is ignored by Git and excluded by the exporter. The full runs listed below were completed during the earlier preparation and were not rerun for this cleanup. No full sensitivity sweep was rerun for this revision.

## Earlier preparation checks

| Command / check | Coverage | Result |
|---|---|---|
| Dataset, normalization, candidate and dependency checks | Four dataset SHA-256/schema checks; 45 normalization cases; candidate/predictability/W=5000 tests; one normalized repetition per dataset; dependency closure | Passed |
| `run_benchmarks()` | **Full benchmark:** 100 repetitions on each of four datasets, 400 total | Passed |
| `run_benchmarks(..., ReplayFile=...)` | End-to-end saved-stream replay and exact parity with existing normalized implementation | Passed |
| `run_full_figure_sweeps` with 16 workers, Protein, K=[1 50], lambda=1e-3, Reps=1 | **Bounded sweep:** four jobs covering both NCMs; partial-run resume; completed-run reuse; incomplete figure-grid rejection | Passed |
| `reconstruct_copa_figures(..., ExportPDF=true)` | Full 2,400-row plotting input reused from the previously completed normalized experiment; PNG/FIG/PDF rendered and visually checked | Passed |
| `run_real(DataDir=localAirDir)` | Full supplied local Air Quality stream: 3,941 monitored observations; exact p-value/NCM/week/martingale/diagnostic parity with the verified normalized epsilon-100 result | Passed |

The maximum independent numerical normalization discrepancy across 45 configurations was 2.4506e-10. Tests include repeated/degenerate histories and uniform fallback cases. These are numerical checks, not a proof of the full anytime-validity assumptions.

## Earlier full benchmark results

| Dataset | Repetitions | Detection | False alarm | Mean delay | Median delay | SD delay |
|---|---:|---:|---:|---:|---:|---:|
| Protein | 100 | 0.84 | 0.00 | 70.845 | 61.500 | 33.837 |
| Particle | 100 | 1.00 | 0.00 | 14.480 | 15.000 | 3.057 |
| Sensorless | 100 | 1.00 | 0.00 | 52.700 | 45.000 | 27.346 |
| Credit | 100 | 1.00 | 0.00 | 73.880 | 67.000 | 26.442 |

Delay summaries are conditional on detection. These are normalized epsilon-100 results; they are not the historical unnormalized Table 1 values. Zero false alarms here describes these 100 repetitions per dataset, not a population guarantee.

## Output locations and limits

Experiments write new subfolders under `results/benchmarks`, `results/figure_sweeps` and `results/airquality`. Short-check logs and reports are kept in separate run folders under `results`. Detailed build-time logs, parity reports and generated figures are retained in the preparation workspace but excluded from the clean upload export. The batch exporter copies only `results/README.md` from the results tree; rerunning the documented commands generates local outputs in new folders.

The complete 2,400-setting, 100-repetition sensitivity experiment was **not recomputed during this packaging check**. The runner was tested on four bounded jobs with 16 workers; full-figure rendering reused the existing completed normalized performance rows. Air Quality input was externally supplied for verification and is excluded from the upload.

Two packaging issues were found and corrected: the normalized-only sweep no longer expects a removed paired-comparison table, and the local cluster object can honor 16 workers when the saved profile has a lower default limit. The saved MATLAB profile was not changed. The initial unsuccessful sweep output remains preserved in the preparation workspace.

The complete published baseline/calibration pipelines and exact Figure 4 provenance remain unavailable. See [reproduction scope](REPRODUCTION_SCOPE.md). All dependency checks resolve to this package and documented MATLAB/toolbox files; no original or instrumented unnormalized betting implementation is required.

Credit results contain information from [Credit Card Fraud Detection](https://www.kaggle.com/datasets/mlg-ulb/creditcardfraud), available under the [Open Database License](https://opendatacommons.org/licenses/odbl/1-0/). Retain the [dataset notices and attribution](DATA_ATTRIBUTION.md).

## Filename cleanup check

The runtime entry points use the shorter names run_full_figure_sweeps, reconstruct_copa_figures and summarize_results. A clean MATLAB smoke check passed one repetition on each of the four datasets, summary parity, renamed entry-point/dependency resolution, and incomplete-grid rejection. The eight-column-only benchmark console and epsilon = 100, W = 5000 and alarm threshold = 100 were retained. No full benchmark or sensitivity sweep was rerun for this filename cleanup. Existing results and original research files were preserved.

## Future output presentation

Future sweeps use full_<timestamp> folders, figure exports use figures_<timestamp>, and the figure files are figure1.*, figure2.* and figure3.*. The exact plotted rows are saved as PerformanceUsed.csv. Air Quality runs use run_<timestamp>, airquality_results.mat and airquality_results.md. Figure captions, titles and legends use cautious betting wording. Internal method identifiers and data fields remain unchanged for compatibility.

A presentation-only check rendered all three figures from previously computed results. Every plotted x/y value exactly matched the saved reference figures, and saved figure text and new artifact names were checked. No experiments were recomputed, and earlier results were left untouched.
