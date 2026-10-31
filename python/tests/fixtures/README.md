# MATLAB reference fixtures

These small numeric fixtures were independently exported from the current sibling MATLAB package on 2026-10-31. Python tests consume them without invoking MATLAB.

- `betting_*.npz`: boundaries, repeated values, alternating values, a 5,010-position rolling window, and extremely small distances. These preserve complete candidate histories and decisions.
- `experiment_*.npz`: independently prepared reference/query matrices, tie-breaking draws, nearest neighbours, scores, p-values and betting outputs for all four benchmark datasets.
- `ellipsoid_k1.npz`: the source scalar-covariance convention at K = 1.
- `neighbor_ties.npz`: deterministic nearest-neighbour ordering at equal distances.
- `benchmark_matlab_v7.mat` and `benchmark_matlab_v73.mat`: minimal saved-p-value archives from the MATLAB benchmark, used to verify both replay readers.

Array indices in the numeric neighbour fixtures are zero-based, except `idx` in `ellipsoid_k1.npz`, which retains MATLAB's one-based values and is converted by its test. The fixtures retain the [dataset attribution](../../DATA_ATTRIBUTION.md) and [license notices](../../LICENSE-NOTICE.md), including Credit result attribution. No Air Quality source data are included.
