# Included benchmark inputs

The four archives contain the same numeric feature matrices as the MATLAB package. No MATLAB installation or conversion step is required to load them in Python. Each NPZ file contains one finite, two-dimensional float64 array named X and is opened with allow_pickle=False.

| File | Dataset identifier | Rows | Features | Dataset terms |
|---|---|---:|---:|---|
| CASP.npz | protein | 45,730 | 9 | CC BY 4.0 |
| MiniBooNEPID.npz | particle | 130,064 | 50 | CC BY 4.0 |
| Sensorlessdrivediagnosis.npz | sensorless | 58,509 | 48 | CC BY 4.0 |
| creditcard_X.npz | credit | 284,807 | 29 | ODbL 1.0 database; DbCL 1.0 contents |

Complete numeric value digests were checked against the current MATLAB inputs on 2026-10-31. Numeric values, row order and feature order agree exactly. Only the storage representation differs. The [manifest](manifest.json) records NPZ file hashes, shapes, the originating MATLAB file hashes, and SHA-256 digests of little-endian float64 values in MATLAB column-major order.

Stored inputs contain no experiment jitter, standardization, artificial change, or sampled row permutation. The benchmark and sweep runners apply experiment preprocessing after loading.

The default directory is this folder. An explicit data_dir argument or the COPA_DATA_DIR environment variable can select another directory containing these NPZ files. The explicit argument takes precedence.

Dataset licenses remain separate from the [MIT code license](../LICENSE). Retain the [sources and notices](SOURCES_AND_LICENSES.md), [preprocessing documentation](PREPROCESSING.md), [citations](CITATIONS.bib), and [file checksums](SHA256SUMS.txt).

Air Quality input is excluded. Its original download and local preparation instructions are in [PREPROCESSING.md](PREPROCESSING.md#air-quality-local-preparation-only).
