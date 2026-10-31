# Included benchmark datasets

The publication includes exactly these four prepared inputs. No download or conversion is needed to run the four benchmarks.

| File | MATLAB variable and type | Rows | Features | Dataset license |
|---|---|---:|---:|---|
| `CASP.mat` | `CASP`, table | 45,730 | 9 | CC BY 4.0 |
| `MiniBooNEPID.mat` | `MiniBooNEPID`, table | 130,064 | 50 | CC BY 4.0 |
| `Sensorlessdrivediagnosis.mat` | `Sensorlessdrivediagnosis`, table | 58,509 | 48 | CC BY 4.0 |
| `creditcard_X.mat` | `creditcard_X`, double matrix | 284,807 | 29 | ODbL 1.0 database; DbCL 1.0 contents |

These dataset terms are separate from the MIT software license. [Sources, original downloads, citations and license notices](SOURCES_AND_LICENSES.md) accompany the files. [PREPROCESSING.md](PREPROCESSING.md) documents the verified feature selections and experiment preprocessing. [CITATIONS.bib](CITATIONS.bib) provides references; [manifest.json](manifest.json) and [SHA256SUMS.txt](SHA256SUMS.txt) identify the exact MAT files.

The drivers use this folder by default. An explicit `DataDir` argument, where supported, or a nonempty `COPA_DATA_DIR` environment variable overrides that default. Clear any previous override to use the bundled inputs:

```matlab
setenv('COPA_DATA_DIR','');
[All, Tables, smokeDir] = run_benchmarks('Reps',1);
```

This short check runs one normalized-betting repetition per dataset and saves detailed outputs in a new `results/benchmarks` folder. It requires Statistics and Machine Learning Toolbox and does not require Air Quality. The console shows only the final eight-column performance table. Exact file hashes, types and dimensions remain documented in `manifest.json` and `SHA256SUMS.txt`; earlier checks are summarized in [the publication verification history](../PUBLICATION_VERIFICATION.md).

Air Quality is excluded because its source page contains conflicting licensing statements. Obtain it independently from the [original UCI download](https://archive.ics.uci.edu/static/public/360/air%2Bquality.zip) for authorized use and follow [local preparation instructions](PREPROCESSING.md#air-quality-local-preparation-only). Its files and all unapproved input files are excluded by the dataset allowlist. Generated results and prepared caches are ignored by Git. No old experimental outputs are required by these entry points.
