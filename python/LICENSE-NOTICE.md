# Separate software and dataset licenses

The Python software and author-written software documentation use the [MIT License](LICENSE). This software license does not grant rights to third-party datasets or replace their terms. No dataset or derived data fixture is relicensed under MIT.

The three UCI inputs (`CASP.npz`, `MiniBooNEPID.npz`, `Sensorlessdrivediagnosis.npz`) retain CC BY 4.0 terms and source attribution. The derived Credit database (`creditcard_X.npz`) retains ODbL 1.0, with its individual contents under Database Contents License (DbCL) 1.0. Retain the [sources, citations, and notices](datasets/SOURCES_AND_LICENSES.md), [documented modifications](datasets/PREPROCESSING.md), and linked complete license terms with redistributed data.

The complete machine-readable Credit derivative is supplied as `datasets/creditcard_X.npz`. Its 29 selected feature columns are documented alongside it. ODbL's derivative-database obligations remain separate from the MIT software license; placing independent software and other datasets in the same repository does not relicense them under ODbL. The complete DbCL warranty and liability provisions remain applicable through its linked terms.

Dataset-derived validation fixtures under tests/fixtures retain the corresponding dataset attribution and modification notices. Credit-related fixtures and public results retain the Credit source and database notices; the software license does not override those terms.

Use this attribution with publicly shared Credit tables, figures, and related outputs:

> Credit results contain information from [Credit Card Fraud Detection](https://www.kaggle.com/datasets/mlg-ulb/creditcardfraud), made available under the [Open Database License 1.0](https://opendatacommons.org/licenses/odbl/1-0/). The data originate from the Worldline and Machine Learning Group - ULB research collaboration; the complete derived dataset and its preprocessing are supplied in this repository's `datasets` directory. Individual contents retain [Database Contents License 1.0](https://opendatacommons.org/licenses/dbcl/1-0/) terms.

Protein, Particle, and Sensorless results use the UCI datasets credited in the [dataset documentation](datasets/SOURCES_AND_LICENSES.md), under CC BY 4.0. Keep appropriate creator/source and modification attribution when sharing relevant data, fixtures, or adaptations.

New experiment outputs are saved under `results` and excluded from Git by default. Apply these attribution notices if publishing them separately. Source database licenses remain distinct from the MIT software license.

Air Quality data and saved Air Quality experiment outputs are not bundled because the original UCI page contains conflicting licensing statements. The download and preparation instructions identify that unresolved conflict; this repository does not resolve it or grant permission to use those data.
