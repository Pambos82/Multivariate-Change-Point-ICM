# Dataset sources, attribution, and licenses

Source and license pages were checked on **2026-10-29**. This publication includes only the four benchmark MAT files listed below. Each is a converted or reduced representation of a third-party dataset. The repository's software uses the [MIT License](../LICENSE); the dataset licenses remain separate. The software license does not override or replace the dataset terms.

The detailed transformations, matrix contents, and preparation commands are documented in [PREPROCESSING.md](PREPROCESSING.md). Machine-readable references are in [CITATIONS.bib](CITATIONS.bib).

| Included file | Original dataset | Dataset terms |
| --- | --- | --- |
| `CASP.mat` | Physicochemical Properties of Protein Tertiary Structure | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |
| `MiniBooNEPID.mat` | MiniBooNE particle identification | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |
| `Sensorlessdrivediagnosis.mat` | Dataset for Sensorless Drive Diagnosis | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |
| `creditcard_X.mat` | Credit Card Fraud Detection | Database: [ODbL 1.0](https://opendatacommons.org/licenses/odbl/1-0/); individual contents: [DbCL 1.0](https://opendatacommons.org/licenses/dbcl/1-0/) |

## Protein: `CASP.mat`

- Creator: **Prashant Rana**.
- Original source: [UCI Physicochemical Properties of Protein Tertiary Structure](https://archive.ics.uci.edu/dataset/265/physicochemical+properties+of+protein+tertiary+structure).
- Original download: [UCI ZIP containing `CASP.csv`](https://archive.ics.uci.edu/static/public/265/physicochemical%2Bproperties%2Bof%2Bprotein%2Btertiary%2Bstructure.zip).
- Dataset citation: Rana, P. (2013). *Physicochemical Properties of Protein Tertiary Structure* [Dataset]. UCI Machine Learning Repository. [doi:10.24432/C5QW3H](https://doi.org/10.24432/C5QW3H).
- License notice: distributed under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/); [complete legal terms, including the warranty disclaimer](https://creativecommons.org/licenses/by/4.0/legalcode.en).
- Conversion notice: the source `RMSD` column is removed, retaining features `F1` through `F9` in their original order, and the numeric matrix is stored in MATLAB format. See [PREPROCESSING.md](PREPROCESSING.md) for the distinction between this stored matrix and subsequent experiment preprocessing.

## Particle: `MiniBooNEPID.mat`

- Creator: **Byron Roe**.
- Original source: [UCI MiniBooNE particle identification](https://archive.ics.uci.edu/dataset/199/miniboone+particle+identification).
- Original download: [UCI ZIP containing `MiniBooNE_PID.txt`](https://archive.ics.uci.edu/static/public/199/miniboone%2Bparticle%2Bidentification.zip).
- Dataset citation: Roe, B. (2005). *MiniBooNE particle identification* [Dataset]. UCI Machine Learning Repository. [doi:10.24432/C5QC87](https://doi.org/10.24432/C5QC87).
- License notice: distributed under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/); [complete legal terms, including the warranty disclaimer](https://creativecommons.org/licenses/by/4.0/legalcode.en).
- Conversion notice: the first source line contains two event counts and is excluded from the numeric feature matrix. All 50 features and the original event order are retained in the MAT file. See [PREPROCESSING.md](PREPROCESSING.md) for subsequent sampling, standardization, and jitter performed by the experiments.

The citation year above follows UCI's recommended citation, which differs from its donation date.

## Sensorless: `Sensorlessdrivediagnosis.mat`

- Creator: **Martyna Bator**.
- Original source: [UCI Dataset for Sensorless Drive Diagnosis](https://archive.ics.uci.edu/dataset/325/dataset+for+sensorless+drive+diagnosis).
- Original download: [UCI ZIP containing `Sensorless_drive_diagnosis.txt`](https://archive.ics.uci.edu/static/public/325/dataset%2Bfor%2Bsensorless%2Bdrive%2Bdiagnosis.zip).
- Dataset citation: Bator, M. (2013). *Dataset for Sensorless Drive Diagnosis* [Dataset]. UCI Machine Learning Repository. [doi:10.24432/C5VP5F](https://doi.org/10.24432/C5VP5F).
- License notice: distributed under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/); [complete legal terms, including the warranty disclaimer](https://creativecommons.org/licenses/by/4.0/legalcode.en).
- Conversion notice: the 49th source column, containing the class label, is removed. The first 48 feature columns are retained in their original order and stored in MATLAB format. See [PREPROCESSING.md](PREPROCESSING.md) for subsequent experiment preprocessing.

The citation year above follows UCI's recommended citation, which differs from its donation date.

For these three CC BY datasets, retain creator/source attribution and supplied notices, link the license, and identify modifications when redistributing. Do not imply creator endorsement or impose restrictions inconsistent with the license. The legal terms linked above govern.

## Credit: `creditcard_X.mat`

- Original source and download: [Credit Card Fraud Detection, Machine Learning Group - ULB, Kaggle](https://www.kaggle.com/datasets/mlg-ulb/creditcardfraud). Use the source page's Download control to obtain `creditcard.csv`.
- Attribution: the dataset was collected and analyzed through a research collaboration between **Worldline** and the **Machine Learning Group of Universite Libre de Bruxelles (ULB)**. Kaggle lists the dataset under Machine Learning Group - ULB.
- Source license label: database under the Open Database License; contents under the Database Contents License.
- Conversion notice: the local MAT file derives from the supplied R data and contains the 29 columns `V1` through `V28`, followed by `Amount`; it excludes `Time` and `Class`. The already anonymized PCA features originate with the dataset providers. The repository does not reconstruct the confidential original features. See [PREPROCESSING.md](PREPROCESSING.md) for the verified conversion and any provenance limitations.

**Dataset license notice.** The database represented by `creditcard_X.mat`, including this derived representation, is made available under [Open Data Commons Open Database License (ODbL) 1.0](https://opendatacommons.org/licenses/odbl/1-0/). Any rights in its individual contents are licensed under [Database Contents License (DbCL) 1.0](https://opendatacommons.org/licenses/dbcl/1-0/). Their complete terms and warranty/liability disclaimers remain applicable. These are dataset notices, not a license for the MATLAB code.

**Access to the derivative.** The complete machine-readable derivative is supplied, without a separate charge, as the adjacent [`creditcard_X.mat`](creditcard_X.mat). Keep this file and its notices available with any public distribution of the derivative or works produced from it. The documented transformations provide additional provenance.

**Redistribution and public outputs.** Preserve license and attribution notices. Publicly used derivative databases are subject to ODbL share-alike and machine-readable access requirements; do not add incompatible restrictions. DbCL also requires ODbL compliance. Independent code and other datasets in this collection retain their separate licensing. See ODbL sections 2.3, 4.2-4.7 and DbCL sections 2-3 for the governing terms.

Use this attribution alongside publicly shared Credit result tables, figures, and related outputs:

> Credit results contain information from [Credit Card Fraud Detection](https://www.kaggle.com/datasets/mlg-ulb/creditcardfraud), made available under the [Open Database License](https://opendatacommons.org/licenses/odbl/1-0/). The data originate from the Worldline and ULB research collaboration; the complete derived dataset and its preprocessing are supplied in this repository's `datasets` directory.

The source requests research citations. Core references are recorded in [CITATIONS.bib](CITATIONS.bib):

- Dal Pozzolo, Caelen, Johnson, and Bontempi (2015), *Calibrating Probability with Undersampling for Unbalanced Classification*. [doi:10.1109/SSCI.2015.33](https://doi.org/10.1109/SSCI.2015.33).
- Dal Pozzolo, Caelen, Le Borgne, Waterschoot, and Bontempi (2014), *Learned lessons in credit card fraud detection from a practitioner perspective*. [doi:10.1016/j.eswa.2014.02.026](https://doi.org/10.1016/j.eswa.2014.02.026).
- Dal Pozzolo, Boracchi, Caelen, Alippi, and Bontempi (2018), *Credit Card Fraud Detection: A Realistic Modeling and a Novel Learning Strategy*. [doi:10.1109/TNNLS.2017.2736643](https://doi.org/10.1109/TNNLS.2017.2736643).

## Air Quality: excluded from the published datasets

Air Quality files are not included in the publication allowlist. Any local copies retained for reproducing past experiments must remain excluded from the GitHub publication.

- Original source: [UCI Air Quality](https://archive.ics.uci.edu/dataset/360/air+quality).
- Original download: [UCI ZIP containing `AirQualityUCI.csv` and `AirQualityUCI.xlsx`](https://archive.ics.uci.edu/static/public/360/air%2Bquality.zip).
- Dataset citation: Vito, S. (2008). *Air Quality* [Dataset]. UCI Machine Learning Repository. [doi:10.24432/C59K5F](https://doi.org/10.24432/C59K5F).
- Source-requested paper: De Vito, Massera, Piga, Martinotto, and Di Francia (2008), *On field calibration of an electronic nose for benzene estimation in an urban pollution monitoring scenario*. Sensors and Actuators B: Chemical 129(2), 750-757. [doi:10.1016/j.snb.2007.09.060](https://doi.org/10.1016/j.snb.2007.09.060).

**Unresolved source licensing conflict, checked 2026-10-29:** the source page's Additional Information restricts use to research and excludes commercial purposes, while its License section states CC BY 4.0 and describes reuse for any purpose. This repository does not resolve that contradiction or offer a license for Air Quality data. Download it from the original source for your own authorized use and follow [PREPROCESSING.md](PREPROCESSING.md) to prepare it locally. The source identifies `-200` as its missing-value sentinel.
