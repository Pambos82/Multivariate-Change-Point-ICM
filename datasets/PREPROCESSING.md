# Dataset preprocessing and provenance

The four distributed files are byte-identical to the author's supplied research inputs. Their exact SHA-256 values, variable names, classes and dimensions are in [manifest.json](manifest.json). Sources, downloads and separate dataset licenses are in [SOURCES_AND_LICENSES.md](SOURCES_AND_LICENSES.md).

## Stored benchmark inputs

No random jitter, experiment standardization, artificial change or row permutation is stored in these files. Labels and regression targets do not construct the injected benchmark changes.

| Input | Verified source-to-matrix transformation | Stored matrix |
|---|---|---|
| Protein | Read `CASP.csv`; discard first column `RMSD`; keep `F1` through `F9`, preserving every row and column order. | 45,730 x 9, table `CASP` |
| Particle | Read `MiniBooNE_PID.txt`; skip the first line containing counts `36499 93565`; retain all 130,064 events and all 50 measurements in source order. No class column is added. Retain finite source sentinel/outlier values. | 130,064 x 50, table `MiniBooNEPID` |
| Sensorless | Read `Sensorless_drive_diagnosis.txt`; retain columns 1:48; omit column 49, the class label. Preserve all rows and their order. | 58,509 x 48, table `Sensorlessdrivediagnosis` |
| Credit | The supplied `creditcard.Rdata` contains `Time`, `V1` through `V28`, `Amount`, `Class`. Retain exactly `V1` through `V28` followed by `Amount` (columns 2:30). Omit `Time` and `Class`, without filtering rows by class. | 284,807 x 29, double matrix `creditcard_X` |

The Credit PCA features were already anonymized by the data providers; this repository did not perform that PCA. Amount remains in its supplied units. The Credit derivative database and contents retain ODbL 1.0 and DbCL 1.0 respectively.

The publication audit compared every value of the three UCI tables directly with the supplied raw text/CSV files in MATLAB; values and row orders matched exactly. A separate read-only decode compared all 8,259,403 Credit values directly with the local R data columns and found an exact match. All four matrices contain finite numbers. The local `creditcard.Rdata` SHA-256 is `0f299e518bcd1aac7a5d41d045e22222a4126b97b292194180b1038b31a5f5c0`.

The historical download date, source version and import script for the R file were not recorded. The mapping to the original Kaggle dataset is documented, but a fresh `creditcard.csv` download has not been checked bit-for-bit against that R representation. Decimal-to-binary import can differ at low bits. Use the included MAT file for exact prepared-input reproduction; re-saving an equivalent matrix need not produce identical MAT bytes.

## Optional conversion of independently downloaded benchmark files

These commands document feature selection and MATLAB schemas. Run them in a separate working directory and save into a new directory, leaving the distributed inputs untouched. `rawDir` contains the downloaded, extracted original files. Obtain Credit through its [Kaggle page](https://www.kaggle.com/datasets/mlg-ulb/creditcardfraud); access may require signing in.

```matlab
rawDir = 'YOUR_RAW_DOWNLOAD_DIRECTORY';
outDir = 'YOUR_NEW_PREPARED_DIRECTORY';
assert(~isfolder(outDir), 'Choose a new output directory.');
mkdir(outDir);

T = readtable(fullfile(rawDir,'CASP.csv'));
assert(height(T)==45730 && width(T)==10);
CASP = T(:,2:10);
save(fullfile(outDir,'CASP.mat'),'CASP');

X = readmatrix(fullfile(rawDir,'MiniBooNE_PID.txt'), ...
    'FileType','text','NumHeaderLines',1);
assert(isequal(size(X),[130064 50]));
MiniBooNEPID = array2table(X,'VariableNames', ...
    cellstr(compose('VarName%d',1:50)));
save(fullfile(outDir,'MiniBooNEPID.mat'),'MiniBooNEPID');

X = readmatrix(fullfile(rawDir,'Sensorless_drive_diagnosis.txt'));
assert(isequal(size(X),[58509 49]));
legacyNames = [{'e07','e06','e05','e1','e2','e3'}, ...
    cellstr(compose('VarName%d',7:48))];
Sensorlessdrivediagnosis = array2table(X(:,1:48), ...
    'VariableNames',legacyNames);
save(fullfile(outDir,'Sensorlessdrivediagnosis.mat'),'Sensorlessdrivediagnosis');

T = readtable(fullfile(rawDir,'creditcard.csv'));
features = [cellstr(compose('V%d',1:28)), {'Amount'}];
creditcard_X = T{:,features};
assert(isequal(size(creditcard_X),[284807 29]));
save(fullfile(outDir,'creditcard_X.mat'),'creditcard_X');
```

The unusual first six Sensorless table labels preserve the supplied MATLAB import metadata; the code selects columns by position. They do not mean that source rows were dropped. The included matrix was checked against every raw row.

## Preprocessing performed by the experiments

1. Particle and Sensorless receive Gaussian tie-scrambling noise, `1e-3*randn(size(X))`, after loading. Protein and Credit do not.
2. Each repetition randomly permutes rows and reserves 4,096 reference observations. Two disjoint 256-observation segments are selected from the remaining pool.
3. Standardization uses the reference-set mean and sample standard deviation, `std(TR,0,1)`. Standard deviations below `1e-12` become 1. The same reference statistics transform both monitored segments. Stored datasets are not standardized globally.
4. One independent shift vector, `0.1*randn(1,d)`, is added to every post-change observation. The pre-change segment is unshifted.
5. Existing p-value, betting, alarm and NCM computations are unchanged by publication packaging. This package uses normalized kNN betting only; the original implementation remains in the separate source research folder.

Generator, seeding and preparation order matter. The legacy Table 1 driver starts Twister at 42 and interleaves preparation and experiments. Historical Figure 1 did not seed initial jitter; Figure 2 seeded preparation. The new full normalized sweep prepares all datasets using Twister/Ziggurat at 42 in Protein, Particle, Sensorless, Credit order and resets each experiment to Threefry/Inversion at 42. See [the package instructions](../README.md). Table 1 and sensitivity sweeps do not use identical streams.

## Air Quality: local preparation only

Air Quality inputs are not distributed. The [UCI source](https://archive.ics.uci.edu/dataset/360/air+quality) has a research-only/commercial-exclusion statement and a CC BY 4.0 statement. This conflict was observed on 2026-10-29 and is unresolved here. Download the [original ZIP](https://archive.ics.uci.edu/static/public/360/air%2Bquality.zip) for your own authorized use and extract **`AirQualityUCI.xlsx`** outside the repository.

Use the XLSX, not the CSV, to reproduce the historical numeric input. The publication audit found that the XLSX measurements match the historical MAT exactly, whereas the CSV rounds some measurements, with a maximum difference of 0.5 in this source copy. For example, the fourth PT08.S1 value is 1375.5 in the XLSX/MAT and 1376 in the CSV. The helper rejects CSV input to prevent a silent change. The CSV also uses semicolons and decimal commas.

The spreadsheet has 9,357 observations and 15 columns including Date and Time. This differs from the 9,358 instance count on the landing page. Converting Date to RowTimes leaves the historical timetable's 9,357 rows and 14 variables: `Time`, followed by the 13 measurements.

The historical timetable stores **date-only RowTimes** (midnight, with repeated dates) and keeps time of day separately as a numeric fraction of a day: 18:00 = 0.75. `prepare_airquality_input` explicitly preserves this convention. Combining Date and Time into hourly RowTimes would change the historical week indexing and is a separate methodological change.

```matlab
localAirDir = 'YOUR_LOCAL_AIR_DIRECTORY';
prepare_airquality_input('YOUR_DOWNLOAD/AirQualityUCI.xlsx',localAirDir);
[Stream,Summary,folder] = run_real('DataDir',localAirDir);
```

The helper refuses to overwrite `airquality.mat`. It constructs `AirQualityUCI` with the positional schema expected by existing scripts. It does not download data or grant a license.

Subsequent preparation selects timetable variables 2:14, discards the third selected measurement (`NMHC_GT_`), changes sentinel `-200` to `NaN`, and removes rows missing any retained measurement. This leaves **6,941 rows x 12 measurements**: CO, PT08.S1, C6H6, PT08.S2, NOx, PT08.S3, NO2, PT08.S4, PT08.S5, temperature, relative humidity and absolute humidity. Original temporal order is retained. The first 3,000 rows form the reference set; 3,941 are monitored. Week is `floor(days(tt-tt(1))/7)+1`, using cleaned date-only RowTimes.

The normalized-only `run_real` generates Figure 3 and writes a new timestamped `results/airquality` subfolder. A short check of the four bundled benchmark inputs uses `run_benchmarks('Reps',1)` without requiring Air Quality. Raw or prepared local Air Quality inputs must remain outside the publication files.
