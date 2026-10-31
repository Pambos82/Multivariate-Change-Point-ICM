# Dataset preprocessing and provenance

The bundled NPZ arrays were converted from the author's prepared MATLAB research inputs. Only the numeric storage representation changed: MATLAB tables were extracted to numeric matrices and saved as compressed NumPy `float64` arrays named `X`. Numeric values, row order, and feature order were checked during conversion. NPZ and MAT file bytes and hashes naturally differ. The [manifest](manifest.json) identifies the distributed NPZ files; [sources, citations, and licenses](SOURCES_AND_LICENSES.md) remain separate from the code license.

## Stored raw benchmark inputs

| Input | Original source-to-feature transformation | Stored shape |
|---|---|---|
| Protein | Read `CASP.csv`; discard the first column, `RMSD`; keep `F1` through `F9` in their original order. | 45,730 Ã— 9 |
| Particle | Read `MiniBooNE_PID.txt`; skip its first line of event counts `36499 93565`; retain all rows and all 50 measurements. Retain finite source sentinel/outlier values. | 130,064 Ã— 50 |
| Sensorless | Read `Sensorless_drive_diagnosis.txt`; keep the first 48 feature columns; omit the 49th class-label column. | 58,509 Ã— 48 |
| Credit | The supplied `creditcard.Rdata` representation contains `Time`, `V1` through `V28`, `Amount`, `Class`. Keep exactly `V1` through `V28`, followed by `Amount`; omit `Time` and `Class` without filtering by class. | 284,807 Ã— 29 |

No experiment jitter, standardization, artificial change, or sampled row permutation is stored in these four arrays. Labels and regression targets do not define the injected benchmark changes.

The Credit PCA features were anonymized by the original data providers; this repository did not perform that PCA. Amount retains its supplied units. The derived database and its individual contents retain ODbL 1.0 and DbCL 1.0 respectively. The earlier prepared-input audit compared the supplied Credit matrix with all 8,259,403 values in the local R representation and documented an exact match. The historical download date and source import script for that R file were not recorded. A fresh `creditcard.csv` download is not certified as bit-identical to that historical decimal-to-binary import. Use the included numeric input for prepared-input reproduction.

## Optional preparation from independent downloads

The original download links and citations are in [SOURCES_AND_LICENSES.md](SOURCES_AND_LICENSES.md). Extract them to a separate raw directory. The following illustrates feature selection; a fresh source version or numeric import can differ at low bits from the included prepared inputs.

```python
from pathlib import Path
import numpy as np

raw = Path("YOUR_RAW_DOWNLOAD_DIRECTORY")
out = Path("YOUR_NEW_PREPARED_DIRECTORY")
out.mkdir(exist_ok=False)

arrays = {
    "CASP": np.loadtxt(raw / "CASP.csv", delimiter=",", skiprows=1,
                       usecols=range(1, 10)),
    "MiniBooNEPID": np.loadtxt(raw / "MiniBooNE_PID.txt", skiprows=1),
    "Sensorlessdrivediagnosis": np.loadtxt(
        raw / "Sensorless_drive_diagnosis.txt", usecols=range(48)),
    "creditcard_X": np.loadtxt(raw / "creditcard.csv", delimiter=",",
                               skiprows=1, usecols=range(1, 30)),
}
for name, X in arrays.items():
    with (out / f"{name}.npz").open("xb") as stream:
        np.savez_compressed(stream, X=np.asarray(X, dtype=np.float64))
```

Credit acquisition through its [Kaggle page](https://www.kaggle.com/datasets/mlg-ulb/creditcardfraud) may require signing in. Keep all applicable dataset notices with independent conversions.

## Experiment preprocessing and random inputs

The benchmark uses the same experiment formulas and scientific parameters as the MATLAB package:

1. Particle and Sensorless receive Gaussian tie-scrambling noise with standard deviation 0.001 after loading. Protein and Credit do not.
2. Each repetition reserves 4,096 randomly selected reference rows and two disjoint 256-row monitored segments.
3. Standardization uses reference-only means and sample standard deviations (ddof=1). Standard deviations below 1e-12 become 1; the same reference statistics transform both monitored segments.
4. A single Gaussian shift vector with component standard deviation 0.1 is added to the standardized post-change segment.
5. NCMs feed randomized sequential conformal p-values and cautious betting with epsilon=100, alarm threshold=100 and look-back W=5000.

The benchmark starts a shared NumPy RandomState at seed 42. It interleaves dataset preparation with experiments in Protein, Particle, Sensorless, Credit order. Each dataset's experiment resets that same generator to seed 42; its final state supplies the following dataset's preparation draws. The full sensitivity sweep prepares all four datasets in order before dispatching independent jobs seeded at 42.

NumPy and MATLAB do not generate identical permutation and Gaussian sequences merely because the numeric seed agrees. Shared-input validation and saved p-value replay distinguish formula consistency from random-input differences. A native Python run independently samples its inputs; no saved MATLAB output is substituted for a Python computation.

## Air Quality: local preparation only

Air Quality inputs are not distributed. The [UCI source](https://archive.ics.uci.edu/dataset/360/air+quality) contains a research-only/commercial-exclusion statement and a CC BY 4.0 statement. That conflict was observed in the source package on 2026-10-29 and remains unresolved here. Download the [original ZIP](https://archive.ics.uci.edu/static/public/360/air%2Bquality.zip) for your own authorized use and extract **`AirQualityUCI.xlsx`** outside this repository.

Use the original XLSX. The source CSV rounds measurements and uses semicolons and decimal commas; the earlier audit found differences up to 0.5 relative to the XLSX/historical input. The helper rejects CSV input rather than silently changing the experiment's numeric input.

```powershell
python prepare_airquality_input.py "YOUR_DOWNLOAD/AirQualityUCI.xlsx" "YOUR_LOCAL_AIR_DIRECTORY"
python run_real.py --data-dir "YOUR_LOCAL_AIR_DIRECTORY"
```

The importer requires the original 9,357 observations and 15 columns, including Date and Time. It preserves date-only row dates and validates the separate fractional-day Time column. It retains the 13 measurements, removes NMHC, changes `-200` to missing, and drops rows missing any retained measurement. This yields **6,941 rows Ã— 12 features**, in original temporal order: CO, PT08.S1, C6H6, PT08.S2, NOx, PT08.S3, NO2, PT08.S4, PT08.S5, temperature, relative humidity, absolute humidity.

Week numbers are `floor((clean_date - first_clean_date)/7 days) + 1`, using date-only values. Combining Date and Time into hourly timestamps would change this convention. The prepared local NPZ contains `X`, `week`, `dates`, and `feature_names`; it is not bundled and the helper refuses to overwrite an existing `airquality.npz`.

The real experiment reserves the first 3,000 cleaned rows as the fixed reference and monitors the remaining 3,941 sequentially, with K=50, covariance diagonal adjustment 0.001, cautious epsilon 100, and threshold 100. Its conformal random draws use native NumPy RNG; exact historical MATLAB trajectory equality is not claimed. Source and prepared Air Quality data and saved Air outputs remain excluded from publication. The helper does not download data or grant a license. Synthetic importer checks do not constitute verification using the actual excluded Air Quality dataset.
