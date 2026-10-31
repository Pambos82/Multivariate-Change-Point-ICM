# Multivariate Change-Point Detection

Code associated with *Multivariate Change-Point Detection Using Feature-Based Inductive Conformal Martingales*, by Charalambos Eliades and Harris Papadopoulos.

## Implementations

- [MATLAB](matlab/README.md): benchmark experiments, sensitivity sweeps, and figure generation.
- [Python](python/README.md): the same benchmark experiments, sensitivity sweeps, and figure generation, with cross-language verification tests.

The implementations run independently. MATLAB is not required to run the Python code. Both use cautious epsilon = 100, look-back window W = 5000, and alarm threshold = 100.

## Run the MATLAB code

Set `matlab/` as the current folder in MATLAB. For a short check, run:

```matlab
restoredefaultpath;
addpath(pwd);
setenv('COPA_DATA_DIR','');
[All, Tables, smokeDir] = run_benchmarks('Reps',1);
```

See the [MATLAB README](matlab/README.md) and [experiment guide](matlab/EXPERIMENTS.md) for the full benchmark and figure commands.

## Run the Python code

Set `python/` as the current folder. On Windows, create a local environment and run a short check:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe run_benchmarks.py --reps 1
```

Omit `--reps 1` for the full 100-repetition benchmark. See the [Python README](python/README.md) and [experiment guide](python/EXPERIMENTS.md) for setup on other platforms, full sensitivity sweeps, and the separate figure reconstruction command.

The [Python verification summary](python/VERIFICATION.md) documents the complete 240,000-experiment sensitivity run, Figures 1–3, and 500 experiments checked against MATLAB using identical inputs. Native random sampling differs across languages; shared-input checks verify the calculations separately.

## Datasets and outputs

The four prepared benchmark datasets are supplied in [`matlab/datasets/`](matlab/datasets/README.md) and [`python/datasets/`](python/datasets/README.md). The MAT and NPZ files contain the same numeric feature matrices. Each dataset folder documents its sources, citations, licenses, and preprocessing. Air Quality input is obtained and prepared separately as documented in each implementation.

Generated experiments, figures, logs, checkpoints, and local environments are excluded from version control. Each implementation saves new outputs under its own `results/` folder.

## Licenses

The code is covered by the [MIT License](LICENSE). Dataset licenses remain separate; retain the [MATLAB dataset notices](matlab/LICENSE-NOTICE.md), [Python dataset notices](python/LICENSE-NOTICE.md), and their data attribution documents. Credit data and derived database contents retain the documented ODbL/Database Contents requirements.
