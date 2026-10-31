# Multivariate Change-Point Detection

Code associated with *Multivariate Change-Point Detection Using Feature-Based Inductive Conformal Martingales*, by Charalambos Eliades and Harris Papadopoulos.

## Implementations

- [MATLAB](matlab/README.md): benchmark experiments, sensitivity sweeps, and figure generation.

The Python implementation will be added alongside the MATLAB package.

## Run the MATLAB code

Set `matlab/` as the current folder in MATLAB. For a short check, run:

```matlab
restoredefaultpath;
addpath(pwd);
setenv('COPA_DATA_DIR','');
[All, Tables, smokeDir] = run_benchmarks('Reps',1);
```

See the [MATLAB README](matlab/README.md) and [experiment guide](matlab/EXPERIMENTS.md) for the full benchmark and figure commands.

The four prepared benchmark datasets are in `matlab/datasets/`. See their [sources and licenses](matlab/datasets/SOURCES_AND_LICENSES.md), [preprocessing instructions](matlab/datasets/PREPROCESSING.md), and [citations](matlab/datasets/CITATIONS.bib). Air Quality input is obtained separately as documented in the MATLAB package.

## Licenses

The code is covered by the [MIT License](LICENSE). Dataset licenses remain separate; retain the [dataset notices](matlab/LICENSE-NOTICE.md) and [data attribution](matlab/DATA_ATTRIBUTION.md).
