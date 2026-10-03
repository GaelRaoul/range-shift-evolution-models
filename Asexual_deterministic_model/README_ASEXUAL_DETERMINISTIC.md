# Deterministic asexual simulation

This folder contains the deterministic asexual research simulation.

## Files

- `asexual-deterministic-macro.py`: creates a macro summary file and launches the core simulation for the requested climate speed(s).
- `core-asexual-deterministic.py`: runs one deterministic simulation for one climate-change speed `c`.

The scripts write to a relative `data/` directory and create it automatically if it does not already exist.

## Current parameters in the supplied core

| Code name | Value | Meaning |
|---|---:|---|
| `Vs` | 1 | Stabilizing-selection scale |
| `B` | 1.2 | Spatial slope of the phenotypic optimum |
| `X` | 100 | Width of the moving spatial computation window |
| `dx` | 0.2 | Spatial grid spacing |
| `dy` | 0.2 | Phenotypic grid spacing |
| `dt` | 1 | Time step |
| `T` | 40 | Final time in the supplied core |
| `sigmax` | 1 | Spatial Gaussian dispersal parameter |
| `sigmay` | `0.8**2` | Phenotypic Gaussian dispersal/mutation parameter |
| `c` | command-line | Climate-change speed |

The phenotype-domain width is `Y = 4 + X*B`.

The initial spatial center is

```python
x0ini = 0.7*x[1] + 0.3*x[-1]
```

and the supplied initial density is rescaled so that the maximum phenotype-integrated population density over space is one.

All variables and parameters are in nondimensional model units.

## Running the supplied macro

From this folder:

```bash
python3 asexual-deterministic-macro.py
```

The supplied macro launches one run with `c = 0.65`.

## Running the core directly

```bash
python3 core-asexual-deterministic.py C MACRO_INDEX
```

For example, if `data/data_macro0.csv` has already been created:

```bash
python3 core-asexual-deterministic.py 0.65 0
```

## Numerical method

The spatial and phenotypic Gaussian dispersal steps are evaluated by the direct one-dimensional convolutions present in the source. Growth and stabilizing selection are then applied on the interior grid. The spatial and phenotypic computational windows are shifted when required to follow the population.

The propagation speed is estimated by linear regression of the recorded front position over the second half of the simulation, using `sklearn.linear_model.LinearRegression`.

## Outputs

For one run, the core writes:

- `data/resultsN.csv`, with columns `t,X,Xtheta,I,x1,y1`;
- `data/coefficientsN.csv`, containing run parameters and the fitted propagation speed;
- one appended row in `data/data_macroM.csv`;
- `data/vectorscT40.csv`, with columns `x,In,Zn,Vn` for the final spatial profile of the supplied `T=40` run.

In the final profile, the exported spatial coordinate is shifted by the initial spatial center `x0ini`, and the exported mean phenotype is shifted by the initial optimum `B*x0ini`. The variance is not changed by this centering.

The filename `vectorscT40.csv` is hard-coded in the supplied core. If `T` is changed to reproduce a different historical snapshot, the output filename should also be changed if distinct files are to be retained.

## Figure 1

Figure 1 used `c=0.65`, with the parameter values recorded in `../documentation/SIMULATION_PARAMETERS.md`. The displayed profiles correspond to `t=0`, `t=40`, and `t=80`. Historically these were stored in files named

- `vectorscT00.csv`,
- `vectorscT40.csv`,
- `vectorscT80.csv`.

The result CSV files themselves are not included in the repository.
