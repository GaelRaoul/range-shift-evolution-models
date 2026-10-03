# Stochastic asexual simulations

This folder contains the stochastic asexual individual-based research simulation.

## Files

- `asexual-stochastic-macro.py`: defines a climate-speed sweep, creates a macro summary file, and launches the core script.
- `core-asexual-stochastic.py`: runs one stochastic simulation for one climate speed `c` and one population-scale parameter `K`.

The scripts write to a relative `data/` directory and create it automatically if it does not already exist.

## Core command line

```bash
python3 core-asexual-stochastic.py C K MACRO_INDEX
```

There is **no seed argument** in the archived research core.

## Current parameters in the supplied core

| Variable | Value | Meaning |
|---|---:|---|
| `Vs` | 1 | Stabilizing-selection scale |
| `B` | 1.2 | Spatial slope of the phenotypic optimum |
| `X` | 25 | Width of the moving spatial computation window |
| `dx` | 0.1 | Spatial grid spacing |
| `dy` | 0.1 | Phenotypic grid spacing |
| `dt` | 0.2 | Time step |
| `T` | 20 | Final simulation time |
| `sigmax` | 1 | Spatial Gaussian dispersal parameter |
| `sigmay` | `0.8**2` | Phenotypic Gaussian dispersal/mutation parameter |
| `c` | command-line | Climate-change speed |
| `K` | command-line | Population-scale parameter |

All variables and parameters are in nondimensional model units.

## Randomness and reproducibility

The demographic update uses NumPy binomial draws through

```python
np.random.binomial(...)
```

without calling `np.random.seed(...)` or constructing a seeded generator. Therefore, repeated executions with identical model parameters will generally produce different stochastic realizations.

This matches the historical Figure 8 workflow: each plotted point was obtained from one stochastic simulation. Additional simulations were run to check that the resulting propagation-speed curve was stable, but those repeat simulations were not averaged into the plotted points. No historical seed was recorded, so exact trajectory-by-trajectory reconstruction is not possible from the archive alone.

## Figure 8 parameter sweep

The manuscript Figure 8 used 20 climate speeds from `c=0` to `c=-1` and two population scales,

\[
K=10^4,\qquad K=10^5.
\]

The complete parameter record is in `../documentation/SIMULATION_PARAMETERS.md`.

The supplied macro is configured with `cmax = -1`, `nb = 20`, and calls the included `core-asexual-stochastic.py` file.

## Numerical method

The two fixed Gaussian convolution operators are precomputed as matrices and applied by matrix multiplication. The demographic update is stochastic and uses independent binomial draws for survival and births. The moving spatial and phenotypic windows follow the same logic as in the research implementation.

The right-front position used for the propagation-speed regression is the last spatial grid point where the local population size exceeds `K/100`. The propagation speed is fitted over the second half of the simulation using `sklearn.linear_model.LinearRegression`.

## Outputs

For each run, the core writes:

- `data/resultsN.csv`, with columns `t,X,Xtheta,I,x1,y1`;
- `data/coefficientsN.csv`, containing parameters and the fitted propagation speed;
- one appended row in `data/data_macroM.csv`.

The macro-summary files `data_macroXX.csv` were used to build the two propagation-speed curves in Figure 8. Generated result files are not included in the repository.
