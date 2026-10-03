# Deterministic sexual simulations

This folder contains the deterministic sexual model with infinitesimal-model reproduction.

## Files

- `sexual-macro.py`: defines the climate-change speed(s), creates a macro summary CSV, and launches the core solver.
- `core-sexual-simulations.py`: runs one simulation for one fixed climate-change speed `c`.

The sexual launcher and core create the `data/` directory automatically if it does not already exist.

## Current macro

The supplied macro contains

```python
c0 = 0.25
cmax = 0.25
nb = 1
```

and therefore launches one simulation with `c=0.25`.

## Current core parameters

| Parameter | Value in supplied file |
|---|---:|
| `m` | 600 |
| `n` | 1024 |
| `tfinal` | 80 |
| `L` | 300 |
| `dt` | 0.1 |
| `sigma` | 1 |
| `rmax` | 1 |
| `Vs` | 1 |
| `eta` | 1 |
| `b` | 0.1 |
| `VLE` | 1 |

All variables and parameters are in nondimensional model units.

For Figure S3, `tfinal=60` was used. For Figures S3 and S4, `VLE` was changed between runs as documented in `../documentation/SIMULATION_PARAMETERS.md`.

## Running

From this folder:

```bash
python3 sexual-macro.py
```

A single simulation can be called as

```bash
python3 core-sexual-simulations.py C MACRO_INDEX
```

For example:

```bash
python3 core-sexual-simulations.py 0.25 0
```

requires the corresponding `data/data_macro_sexual0.csv` file if the final macro-summary row is to be appended.

## Derived coefficients

For each run the code computes

\[
\alpha=\frac{V_{LE}}{r_{\max}V_s},\qquad
\beta=b\sqrt{\frac{\sigma^2}{2V_{LE}r_{\max}}},\qquad
\gamma=\eta,
\]

\[
\nu=c\sqrt{\frac{2}{r_{\max}\sigma^2}},\qquad
K=L\beta+40,\qquad dv=K/n.
\]

Changing `VLE` therefore changes both the selection/reproduction scaling and the computational phenotype interval.

## Initial condition

The initial population is Gaussian in space and phenotype around the local optimum. In the supplied code it is additionally rescaled so that the maximum local population size

\[
I_0(x)=\int n_0(x,v)\,dv
\]

is one.

## Numerical method

The infinitesimal-model reproduction convolution is evaluated using FFTs. Spatial diffusion is treated implicitly using a symmetric banded linear solve (`scipy.linalg.solveh_banded`). The spatial calculation window is recentered when the population moves sufficiently far through the domain.

## Outputs

For each value of `c`, the code writes:

- `data/eachspeedVLEi...csv`: time-dependent front and speed quantities;
- `data/RECVLEi...csv`: final and initial spatial profiles;
- one row in `data/data_macro_sexualM.csv`.

The `eachspeed...csv` header is

```text
t,c,vS,vN,vZS,vZN,xS,xN,Ik1,Ik2,lambda,VLE,b,C,tfinal
```

with `vS` and `vN` expressed using the same spatial scaling as the macro-summary speeds, `vZS`/`vZN` the phenotype-front speeds, `xS`/`xN` the left/right front positions, and `C = eta`.

The `REC...csv` header is

```text
x,I,Z,V,Zopt,x0,I0,Z0,V0,Zopt0,xshift
```

The exported spatial and phenotype coordinates are centered according to the transformations implemented in `export_final_profiles`. The phenotype variance `V` is computed around the uncentered local mean in the internal phenotype coordinate and is then multiplied by `VLE`; the additive centering used for the exported mean phenotype `Z` therefore does not alter `V`.

The manuscript figure mapping and the values of `VLE` used in Figures S3 and S4 are recorded in `../documentation/SIMULATION_PARAMETERS.md`.

## Reproducibility

The sexual model is deterministic; no random seed is used.
