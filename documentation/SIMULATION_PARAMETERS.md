# Simulation parameters and manuscript figure mapping

This file records the simulation parameters used for the figures in the
manuscript and the default parameters encoded in the supplied Python research scripts.

The goal is to make it possible to identify, for each simulation-based figure,
the model, numerical parameters, biological parameters, stochastic-replicate
information where relevant, and the code/output files used to generate the figure.

# 1. Parameters used for manuscript figures

## Figure 1

Figure 1 was obtained with the deterministic asexual model.

### Biological/model parameters

| Parameter | Value | Comment |
|---|---:|---|
| \(r_{\max}\) | 1 | Maximum intrinsic growth rate; encoded as the constant `1` in the growth term |
| \(V_s\) | 1 | Strength/scale of stabilizing selection |
| \(\sigma\) | 1 | Spatial dispersal parameter (`sigmax` in the code) |
| \(k\) | 1 | Density-dependence parameter; implicit because the code uses the local density directly |
| \(b\) | 1.2 | Spatial slope of the phenotypic optimum (`B` in the code) |
| \(\mu\) | 0.8 | Trait-diffusion parameter; the code stores `sigmay = 0.8**2` |
| \(c\) | 0.65 | Climate-change speed |

### Numerical parameters

| Parameter | Value |
|---|---:|
| \(X\) | 100 |
| \(dx\) | 0.2 |
| \(dy\) | 0.2 |
| \(dt\) | 1 |
| \(Y\) | \(4+Xb=124\) |

The figure displays simulated population curves at

\[
t=0,\qquad t=40,\qquad t=80.
\]

Figure 1 was generated from two runs of the same deterministic code:
one with `T = 40` and one with `T = 80`. The curve at \(t=0\) is the
initial condition, the curve at \(t=40\) is obtained from the `T = 40`
simulation (and is also available in the longer run), and the curve at
\(t=80\) is obtained from the `T = 80` simulation.

### Initial condition

The deterministic code first constructs

\[
\widetilde n_0(x,y)=
\frac{1}{5\sqrt{2\pi\,0.74}}
\exp\!\left(
-\frac{(x-x_{\mathrm{ini}})^2}{2}
-\frac{(y-bx)^2}{2\times0.74}
\right),
\]

and then rescales it by the maximum local population size so that
\(\max_x I_0(x)=1\), where \(I_0(x)=\int n_0(x,y)\,dy\).

The spatial center is

\[
x_{\mathrm{ini}}=0.7\,x_1+0.3\,x_{\max}.
\]

- Code files: `asexual-deterministic-macro.py` and `core-asexual-deterministic.py`.
- Stochastic seed: N.A. (deterministic simulation).
- The curves used for the figure were read from files named
  `vectorscT00.csv`, `vectorscT40.csv`, and `vectorscT80.csv`. The time label
  was set explicitly in the simulation code for the corresponding run.
- These result files are not included in the archive; they can be regenerated
  from the supplied scripts.
- Units: all variables and parameters are expressed in nondimensional model units.
## Figure 8

Figure 8 was obtained with the stochastic asexual individual-based model.
Two population-scale values were used:

\[
K=10^4 \qquad\text{and}\qquad K=10^5.
\]

Figure 8 contains two curves, one for \(K=10^4\) and one for
\(K=10^5\).

### Fixed biological/model parameters

| Parameter | Value | Comment |
|---|---:|---|
| \(r_{\max}\) | 1 | Implicit in the birth term `exp((1-rate)*dt)` |
| \(V_s\) | 1 | Strength/scale of stabilizing selection |
| \(\sigma\) | 1 | Spatial dispersal parameter (`sigmax` in the code) |
| \(\mu\) | 0.8 | Trait-diffusion parameter; the code stores `sigmay = 0.8**2` |
| \(b\) | 1.2 | Spatial slope of the phenotypic optimum (`B` in the code) |
| \(k\) | 1 | Implicit density-dependence coefficient; the code uses \(I/K\) |
| \(K\) | \(10^4\) or \(10^5\) | Population-scale parameter |

### Numerical parameters

| Parameter | Value |
|---|---:|
| \(X\) | 25 |
| \(dx\) | 0.1 |
| \(dy\) | 0.1 |
| \(dt\) | 0.2 |
| \(T\) | 20 |
| \(Y\) | \(4+Xb=34\) |

### Climate-change speeds

The macro uses `nb = 20` and

```python
c_macro = np.linspace(0, -1, num=20)
```

so that the 20 simulated climate-change speeds are

\[
c \in \{0, -0.0526316, -0.105263, -0.157895, -0.210526, -0.263158, -0.315789, -0.368421, -0.421053, -0.473684, -0.526316, -0.578947, -0.631579, -0.684211, -0.736842, -0.789474, -0.842105, -0.894737, -0.947368, -1\}.
\]

The simulations for Figure 8 were run for both \(K=10^4\) and \(K=10^5\),
with the other coefficients unchanged.

### Initial condition

For each value of \(K\), the initial integer-valued population is constructed
as

\[
n_0(x,y)=
\max\!\left(
0,\,
\left\lfloor
K\,dy\,
\exp\!\left[
-\frac{(x-x_{\mathrm{ini}})^2}{10}
-\frac{(y-bx_{\mathrm{ini}})^2}{10}
\right]
\right\rfloor
-50
\right),
\]

where

\[
x_{\mathrm{ini}}
=\frac{x_1+x_{\max}}{2}.
\]

Thus changing \(K\) changes both the demographic-noise scale and the initial
number of individuals.

### Random seed and stochastic reproducibility

The supplied Figure 8 code uses NumPy binomial draws through
`np.random.binomial(...)`, but it does not initialize NumPy with an explicit
random seed.

For each value of the climate-shift speed and each value of \(K\), the point
shown in Figure 8 was obtained from a single stochastic simulation. Additional
repeat simulations were performed as a robustness check and gave very stable
results; they were not averaged to produce the plotted points.

The code does not set an explicit NumPy random seed. Therefore the exact
historical stochastic realizations are not reproducible from the archived code
alone, and repeated new runs with the same parameters can produce different
realizations.

### Code and outputs

- Macro: `asexual-stochastic-macro.py`
- Core: `core-asexual-stochastic.py`
- The supplied macro currently contains `K_macro = 100000`; the
  \(K=10^4\) set is obtained by running the same code with
  `K_macro = 10000`.
- Per-run time series: `data/resultsN.csv`
- Per-run parameter/speed file: `data/coefficientsN.csv`
- The two curves in Figure 8 were constructed from macro-summary files of the
  form `data_macroXX.csv` (where `XX` is the automatically assigned run index).
- The plotted quantity is the spatial propagation speed.
- These result files are not included in the archive; they can be regenerated
  by running the supplied scripts.
- Units: all variables and parameters are expressed in nondimensional model units.
## Figure S3

Figure S3 was obtained with the sexual model. All panels use
\(t_{\mathrm{final}}=60\). The corrected archived sexual code is used
with the same coefficients in every panel except for \(V_{LE}\).

### Fixed parameters

| Parameter | Value | Comment |
|---|---:|---|
| \(c\) | 0.25 | Climate-change speed, supplied by `sexual-macro.py` |
| \(\sigma\) | 1 | Spatial dispersal parameter |
| \(r_{\max}\) | 1 | Maximum intrinsic growth rate |
| \(V_s\) | 1 | Strength/scale of stabilizing selection |
| \(\eta\) | 1 | Reproduction-rate parameter |
| \(b\) | 0.1 | Spatial slope of the phenotypic optimum |
| \(m\) | 600 | Number of spatial grid points |
| \(n\) | 1024 | Number of phenotypic grid points |
| \(L\) | 300 | Length of the spatial computational interval |
| \(dx=L/m\) | 0.5 | Spatial grid spacing |
| \(dt\) | 0.1 | Time step |
| \(t_{\mathrm{final}}\) | 60 | Figure-specific final simulation time |

The corrected macro contains `c0 = cmax = 0.25` and `nb = 1`, so \(c=0.25\)
is passed directly to the core simulation.

### Values of \(V_{LE}\) used in Figure S3

| Figure S3 panel | \(V_{LE}\) |
|---|---:|
| (a) | 0.001 |
| (b) | 0.01 |
| (c) | 0.03 |
| (d) | 0.1 |
| (e) | 0.3 |
| (f) | 1 |

Thus \(V_{LE}\) is the only independent coefficient changed between the six
panels.

### Derived coefficients

The core computes

\[
\alpha=\frac{V_{LE}}{r_{\max}V_s},
\qquad
\beta=b\sqrt{\frac{\sigma^2}{2V_{LE}r_{\max}}},
\qquad
\gamma=\eta,
\]

\[
\nu=c\sqrt{\frac{2}{r_{\max}\sigma^2}},
\qquad
K=L\beta+40,
\qquad
dv=\frac{K}{n}.
\]

For Figure S3:

| Panel | \(V_{LE}\) | \(\alpha\) | \(\beta\) | \(K\) | \(dv\) |
|---|---:|---:|---:|---:|---:|
| (a) | 0.001 | 0.001 | 2.236068 | 710.820393 | 0.694161 |
| (b) | 0.01 | 0.01 | 0.707107 | 252.132034 | 0.246223 |
| (c) | 0.03 | 0.03 | 0.408248 | 162.474487 | 0.158666 |
| (d) | 0.1 | 0.1 | 0.223607 | 107.082039 | 0.104572 |
| (e) | 0.3 | 0.3 | 0.129099 | 78.729833 | 0.076885 |
| (f) | 1 | 1 | 0.070711 | 61.213203 | 0.059779 |

The initial condition is the one defined in `core-sexual-simulations.py`, with
the population centered at \(L/2\) in space and around the local phenotypic
optimum.

- Code files: `sexual-macro.py` and `core-sexual-simulations.py`.
- Stochastic seed: N.A. (deterministic sexual model).
- Result files are not included in the archive; they are regenerated by running the simulation scripts.
- Units: all variables and parameters are expressed in nondimensional model units.
## Figure S4

Figure S4 was obtained with the corrected sexual-model code. The independent
model and numerical parameters were kept fixed between groups of panels,
except for \(V_{LE}\).

### Fixed parameters

| Parameter | Value | Meaning |
|---|---:|---|
| \(c\) | 0.25 | Climate-change speed, supplied directly by the macro |
| \(\sigma\) | 1 | Spatial dispersal parameter |
| \(r_{\max}\) | 1 | Maximum intrinsic growth rate |
| \(V_s\) | 1 | Strength/scale of stabilizing selection |
| \(\eta\) | 1 | Reproduction-rate parameter |
| \(b\) | 0.1 | Spatial slope of the phenotypic optimum |
| \(m\) | 600 | Number of spatial grid points |
| \(n\) | 1024 | Number of phenotypic grid points |
| \(L\) | 300 | Length of the spatial computational interval |
| \(dx=L/m\) | 0.5 | Spatial grid spacing |
| \(dt\) | 0.1 | Time step |
| \(t_{\mathrm{final}}\) | 80 | Final simulation time |

The curves displayed in Figure S4 correspond to

\[
t=0,\qquad t=40,\qquad t=80.
\]

The corrected code therefore matches the stated snapshot times directly.

Units: all variables and parameters are expressed in nondimensional model units.

The initial population is centered at the middle of the spatial interval
(\(\mathrm{prop}=0.5\)) and around the local phenotypic optimum. The code first constructs

\[
\widetilde n_0(x,v)=\frac{0.05}{dv}
\exp\!\left(-\frac{|x-L/2|^2}{20}\right)
\exp\!\left(-\frac{(v-Z_{\mathrm{opt},0}(x))^2}{2}\right),
\]

and then rescales it so that the maximum local population size \(I_0(x)\) is one.

### Values of \(V_{LE}\) used in Figure S4

| Figure S4 panels | \(V_{LE}\) |
|---|---:|
| (a,b,c) | 0.001 |
| (d,e,f) | 0.01 |
| (g,h,i) | 0.03 |
| (j,k,l) | 0.1 |
| (m,n,o) | 1 |

Thus, among the independent inputs listed above, only \(V_{LE}\) changes
from one group of panels to another.

### Derived coefficients

The code computes

\[
\alpha=\frac{V_{LE}}{r_{\max}V_s},
\qquad
\beta=b\sqrt{\frac{\sigma^2}{2V_{LE}r_{\max}}},
\qquad
\gamma=\eta,
\]

\[
\nu=c\sqrt{\frac{2}{r_{\max}\sigma^2}},
\qquad
r_* = 1-\frac{V_{LE}}{2r_{\max}V_s},
\qquad
K=L\beta+40,
\qquad
dv=\frac{K}{n}.
\]

For the five values used in Figure S4:

| Panels | \(V_{LE}\) | \(\alpha\) | \(\beta\) | \(\gamma\) | \(\nu\) | \(K\) | \(dv\) |
|---|---:|---:|---:|---:|---:|---:|---:|
| (a,b,c) | 0.001 | 0.001 | 2.236068 | 1 | 0.353553 | 710.820393 | 0.694161 |
| (d,e,f) | 0.01 | 0.01 | 0.707107 | 1 | 0.353553 | 252.132034 | 0.246223 |
| (g,h,i) | 0.03 | 0.03 | 0.408248 | 1 | 0.353553 | 162.474487 | 0.158666 |
| (j,k,l) | 0.1 | 0.1 | 0.223607 | 1 | 0.353553 | 107.082039 | 0.104572 |
| (m,n,o) | 1 | 1 | 0.070711 | 1 | 0.353553 | 61.213203 | 0.059779 |

The code also computes

\[
A_{KB}=\frac{V_{LE}}{r_*r_{\max}V_s},
\qquad
B_{KB}=\frac{b\sigma}{2r_*r_{\max}\sqrt{V_s}},
\qquad
C_{\mathrm{climate}}
=c\sqrt{\frac{2}{r_{\max}r_*\sigma^2}}.
\]

Their values are:

| Panels | \(V_{LE}\) | \(r_*\) | \(A_{KB}\) | \(B_{KB}\) | \(C_{\mathrm{climate}}\) |
|---|---:|---:|---:|---:|---:|
| (a,b,c) | 0.001 | 0.999500 | 0.001001 | 0.050025 | 0.353642 |
| (d,e,f) | 0.01 | 0.995000 | 0.010050 | 0.050251 | 0.354441 |
| (g,h,i) | 0.03 | 0.985000 | 0.030457 | 0.050761 | 0.356235 |
| (j,k,l) | 0.1 | 0.950000 | 0.105263 | 0.052632 | 0.362738 |
| (m,n,o) | 1 | 0.500000 | 2.000000 | 0.100000 | 0.500000 |

### Code used

- Macro/launcher: `sexual-macro.py`
- Core simulation: `core-sexual-simulations.py`
- The macro fixes \(c_0=c_{\max}=0.25\) and `nb=1`, and passes \(c=0.25\)
  directly to the core.
- The profiles used for Figure S4 were read from files of the form
  `RECVLEi0.1b0.1tfinal80c0.25.csv`, with the numerical values changing
  according to the parameter set represented in each panel.
- These result files are not included in the archive; they can be regenerated
  by running the supplied scripts.
# 2. Default parameters in the supplied Python archive

These defaults describe the supplied research scripts. The figure-specific
records in Section 1 should be regarded as the authoritative description of
the manuscript simulations.

## 2.1 Deterministic asexual model

| Parameter | Default |
|---|---:|
| \(V_s\) | 1 |
| \(B\) | 1.2 |
| \(X\) | 100 |
| \(dx\) | 0.2 |
| \(dy\) | 0.2 |
| \(dt\) | 1 |
| \(T\) | 40 in the supplied core; Figure 1 also used a second run with \(T=80\) |
| \(\sigma_x\) | 1 |
| \(\sigma_y\) | \(0.8^2=0.64\) |
| \(c\) | supplied by the macro |

The manuscript parameter sets are documented figure by figure in Section 1.

## 2.2 Stochastic asexual model

| Parameter | Default |
|---|---:|
| \(V_s\) | 1 |
| \(B\) | 1.2 |
| \(X\) | 25 |
| \(dx\) | 0.1 |
| \(dy\) | 0.1 |
| \(dt\) | 0.2 |
| \(T\) | 20 |
| \(\sigma_x\) | 1 |
| \(\sigma_y\) | \(0.8^2=0.64\) |
| \(K\) | 100000 in the latest supplied macro |
| \(c\) | supplied by the macro |
| random seed | no explicit seed is set |

The supplied macro is configured for the manuscript sweep with `c0 = 0`, `cmax = -1`, and `nb = 20`, and it calls the included core file `core-asexual-stochastic.py`.

Historical manuscript seeds: no explicit seed was set or recorded for the Figure 8 simulations.

Number of simulations contributing to each plotted Figure 8 point: 1.
Additional simulations were used only to check stability and were not averaged
into the plotted curve.

## 2.3 Sexual model

The supplied sexual core uses:

| Parameter | Default/current value |
|---|---:|
| \(m\) | 600 |
| \(n\) | 1024 |
| \(L\) | 300 |
| \(dx=L/m\) | 0.5 |
| \(dt\) | 0.1 |
| \(t_{\mathrm{final}}\) | 80 |
| \(\sigma\) | 1 |
| \(r_{\max}\) | 1 |
| \(V_s\) | 1 |
| \(\eta\) | 1 |
| \(b\) | 0.1 |
| \(V_{LE}\) | 1 in the supplied core; changed manually for Figures S3/S4 |
| \(c\) | supplied directly by the macro; current macro value 0.25 |

Figure S3 uses the same coefficients except for its explicitly stated
figure-specific final time \(t_{\mathrm{final}}=60\) and its panel-dependent
\(V_{LE}\). Figure S4 uses \(t_{\mathrm{final}}=80\) and the panel-dependent
\(V_{LE}\) values listed in Section 1.
# 3. Software environment

- Operating system: Ubuntu 24.04.5 LTS
- Python version: `3.12.3`
- NumPy version: `1.26.4`
- SciPy version: `1.11.4`
- scikit-learn version: `1.4.1.post1`
- Hardware/cluster information: not required for reproducing these simulations.
- Repository DOI / release / Git commit: to be assigned at deposition.

# 4. Manuscript-to-output mapping

| Manuscript item | Model | Parameter record | Seed/replicate | Output file(s) |
|---|---|---|---|---|
| Figure 1 | Deterministic asexual | Section 1, Figure 1 | N.A. | `vectorscT00.csv`, `vectorscT40.csv`, `vectorscT80.csv`; not archived |
| Figure 8 | Stochastic asexual | Section 1, Figure 8 | One plotted realization per `(c,K)`; additional runs only checked stability; historical seed not explicitly set | `data_macroXX.csv`, regenerated and not archived |
| Figure S3 | Sexual | Section 1, Figure S3 | N.A. | `Propspeeds...csv`; not archived |
| Figure S4 | Sexual | Section 1, Figure S4 | N.A. | `REC...csv`; not archived |
