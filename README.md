# Simulation code for *When do leading and rear edges of the range shift slower or faster than climate? Insights from a mathematical model*
[![DOI](https://zenodo.org/badge/1401622942.svg)](https://doi.org/10.5281/zenodo.23125708)
## Interactive simulator

An interactive browser-based version of the three models is available here:

**[Launch the interactive simulator](https://gaelraoul.github.io/range-shift-evolution-models/)**

The simulations run directly in the browser; no installation is required.

**Authors:** Gaël Raoul, Matthieu Alfaro, Ophélie Ronce  
**Corresponding author / contact:** Ophélie Ronce — ophelie.ronce@umontpellier.fr  
**Publication or preprint:** Manuscript in preparation; DOI/URL to be added after deposition.  
**Software archive DOI / release:** v1.0.0 — [10.5281/zenodo.23125709](https://doi.org/10.5281/zenodo.23125709)  
**License:** MIT License (see `LICENSE`).

**Code development:** Gaël Raoul developed the simulation codes, in scientific discussion with Matthieu Alfaro and Ophélie Ronce.

## 1. Contents

This repository contains the Python simulation codes associated with the study:

1. a deterministic asexual population model structured by space and phenotype;
2. a stochastic asexual individual-based model;
3. a deterministic sexual population model with infinitesimal-model reproduction;
4. a browser-based interactive simulator for exploring the three models side by side.

All variables and parameters are expressed in nondimensional model units.

Repository structure:

```text
Asexual_deterministic_model/
    asexual-deterministic-macro.py
    core-asexual-deterministic.py
    README_ASEXUAL_DETERMINISTIC.md

Asexual_stochastic_model/
    asexual-stochastic-macro.py
    core-asexual-stochastic.py
    README_ASEXUAL_STOCHASTIC.md

sexual_deterministic_model/
    sexual-macro.py
    core-sexual-simulations.py
    README_SEXUAL.md

documentation/
    SIMULATION_PARAMETERS.md

interactive/
    browser interface and browser-adapted model implementations
```

The generated simulation-result CSV files used during figure preparation are not included. The parameter values, historical filename patterns, and manuscript-to-code mapping are recorded in `documentation/SIMULATION_PARAMETERS.md`.

## 2. Software requirements

Recommended setup:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Tested scientific-code environment:

- Python 3.12.3
- NumPy 1.26.4
- SciPy 1.11.4
- scikit-learn 1.4.1.post1
- Ubuntu 24.04.5 LTS

The two asexual research scripts currently import Matplotlib, although plotting is not required for the saved numerical outputs; `matplotlib` is therefore included in `requirements.txt`.

## 3. Running the research simulations

The asexual launchers use paths relative to their own model folders. Run them from the corresponding folder. The scripts create the `data/` output directory automatically if it does not already exist.

### Deterministic asexual model

```bash
cd Asexual_deterministic_model
python3 asexual-deterministic-macro.py
```

The core can also be called directly as

```bash
python3 core-asexual-deterministic.py C MACRO_INDEX
```

The supplied macro uses `c = 0.65`.

### Stochastic asexual model

```bash
cd Asexual_stochastic_model
python3 asexual-stochastic-macro.py
```

The core interface is

```bash
python3 core-asexual-stochastic.py C K MACRO_INDEX
```

The research stochastic code does **not** set an explicit random seed. Repeating a run with the same input parameters can therefore give a different stochastic realization. This matches the historical Figure 8 workflow: one realization was used for each plotted point, with additional realizations used only to check stability.

The parameter values used for Figure 8, including the 20 climate speeds from `0` to `-1` and the two values `K=10^4` and `K=10^5`, are recorded in `documentation/SIMULATION_PARAMETERS.md`.

### Sexual deterministic model

```bash
cd sexual_deterministic_model
python3 sexual-macro.py
```

The sexual launcher/core creates its `data/` directory automatically if needed. The core interface is

```bash
python3 core-sexual-simulations.py C MACRO_INDEX
```

The supplied macro launches `c = 0.25`. The supplied core has `VLE = 1` and `tfinal = 80`; figure-specific changes of `VLE` and the `tfinal = 60` setting used for Figure S3 are documented in `documentation/SIMULATION_PARAMETERS.md`.

## 4. Output files

The current research scripts retain their historical output organization rather than using a single standardized schema.

### Deterministic asexual

- `data/resultsN.csv`: time series;
- `data/coefficientsN.csv`: parameters and estimated propagation speed;
- `data/data_macroM.csv`: macro summary;
- `data/vectorscT40.csv`: final spatial profile exported by the supplied `T=40` core. For the historical Figure 1 workflow, analogous files were saved as `vectorscT00.csv`, `vectorscT40.csv`, and `vectorscT80.csv`.

### Stochastic asexual

- `data/resultsN.csv`: time series;
- `data/coefficientsN.csv`: parameters and estimated propagation speed;
- `data/data_macroM.csv`: macro summary used to assemble the speed-versus-climate curves.

### Sexual

- `data/eachspeedVLEi...csv`: time-dependent front/speed quantities;
- `data/RECVLEi...csv`: final spatial profile;
- `data/data_macro_sexualM.csv`: macro summary.

See the model-specific README files and `documentation/SIMULATION_PARAMETERS.md` for details.

## 5. Stochastic simulations and seeds

The archived **research** stochastic asexual code uses NumPy binomial draws without explicitly setting a seed. Consequently, exact historical stochastic trajectories cannot be reconstructed from the code alone. This is intentional and reflects the historical simulation workflow.

The browser-based interactive simulator is a separate exploratory implementation. At present it exposes a user-controlled seed so that an interactive stochastic example can be repeated exactly; this browser seed is not a manuscript seed and is not used by the archived research stochastic script.

## 6. Interactive simulator

The `interactive/` folder contains a browser-only Pyodide/WebAssembly interface in which simulations A and B can independently use any of the three models. It displays the two-dimensional population distribution together with the spatial population size, mean phenotype `Z`, phenotypic variance `V`, and the environmental optimum `Zopt` on the `Z` graphs.

The vertical dashed range markers are defined by the outer crossings of `N(x)=0.1 k` for the deterministic asexual and sexual browser models and `N(x)=0.1 K` for the stochastic browser model. The same marker positions are shown on the heatmap and on the `N`, `Z`, and `V` panels.

See `interactive/README.md` and `interactive/DEPLOY_GITHUB_PAGES.md`.

## 7. Manuscript figure mapping

`documentation/SIMULATION_PARAMETERS.md` records the parameter sets and historical output filenames for:

- Figure 1 — deterministic asexual model;
- Figure 8 — stochastic asexual model;
- Figure S3 — sexual model, varying `VLE`, `tfinal = 60`;
- Figure S4 — sexual model, varying `VLE`, snapshots at `t=0,40,80`.

## 8. Citation and reuse

Please cite the associated article once available and the archived software release/DOI once assigned.
