# Browser-only population range-shift simulator

This folder contains the browser version of the side-by-side three-model simulator. It runs locally in the browser through Pyodide/WebAssembly.

Simulation A and Simulation B can independently use:

- the deterministic asexual model;
- the stochastic asexual model;
- the sexual model.

The interface exposes model-specific and shared controls including climate speed, environmental gradient, asexual phenotype diffusion, sexual `VLE`, stochastic population scale, selection/dispersal/growth parameters, final time, resolution, and snapshot count. The current browser stochastic implementation also exposes a user-controlled seed so that an interactive realization can be repeated exactly. This browser seed is separate from the archived research stochastic code, which does not set an explicit seed.

Very small `VLE` values enlarge the computational phenotype interval, so browser runs can become substantially slower because the trait grid is automatically refined to avoid losing resolution.

For each simulation and snapshot, the page displays:

1. the population heatmap over space and phenotype;
2. the spatial population size `N(x)`;
3. the local mean phenotype `Z(x)`, together with the environmental optimum `Zopt(x)`;
4. the local phenotypic variance `V(x)`.

The same two vertical dashed range markers are drawn on the heatmap and on the `N(x)`, `Z(x)` and `V(x)` panels. Their positions are defined from the local population size as the outer crossings of `N(x)=0.1 k` for the deterministic asexual and sexual browser models and `N(x)=0.1 K` for the stochastic asexual browser model. Crossing positions are linearly interpolated between neighboring spatial grid points. These markers are display diagnostics only and do not modify the simulations.

For the sexual browser model, the displayed phenotype moments use the physical scaling

- `Z(x) = sqrt(VLE) E[v|x]` followed by the display centering used by the interface;
- `V(x) = VLE Var[v|x]`.

Thus the displayed sexual phenotypic variance is naturally of order `VLE` when the variance in the computational phenotype variable is of order one. `VLE` is an editable model parameter and enters the sexual simulation itself, not only the output display.

## Local test

From this directory:

```bash
python3 -m http.server 8000
```

then open

```text
http://localhost:8000
```

Do not open `index.html` directly with a `file://` URL because browser workers and module loading require HTTP(S).

## GitHub Pages

The folder is static and can be published with GitHub Pages. See `DEPLOY_GITHUB_PAGES.md`.

## Scientific note

This browser interface is intended for interactive exploration and uses browser-adapted implementations at reduced/selectable numerical resolution. The research scripts in the three model folders remain the reference code for manuscript reproducibility.
