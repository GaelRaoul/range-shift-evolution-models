from __future__ import annotations

import json
import numpy as np

from unified_models import simulate_model

PRESETS = {
    "asexual_deterministic": {
        "fast": {"nx": 140, "ny": 190, "dt": 1.0},
        "balanced": {"nx": 200, "ny": 270, "dt": 1.0},
        "fine": {"nx": 280, "ny": 380, "dt": 1.0},
    },
    "asexual_stochastic": {
        "fast": {"nx": 140, "ny": 190, "dt": 0.2},
        "balanced": {"nx": 200, "ny": 270, "dt": 0.2},
        "fine": {"nx": 280, "ny": 380, "dt": 0.2},
    },
    "sexual": {
        "fast": {"m": 320, "n": 240, "dt": 0.20},
        "balanced": {"m": 480, "n": 320, "dt": 0.15},
        "fine": {"m": 800, "n": 420, "dt": 0.10},
    },
}


def _float_arg(data, name, default, lo=None, hi=None):
    try:
        value = float(data.get(name, default))
    except (TypeError, ValueError):
        value = float(default)
    if lo is not None:
        value = max(lo, value)
    if hi is not None:
        value = min(hi, value)
    return value


def _int_arg(data, name, default, lo=None, hi=None):
    try:
        value = int(data.get(name, default))
    except (TypeError, ValueError):
        value = int(default)
    if lo is not None:
        value = max(lo, value)
    if hi is not None:
        value = min(hi, value)
    return value


def _density_threshold_crossings(x, density, threshold):
    """Return the outer left/right x-locations where density crosses threshold.

    Crossing positions are linearly interpolated between neighboring grid points.
    If the threshold is never reached, both values are None.
    """
    x = np.asarray(x, dtype=float)
    density = np.asarray(density, dtype=float)
    threshold = float(threshold)
    if x.size == 0 or density.size != x.size or not np.isfinite(threshold):
        return None, None

    finite = np.isfinite(x) & np.isfinite(density)
    above = finite & (density >= threshold)
    inds = np.flatnonzero(above)
    if inds.size == 0:
        return None, None

    def interpolate(i0, i1):
        x0, x1 = float(x[i0]), float(x[i1])
        y0, y1 = float(density[i0]), float(density[i1])
        if not (np.isfinite(x0) and np.isfinite(x1) and np.isfinite(y0) and np.isfinite(y1)):
            return x1
        if abs(y1 - y0) < 1e-15:
            return 0.5 * (x0 + x1)
        a = (threshold - y0) / (y1 - y0)
        a = min(1.0, max(0.0, float(a)))
        return x0 + a * (x1 - x0)

    ileft = int(inds[0])
    iright = int(inds[-1])
    left = float(x[ileft]) if ileft == 0 else interpolate(ileft - 1, ileft)
    right = float(x[iright]) if iright == x.size - 1 else interpolate(iright, iright + 1)
    return left, right


def _attach_range_markers(result, model):
    """Attach common population-range markers to each displayed snapshot.

    Deterministic asexual and sexual models use N(x)=0.1*k.
    The stochastic asexual model uses N(x)=0.1*K, exactly as requested.
    These markers are display diagnostics only and do not alter the model or
    the front-speed calculations.
    """
    params = result.get("parameters", {})
    if model == "asexual_stochastic":
        threshold = 0.1 * float(params.get("K", 0.0))
        threshold_label = "0.1 K"
    else:
        threshold = 0.1 * float(params.get("k", 1.0))
        threshold_label = "0.1 k"

    for snap in result.get("snapshots", []):
        left, right = _density_threshold_crossings(snap.get("x", []), snap.get("density", []), threshold)
        snap["range_left"] = left
        snap["range_right"] = right
        snap["range_threshold"] = float(threshold)
        snap["range_threshold_label"] = threshold_label
    return result


def _downsample_snapshot(snapshot, *, max_x=140, max_trait=180, x_domain=None, trait_domain=None):
    pop = np.asarray(snapshot["population"], dtype=float)
    x = np.asarray(snapshot["x"], dtype=float)
    trait = np.asarray(snapshot["trait"], dtype=float)
    density = np.asarray(snapshot["density"], dtype=float)
    mean_trait = np.asarray(snapshot["mean_trait"], dtype=float)
    trait_variance = np.asarray(snapshot["trait_variance"], dtype=float)
    optimum = np.asarray(snapshot["optimum"], dtype=float)
    mean_trait_heatmap = np.asarray(snapshot.get("mean_trait_heatmap", snapshot["mean_trait"]), dtype=float)
    optimum_heatmap = np.asarray(snapshot.get("optimum_heatmap", snapshot["optimum"]), dtype=float)

    def contiguous_slice(values, domain):
        if domain is None or values.size == 0:
            return slice(0, len(values))
        lo, hi = map(float, domain)
        idx = np.flatnonzero((values >= lo) & (values <= hi))
        if idx.size == 0:
            mid = 0.5 * (lo + hi)
            k = int(np.argmin(np.abs(values - mid)))
            return slice(k, k + 1)
        i0 = max(0, int(idx[0]) - 1)
        i1 = min(len(values), int(idx[-1]) + 2)
        return slice(i0, i1)

    xs = contiguous_slice(x, x_domain)
    zs = contiguous_slice(trait, trait_domain)
    x = x[xs]
    trait = trait[zs]
    pop = pop[xs, zs]
    density = density[xs]
    mean_trait = mean_trait[xs]
    trait_variance = trait_variance[xs]
    optimum = optimum[xs]
    mean_trait_heatmap = mean_trait_heatmap[xs]
    optimum_heatmap = optimum_heatmap[xs]

    sx = max(1, int(np.ceil(len(x) / max_x)))
    sz = max(1, int(np.ceil(len(trait) / max_trait)))
    return {
        "time": float(snapshot["time"]),
        "x": x[::sx].tolist(),
        "trait": trait[::sz].tolist(),
        "population": pop[::sx, ::sz].tolist(),
        "density": density[::sx].tolist(),
        "mean_trait": mean_trait[::sx].tolist(),
        "trait_variance": trait_variance[::sx].tolist(),
        "optimum": optimum[::sx].tolist(),
        "mean_trait_heatmap": mean_trait_heatmap[::sx].tolist(),
        "optimum_heatmap": optimum_heatmap[::sx].tolist(),
        "front_left": snapshot.get("front_left"),
        "front_right": snapshot.get("front_right"),
        "range_left": snapshot.get("range_left"),
        "range_right": snapshot.get("range_right"),
        "range_threshold": snapshot.get("range_threshold"),
        "range_threshold_label": snapshot.get("range_threshold_label"),
        "plot_x_min": snapshot.get("plot_x_min"),
        "plot_x_max": snapshot.get("plot_x_max"),
        "grid_x_min": snapshot.get("grid_x_min"),
        "grid_x_max": snapshot.get("grid_x_max"),
    }


def _serialize_simulation(sim, plot_domain):
    model = sim["model"]
    if model == "sexual":
        max_x, max_trait = 260, 420
    else:
        max_x, max_trait = 140, 180
    x_domain = (plot_domain["x_min"], plot_domain["x_max"])
    trait_domain = (plot_domain["trait_min"], plot_domain["trait_max"])
    return {
        "model": model,
        "c": float(sim["c"]),
        "speed_left": None if sim.get("speed_left") is None else float(sim["speed_left"]),
        "speed_right": float(sim.get("speed_right", 0.0)),
        "snapshots": [
            _downsample_snapshot(
                snap, max_x=max_x, max_trait=max_trait,
                x_domain=x_domain, trait_domain=trait_domain,
            )
            for snap in sim.get("snapshots", [])
        ],
        "parameters": sim.get("parameters", {}),
    }


def _compute_plot_domain(*simulations):
    y_values = []
    for sim in simulations:
        for snap in sim.get("snapshots", []):
            x = np.asarray(snap.get("x", []), dtype=float)
            pop = np.asarray(snap.get("population", []), dtype=float)
            trait = np.asarray(snap.get("trait", []), dtype=float)
            optimum = np.asarray(snap.get("optimum_heatmap", snap.get("optimum", [])), dtype=float)
            density = np.asarray(snap.get("density", []), dtype=float)
            mean_trait = np.asarray(snap.get("mean_trait_heatmap", snap.get("mean_trait", [])), dtype=float)
            xmask = (x >= 0.0) & (x <= 60.0) if x.size else np.array([], dtype=bool)

            if pop.ndim == 2 and pop.size and trait.size == pop.shape[1] and xmask.size == pop.shape[0]:
                pop_view = pop[xmask, :]
                if pop_view.size:
                    col = np.sum(np.maximum(pop_view, 0.0), axis=0)
                    if np.max(col) > 0:
                        mask = col > np.max(col) * 1e-4
                        if np.any(mask):
                            y_values.extend([float(np.min(trait[mask])), float(np.max(trait[mask]))])
            if optimum.size and xmask.size == optimum.size:
                finite = optimum[xmask & np.isfinite(optimum)]
                if finite.size:
                    y_values.extend([float(np.min(finite)), float(np.max(finite))])
            if density.size and mean_trait.size == density.size and xmask.size == density.size:
                density_view = density[xmask]
                if density_view.size and np.max(density_view) > 0:
                    mask = xmask & (density > np.max(density_view) * 1e-3) & np.isfinite(mean_trait)
                    vals = mean_trait[mask]
                    if vals.size:
                        y_values.extend([float(np.min(vals)), float(np.max(vals))])

    if not y_values:
        ymin, ymax = -5.0, 35.0
    else:
        ymin, ymax = min(y_values), max(y_values)
        span = max(ymax - ymin, 10.0)
        margin = max(2.0, 0.08 * span)
        ymin -= margin
        ymax += margin
    return {"x_min": 0.0, "x_max": 60.0, "trait_min": float(ymin), "trait_max": float(ymax)}


def _transform_for_common_axes(result, model, c, b):
    out = result
    params = out.get("parameters", {})

    if model in ("asexual_deterministic", "asexual_stochastic"):
        for snap in out.get("snapshots", []):
            x = np.asarray(snap["x"], dtype=float)
            trait = np.asarray(snap["trait"], dtype=float)
            origin_x = float(snap.get("plot_x_min", x[0] if x.size else 0.0))
            origin_y = b * origin_x
            snap["x"] = (x - origin_x).astype(np.float32)
            snap["trait"] = (trait - origin_y).astype(np.float32)
            snap["mean_trait"] = (np.asarray(snap["mean_trait"], dtype=float) - origin_y).astype(np.float32)
            # Variance is unchanged by a translation of phenotype coordinates.
            snap["trait_variance"] = np.asarray(snap["trait_variance"], dtype=np.float32)
            snap["optimum"] = (np.asarray(snap["optimum"], dtype=float) - origin_y).astype(np.float32)
            if snap.get("front_left") is not None:
                snap["front_left"] = float(snap["front_left"] - origin_x)
            if snap.get("front_right") is not None:
                snap["front_right"] = float(snap["front_right"] - origin_x)
            snap["grid_x_min"] = float(snap.get("grid_x_min", origin_x) - origin_x)
            snap["grid_x_max"] = float(snap.get("grid_x_max", origin_x + 50.0) - origin_x)
            snap["plot_x_min"] = 0.0
            snap["plot_x_max"] = 60.0
        params["display_coordinates"] = "common physical x/phenotype"
        params["b"] = float(b)
        return out

    sigma = float(params.get("sigma", 1.0))
    rmax = float(params.get("rmax", 1.0))
    VLE = float(params.get("VLE", 0.3))
    L = float(params.get("L", 350.0))
    m = int(params.get("m", 120))
    xscale = np.sqrt(sigma**2 / (2.0 * rmax))
    yscale = np.sqrt(VLE)
    beta = np.sqrt(sigma**2 / (2.0 * VLE * rmax)) * b
    Ktrait = L * beta + 40.0
    dx = L / max(m, 1)

    for snap in out.get("snapshots", []):
        xraw = np.asarray(snap["x"], dtype=float)
        traw = np.asarray(snap["trait"], dtype=float)
        xshift = float(xraw[0] - dx) if xraw.size else 0.0
        raw_left = L / 2.0 + xshift - (60.0 / 3.0) / xscale
        xdisp = xscale * (xraw - raw_left)

        # Preserve the phenotype origin of the previous heatmap.  The Z(x)
        # profile itself is *not* based on this display origin: sexual_model.py
        # already computes it in the exact Zk_export convention of the
        # archived research output.
        raw_opt_left_at_t0 = Ktrait / 2.0 + beta * (raw_left - L / 2.0)
        yorigin = yscale * raw_opt_left_at_t0
        traitdisp = yscale * traw - yorigin

        Zk = np.asarray(snap.get("mean_trait_physical", []), dtype=float)
        optimum_raw = np.asarray(snap.get("optimum_raw", []), dtype=float)
        if Zk.size:
            snap["mean_trait_heatmap"] = (Zk - yorigin).astype(np.float32)
        else:
            # Backward-compatible fallback for old snapshots.
            snap["mean_trait_heatmap"] = np.asarray(snap["mean_trait"], dtype=np.float32)
        if optimum_raw.size:
            snap["optimum_heatmap"] = (yscale * optimum_raw - yorigin).astype(np.float32)
        else:
            snap["optimum_heatmap"] = np.asarray(snap.get("optimum", []), dtype=np.float32)

        snap["x"] = xdisp.astype(np.float32)
        snap["trait"] = traitdisp.astype(np.float32)
        # mean_trait is already Zk_export and trait_variance is already Vk.
        snap["mean_trait"] = np.asarray(snap["mean_trait"], dtype=np.float32)
        snap["trait_variance"] = np.asarray(snap["trait_variance"], dtype=np.float32)
        # Keep optimum in the same centered research-output convention as Z.
        snap["optimum"] = np.asarray(snap.get("optimum_export", snap["optimum_heatmap"]), dtype=np.float32)
        if snap.get("front_left") is not None:
            snap["front_left"] = float(xscale * (snap["front_left"] - raw_left))
        if snap.get("front_right") is not None:
            snap["front_right"] = float(xscale * (snap["front_right"] - raw_left))
        snap["grid_x_min"] = float(xscale * (xraw[0] - raw_left)) if xraw.size else 0.0
        snap["grid_x_max"] = float(xscale * (xraw[-1] - raw_left)) if xraw.size else 60.0
        snap["plot_x_min"] = 0.0
        snap["plot_x_max"] = 60.0

    params.update({
        "b": float(b), "sigma": sigma, "Vs": float(params.get("Vs", 1.0)),
        "VLE": VLE, "rmax": rmax, "display_x_scale": float(xscale),
        "display_trait_scale": float(yscale),
        "Z_convention": "sqrt(VLE)*E[v|x] - Zrefk (same as REC Z column)",
        "V_convention": "VLE*Var[v|x] (same as REC V column)",
        "display_coordinates": "common physical x/phenotype",
    })
    return out


def _simulate_one(model, c, *, preset, tfinal, snapshot_count, K, seed, Vs, b, sigmax, sigmay, VLE, rmax, competition_k):
    if model not in PRESETS:
        raise ValueError(f"Unknown model: {model}")
    if preset not in PRESETS[model]:
        preset = "fast"

    kwargs = dict(PRESETS[model][preset])
    kwargs["tfinal"] = tfinal
    kwargs["snapshot_count"] = snapshot_count

    if model == "sexual":
        base_n = int(kwargs["n"])
        L = 350.0
        # Keep roughly the same dv as the default VLE=0.3, b=0.1 case when
        # changing VLE or b.  This matters particularly for small VLE, for
        # which the computational trait interval becomes much wider.
        beta0 = np.sqrt(1.0 / (2.0 * 0.3)) * 0.1
        beta_current = np.sqrt(1.0 / (2.0 * VLE * rmax)) * b
        K0 = L * beta0 + 40.0
        Kcurrent = L * beta_current + 40.0
        kwargs["n"] = min(2600, max(base_n, int(np.ceil(base_n * Kcurrent / K0))))

    if model in ("asexual_deterministic", "asexual_stochastic"):
        kwargs.update(Vs=Vs, B=b, sigmax=sigmax, sigmay=sigmay, rmax=rmax, competition_k=competition_k)
        if model == "asexual_stochastic":
            kwargs.update(K=K, seed=seed)
    else:
        # VLE is a genuine sexual-model parameter: it enters the dynamics
        # (alpha, beta, trait-domain scaling) as well as the Z/V output scaling.
        kwargs.update(b=b, VLE=VLE, rmax=rmax, competition_k=competition_k)

    result = simulate_model(model, c, **kwargs)
    result = _transform_for_common_axes(result, model, c, b)
    result = _attach_range_markers(result, model)
    return {
        "model": model,
        "c": float(result["c"]),
        "speed_left": result.get("speed_left"),
        "speed_right": result.get("speed_right", 0.0),
        "snapshots": result["snapshots"],
        "parameters": result["parameters"],
    }


def run_compare(data):
    model_a = str(data.get("model_a", "asexual_deterministic"))
    model_b = str(data.get("model_b", "sexual"))
    if model_a not in PRESETS or model_b not in PRESETS:
        raise ValueError("Unknown model")

    c_a = _float_arg(data, "c_a", 0.2, -4.0, 8.0)
    c_b = _float_arg(data, "c_b", 0.2, -4.0, 8.0)
    tfinal = _float_arg(data, "tfinal", 20.0, 2.0, 100.0)
    snapshot_count = _int_arg(data, "snapshot_count", 16, 4, 32)
    preset = str(data.get("preset", "fast"))
    if preset not in ("fast", "balanced", "fine"):
        preset = "fast"

    K_a = _int_arg(data, "K_a", 5000, 1000, 1000000)
    K_b = _int_arg(data, "K_b", 5000, 1000, 1000000)
    seed_a = _int_arg(data, "seed_a", 1, 0, 2**31 - 1)
    seed_b = _int_arg(data, "seed_b", 1, 0, 2**31 - 1)
    Vs = _float_arg(data, "Vs", 1.0, 0.05, 20.0)
    sigmax = _float_arg(data, "sigmax", 1.0, 0.01, 10.0)
    rmax = _float_arg(data, "rmax", 1.0, 0.05, 3.0)
    competition_k = _float_arg(data, "competition_k", 1.0, 0.05, 20.0)
    b_a = _float_arg(data, "b_a", 1.0, 0.0, 3.0)
    b_b = _float_arg(data, "b_b", 1.0, 0.0, 3.0)
    sigmay_a = _float_arg(data, "sigmay_a", 0.64, 0.001, 10.0)
    sigmay_b = _float_arg(data, "sigmay_b", 0.64, 0.001, 10.0)
    VLE_a = _float_arg(data, "VLE_a", 0.3, 0.001, 5.0)
    VLE_b = _float_arg(data, "VLE_b", 0.3, 0.001, 5.0)

    # Sequential by design in Pyodide/WebAssembly.
    sim_a = _simulate_one(
        model_a, c_a, preset=preset, tfinal=tfinal, snapshot_count=snapshot_count,
        K=K_a, seed=seed_a, Vs=Vs, b=b_a, sigmax=sigmax, sigmay=sigmay_a,
        VLE=VLE_a, rmax=rmax, competition_k=competition_k,
    )
    sim_b = _simulate_one(
        model_b, c_b, preset=preset, tfinal=tfinal, snapshot_count=snapshot_count,
        K=K_b, seed=seed_b, Vs=Vs, b=b_b, sigmax=sigmax, sigmay=sigmay_b,
        VLE=VLE_b, rmax=rmax, competition_k=competition_k,
    )

    plot_domain = _compute_plot_domain(sim_a, sim_b)
    return {
        "simulation_a": _serialize_simulation(sim_a, plot_domain),
        "simulation_b": _serialize_simulation(sim_b, plot_domain),
        "plot_domain": plot_domain,
    }


def run_compare_json(payload_json: str) -> str:
    data = json.loads(payload_json)
    return json.dumps(run_compare(data), separators=(",", ":"), allow_nan=False)
