from __future__ import annotations

import numpy as np
from scipy.signal import fftconvolve

from sexual_model import simulate as simulate_sexual





INITIAL_X_FRACTION = 1.0 / 3.0
INITIAL_X_VARIANCE = 10.0

# All three web panels display the same horizontal width.
COMMON_PLOT_X_WIDTH = 60.0

# The asexual computation window is deliberately a little smaller than
# the displayed interval, so the visualisation can show zero population
# outside the numerical grid instead of stretching the grid to the panel.
ASEXUAL_X_DOMAIN = 50.0
ASEXUAL_PLOT_RIGHT_PADDING = COMMON_PLOT_X_WIDTH - ASEXUAL_X_DOMAIN


def _initial_x_center(x):
    """Place the asexual initial population at one third of the plot window."""
    # The 50-unit numerical grid starts at the left edge of the 60-unit
    # displayed interval; the remaining 10 units ahead of the population
    # are visualized with n=0 outside the computation grid.
    return float(x[0] + INITIAL_X_FRACTION * COMMON_PLOT_X_WIDTH)


def _asexual_initial_profile(x, y, B):
    """Shared initial profile for the two asexual demo models.

    Both models start from the same continuous Gaussian density.  The
    stochastic model later converts this density into individual counts,
    so that In/K approximates the deterministic local density.
    """
    x0 = _initial_x_center(x)
    spatial_variance = INITIAL_X_VARIANCE
    trait_variance = 0.74

    Xgrid = x[:, None]
    Ygrid = y[None, :]
    dy = float(y[1] - y[0]) if len(y) > 1 else 1.0

    # Prescribe the x-marginal exactly, then distribute it in trait.
    # This makes the deterministic and stochastic asexual demos start
    # from the same spatial Gaussian even though their computational
    # domains and trait grids are different.
    spatial_profile = np.exp(
        -(x - x0) ** 2 / (2.0 * spatial_variance)
    ) / 5.0

    trait_profile = np.exp(
        -(Ygrid - B * Xgrid) ** 2 / (2.0 * trait_variance)
    )
    trait_mass = np.sum(trait_profile, axis=1, keepdims=True) * dy
    trait_profile /= np.maximum(trait_mass, 1e-300)

    return spatial_profile[:, None] * trait_profile

def _recent_regression_speed(times, positions):
    times = np.asarray(times, dtype=float)
    positions = np.asarray(positions, dtype=float)
    if len(times) < 4:
        return 0.0
    start = len(times) // 2
    if np.allclose(times[start:], times[start]):
        return 0.0
    return float(np.polyfit(times[start:], positions[start:], 1)[0])


def _separable_dispersal(n, gamma_x, gamma_y):
    """Apply the same separable discrete convolution as the research code."""
    temp = fftconvolve(n, gamma_y[None, :], mode="same", axes=1)
    return fftconvolve(temp, gamma_x[:, None], mode="same", axes=0)


def _shift_asexual_window(n, x, y, front, t, c, dx, dy, B):
    """Moving x/y windows from the asexual research programs."""
    Nx, Ny = n.shape
    shift = 5

    if front > 0.3 * x[1] + 0.7 * x[-1]:
        x = x + shift * dx
        old = n
        n = np.zeros_like(old)
        n[:-shift, :] = old[shift:, :]

    if front < 0.5 * x[1] + 0.5 * x[-1]:
        x = x - shift * dx
        old = n
        n = np.zeros_like(old)
        n[shift:, :] = old[:-shift, :]

        # Preserve the special refill rule on the left boundary.
        for i in range(0, shift - 1):
            trait_shift = int(i * B * dx / dy)
            stop = Ny - trait_shift
            if stop <= 0:
                continue
            if trait_shift == 0:
                n[i, :stop] = old[shift, :stop]
            else:
                n[i, :stop] = np.maximum(
                    old[shift, trait_shift:],
                    old[shift, :stop],
                )

    deltay = B * ((x[1] + x[-1]) / 2.0 - c * t) - (y[-1] + y[1]) / 2.0

    if deltay > 1:
        yshift = 5
        y = y + yshift * dy
        old = n
        n = np.zeros_like(old)
        n[:, :-yshift] = old[:, yshift:]

    if deltay < -1:
        yshift = -5
        y = y + yshift * dy
        old = n
        n = np.zeros_like(old)
        n[:, -yshift:] = old[:, :yshift]

    return n, x, y


def _snapshot_asexual(n, x, y, time, c, density_kind, B, K=None):
    """Create a snapshot, padding the x-grid with zeros for plotting.

    The numerical asexual grid has width 50, while every web panel displays
    60 x-units.  We explicitly extend the returned population with zeros on
    both sides.  This is especially useful for the stochastic model: outside
    the computational grid there are simply no individuals, rather than a
    stretched image.
    """
    if density_kind == "deterministic":
        dy = float(y[1] - y[0]) if len(y) > 1 else 1.0
        density = np.sum(n, axis=1) * dy
    else:
        density = np.sum(n, axis=1)

    den = np.sum(n, axis=1)
    mean_trait = np.zeros_like(den, dtype=float)
    trait_variance = np.zeros_like(den, dtype=float)
    mask = den > 1e-14
    mean_trait[mask] = np.sum(n[mask, :] * y[None, :], axis=1) / den[mask]
    trait_variance[mask] = (
        np.sum(n[mask, :] * (y[None, :] - mean_trait[mask, None]) ** 2, axis=1)
        / den[mask]
    )

    if density_kind == "deterministic":
        inds = np.flatnonzero(density > np.max(density) / 2.0) if np.max(density) > 0 else np.array([], dtype=int)
    else:
        inds = np.flatnonzero(density > (K / 100.0)) if K is not None else np.array([], dtype=int)
    front = float(x[inds[-1]]) if inds.size else None

    # Exact displayed interval.  Its width is common to all models.
    plot_x_min = float(x[0])
    plot_x_max = plot_x_min + COMMON_PLOT_X_WIDTH

    # Extend the sampled x-grid by multiples of dx and fill with zeros.
    dx = float(x[1] - x[0]) if len(x) > 1 else 1.0
    nleft = int(np.ceil(max(0.0, x[0] - plot_x_min) / dx))
    nright = int(np.ceil(max(0.0, plot_x_max - x[-1]) / dx))
    left_x = x[0] - dx * np.arange(nleft, 0, -1)
    right_x = x[-1] + dx * np.arange(1, nright + 1)
    x_plot = np.concatenate((left_x, x, right_x))

    population_plot = np.pad(np.asarray(n), ((nleft, nright), (0, 0)), mode="constant")
    density_plot = np.pad(np.asarray(density), (nleft, nright), mode="constant")
    mean_trait_plot = np.pad(np.asarray(mean_trait), (nleft, nright), mode="constant")
    trait_variance_plot = np.pad(np.asarray(trait_variance), (nleft, nright), mode="constant")
    optimum_plot = B * (x_plot - c * time)

    return {
        "time": float(time),
        "x": np.asarray(x_plot, dtype=np.float32),
        "trait": np.asarray(y, dtype=np.float32).copy(),
        "population": np.asarray(population_plot, dtype=np.float32),
        "density": np.asarray(density_plot, dtype=np.float32),
        "mean_trait": np.asarray(mean_trait_plot, dtype=np.float32),
        "trait_variance": np.asarray(trait_variance_plot, dtype=np.float32),
        "optimum": np.asarray(optimum_plot, dtype=np.float32),
        "front_right": front,
        "front_left": None,
        "plot_x_min": plot_x_min,
        "plot_x_max": plot_x_max,
        "grid_x_min": float(x[0]),
        "grid_x_max": float(x[-1]),
    }


def simulate_asexual_deterministic(
    c: float,
    *,
    nx: int = 120,
    ny: int = 160,
    tfinal: float = 40.0,
    dt: float = 1.0,
    snapshot_count: int = 16,
    Vs: float = 1.0,
    B: float = 1.2,
    sigmax: float = 1.0,
    sigmay: float = 0.8**2,
    rmax: float = 1.0,
    competition_k: float = 1.0,
):
    A = 1.0 / (2.0 * Vs)
    # Use the same computation window as the stochastic model.
    # It is slightly larger than before, which improves visualization
    # while keeping the two asexual demos directly comparable.
    X = ASEXUAL_X_DOMAIN
    Y = 4.0 + X * B

    x = np.linspace(0.0, X, nx)
    y = np.linspace(0.0, Y, ny)
    dx = X / (nx - 1)
    dy = Y / (ny - 1)
    y = y - ((y[1] + y[-1]) / 2.0 - B * (x[1] + x[-1]) / 2.0)

    # Same continuous initial density as in the stochastic asexual demo.
    n = _asexual_initial_profile(x, y, B)

    xconv = np.linspace(-X, X, 2 * nx - 1)
    gamma_x = np.exp(-xconv**2 / (2.0 * sigmax * dt))
    gamma_x /= np.sum(gamma_x)
    yconv = np.linspace(-Y, Y, 2 * ny - 1)
    gamma_y = np.exp(-yconv**2 / (2.0 * sigmay * dt))
    gamma_y /= np.sum(gamma_y)

    snapshot_targets = np.linspace(0.0, tfinal, max(2, snapshot_count))
    snapshots = [_snapshot_asexual(n, x, y, 0.0, c, "deterministic", B)]
    next_snapshot = 1

    times = []
    fronts = []
    total_mass = []

    t = 0.0
    while t < tfinal - 1e-12:
        dt_step = min(dt, tfinal - t)
        ntemp = _separable_dispersal(n, gamma_x, gamma_y)

        In = np.sum(ntemp, axis=1) * dy
        maladaptation = y[None, 1:-1] - B * (x[1:-1, None] - c * t)
        death_rate = A * maladaptation**2 + In[1:-1, None] / competition_k
        n[1:-1, 1:-1] = ntemp[1:-1, 1:-1] * np.exp((rmax - death_rate) * dt_step)

        # Boundary conditions from the research code.
        n[0, 1:-1] = 0.0
        n[-1, 1:-1] = 0.0

        half_inds = np.flatnonzero(In > np.max(In) / 2.0) if np.max(In) > 0 else np.array([], dtype=int)
        front = float(x[half_inds[-1]]) if half_inds.size else float(x[0])

        times.append(t)
        fronts.append(front)
        total_mass.append(float(np.sum(In) * dx))

        n, x, y = _shift_asexual_window(n, x, y, front, t, c, dx, dy, B)

        t += dt_step
        while next_snapshot < len(snapshot_targets) and t + 1e-12 >= snapshot_targets[next_snapshot]:
            snapshots.append(_snapshot_asexual(n, x, y, t, c, "deterministic", B))
            next_snapshot += 1

    speed = _recent_regression_speed(times, fronts)

    return {
        "model": "asexual_deterministic",
        "c": float(c),
        "speed_right": speed,
        "speed_left": None,
        "times": np.asarray(times, dtype=float),
        "right_front": np.asarray(fronts, dtype=float),
        "left_front": None,
        "total_mass": np.asarray(total_mass, dtype=float),
        "snapshots": snapshots,
        "parameters": {"nx": nx, "ny": ny, "tfinal": tfinal, "dt": dt, "Vs": Vs, "B": B, "sigmax": sigmax, "sigmay": sigmay, "rmax": rmax, "k": competition_k, "initial_x_center": _initial_x_center(np.linspace(0.0, X, nx)), "initial_x_fraction": INITIAL_X_FRACTION, "initial_x_variance": INITIAL_X_VARIANCE, "plot_x_width": COMMON_PLOT_X_WIDTH, "computation_x_width": ASEXUAL_X_DOMAIN},
    }


def simulate_asexual_stochastic(
    c: float,
    *,
    K: int = 5000,
    seed: int = 1,
    nx: int = 120,
    ny: int = 160,
    tfinal: float = 20.0,
    dt: float = 0.2,
    snapshot_count: int = 16,
    Vs: float = 1.0,
    B: float = 1.2,
    sigmax: float = 1.0,
    sigmay: float = 0.8**2,
    rmax: float = 1.0,
    competition_k: float = 1.0,
):
    A = 1.0 / (2.0 * Vs)
    X = ASEXUAL_X_DOMAIN
    Y = 4.0 + X * B

    x = np.linspace(0.0, X, nx)
    y = np.linspace(0.0, Y, ny)
    dx = X / (nx - 1)
    dy = Y / (ny - 1)
    y = y - ((y[1] + y[-1]) / 2.0 - B * (x[1] + x[-1]) / 2.0)

    # Start from exactly the same continuous profile as the deterministic
    # asexual demo, then convert density to integer individual counts.
    # With this scaling, In/K approximates the deterministic integral
    # sum_y n(x,y) dy.
    initial_density = _asexual_initial_profile(x, y, B)

    # Preserve the deterministic x-marginal as closely as possible in
    # integer counts.  We first choose the exact row total, then distribute
    # those individuals over trait with the largest-remainder method.
    row_density = np.sum(initial_density, axis=1) * dy
    row_counts = np.rint(K * row_density).astype(np.int64)
    trait_prob = initial_density * dy / np.maximum(row_density[:, None], 1e-300)
    n = np.zeros_like(initial_density, dtype=np.int64)
    for i in range(nx):
        raw = row_counts[i] * trait_prob[i]
        base = np.floor(raw).astype(np.int64)
        remainder = int(row_counts[i] - np.sum(base))
        if remainder > 0:
            frac = raw - base
            add = np.argpartition(frac, -remainder)[-remainder:]
            base[add] += 1
        n[i] = base

    xconv = np.linspace(-X, X, 2 * nx - 1)
    gamma_x = np.exp(-xconv**2 / (2.0 * sigmax * dt))
    gamma_x /= np.sum(gamma_x)
    yconv = np.linspace(-Y, Y, 2 * ny - 1)
    gamma_y = np.exp(-yconv**2 / (2.0 * sigmay * dt))
    gamma_y /= np.sum(gamma_y)

    rng = np.random.default_rng(seed)

    snapshot_targets = np.linspace(0.0, tfinal, max(2, snapshot_count))
    snapshots = [_snapshot_asexual(n, x, y, 0.0, c, "stochastic", B, K=K)]
    next_snapshot = 1

    times = []
    fronts = []
    total_mass = []

    t = 0.0
    while t < tfinal - 1e-12:
        dt_step = min(dt, tfinal - t)
        ntemp = _separable_dispersal(n.astype(float), gamma_x, gamma_y)

        In = np.sum(ntemp, axis=1)
        maladaptation = y[None, 1:-1] - B * (x[1:-1, None] - c * t)
        rate = A * maladaptation**2 + In[1:-1, None] / (competition_k * K)
        p_survival = np.exp(-rate * dt_step)
        p_birth = np.exp((rmax - rate) * dt_step) - p_survival

        if np.max(p_birth) > 1.0 + 1e-12:
            raise ValueError("rmax and dt give a birth probability above 1; choose a smaller rmax or a finer time step")
        p_birth = np.clip(p_birth, 0.0, 1.0)

        trials = np.maximum(ntemp[1:-1, 1:-1].astype(np.int64), 0)
        n[1:-1, 1:-1] = (
            rng.binomial(trials, p_survival)
            + rng.binomial(trials, p_birth)
        )

        n[0, 1:-1] = 0
        n[-1, 1:-1] = 0

        front_inds = np.flatnonzero(In > K / 100.0)
        front = float(x[front_inds[-1]]) if front_inds.size else float(x[0])

        times.append(t)
        fronts.append(front)
        total_mass.append(float(np.sum(In) * dx))

        n, x, y = _shift_asexual_window(n, x, y, front, t, c, dx, dy, B)

        t += dt_step
        while next_snapshot < len(snapshot_targets) and t + 1e-12 >= snapshot_targets[next_snapshot]:
            snapshots.append(_snapshot_asexual(n, x, y, t, c, "stochastic", B, K=K))
            next_snapshot += 1

    speed = _recent_regression_speed(times, fronts)

    return {
        "model": "asexual_stochastic",
        "c": float(c),
        "speed_right": speed,
        "speed_left": None,
        "times": np.asarray(times, dtype=float),
        "right_front": np.asarray(fronts, dtype=float),
        "left_front": None,
        "total_mass": np.asarray(total_mass, dtype=float),
        "snapshots": snapshots,
        "parameters": {
            "K": int(K), "seed": int(seed), "nx": nx, "ny": ny,
            "tfinal": tfinal, "dt": dt, "Vs": Vs, "B": B,
            "sigmax": sigmax, "sigmay": sigmay,
            "rmax": rmax, "k": competition_k,
            "initial_x_center": _initial_x_center(np.linspace(0.0, X, nx)),
            "initial_x_fraction": INITIAL_X_FRACTION,
            "initial_x_variance": INITIAL_X_VARIANCE,
            "plot_x_width": COMMON_PLOT_X_WIDTH,
            "computation_x_width": ASEXUAL_X_DOMAIN,
        },
    }


def simulate_model(model: str, c: float, **kwargs):
    if model == "asexual_deterministic":
        return simulate_asexual_deterministic(c, **kwargs)
    if model == "asexual_stochastic":
        return simulate_asexual_stochastic(c, **kwargs)
    if model == "sexual":
        # Run the browser-adapted sexual solver.  Its scientific equations
        # follow the research code; the browser layer adds presentation metadata.
        result = simulate_sexual(c, keep_snapshots=True, **kwargs)
        result["model"] = "sexual"
        left = np.asarray(result["left_front"], dtype=float)
        right = np.asarray(result["right_front"], dtype=float)
        times = np.asarray(result["times"], dtype=float)

        L = float(result.get("parameters", {}).get("L", 350.0))
        m = int(result.get("parameters", {}).get("m", 120))
        dx = L / max(m, 1)
        plot_width = 60.0

        for snap in result["snapshots"]:
            st = float(snap["time"])
            if len(times):
                idx = int(np.argmin(np.abs(times - st)))
                snap["front_left"] = float(left[idx]) if len(left) else None
                snap["front_right"] = float(right[idx]) if len(right) else None
            else:
                snap["front_left"] = None
                snap["front_right"] = None

            # The original solver returns x+xshift.  Infer xshift from the first
            # grid point and choose a display-only 60-unit window.  The original
            # Gaussian center L/2 is at one third of this displayed interval.
            xarr = np.asarray(snap["x"], dtype=float)
            shift = float(xarr[0] - dx) if xarr.size else 0.0
            center = L / 2.0 + shift
            snap["plot_x_min"] = center - plot_width / 3.0
            snap["plot_x_max"] = center + 2.0 * plot_width / 3.0
            snap["grid_x_min"] = float(xarr[0]) if xarr.size else 0.0
            snap["grid_x_max"] = float(xarr[-1]) if xarr.size else L

        params = result.setdefault("parameters", {})
        params.update({
            "Vs": float(params.get("Vs", 1.0)),
            "sigma": float(params.get("sigma", 1.0)),
            "b": float(params.get("b", 0.1)),
            "VLE": float(params.get("VLE", 0.3)),
            "initial_x_center": L / 2.0,
            "initial_x_variance": 10.0,
            "plot_x_width": plot_width,
            "sexual_solver": "browser_adaptation_of_research_model",
        })
        return result
    raise ValueError(f"Unknown model: {model}")
