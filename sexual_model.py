from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
from scipy.fft import irfft, next_fast_len, rfft
from scipy.linalg import solveh_banded


def octave_round(x):
    """MATLAB/Octave rounding: half-integers are rounded away from zero."""
    x = np.asarray(x)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def _reproduction_step_fft(
    nd: np.ndarray,
    Ik: np.ndarray,
    im: np.ndarray,
    ip: np.ndarray,
    di_local: np.ndarray,
    Q: np.ndarray,
    q_center: int,
    dv: float,
    gamma: float,
    cache: dict,
) -> np.ndarray:
    """Compute the birth term using batched FFTs for rows with equal support size."""
    birth = np.zeros_like(nd)
    active = np.flatnonzero(Ik > 1e-5)
    if active.size == 0:
        return birth

    birth[active, :] = nd[active, :]

    for di_i in np.unique(di_local[active]):
        di_i = int(di_i)
        rows = active[di_local[active] == di_i]
        im0 = im[rows] - 1

        cols = im0[:, None] + np.arange(di_i)[None, :]
        f_local = nd[rows[:, None], cols]

        cached = cache.get(di_i)
        if cached is None:
            Qtemp = Q[q_center - di_i : q_center + di_i + 1]
            fft_length = next_fast_len(4 * di_i - 1)
            Qhat = rfft(Qtemp, n=fft_length)
            cached = (fft_length, Qhat)
            cache[di_i] = cached

        fft_length, Qhat = cached
        Fhat = rfft(f_local, n=fft_length, axis=1)
        tric = irfft((Fhat * Fhat) * Qhat[None, :], n=fft_length, axis=1)
        tric = tric[:, : 4 * di_i - 1] * dv**2

        # MATLAB: truc(j) = tric(2*j + di - 1), j=1,...,di
        truc = tric[:, di_i : 3 * di_i - 1 : 2]
        birth[rows[:, None], cols] = gamma / Ik[rows, None] * truc

    return birth


def simulate(
    c: float,
    *,
    m: int = 120,
    n: int = 256,
    tfinal: float = 40.0,
    L: float = 350.0,
    dt: float = 0.2,
    b: float = 0.1,
    VLE: float = 0.3,
    rmax: float = 1.0,
    competition_k: float = 1.0,
    snapshot_count: int = 24,
    keep_snapshots: bool = True,
    progress_callback=None,
):
    """Run one simulation of the sexual model.

    The equations and splitting follow the research code, but the implementation is
    adapted for interactive use: vectorized selection, batched FFT reproduction,
    and a banded implicit diffusion solve.
    """

    # Biological parameters (same default values as the research script)
    sigma = 1.0
    Vs = 1.0
    eta = 1.0
    b = float(b)
    VLE = float(VLE)
    rmax = float(rmax)
    competition_k = float(competition_k)
    if VLE <= 0:
        raise ValueError("VLE must be positive")
    if rmax <= 0:
        raise ValueError("rmax must be positive")
    if competition_k <= 0:
        raise ValueError("k must be positive")

    alpha = VLE / (rmax * Vs)
    beta = np.sqrt(sigma**2 / (2.0 * VLE * rmax)) * b
    gamma = eta
    birth_gamma = gamma + rmax - 1.0
    nu = np.sqrt(2.0 / (rmax * sigma**2)) * c

    K = L * beta + 40.0

    dx = L / m
    x = np.arange(1, m + 1, dtype=float) * dx
    dv = K / n
    v = np.arange(1, n + 1, dtype=float) * dv

    # Initial condition
    Zopt0 = K / 2.0 + beta * (x - L / 2.0)
    prop = 0.5
    spatial_profile = np.exp(-(np.abs(x - L * prop)) ** 2 / 20.0)[:, None]
    trait_profile = np.exp(-(v[None, :] - Zopt0[:, None]) ** 2 / 2.0)
    nd = (0.05 / dv) * spatial_profile * trait_profile

    # Fixed output references, following export_final_profiles() in the
    # research code.  They are used only to express Z(x) in the same
    # centered phenotype coordinate as the archived REC*.csv files.
    I0_output = np.sum(np.abs(nd), axis=1) * dv
    I0_output_max = float(np.max(I0_output))
    x0_output = float(np.sum(x * I0_output) / np.sum(I0_output))
    Zref_output = float(K / 2.0 + beta * (x0_output - L / 2.0))
    sqrt_VLE = float(np.sqrt(VLE))
    Zrefk = sqrt_VLE * Zref_output

    # Reproduction kernel
    v2 = np.concatenate((-v[::-1], [0.0], v))
    Q = np.exp(-(v2**2) / 4.0) / np.sqrt(np.pi)
    q_center = len(Q) // 2
    reproduction_fft_cache = {}

    # Implicit spatial diffusion matrix in symmetric banded form
    rdiff = dt / (dx * dx)
    Mdiff = m - 2
    diffusion_ab = np.zeros((2, Mdiff))
    diffusion_ab[0, 1:] = -rdiff
    diffusion_ab[1, :] = 1.0 + 2.0 * rdiff

    interp = dx * beta / dv
    deltaj = int(np.floor(interp) + 1)
    thresholdI = 0.03
    spatial_shift_size = int(octave_round(len(x) / 5.0))
    support_half_width = (10.0 + np.sqrt(2.0 / alpha)) / dv

    # Moving window and time records
    xshift = 0.0
    yshift = 0.0
    temps = 0.0
    kpt = 0
    kplot = max(1, round(tfinal))

    tx = []
    xt = []
    xtm = []
    Zt = []
    Ztm = []
    total_mass_history = []

    snapshots = []
    if keep_snapshots and snapshot_count > 0:
        snapshot_targets = np.linspace(0.0, tfinal, snapshot_count)
        next_snapshot = 0
    else:
        snapshot_targets = np.array([])
        next_snapshot = 0

    def save_snapshot(time_value):
        nonlocal next_snapshot
        abs_nd = np.abs(nd)
        den = np.sum(abs_nd, axis=1)
        density = den * dv

        # Compute the phenotype moments exactly in the convention used by
        # export_final_profiles() of the research code.  The computational
        # trait coordinate v is scaled by sqrt(VLE), hence its variance by VLE.
        # In particular, Vk is expected to be of order VLE when the variance
        # in v is of order one.
        raw_mean = np.zeros_like(density)
        Zk = np.zeros_like(density)
        Vk = np.zeros_like(density)
        mask = density > 1e-6
        weighted_trait = np.sum(abs_nd * v[None, :], axis=1)
        raw_mean[mask] = weighted_trait[mask] / den[mask]
        Zk[mask] = sqrt_VLE * raw_mean[mask]

        centered_trait = v[None, :] - Zk[:, None] / sqrt_VLE
        variance_numerator = np.sum(abs_nd * centered_trait**2, axis=1)
        Vk[mask] = VLE * variance_numerator[mask] / den[mask]

        # Z exported by the research code is centered on the initial optimum
        # at the initial population center of mass.  Keep the same zero
        # sentinel where the population is negligible.
        Z_export = np.zeros_like(Zk)
        Z_export[mask] = Zk[mask] - Zrefk

        optimum_raw = K / 2.0 + beta * (x + xshift - L / 2.0 - nu * time_value)
        optimum_export = sqrt_VLE * optimum_raw - Zrefk

        snapshots.append(
            {
                "time": float(time_value),
                "x": (x + xshift).astype(np.float32),
                "trait": (v - yshift).astype(np.float32),
                "population": nd.astype(np.float32).copy(),
                "density": density.astype(np.float32),
                # Values used by the Z and V profile plots: same convention as
                # Zk_export and Vk in the archived REC*.csv output.
                "mean_trait": Z_export.astype(np.float32),
                "trait_variance": Vk.astype(np.float32),
                # Physical (scaled but not centered) mean, retained so the
                # heatmap overlay can preserve the previous visual origin.
                "mean_trait_physical": Zk.astype(np.float32),
                "optimum_raw": optimum_raw.astype(np.float32),
                "optimum_export": optimum_export.astype(np.float32),
            }
        )
        next_snapshot += 1

    if keep_snapshots and snapshot_targets.size:
        save_snapshot(0.0)

    max_steps = int(np.ceil(tfinal / dt)) + 2
    step = 0

    while temps < tfinal and step < max_steps:
        Zopt = K / 2.0 + beta * (x + xshift - L / 2.0 - nu * temps)
        dt_step = min(dt, tfinal - temps)

        Ik = np.sum(np.abs(nd), axis=1) * dv

        iopt = octave_round((Zopt - v[0]) / (v[-1] - v[0]) * n).astype(int)
        im = np.maximum(octave_round(iopt - support_half_width).astype(int), 1)
        ip = np.minimum(octave_round(iopt + support_half_width).astype(int), n)
        di_local = ip - im + 1

        trait_difference = v[None, :] - yshift - Zopt[:, None]
        selection_cost = np.minimum(alpha / 2.0 * trait_difference**2, 20000.0)
        total_rate = selection_cost + (gamma - 1.0) + Ik[:, None] / competition_k

        ndd = nd * np.exp(-dt_step * total_rate)
        birth = _reproduction_step_fft(
            nd, Ik, im, ip, di_local, Q, q_center, dv, birth_gamma, reproduction_fft_cache
        )
        ndd += dt_step * birth

        # Implicit diffusion in x for all trait columns at once
        nd[1:-1, 1:-1] = solveh_banded(
            diffusion_ab,
            ndd[1:-1, 1:-1],
            lower=False,
            overwrite_ab=False,
            overwrite_b=False,
            check_finite=False,
        )

        # Boundary conditions
        nd[-1, 1:-1] = nd[-2, 1:-1]
        jboundary = np.arange(1, n - 10)
        nd[0, jboundary] = (
            interp / deltaj * nd[1, jboundary + deltaj]
            + (deltaj - interp) / deltaj * nd[1, jboundary]
        )
        nd[0, n - 11 : n - 1] = nd[1, n - 11 : n - 1]
        nd[1:-1, 0] = 0.0
        nd[1:-1, -1] = 0.0
        nd[0, 0] = 0.5 * nd[1, 0] + 0.5 * nd[0, 1]
        nd[0, -1] = 0.5 * nd[1, -1] + 0.5 * nd[0, -2]
        nd[-1, 0] = 0.0
        nd[-1, -1] = 0.5 * nd[-1, -2] + 0.5 * nd[-2, -1]

        # Front positions and speeds, as in the research script
        tx.append(temps)
        total_mass_history.append(float(np.linalg.norm(Ik, ord=1)))

        k_idx = m - 1
        while Ik[k_idx] < thresholdI and k_idx > 0:
            k_idx -= 1
        k_idx = min(k_idx, len(Ik) - 2)
        k1 = k_idx + 1
        denom = Ik[k_idx] - Ik[k_idx + 1]
        thetaI = 0.0 if abs(denom) < 1e-15 else (thresholdI - Ik[k_idx + 1]) / denom
        thetaI = float(np.clip(thetaI, 0.0, 1.0))
        XtI = thetaI * k1 * dx + (1.0 - thetaI) * (k1 + 1) * dx + xshift
        xt.append(XtI)

        zk = np.sum(nd[k_idx, :] * v) / max(np.sum(nd[k_idx, :]), 1e-300)
        zk1 = np.sum(nd[k_idx + 1, :] * v) / max(np.sum(nd[k_idx + 1, :]), 1e-300)
        Zt.append(thetaI * zk + (1.0 - thetaI) * zk1 - yshift)

        km_idx = 1
        while Ik[km_idx] < thresholdI and km_idx < m - 1:
            km_idx += 1
        km1 = km_idx + 1
        denom = Ik[km_idx] - Ik[km_idx - 1]
        thetaI = 0.0 if abs(denom) < 1e-15 else (thresholdI - Ik[km_idx - 1]) / denom
        thetaI = float(np.clip(thetaI, 0.0, 1.0))
        XtIm = thetaI * km1 * dx + (1.0 - thetaI) * (km1 - 1) * dx + xshift
        xtm.append(XtIm)

        zk = np.sum(nd[km_idx, :] * v) / max(np.sum(nd[km_idx, :]), 1e-300)
        zk1 = np.sum(nd[km_idx - 1, :] * v) / max(np.sum(nd[km_idx - 1, :]), 1e-300)
        Ztm.append(thetaI * zk + (1.0 - thetaI) * zk1 - yshift)

        # Moving spatial window
        if temps > tfinal * (kpt / kplot) or temps > tfinal - 2.0 * dt_step:
            kpt += 1
            Mx = np.sum(nd, axis=1)
            total = np.sum(Mx)
            if total > 0:
                xcm = np.dot(x, Mx) / total
                if xcm > 0.6 * x[-1]:
                    xshift += x[spatial_shift_size - 1]
                    shift = spatial_shift_size
                    nd[:-shift, :] = nd[shift:, :].copy()
                    nd[-shift:, :] = 0.0
                    nd[:, :shift] = 0.0

        temps += dt_step
        step += 1

        if keep_snapshots and next_snapshot < len(snapshot_targets):
            while next_snapshot < len(snapshot_targets) and temps + 1e-12 >= snapshot_targets[next_snapshot]:
                save_snapshot(temps)

        if progress_callback is not None:
            progress_callback(min(1.0, temps / tfinal))

    tx = np.asarray(tx)
    xt = np.asarray(xt)
    xtm = np.asarray(xtm)
    Zt = np.asarray(Zt)
    Ztm = np.asarray(Ztm)

    def recent_speed(position):
        if len(position) < 2:
            return 0.0
        idx = int(np.floor(len(position) * 0.9))
        idx = min(idx, len(position) - 1)
        delta_t = tx[-1] - tx[idx]
        if delta_t <= 0:
            return 0.0
        return float((position[-1] - position[idx]) / delta_t)

    speed_right = recent_speed(xt)
    speed_left = recent_speed(xtm)
    trait_speed_right = recent_speed(Zt)
    trait_speed_left = recent_speed(Ztm)

    return {
        "c": float(c),
        "speed_right": speed_right * np.sqrt(sigma**2 / (2.0 * rmax)),
        "speed_left": speed_left * np.sqrt(sigma**2 / (2.0 * rmax)),
        "trait_speed_right": trait_speed_right * np.sqrt(VLE * Vs),
        "trait_speed_left": trait_speed_left * np.sqrt(VLE * Vs),
        "times": tx,
        "right_front": xt,
        "left_front": xtm,
        "total_mass": np.asarray(total_mass_history),
        "snapshots": snapshots,
        "parameters": {
            "m": m,
            "n": n,
            "tfinal": tfinal,
            "dt": dt,
            "L": L,
            "VLE": VLE,
            "b": b,
            "sigma": sigma,
            "Vs": Vs,
            "rmax": rmax,
            "k": competition_k,
            "eta": eta,
            "I0_output_max": I0_output_max,
            "x0_output": x0_output,
            "Zref_output": Zref_output,
            "Zrefk": Zrefk,
        },
    }


def speed_sweep(
    c_values,
    *,
    m=100,
    n=192,
    tfinal=30.0,
    L=350.0,
    dt=0.2,
    progress_callback=None,
):
    c_values = np.asarray(c_values, dtype=float)
    speed_left = np.zeros_like(c_values)
    speed_right = np.zeros_like(c_values)

    for k, c in enumerate(c_values):
        result = simulate(
            float(c),
            m=m,
            n=n,
            tfinal=tfinal,
            L=L,
            dt=dt,
            snapshot_count=0,
            keep_snapshots=False,
        )
        speed_left[k] = result["speed_left"]
        speed_right[k] = result["speed_right"]
        if progress_callback is not None:
            progress_callback((k + 1) / len(c_values))

    return {
        "c": c_values,
        "speed_left": speed_left,
        "speed_right": speed_right,
    }
