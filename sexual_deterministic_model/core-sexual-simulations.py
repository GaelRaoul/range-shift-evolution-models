#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np
from scipy.fft import irfft, next_fast_len, rfft
from scipy.linalg import solveh_banded


# ============================================================
# Utilities
# ============================================================

def octave_round(x):
    """
    Reproduce MATLAB/Octave round() semantics:
    half-integers are rounded away from zero.
    """
    x = np.asarray(x)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def write_csv_row(file_handle, values):
    """Write one numerical CSV row."""
    file_handle.write(",".join(f"{float(v):.16g}" for v in values) + "\n")


# ============================================================
# One simulation for a fixed climate-change speed c
# ============================================================

def run_simulation(
    c,
    m=600,
    n=2**10,
    tfinal=80,
    L=300.0,
    dt=0.1,
    data_dir=Path("data"),
):

    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # Biological parameters
    # --------------------------------------------------------

    sigma = 1.0
    rmax = 1.0
    Vs = 1.0
    eta = 1.0

    b = 0.1
    VLE = 1
    c = float(c)

    alpha = VLE / (rmax * Vs)
    beta = np.sqrt(sigma**2 / (2.0 * VLE * rmax)) * b
    gamma = eta
    nu = np.sqrt(2.0 / (rmax * sigma**2)) * c

    r_ast = 1.0 - VLE / (2.0 * rmax * Vs)
    AKB = VLE / (r_ast * rmax * Vs)
    BKB = b * sigma / (2.0 * r_ast * rmax * np.sqrt(Vs))
    climatespeed = np.sqrt(2.0 / (rmax * r_ast * sigma**2)) * c

    K = L * beta + 40.0


    # --------------------------------------------------------
    # Spatial and trait grids
    # --------------------------------------------------------

    dx = L / m
    x = np.arange(1, m + 1, dtype=float) * dx

    dv = K / n
    v = np.arange(1, n + 1, dtype=float) * dv


    # --------------------------------------------------------
    # Initial population
    # --------------------------------------------------------

    Zopt0 = K / 2.0 + beta * (x - L / 2.0)
    Zopt = Zopt0.copy()
    prop = 0.5

    x0ini = L * prop
    y0ini = K / 2.0 + beta * (x0ini - L / 2.0)

    spatial_profile = np.exp(
        -(np.abs(x - x0ini))**2 / 20.0
    )[:, None]

    trait_profile = np.exp(
        -(v[None, :] - Zopt[:, None])**2 / 2.0
    )

    nd = (0.05 / dv) * spatial_profile * trait_profile
    nd=nd/max(np.sum(np.abs(nd), axis=1) * dv)
    nd0 = nd.copy()


    # --------------------------------------------------------
    # Reproduction kernel
    # --------------------------------------------------------

    v2 = np.concatenate((-v[::-1], [0.0], v))
    Q = np.exp(-(v2**2) / 4.0) / np.sqrt(np.pi)
    q_center = len(Q) // 2

    # Cache FFT information for each support size di:
    #     di -> (fft_length, FFT(Qtemp))
    reproduction_fft_cache = {}


    # --------------------------------------------------------
    # Implicit spatial diffusion
    # --------------------------------------------------------

    # A is the symmetric tridiagonal matrix
    #
    # diag(A)     = 1 + 2 rdiff
    # offdiag(A)  = -rdiff
    #
    # scipy.linalg.solveh_banded stores the upper diagonal in
    # row 0 and the main diagonal in row 1.

    rdiff = dt / (dx * dx)
    Mdiff = m - 2

    diffusion_ab = np.zeros((2, Mdiff))
    diffusion_ab[0, 1:] = -rdiff
    diffusion_ab[1, :] = 1.0 + 2.0 * rdiff


    # --------------------------------------------------------
    # Constants reused in the time loop
    # --------------------------------------------------------

    interp = dx * beta / dv
    deltaj = int(np.floor(interp) + 1)

    thresholdI = 0.03
    spatial_shift_size = int(octave_round(len(x) / 5.0))

    support_half_width = (10.0 + np.sqrt(2.0 / alpha)) / dv

    kplot = max(1, round(tfinal))


    # --------------------------------------------------------
    # Recording arrays
    # --------------------------------------------------------

    max_steps = int(np.ceil(tfinal / dt)) + 2

    tx = np.zeros(max_steps)
    xt = np.zeros(max_steps)
    xtm = np.zeros(max_steps)
    Zt = np.zeros(max_steps)
    Ztm = np.zeros(max_steps)
    Ikrec = np.zeros(max_steps)

    speedrec = np.zeros(max_steps)
    speedrecm = np.zeros(max_steps)
    speedrecz = np.zeros(max_steps)
    speedreczm = np.zeros(max_steps)
    maxrec = np.zeros(max_steps)

    Ik = np.zeros(m)

    step = 0
    kpt = 0
    temps = 0.0

    xshift = 0.0
    yshift = 0.0


    # --------------------------------------------------------
    # Output file for this climate-change speed
    # --------------------------------------------------------

    tfinal_label = int(tfinal) if float(tfinal).is_integer() else tfinal

    filename3 = data_dir / (
        f"eachspeedVLEi{VLE}"
        f"b{b}"
        f"tfinal{tfinal_label}"
        f"c{c}.csv"
    )

    with filename3.open("w") as fid3:

        fid3.write(
            "t,c,vS,vN,vZS,vZN,xS,xN,Ik1,Ik2,"
            "lambda,VLE,b,C,tfinal\n"
        )


        # ====================================================
        # Time loop
        # ====================================================

        while temps < tfinal:

            s = step

            Zopt = (
                K / 2.0
                + beta * (x + xshift - L / 2.0 - nu * temps)
            )

            dt_step = min(dt, tfinal - temps)


            # ------------------------------------------------
            # Local population sizes
            # ------------------------------------------------

            # Equivalent to MATLAB norm(nd(i,:),1)*dv row by row.
            Ik[:] = np.sum(np.abs(nd), axis=1) * dv


            # ------------------------------------------------
            # Reproduction support indices
            # ------------------------------------------------

            iopt = octave_round(
                (Zopt - v[0]) / (v[-1] - v[0]) * n
            ).astype(int)

            im = np.maximum(
                octave_round(iopt - support_half_width).astype(int),
                1,
            )

            ip = np.minimum(
                octave_round(iopt + support_half_width).astype(int),
                n,
            )

            di_local = ip - im + 1


            # ------------------------------------------------
            # Selection / mortality term
            # ------------------------------------------------

            trait_difference = (
                v[None, :]
                - yshift
                - Zopt[:, None]
            )

            selection_cost = np.minimum(
                alpha / 2.0 * trait_difference**2,
                20000.0,
            )

            total_rate = (
                selection_cost
                + (gamma - 1.0)
                + Ik[:, None]
            )

            ndd = nd * np.exp(-dt_step * total_rate)


            # ------------------------------------------------
            # Reproduction term
            # ------------------------------------------------

            # MATLAB computes, row by row,
            #
            #   tric = conv(conv(Qtemp,f)*dv,f)*dv
            #
            # By associativity,
            #
            #   FFT(tric) = FFT(Qtemp) * FFT(f)^2 * dv^2.
            #
            # Rows with the same support size di are therefore
            # evaluated together with a batched FFT.

            active = np.flatnonzero(Ik > 1e-5)

            if active.size:

                # For active rows, the birth term equals nd
                # outside [im,ip].
                birth = np.zeros_like(nd)
                birth[active, :] = nd[active, :]

                for di_i in np.unique(di_local[active]):

                    di_i = int(di_i)

                    rows = active[di_local[active] == di_i]
                    im0 = im[rows] - 1

                    cols = im0[:, None] + np.arange(di_i)[None, :]
                    f_local = nd[rows[:, None], cols]

                    cached = reproduction_fft_cache.get(di_i)

                    if cached is None:

                        Qtemp = Q[
                            q_center - di_i :
                            q_center + di_i + 1
                        ]

                        # Qtemp * f * f has length 4*di_i - 1.
                        fft_length = next_fast_len(4 * di_i - 1)

                        Qhat = rfft(Qtemp, n=fft_length)

                        cached = (fft_length, Qhat)
                        reproduction_fft_cache[di_i] = cached

                    fft_length, Qhat = cached

                    Fhat = rfft(
                        f_local,
                        n=fft_length,
                        axis=1,
                    )

                    tric = irfft(
                        (Fhat * Fhat) * Qhat[None, :],
                        n=fft_length,
                        axis=1,
                    )

                    tric = (
                        tric[:, :4 * di_i - 1]
                        * dv**2
                    )

                    # MATLAB:
                    # truc(j)=tric(2*j+di-1), j=1,...,di
                    truc = tric[
                        :,
                        di_i : 3 * di_i - 1 : 2,
                    ]

                    birth[
                        rows[:, None],
                        cols,
                    ] = (
                        gamma
                        / Ik[rows, None]
                        * truc
                    )

                ndd += dt_step * birth


            # ------------------------------------------------
            # Implicit spatial diffusion
            # ------------------------------------------------

            # Solve all trait columns simultaneously.
            nd[1:-1, 1:-1] = solveh_banded(
                diffusion_ab,
                ndd[1:-1, 1:-1],
                lower=False,
                overwrite_ab=False,
                overwrite_b=False,
                check_finite=False,
            )


            # ------------------------------------------------
            # Boundary conditions
            # ------------------------------------------------

            nd[-1, 1:-1] = nd[-2, 1:-1]

            # MATLAB: j = 2:n-10
            jboundary = np.arange(1, n - 10)

            nd[0, jboundary] = (
                interp / deltaj
                * nd[1, jboundary + deltaj]
                + (deltaj - interp) / deltaj
                * nd[1, jboundary]
            )

            # MATLAB: nd(1,n-10:n-1)=...
            nd[0, n - 11:n - 1] = nd[1, n - 11:n - 1]

            nd[1:-1, 0] = 0.0
            nd[1:-1, -1] = 0.0

            nd[0, 0] = 0.5 * nd[1, 0] + 0.5 * nd[0, 1]
            nd[0, -1] = 0.5 * nd[1, -1] + 0.5 * nd[0, -2]

            nd[-1, 0] = 0.0
            nd[-1, -1] = 0.5 * nd[-1, -2] + 0.5 * nd[-2, -1]


            # ------------------------------------------------
            # Recording
            # ------------------------------------------------

            tx[s] = temps
            Ikrec[s] = np.linalg.norm(Ik, ord=1)


            # ------------------------------------------------
            # Rightmost occupied position
            # ------------------------------------------------

            k_idx = m - 1

            while Ik[k_idx] < thresholdI and k_idx > 0:
                k_idx -= 1

            # MATLAB: k=min(k,length(Ik)-1)
            k_idx = min(k_idx, len(Ik) - 2)
            k1 = k_idx + 1

            thetaI = (
                (thresholdI - Ik[k_idx + 1])
                / (Ik[k_idx] - Ik[k_idx + 1])
            )

            thetaI = min(1.0, thetaI)
            thetaI = max(0.0, thetaI)

            XtI = (
                thetaI * k1 * dx
                + (1.0 - thetaI) * (k1 + 1) * dx
                + xshift
            )

            xt[s] = XtI

            # MATLAB index:
            # floor(length(tx)*0.9)+1
            # -> Python zero-based floor(step_count*0.9)
            index90 = int(np.floor((s + 1) * 0.9))

            speed = (
                xt[s] - xt[index90]
            ) / max(
                dt_step,
                tx[s] - tx[index90],
            )

            speedrec[s] = speed

            zk = np.sum(nd[k_idx, :] * v) / np.sum(nd[k_idx, :])
            zk1 = (
                np.sum(nd[k_idx + 1, :] * v)
                / np.sum(nd[k_idx + 1, :])
            )

            ZtI = thetaI * zk + (1.0 - thetaI) * zk1 - yshift
            Zt[s] = ZtI

            speedz = (
                Zt[s] - Zt[index90]
            ) / max(
                dt_step,
                tx[s] - tx[index90],
            )

            speedrecz[s] = speedz


            # ------------------------------------------------
            # Leftmost occupied position
            # ------------------------------------------------

            km_idx = 1

            while Ik[km_idx] < thresholdI and km_idx < m - 1:
                km_idx += 1

            km1 = km_idx + 1

            thetaI = (
                (thresholdI - Ik[km_idx - 1])
                / (Ik[km_idx] - Ik[km_idx - 1])
            )

            thetaI = max(0.0, thetaI)
            thetaI = min(1.0, thetaI)

            XtIm = (
                thetaI * km1 * dx
                + (1.0 - thetaI) * (km1 - 1) * dx
                + xshift
            )

            xtm[s] = XtIm

            speedm = (
                xtm[s] - xtm[index90]
            ) / max(
                dt_step,
                tx[s] - tx[index90],
            )

            speedrecm[s] = speedm

            zk = np.sum(nd[km_idx, :] * v) / np.sum(nd[km_idx, :])
            zk1 = (
                np.sum(nd[km_idx - 1, :] * v)
                / np.sum(nd[km_idx - 1, :])
            )

            ZtIm = thetaI * zk + (1.0 - thetaI) * zk1 - yshift
            Ztm[s] = ZtIm

            speedzm = (
                Ztm[s] - Ztm[index90]
            ) / max(
                dt_step,
                tx[s] - tx[index90],
            )

            speedreczm[s] = speedzm
            maxrec[s] = np.linalg.norm(Ik, ord=np.inf)


            # ------------------------------------------------
            # Population-size growth rate
            # ------------------------------------------------

            # MATLAB index:
            # floor(length(tx)/2)+1
            # -> Python zero-based floor(step_count/2)
            index50 = int(np.floor((s + 1) / 2.0))

            lambda_value = np.log(
                Ikrec[s] / Ikrec[index50]
            ) / max(
                dt_step,
                tx[s] - tx[index50],
            )


            # ------------------------------------------------
            # Recenter the spatial calculation window
            # ------------------------------------------------

            if (
                temps > tfinal * (kpt / kplot)
                or temps > tfinal - 2.0 * dt_step
            ):

                kpt += 1

                Mx = np.sum(nd, axis=1)
                xcm = np.dot(x, Mx) / np.sum(Mx)

                if xcm > 0.6 * x[-1]:

                    xshift += x[spatial_shift_size - 1]

                    shift = spatial_shift_size

                    # Equivalent to the two matrix operations
                    # in the active Octave code.
                    nd[:-shift, :] = nd[shift:, :].copy()
                    nd[-shift:, :] = 0.0
                    nd[:, :shift] = 0.0


            # ------------------------------------------------
            # Advance time and write current values
            # ------------------------------------------------

            temps += dt_step

            write_csv_row(
                fid3,
                [
                    temps,
                    nu * np.sqrt(sigma**2 / 2.0),
                    speedm * np.sqrt(sigma**2 / (2.0 * rmax)),
                    speed * np.sqrt(sigma**2 / (2.0 * rmax)),
                    speedzm * np.sqrt(VLE * Vs),
                    speedz * np.sqrt(VLE * Vs),
                    XtIm,
                    XtI,
                    np.linalg.norm(Ik, ord=np.inf),
                    np.linalg.norm(Ik, ord=1),
                    lambda_value,
                    VLE,
                    b,
                    gamma,
                    tfinal,
                ],
            )

            step += 1


    # --------------------------------------------------------
    # Trim recording arrays
    # --------------------------------------------------------

    tx = tx[:step]
    xt = xt[:step]
    xtm = xtm[:step]
    Zt = Zt[:step]
    Ztm = Ztm[:step]
    Ikrec = Ikrec[:step]

    speedrec = speedrec[:step]
    speedrecm = speedrecm[:step]
    speedrecz = speedrecz[:step]
    speedreczm = speedreczm[:step]
    maxrec = maxrec[:step]


    # --------------------------------------------------------
    # Summary for this value of c
    # --------------------------------------------------------

    summary = {
        "c": nu * np.sqrt(sigma**2 / 2.0),
        "vS": speedm * np.sqrt(sigma**2 / (2.0 * rmax)),
        "vN": speed * np.sqrt(sigma**2 / (2.0 * rmax)),
        "vZS": speedzm * np.sqrt(VLE * Vs),
        "vZN": speedz * np.sqrt(VLE * Vs),
        "Ik1": Ikrec[int(np.floor(len(Ikrec) * 0.8))],
        "Ik2": Ikrec[-1],
        "lambda": lambda_value,
        "A": alpha,
        "B": beta,
        "C": gamma,
    }


    # --------------------------------------------------------
    # State needed for the final-profile export
    # --------------------------------------------------------

    final_state = {
        "m": m,
        "n": n,
        "x": x,
        "v": v,
        "nd": nd,
        "nd0": nd0,
        "Zopt": Zopt,
        "Zopt0": Zopt0,
        "xshift": xshift,
        "dv": dv,
        "VLE": VLE,
        "Vs": Vs,
        "sigma": sigma,
        "rmax": rmax,
        "b": b,
        "tfinal": tfinal,
        "c": c,
        "x0ini": x0ini,
        "y0ini": y0ini,
    }

    return summary, final_state


# ============================================================
# Final profile export
# ============================================================

def export_final_profiles(state, data_dir):

    x = state["x"]
    v = state["v"]
    nd = state["nd"]
    nd0 = state["nd0"]

    x0ini = state["x0ini"]
    y0ini = state["y0ini"]

    Zopt = state["Zopt"]
    Zopt0 = state["Zopt0"]

    xshift = state["xshift"]
    dv = state["dv"]

    VLE = state["VLE"]
    sigma = state["sigma"]
    rmax = state["rmax"]

    m = state["m"]

    sqrt_VLE = np.sqrt(VLE)

    xk0 = (x- x0ini) * np.sqrt(sigma**2 / (2.0 * rmax))
    xk = (x -x0ini+ xshift) * np.sqrt(sigma**2 / (2.0 * rmax))

    Zoptk = sqrt_VLE * (Zopt-y0ini)
    Zoptk0 = sqrt_VLE * (Zopt0-y0ini)

    xshiftrec = np.full(m, xshift)

    abs_nd = np.abs(nd)
    abs_nd0 = np.abs(nd0)

    den = np.sum(abs_nd, axis=1)
    den0 = np.sum(abs_nd0, axis=1)

    Ik = den * dv
    Ik0 = den0 * dv

    Zk = np.zeros(m)
    Vk = np.zeros(m)
    Zk0 = np.zeros(m)
    Vk0 = np.zeros(m)

    mask = Ik > 1e-6
    mask0 = Ik0 > 1e-6

    weighted_trait = np.sum(abs_nd * v[None, :], axis=1)
    weighted_trait0 = np.sum(abs_nd0 * v[None, :], axis=1)

    # Mean phenotype in the internal computational v-coordinate.  This
    # uncentered mean must be used when computing the variance.
    mean_trait = np.divide(
        weighted_trait,
        den,
        out=np.zeros_like(weighted_trait),
        where=den > 0,
    )
    mean_trait0 = np.divide(
        weighted_trait0,
        den0,
        out=np.zeros_like(weighted_trait0),
        where=den0 > 0,
    )

    # Exported mean phenotype: scale by sqrt(VLE) and center at the
    # initial optimum y0ini.
    Zk[mask] = sqrt_VLE * (mean_trait[mask] - y0ini)
    Zk0[mask0] = sqrt_VLE * (mean_trait0[mask0] - y0ini)

    # Variance is translation invariant: center around the raw mean in
    # computational v-coordinates, not around the centered exported Zk.
    centered_trait = v[None, :] - mean_trait[:, None]
    centered_trait0 = v[None, :] - mean_trait0[:, None]

    variance_numerator = np.sum(
        abs_nd * centered_trait**2,
        axis=1,
    )

    variance_numerator0 = np.sum(
        abs_nd0 * centered_trait0**2,
        axis=1,
    )

    Vk[mask] = (
        VLE
        * variance_numerator[mask]
        / den[mask]
    )

    Vk0[mask0] = (
        VLE
        * variance_numerator0[mask0]
        / den0[mask0]
    )

    M = np.column_stack(
        [
            xk,
            Ik,
            Zk,
            Vk,
            Zoptk,
            xk0,
            Ik0,
            Zk0,
            Vk0,
            Zoptk0,
            xshiftrec,
        ]
    )

    tfinal_label = (
        int(state["tfinal"])
        if float(state["tfinal"]).is_integer()
        else state["tfinal"]
    )

    filename = data_dir / (
        f"RECVLEi{state['VLE']}"
        f"b{state['b']}"
        f"tfinal{tfinal_label}"
        f"c{state['c']}.csv"
    )

    np.savetxt(
        filename,
        M,
        delimiter=",",
        header="x,I,Z,V,Zopt,x0,I0,Z0,V0,Zopt0,xshift",
        comments="",
    )

# ============================================================
# Command-line interface: one simulation only
# ============================================================

def main():
    """Run one sexual simulation and append its summary to the macro file.

    Usage
    -----
    python3 core-sexual-simulations.py c idata_macrofile

    c
        Climate-change speed for this simulation.
    idata_macrofile
        Index of the macro summary file created by sexual-macro.py.
    """
    if len(sys.argv) < 3:
        raise SystemExit(
            "Usage: python3 core-sexual-simulations.py "
            "c idata_macrofile"
        )

    c = float(sys.argv[1])
    idata_macrofile = int(float(sys.argv[2]))

    data_dir = Path("data")
    data_dir.mkdir(parents=True, exist_ok=True)

    print(f"c = {c}")

    summary, final_state = run_simulation(
        c=c,
        data_dir=data_dir,
    )

    # Export the final profile for this simulation.
    export_final_profiles(
        state=final_state,
        data_dir=data_dir,
    )

    # Append one line to the summary file created by sexual-macro.py.
    macro_file = data_dir / f"data_macro_sexual{idata_macrofile}.csv"
    with macro_file.open("a") as fmacro:
        write_csv_row(
            fmacro,
            [
                summary["c"],
                summary["vS"],
                summary["vN"],
                summary["vZS"],
                summary["vZN"],
                summary["Ik1"],
                summary["Ik2"],
                summary["lambda"],
                summary["A"],
                summary["B"],
                summary["C"],
                final_state["tfinal"],
            ],
        )

    print(
        "c,vS,vN,vZS,vZN,Ik1,Ik2,lambda,A,B,C,tfinal = "
        f"{summary['c']},{summary['vS']},{summary['vN']},"
        f"{summary['vZS']},{summary['vZN']},{summary['Ik1']},"
        f"{summary['Ik2']},{summary['lambda']},{summary['A']},"
        f"{summary['B']},{summary['C']},{final_state['tfinal']}"
    )


if __name__ == "__main__":
    main()
