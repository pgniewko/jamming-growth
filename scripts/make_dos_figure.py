#!/usr/bin/env python3
"""Make figure  -- vibrational density of states of arrested packings.

Builds the dynamical matrix of the mechanically arrested budding-cell packings
(``scripts/hessian_dos.py``) and tests the three spectral predictions of the
post-jamming theory:

(A) D(omega) at several packing fractions from phi_J to phi_2 at fixed strong
    feedback: the low-frequency band (the unconstrained-bud rotational modes)
    drains as buds are completed.
(B) The band-edge frequency omega* collapses onto a single near-linear branch
    when plotted against the *modified* excess coordination Delta Z = Z-(6-2u);
    against the *naive* count Z-6 the same points do not collapse (and lie at
    negative excess for more than half the states).
(C) The count of quasi-zero modes N_zero (lambda < lambda_tol) equals the
    independently measured unconstrained-bud number N_u.

Usage:
    python3 scripts/make_dos_figure.py [--size 15] [--seeds 1201-1212]
"""

from __future__ import annotations

import argparse
import re
import warnings
from pathlib import Path

import numpy as np

from make_paper_figures import (
    P0_COLORS,
    P0_LABELS,
    P0_ORDER,
    add_panel_label,
    configure_style,
    save_figure,
)
from hessian_dos import analyze_packing, count_zero_modes, omega_star

warnings.filterwarnings("ignore")

REPO_ROOT = Path(__file__).resolve().parents[1]
GROWTH = REPO_ROOT / "data" / "output" / "growth"

LAMBDA_TOL = 1e-6
N_EDGE = 3
PANEL_A_P0 = "1e-4"  # one strong-feedback value for the draining panel
PANEL_A_DPHIS = ["1e-3", "3e-3", "1e-2"]  # plus the jamming state


def _files(lx, p0, seed=None, tag="LF_DPHI"):
    pat = f"{tag}_*seed_{'*' if seed is None else seed}_Lx{lx}_Ly{lx}_*_P{p0}.dat.gz"
    return sorted((GROWTH).glob(pat))


def _dphi_of(path):
    m = re.search(r"dphi([0-9.e-]+)_P", str(path))
    return m.group(1) if m else None


# ----------------------------------------------------------------------------
# Panel A: averaged D(omega) at several phi for one strong P0
# ----------------------------------------------------------------------------
def panel_a_data(lx, seeds, n_bins=42, omega_max=2.6):
    """Averaged extended-band D(omega) at several phi for one strong P0.

    The trivial quasi-zero modes (lambda < LAMBDA_TOL, the unconstrained-bud
    rotations counted in Panel C) are excluded from the histogram so that
    D(omega) shows the genuine vibrational band; its low-frequency weight
    drains and its onset omega* shifts up as buds are completed.
    """
    edges = np.linspace(0.0, omega_max, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    states = [("jamm", "LF_JAMM")] + [(d, "LF_DPHI") for d in PANEL_A_DPHIS]
    curves = []
    for key, tag in states:
        hist_stack = []
        phis, nus, wstars = [], [], []
        for seed in seeds:
            if tag == "LF_JAMM":
                cands = _files(lx, PANEL_A_P0, seed, tag="LF_JAMM")
            else:
                cands = [
                    p for p in _files(lx, PANEL_A_P0, seed, tag="LF_DPHI")
                    if _dphi_of(p) == key
                ]
            if not cands:
                continue
            try:
                spec = analyze_packing(cands[0], float(lx))
            except Exception:
                continue
            omega_ext = spec.omega[spec.eigvals >= LAMBDA_TOL]
            counts, _ = np.histogram(omega_ext, bins=edges)
            hist_stack.append(counts / spec.ndof)
            phis.append(spec.phi)
            nus.append(spec.n_u)
            wstars.append(omega_star(spec, LAMBDA_TOL, N_EDGE))
        if hist_stack:
            curves.append(
                {
                    "key": key,
                    "D": np.mean(hist_stack, axis=0),
                    "phi": float(np.mean(phis)),
                    "n_u": float(np.mean(nus)),
                    "omega_star": float(np.mean(wstars)),
                    "n": len(hist_stack),
                }
            )
    return centers, curves


# ----------------------------------------------------------------------------
# Panels B/C: omega*, Delta Z, N_zero, N_u over the (P0, dphi, seed) grid
# ----------------------------------------------------------------------------
def grid_data(lx, p0s, seeds):
    rows = []
    for p0 in p0s:
        for seed in seeds:
            for path in _files(lx, p0, seed, tag="LF_DPHI"):
                try:
                    spec = analyze_packing(path, float(lx))
                except Exception:
                    continue
                ws = omega_star(spec, LAMBDA_TOL, N_EDGE)
                if not np.isfinite(ws):
                    continue
                rows.append(
                    {
                        "p0": p0,
                        "seed": seed,
                        "dphi": _dphi_of(path),
                        "dZ_mod": spec.dZ_mod,
                        "dZ_naive": spec.dZ_naive,
                        "omega_star": ws,
                        "n_zero": count_zero_modes(spec, LAMBDA_TOL),
                        "n_u": spec.n_u,
                    }
                )
    return rows


def compute_data(lx, seeds, cache=None):
    """Compute (or load cached) the expensive (B/C) grid; Panel A is cheap."""
    import pickle

    seeds_a = seeds[: min(10, len(seeds))]
    centers, a_curves = panel_a_data(lx, seeds_a)
    if cache and Path(cache).exists():
        with open(cache, "rb") as fh:
            rows = pickle.load(fh)
    else:
        rows = grid_data(lx, P0_ORDER, seeds)
        if cache:
            Path(cache).parent.mkdir(parents=True, exist_ok=True)
            with open(cache, "wb") as fh:
                pickle.dump(rows, fh)
    return centers, a_curves, rows


def make_figure(lx, seeds, output_dir, cache=None):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    centers, a_curves, rows = compute_data(lx, seeds, cache=cache)

    dZm = np.array([r["dZ_mod"] for r in rows])
    dZn = np.array([r["dZ_naive"] for r in rows])
    ws = np.array([r["omega_star"] for r in rows])
    nz = np.array([r["n_zero"] for r in rows])
    nu = np.array([r["n_u"] for r in rows])
    colors = [P0_COLORS[r["p0"]] for r in rows]

    # fits / diagnostics
    slope_origin = float(np.sum(ws * dZm) / np.sum(dZm * dZm))
    A = np.polyfit(dZm, ws, 1)
    R2_mod = 1 - np.var(ws - np.polyval(A, dZm)) / np.var(ws)
    An = np.polyfit(dZn, ws, 1)
    R2_naive = 1 - np.var(ws - np.polyval(An, dZn)) / np.var(ws)
    frac_neg_naive = float(np.mean(dZn < 0))
    exact = int(np.sum(nz == nu))
    czfit = float(np.sum(nz * nu) / np.sum(nu * nu))  # slope of N_zero vs N_u thru 0

    print(f"[Panel B] omega* = {slope_origin:.4f} * dZ_mod (through origin); "
          f"linear fit slope={A[0]:.4f} intercept={A[1]:+.4f} R^2(mod)={R2_mod:.3f} "
          f"R^2(naive)={R2_naive:.3f}; frac(dZ_naive<0)={frac_neg_naive:.2f}; N={len(rows)}")
    print(f"[Panel C] N_zero == N_u in {exact}/{len(rows)} states; "
          f"slope(N_zero vs N_u thru 0)={czfit:.4f}; "
          f"max|N_zero-N_u|={int(np.max(np.abs(nz - nu)))}")

    configure_style()
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5))
    axA, axB, axC = axes

    # ---- Panel A: extended-band D(omega) draining; omega* edge shifts up ----
    cmap = plt.cm.viridis
    n_curves = len(a_curves)
    ymax = 0.0
    for i, cur in enumerate(a_curves):
        frac = i / max(n_curves - 1, 1)
        color = cmap(0.10 + 0.80 * frac)
        if cur["key"] == "jamm":
            label = rf"$\phi_J$"
        else:
            label = rf"$\Delta\phi\!=\!{cur['key']}$"
        axA.plot(centers, cur["D"], "-", color=color, lw=1.3, label=label)
        ymax = max(ymax, cur["D"].max())
        # mark the band edge omega* with a tick along the baseline
        axA.plot([cur["omega_star"]], [0.0], marker="^", color=color, ms=4.5,
                 clip_on=False, zorder=5)
    axA.set_xlabel(r"$\omega$")
    axA.set_ylabel(r"$D(\omega)$")
    axA.set_xlim(0, 1.6)
    axA.set_ylim(0, 1.55 * ymax)  # headroom so the legend clears the curves
    axA.annotate(r"$\omega^\ast$ shifts up", xy=(0.10, 0.42 * ymax),
                 fontsize=6.3, color="0.3")
    axA.legend(fontsize=6.4, loc="best", handlelength=1.0, labelspacing=0.2,
               borderpad=0.3, handletextpad=0.4, ncol=2, columnspacing=0.9,
               title=rf"$P_0=10^{{{int(np.log10(float(PANEL_A_P0)))}}}$ (curves: $\phi_J\!\to\!\phi_2$)",
               title_fontsize=6.4)
    add_panel_label(axA, "A")

    # ---- Panel B: omega* vs Delta Z (modified collapse + naive contrast) ----
    axB.axvline(0.0, color="0.8", lw=0.8, zorder=0)
    axB.scatter(dZn, ws, s=10, facecolors="none", edgecolors="0.62",
                linewidths=0.6, label=r"naive $Z\!-\!6$", zorder=2)
    axB.scatter(dZm, ws, s=11, c=colors, edgecolors="white", linewidths=0.3,
                label=r"mod. $Z\!-\!(6\!-\!2u)$", zorder=3)
    xs = np.linspace(0, dZm.max(), 50)
    axB.plot(xs, slope_origin * xs, "-", color="0.15", lw=1.1, zorder=4)
    axB.set_xlabel(r"$\Delta Z$")
    axB.set_ylabel(r"$\omega^\ast$")
    axB.legend(fontsize=6.3, loc="upper left", handlelength=1.1, labelspacing=0.25,
               borderpad=0.3, handletextpad=0.4)
    add_panel_label(axB, "B")

    # ---- Panel C: N_zero vs N_u ----
    nmax = max(nz.max(), nu.max())
    axC.plot([0, nmax], [0, nmax], "--", color="0.5", lw=1.0, zorder=1,
             label=r"$N_{\rm zero}\!=\!N_u$")
    axC.scatter(nu, nz, s=12, c=colors, edgecolors="white", linewidths=0.3, zorder=3)
    axC.set_xlabel(r"$N_u$")
    axC.set_ylabel(r"$N_{\rm zero}$")
    axC.set_xlim(left=-1)
    axC.set_ylim(bottom=-1)
    add_panel_label(axC, "C")

    # shared P0 legend on panel C
    handles = [Line2D([0], [0], ls="--", color="0.5", label=r"$N_{\rm zero}=N_u$")]
    handles += [
        Line2D([0], [0], marker="o", ls="", mfc=P0_COLORS[p0], mec="white",
               ms=4.5, label=P0_LABELS[p0])
        for p0 in P0_ORDER
    ]
    axC.legend(handles=handles, title=r"$P_0$", fontsize=6.0, title_fontsize=6.8,
               loc="lower right", ncol=2, handlelength=1.0, labelspacing=0.22,
               columnspacing=0.7, borderpad=0.3)

    fig.tight_layout(w_pad=1.1)
    save_figure(fig, output_dir, f"figureS4_dos_l{lx}", dpi=300)
    plt.close(fig)
    return {
        "slope_origin": slope_origin,
        "slope_lin": float(A[0]),
        "intercept_lin": float(A[1]),
        "R2_mod": float(R2_mod),
        "R2_naive": float(R2_naive),
        "frac_neg_naive": frac_neg_naive,
        "n_states": len(rows),
        "nzero_eq_nu": exact,
        "cz_slope": czfit,
    }


def parse_seeds(spec):
    if "-" in spec:
        a, b = spec.split("-")
        return [str(s) for s in range(int(a), int(b) + 1)]
    return [s.strip() for s in spec.split(",") if s.strip()]


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--size", default="15")
    ap.add_argument("--seeds", default="1201-1212")
    ap.add_argument("--output-dir", default=str(REPO_ROOT / "paper" / "figures"))
    ap.add_argument("--cache", default=None,
                    help="optional pickle cache of computed spectra (speeds re-plotting)")
    return ap.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    seeds = parse_seeds(args.seeds)
    stats = make_figure(args.size, seeds, output_dir, cache=args.cache)
    print(f"Wrote figure to: {output_dir}")
    return stats


if __name__ == "__main__":
    main()
