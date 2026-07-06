#!/usr/bin/env python3
"""Figure 4 -- post-jamming state diagram in the (P/P_max, G/P_max) plane.

For each feedback strength ``P0`` the arrested packing is traced from the jamming
point ``phi_J`` to the depletion density ``phi_2`` by sweeping the post-jamming
compaction ``dphi``. Each ``(P0, seed, dphi)`` run contributes one point
``(phi, P, G)`` where ``phi`` and ``P`` are read from the last
``POSTJAMM_SUMMARY`` row and ``G`` is the shear modulus (linear fit of the
``G_data_LF_DPHI`` shear-stress response). The arrest endpoint ``(P_2, G_2)``
uses the depletion-state shear ``G_data_LF_PHI2`` where available. Trajectories
are averaged across runs on a normalized ``phi`` grid; pressure is plotted on a
logarithmic axis and the shear modulus on a linear axis, with the rigidity
plateau ``G_2(0+) = n_J k l_c^2 M_4 u_J`` drawn as a dashed line.

Reuses the palette/parsers/styling of ``make_paper_figures.py`` so the figure is
visually consistent with Figs. 2-3.

Usage:
    python3 scripts/make_state_diagram_figure.py [--size 15] [--output-dir paper/figures]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from make_paper_figures import (
    P0_COLORS,
    P0_LABELS,
    P0_ORDER,
    P_MAX,
    configure_style,
    estimate_shear_modulus,
    parse_post_rows,
    save_figure,
)
from hessian_dos import analyze_packing

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPO_ROOT / "data" / "output"
GROWTH = DATA_ROOT / "growth"
SHEAR = DATA_ROOT / "shear"

SEEDS = [str(s) for s in range(1201, 1221)]
DPHIS = ["1e-4", "3e-4", "1e-3", "3e-3", "1e-2", "3e-2", "6e-2", "1e-1", "1.2e-1"]
DPHI_VAL = {d: float(d) for d in DPHIS}
ACTIVE_CUTOFF = 1  # SECONDARY_ARREST_ACTIVE_CUTOFF
M4 = 1.0 / 8.0
K = 1.0
N_GRID = 24  # points on the normalized phi grid


def basename(seed, lx, dphi, p0):
    return (
        f"v1.0_ar1.01_div_4_desync0.4_seed_{seed}_Lx{lx}_Ly{lx}"
        f"_att0.0_dphi{dphi}_P{p0}.dat"
    )


def _arrest_row(rows):
    """First post-jamming row with the initial-free reservoir essentially gone."""
    for r in rows:
        if r.initial_free_active <= ACTIVE_CUTOFF:
            return r
    return None


def collect_runs(lx):
    """Return per-(p0, seed) trajectories and arrest endpoints."""
    out = {}  # (p0, seed) -> dict
    for p0 in P0_ORDER:
        for seed in SEEDS:
            traj = []  # (phi, P, G) at dphi endpoints
            phi_j = None
            arrest = None  # (phi2, P2, G2)
            for dphi in DPHIS:
                name = basename(seed, lx, dphi, p0)
                post = GROWTH / f"POSTJAMM_SUMMARY_{name}"
                gdphi = SHEAR / f"G_data_LF_DPHI_{name}"
                if not post.exists():
                    continue
                rows = parse_post_rows(post)
                if not rows:
                    continue
                if phi_j is None:
                    phi_j = rows[0].phi
                # endpoint of this run
                if gdphi.exists():
                    g = estimate_shear_modulus(gdphi)
                    if np.isfinite(g) and g > 0:
                        traj.append((rows[-1].phi, rows[-1].pressure, g))
                # arrest endpoint (depletion state) from this run, if reached
                row2 = _arrest_row(rows)
                if row2 is not None:
                    gphi2 = SHEAR / f"G_data_LF_PHI2_{name}"
                    g2 = estimate_shear_modulus(gphi2) if gphi2.exists() else np.nan
                    cand = (row2.phi, row2.pressure, g2, DPHI_VAL[dphi])
                    # keep the smallest-dphi run that reaches arrest
                    if arrest is None or cand[3] < arrest[3]:
                        arrest = cand
            if phi_j is None or len(traj) < 2:
                continue
            out[(p0, seed)] = {
                "phi_j": phi_j,
                "traj": sorted(traj, key=lambda t: t[0]),
                "arrest": arrest,
            }
    return out


def build_mean_curves(runs):
    """Average each P0's (P/P_max, G) trajectory on a normalized phi grid."""
    s_grid = np.linspace(0.0, 1.0, N_GRID)
    curves = {}
    endpoints = {}
    run_endpoints = {}
    for p0 in P0_ORDER:
        logP_stack = []
        G_stack = []
        end_P = []
        end_G = []
        for (q0, seed), rec in runs.items():
            if q0 != p0:
                continue
            arrest = rec["arrest"]
            phi_j = rec["phi_j"]
            # define phi_2: arrest phi if available else last trajectory phi
            if arrest is not None:
                phi_2 = arrest[0]
            else:
                phi_2 = rec["traj"][-1][0]
            if phi_2 <= phi_j:
                continue
            pts = [(p, P, G) for (p, P, G) in rec["traj"] if p <= phi_2 + 1e-9]
            # append arrest endpoint with a valid G2
            if arrest is not None and np.isfinite(arrest[2]) and arrest[2] > 0:
                pts.append((arrest[0], arrest[1], arrest[2]))
                end_P.append(arrest[1])
                end_G.append(arrest[2])
            elif pts:
                end_P.append(pts[-1][1])
                end_G.append(pts[-1][2])
            pts = sorted(set(pts), key=lambda t: t[0])
            if len(pts) < 2:
                continue
            phi = np.array([t[0] for t in pts])
            P = np.array([t[1] for t in pts])
            G = np.array([t[2] for t in pts])
            s = (phi - phi_j) / (phi_2 - phi_j)
            s = np.clip(s, 0.0, 1.0)
            # ensure strictly increasing s for interpolation
            order = np.argsort(s)
            s, P, G = s[order], P[order], G[order]
            s_u, idx = np.unique(s, return_index=True)
            if s_u.size < 2:
                continue
            logP = np.interp(s_grid, s_u, np.log10(np.maximum(P[idx], 1e-12) / P_MAX))
            Gi = np.interp(s_grid, s_u, G[idx])
            logP_stack.append(logP)
            G_stack.append(Gi)
        if logP_stack:
            curves[p0] = {
                "s": s_grid,
                "logP": np.mean(logP_stack, axis=0),
                "G": np.mean(G_stack, axis=0),
                "G_se": (np.std(G_stack, axis=0, ddof=1) / np.sqrt(len(G_stack)))
                if len(G_stack) > 1
                else np.zeros_like(s_grid),
                "n": len(logP_stack),
            }
        if end_P:
            endpoints[p0] = {
                "P2": float(np.median(end_P)),
                "G2": float(np.median(end_G)),
                "n": len(end_P),
            }
            run_endpoints[p0] = list(zip(end_P, end_G))
    return curves, endpoints, run_endpoints


def plateau_inputs(lx):
    """Recompute n_J, u_J, l_c from the strong-feedback (P0=1e-5) jammed packings."""
    files = sorted(GROWTH.glob(f"LF_JAMM_*Lx{lx}_Ly{lx}_*_P1e-5.dat.gz"))
    lxi = int(lx)
    n_js, u_js, l_cs = [], [], []
    for f in files:
        try:
            spec = analyze_packing(f, float(lxi))
        except Exception:
            continue
        n_js.append(spec.n_nonrattler / (lxi * lxi))
        u_js.append(spec.u)
        l_cs.append(spec.l_c)
    n_j = float(np.mean(n_js))
    u_j = float(np.mean(u_js))
    l_c = float(np.mean(l_cs))
    g_plateau = n_j * K * l_c**2 * M4 * u_j
    return n_j, u_j, l_c, g_plateau, len(n_js)


def make_figure(lx, output_dir):
    import matplotlib.pyplot as plt

    runs = collect_runs(lx)
    curves, endpoints, run_endpoints = build_mean_curves(runs)
    n_j, u_j, l_c, g_plateau, n_jam = plateau_inputs(lx)

    print(f"plateau inputs (P0=1e-5, {n_jam} jammed packings): "
          f"n_J={n_j:.3f} u_J={u_j:.3f} l_c={l_c:.3f} -> G2(0+)={g_plateau:.4f}")
    med_P = {p0: endpoints[p0]["P2"] for p0 in endpoints}
    if med_P:
        span = max(med_P.values()) / min(med_P.values())
        print(f"median-endpoint pressure fan-out: {np.log10(span):.2f} decades "
              f"(P2/Pmax {min(med_P.values())/P_MAX:.2e} .. {max(med_P.values())/P_MAX:.2e})")
    allP = [p for pts in run_endpoints.values() for (p, g) in pts]
    if allP:
        print(f"per-run pressure fan-out: {np.log10(max(allP)/min(allP)):.2f} decades "
              f"(P2/Pmax {min(allP)/P_MAX:.2e} .. {max(allP)/P_MAX:.2e})")
    for p0 in P0_ORDER:
        if p0 in endpoints:
            e = endpoints[p0]
            print(f"  P0={p0:>4}: endpoint (P2/Pmax={e['P2']/P_MAX:.3e}, G2={e['G2']:.4f}) "
                  f"n_runs={e['n']}, traj n={curves.get(p0, {}).get('n', 0)}")

    configure_style()
    fig, ax = plt.subplots(figsize=(3.7, 3.3))

    # rigidity plateau (the P0 -> 0 arrest limit of G_2), normalized by P_max
    gp = g_plateau / P_MAX
    ax.axhline(gp, ls="--", color="0.35", lw=1.0, zorder=1)
    ax.text(
        0.985, gp, r"rigidity plateau $G_2(0^+)/P_{\max}$",
        transform=ax.get_yaxis_transform(),
        ha="right", va="bottom", color="0.35", fontsize=7,
    )

    for p0 in P0_ORDER:
        color = P0_COLORS[p0]
        label = P0_LABELS[p0]
        if p0 in curves:
            c = curves[p0]
            x = 10.0 ** c["logP"]
            ax.plot(x, c["G"] / P_MAX, "-", color=color, lw=1.4, zorder=3, label=label)
        # individual-run arrest endpoints (faint), showing the full fan-out
        if p0 in run_endpoints:
            pts = np.array(run_endpoints[p0])
            ax.scatter(pts[:, 0] / P_MAX, pts[:, 1] / P_MAX, s=7, color=color,
                       alpha=0.18, edgecolors="none", zorder=2)
        if p0 in endpoints:
            e = endpoints[p0]
            ax.plot(
                e["P2"] / P_MAX, e["G2"] / P_MAX, "o", color=color, ms=7,
                mec="white", mew=0.7, zorder=5,
            )

    # directional cue placed in the empty upper-left region (no datapoints there):
    # as feedback strengthens, arrest endpoints move to low pressure / the plateau.
    ax.annotate(
        "", xy=(3.16e-3, 0.82), xytext=(1.8e-2, 1.05),
        arrowprops=dict(arrowstyle="->", color="0.45", lw=1.3,
                        connectionstyle="arc3,rad=0.18"), zorder=4,
    )
    ax.text(7.5e-3, 1.16, "stronger feedback", fontsize=7.5, color="0.45",
            ha="center", va="center")

    ax.set_xscale("log")
    ax.set_xlabel(r"$P/P_{\max}$")
    ax.set_ylabel(r"$G/P_{\max}$")
    ax.set_ylim(0.0, 1.7)
    leg = ax.legend(title=r"$P_0$", fontsize=7.5, title_fontsize=8,
                    loc="upper left", handlelength=1.4, labelspacing=0.25)
    leg._legend_box.align = "left"

    save_figure(fig, output_dir, f"figure4_state_diagram_l{lx}", dpi=300)
    plt.close(fig)
    return curves, endpoints, (n_j, u_j, l_c, g_plateau)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--size", default="15")
    ap.add_argument("--output-dir", default=str(REPO_ROOT / "paper" / "figures"))
    return ap.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    make_figure(args.size, output_dir)
    print(f"Wrote figure to: {output_dir}")


if __name__ == "__main__":
    main()
