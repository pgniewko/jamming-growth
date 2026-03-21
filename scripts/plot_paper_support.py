#!/usr/bin/env python3

import argparse
import gzip
import math
import os
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PLOT_CACHE_DIR = REPO_ROOT / "output" / ".plot_cache"
PLOT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(PLOT_CACHE_DIR / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(PLOT_CACHE_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from run_growth_shear import DPHIS as SWEEP_DPHIS, P0S as SWEEP_P0S, SIZES as SWEEP_SIZES

DEFAULT_OUTPUT_DIR = REPO_ROOT / "output" / "plots"
SWEEP_SIZE_VALUES = [int(value) for value in SWEEP_SIZES]
SWEEP_P0_VALUES = [float(value) for value in SWEEP_P0S]
SWEEP_DPHI_VALUES = [float(value) for value in SWEEP_DPHIS]
DEFAULT_MAIN_L = SWEEP_SIZE_VALUES[0]
DEFAULT_MAIN_P0S = list(SWEEP_P0_VALUES)
DEFAULT_HIST_P0S = [value for value in DEFAULT_MAIN_P0S if value > 0.0]
SHEAR_FIT_MAX_STRAIN = 5e-5
COMPRESSED_PRESSURE_TOL = 1e-12
MIN_RATIO_DELTA_Z = 0.05

NC_FILENAME_RE = re.compile(
    r"^NC_v(?P<version>[^_]+)_ar(?P<ar>[^_]+)_div_(?P<divtype>[^_]+)"
    r"_desync(?P<desync>[^_]+)_seed_(?P<seed>[^_]+)"
    r"_Lx(?P<L>[^_]+)_Ly(?P<Ly>[^_]+)_att(?P<att>[^_]+)"
    r"_dphi(?P<dphi>[^_]+)_P(?P<p0>.+)\.dat$"
)
BEXT_FILENAME_RE = re.compile(
    r"^B_ext_data_dphiprobe(?P<probe_tag>[^_]+)_LF_DPHI_"
    r"v(?P<version>[^_]+)_ar(?P<ar>[^_]+)_div_(?P<divtype>[^_]+)"
    r"_desync(?P<desync>[^_]+)_seed_(?P<seed>[^_]+)"
    r"_Lx(?P<L>[^_]+)_Ly(?P<Ly>[^_]+)_att(?P<att>[^_]+)"
    r"_dphi(?P<dphi>[^_]+)_P(?P<p0>.+)\.dat$"
)
G_FILENAME_RE = re.compile(
    r"^G_data_LF_DPHI_v(?P<version>[^_]+)_ar(?P<ar>[^_]+)_div_(?P<divtype>[^_]+)"
    r"_desync(?P<desync>[^_]+)_seed_(?P<seed>[^_]+)"
    r"_Lx(?P<L>[^_]+)_Ly(?P<Ly>[^_]+)_att(?P<att>[^_]+)"
    r"_dphi(?P<dphi>[^_]+)_P(?P<p0>.+)\.dat$"
)
STATS_FILENAME_RE = re.compile(
    r"^STATS_LF_DPHI_v(?P<version>[^_]+)_ar(?P<ar>[^_]+)_div_(?P<divtype>[^_]+)"
    r"_desync(?P<desync>[^_]+)_seed_(?P<seed>[^_]+)"
    r"_Lx(?P<L>[^_]+)_Ly(?P<Ly>[^_]+)_att(?P<att>[^_]+)"
    r"_dphi(?P<dphi>[^_]+)_P(?P<p0>.+)\.dat$"
)
STEPLOG_FILENAME_RE = re.compile(
    r"^STEPLOG_v(?P<version>[^_]+)_ar(?P<ar>[^_]+)_div_(?P<divtype>[^_]+)"
    r"_desync(?P<desync>[^_]+)_seed_(?P<seed>[^_]+)"
    r"_Lx(?P<L>[^_]+)_Ly(?P<Ly>[^_]+)_att(?P<att>[^_]+)"
    r"_dphi(?P<dphi>[^_]+)_P(?P<p0>.+)\.dat$"
)

def nearest_log_value(target, values):
    return min(values, key=lambda value: abs(math.log10(value) - math.log10(target)))


def build_p0_colors(values):
    colors = {}
    ordered_values = sorted(set(values), reverse=True)
    positive = [value for value in ordered_values if value > 0.0]
    if -1.0 in ordered_values:
        colors[-1.0] = "#4c566a"
    if positive:
        cmap = matplotlib.colormaps["viridis"]
        for value, position in zip(positive, np.linspace(0.18, 0.9, len(positive))):
            colors[value] = cmap(position)
    return colors


DEFAULT_HIST_DPHI = nearest_log_value(2e-2, SWEEP_DPHI_VALUES)
P0_COLORS = build_p0_colors(DEFAULT_MAIN_P0S)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate manuscript-support plots from existing growth, shear, and B_ext outputs."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for figures and CSV summaries. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--main-L",
        type=int,
        default=DEFAULT_MAIN_L,
        help=f"Primary system size used for replicated manuscript plots. Default: {DEFAULT_MAIN_L}",
    )
    parser.add_argument(
        "--hist-dphi",
        type=float,
        default=DEFAULT_HIST_DPHI,
        help=f"dphi slice used for growth-rate histograms. Default: {DEFAULT_HIST_DPHI:.3g}",
    )
    parser.add_argument(
        "--formats",
        default="png,pdf",
        help="Comma-separated image formats to write. Default: png,pdf",
    )
    return parser.parse_args()


def plot_formats(value):
    formats = [item.strip().lower() for item in value.split(",") if item.strip()]
    if not formats:
        raise ValueError("at least one output format is required")
    return formats


def p0_order(values):
    ordered = [value for value in DEFAULT_MAIN_P0S if value in values]
    leftovers = sorted(set(values) - set(ordered))
    return ordered + leftovers


def p0_label(value):
    if np.isclose(value, -1.0):
        return "No feedback"
    return f"$P_0 = {value:.0e}$"


def save_dataframe(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def save_figure(fig, output_dir, stem, formats):
    for ext in formats:
        fig.savefig(output_dir / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def cast_metadata(frame):
    if frame.empty:
        return frame
    for column in ("seed", "L", "Ly"):
        if column in frame.columns:
            frame[column] = frame[column].astype(int)
    for column in ("ar", "att", "desync", "dphi", "p0"):
        if column in frame.columns:
            frame[column] = frame[column].astype(float)
    return frame


def add_contact_observables(frame):
    frame = frame.copy()
    frame["n_nonfloat"] = frame["N"] - frame["Nf"]
    frame["u"] = frame["Nu"] / frame["n_nonfloat"]
    frame["Z"] = frame["Nc"] / frame["n_nonfloat"]
    frame["Ziso"] = frame["Ziso_total"] / frame["n_nonfloat"]
    frame["DeltaZ"] = (frame["Nc"] - frame["Ziso_total"]) / frame["n_nonfloat"]
    return frame


def load_growth_endpoints():
    rows = []
    bad_files = []
    for path in sorted((REPO_ROOT / "output" / "growth").glob("NC_*.dat")):
        match = NC_FILENAME_RE.match(path.name)
        if match is None:
            continue
        try:
            fields = path.read_text(encoding="utf-8").split()
            if len(fields) != 26:
                raise ValueError("unexpected NC value count")
            values = list(map(float, fields))
            for stage, offset in (("jamm", 0), ("final", 13)):
                N, Ziso_total, Nc, Nf, Nu, Nmm, Nbb, Nmb = map(int, values[offset : offset + 8])
                phi, pressure, fret, p0_file, total_growthrate = values[offset + 8 : offset + 13]
                rows.append(
                    {
                        **match.groupdict(),
                        "stage": stage,
                        "N": N,
                        "Ziso_total": Ziso_total,
                        "Nc": Nc,
                        "Nf": Nf,
                        "Nu": Nu,
                        "Nmm": Nmm,
                        "Nbb": Nbb,
                        "Nmb": Nmb,
                        "phi": phi,
                        "P": pressure,
                        "fret": fret,
                        "p0_file": p0_file,
                        "total_growthrate": total_growthrate,
                        "source_file": path.name,
                    }
                )
        except (OSError, ValueError):
            bad_files.append(path.name)
            continue
    if not rows:
        raise SystemExit("No complete NC endpoint files found under output/growth")
    frame = pd.DataFrame(rows)
    cast_metadata(frame)
    for column in ("phi", "P", "fret", "p0_file", "total_growthrate"):
        frame[column] = frame[column].astype(float)
    frame = add_contact_observables(frame)
    return frame, bad_files


def load_bext_endpoints():
    rows = []
    for path in sorted((REPO_ROOT / "output" / "bext").glob("B_ext_data_dphiprobe*.dat")):
        match = BEXT_FILENAME_RE.match(path.name)
        if match is None:
            continue
        try:
            lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
            if len(lines) != 2:
                raise ValueError("unexpected B_ext line count")
            values = lines[1].split()
            if len(values) != 20:
                raise ValueError("unexpected B_ext column count")
            rows.append(
                {
                    **match.groupdict(),
                    "phi0": float(values[0]),
                    "P_measured": float(values[1]),
                    "phi1": float(values[2]),
                    "B_ext": float(values[5]),
                    "N": int(values[10]),
                    "Nc": int(values[11]),
                    "Nf": int(values[12]),
                    "Nu": int(values[13]),
                    "Ziso_total": int(values[14]),
                    "source_file": path.name,
                }
            )
        except (OSError, ValueError):
            continue
    if not rows:
        raise SystemExit("No B_ext rows found under output/bext")
    frame = pd.DataFrame(rows)
    cast_metadata(frame)
    for column in ("phi0", "phi1", "P_measured", "B_ext"):
        frame[column] = frame[column].astype(float)
    frame = add_contact_observables(frame)
    return frame


def estimate_shear_modulus(data):
    if data.ndim == 1:
        data = data[None, :]
    fit_window = data[data[:, 0] <= SHEAR_FIT_MAX_STRAIN]
    if len(fit_window) < 5:
        fit_window = data[:10]
    if len(fit_window) < 2:
        return float("nan")
    return float(np.polyfit(fit_window[:, 0], fit_window[:, 2], 1)[0])


def load_shear_endpoints():
    rows = []
    for path in sorted((REPO_ROOT / "output" / "shear").glob("G_data_*.dat")):
        match = G_FILENAME_RE.match(path.name)
        if match is None:
            continue
        try:
            data = np.loadtxt(path, comments="#")
        except (OSError, ValueError):
            continue
        if data.size == 0:
            continue
        if data.ndim == 1:
            data = data[None, :]
        if data.ndim != 2 or data.shape[1] != 12 or len(data) < 2 or not np.all(np.isfinite(data)):
            continue
        G = estimate_shear_modulus(data)
        if not np.isfinite(G):
            continue
        try:
            N_total = int(round(data[0, 3]))
            Nf = int(round(data[0, 5]))
            Nu = int(round(data[0, 6]))
            Ziso_total = int(round(data[0, 7]))
            Nc = int(round(data[0, 4]))
        except (TypeError, ValueError, OverflowError):
            continue
        rows.append(
            {
                **match.groupdict(),
                "phi": float(data[0, 8]),
                "Nc": Nc,
                "Nf": Nf,
                "Nu": Nu,
                "N": N_total,
                "Ziso_total": Ziso_total,
                "G": G,
                "source_file": path.name,
            }
        )
    if not rows:
        raise SystemExit("No shear trajectories found under output/shear")
    frame = pd.DataFrame(rows)
    cast_metadata(frame)
    frame["phi"] = frame["phi"].astype(float)
    frame["G"] = frame["G"].astype(float)
    frame = add_contact_observables(frame)
    return frame


def open_text_maybe_gzip(path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def wrap_periodic(delta, length):
    return delta - np.rint(delta / length) * length


def build_lobe_geometry(cell_data):
    x = cell_data[:, 0]
    y = cell_data[:, 1]
    diam = cell_data[:, 2]
    alpha = cell_data[:, 3]
    theta = cell_data[:, 4]
    dd = alpha - 1.0
    prefactor = (1.0 + dd) / (1.0 + dd**2)
    dr1 = prefactor * dd**2 * diam / 2.0
    dr2 = -prefactor * diam / 2.0
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    xa = np.column_stack((x + dr1 * cos_theta, x + dr2 * cos_theta))
    ya = np.column_stack((y + dr1 * sin_theta, y + dr2 * sin_theta))
    dk = np.column_stack((diam, dd * diam))
    return xa, ya, dk


def estimate_contact_m4_from_packing(path, lx, ly):
    with open_text_maybe_gzip(path) as handle:
        header = handle.readline().split()
        if len(header) < 2:
            raise ValueError(f"Malformed packing header in {path.name}")
        n_cells = int(header[0])
        cell_data = np.loadtxt(handle)
    if cell_data.ndim == 1:
        cell_data = cell_data[None, :]
    if len(cell_data) != n_cells or cell_data.shape[1] < 5:
        raise ValueError(f"Unexpected packing body in {path.name}")

    xa, ya, dk = build_lobe_geometry(cell_data[:, :5])
    terms = []
    contact_pairs = 0
    for i in range(n_cells - 1):
        for ki in range(2):
            xi = xa[i, ki]
            yi = ya[i, ki]
            di = dk[i, ki]
            for kj in range(2):
                dx = wrap_periodic(xi - xa[i + 1 :, kj], lx)
                dy = wrap_periodic(yi - ya[i + 1 :, kj], ly)
                dij = 0.5 * (di + dk[i + 1 :, kj])
                rsq = dx * dx + dy * dy
                mask = rsq < (dij * dij)
                if not np.any(mask):
                    continue
                dx = dx[mask]
                dy = dy[mask]
                rsq = rsq[mask]
                nonzero = rsq > 0.0
                if not np.any(nonzero):
                    continue
                dx = dx[nonzero]
                dy = dy[nonzero]
                rsq = rsq[nonzero]
                terms.append((dx * dx * dy * dy) / (rsq * rsq))
                contact_pairs += len(rsq)
    if not terms:
        return {"M4": np.nan, "contact_pairs_geom": 0, "packing_file": path.name}
    return {
        "M4": float(np.mean(np.concatenate(terms))),
        "contact_pairs_geom": int(contact_pairs),
        "packing_file": path.name,
    }


def estimate_contact_m4(frame):
    rows = []
    growth_dir = REPO_ROOT / "output" / "growth"
    for row in frame.itertuples(index=False):
        try:
            source_name = getattr(row, "packing_source_file", None)
            if source_name is None:
                source_name = row.source_file
            if source_name.startswith("G_data_"):
                packing_name = source_name.removeprefix("G_data_")
            elif source_name.startswith("NC_"):
                packing_name = "LF_DPHI_" + source_name.removeprefix("NC_")
            elif source_name.startswith("LF_DPHI_"):
                packing_name = source_name
            else:
                continue
            packing_path = growth_dir / packing_name
            if not packing_path.exists():
                gz_path = Path(str(packing_path) + ".gz")
                if gz_path.exists():
                    packing_path = gz_path
                else:
                    continue
            estimate = estimate_contact_m4_from_packing(packing_path, row.L, row.Ly)
            expected_pairs = int(round(row.Nc / 2))
            if estimate["contact_pairs_geom"] != expected_pairs:
                continue
            rows.append(
                {
                    "seed": row.seed,
                    "L": row.L,
                    "Ly": row.Ly,
                    "dphi": row.dphi,
                    "p0": row.p0,
                    "M4": estimate["M4"],
                    "contact_pairs_geom": estimate["contact_pairs_geom"],
                    "contact_pairs_saved": expected_pairs,
                    "packing_file": estimate["packing_file"],
                }
            )
        except (OSError, ValueError):
            continue
    if not rows:
        raise SystemExit("No valid shear packings with matching saved growth packings found for contact reconstruction")
    m4_frame = pd.DataFrame(rows)
    cast_metadata(m4_frame)
    m4_frame["M4"] = m4_frame["M4"].astype(float)
    return m4_frame


def load_growth_rate_distributions():
    rows = []
    histogram_rows = []
    for path in sorted((REPO_ROOT / "output" / "growth").glob("STATS_LF_DPHI_*.dat")):
        match = STATS_FILENAME_RE.match(path.name)
        if match is None:
            continue
        try:
            with path.open("r", encoding="utf-8") as handle:
                header = handle.readline().split()
                if len(header) != 3:
                    raise ValueError("unexpected STATS header")
                phi = float(header[1])
                bud_growth = []
                bud_pressure = []
                compressed_growth = []
                for line in handle:
                    values = line.split()
                    if not values:
                        continue
                    if len(values) != 8:
                        raise ValueError("unexpected STATS row width")
                    lobe_type = int(values[3])
                    if lobe_type != 1:
                        continue
                    pressure = float(values[4])
                    rate_eff = float(values[5])
                    rate_raw = float(values[6])
                    growth = rate_eff / rate_raw if rate_raw else np.nan
                    if not np.isfinite(growth):
                        continue
                    bud_pressure.append(pressure)
                    bud_growth.append(growth)
                    if pressure > COMPRESSED_PRESSURE_TOL:
                        compressed_growth.append(growth)
        except (OSError, ValueError):
            continue
        if not bud_growth:
            continue
        bud_growth = np.asarray(bud_growth)
        bud_pressure = np.asarray(bud_pressure)
        compressed_growth = np.asarray(compressed_growth)
        rows.append(
            {
                **match.groupdict(),
                "phi": phi,
                "n_buds": len(bud_growth),
                "n_compressed": int(len(compressed_growth)),
                "chi_c": float(np.nanmean(compressed_growth)) if len(compressed_growth) else 1.0,
                "g_mean": float(np.nanmean(bud_growth)),
                "g_var": float(np.nanvar(bud_growth)),
                "h_gamma": float(np.nanvar(bud_growth) / (np.nanmean(bud_growth) ** 2))
                if np.nanmean(bud_growth) > 0
                else np.nan,
                "frac_near1": float(np.mean(bud_growth > 0.98)),
                "g_p10": float(np.nanquantile(bud_growth, 0.10)),
                "g_p50": float(np.nanquantile(bud_growth, 0.50)),
                "g_p90": float(np.nanquantile(bud_growth, 0.90)),
                "p_p50": float(np.nanquantile(bud_pressure, 0.50)),
                "source_file": path.name,
            }
        )
        for value in bud_growth:
            histogram_rows.append({**match.groupdict(), "g": float(value), "phi": phi})
    if not rows:
        raise SystemExit("No bud-level STATS_LF_DPHI files found under output/growth")
    summary = pd.DataFrame(rows)
    cast_metadata(summary)
    for column in (
        "phi",
        "chi_c",
        "g_mean",
        "g_var",
        "h_gamma",
        "frac_near1",
        "g_p10",
        "g_p50",
        "g_p90",
        "p_p50",
    ):
        summary[column] = summary[column].astype(float)
    histogram = pd.DataFrame(histogram_rows)
    cast_metadata(histogram)
    histogram["g"] = histogram["g"].astype(float)
    histogram["phi"] = histogram["phi"].astype(float)
    return summary, histogram


def load_bgrow_local_slopes():
    rows = []
    for path in sorted((REPO_ROOT / "output" / "growth").glob("STEPLOG_*.dat")):
        match = STEPLOG_FILENAME_RE.match(path.name)
        if match is None:
            continue
        try:
            data = np.loadtxt(path, comments="#")
        except (OSError, ValueError):
            continue
        if data.size == 0:
            continue
        if data.ndim == 1:
            data = data[None, :]
        if data.ndim != 2 or data.shape[1] != 10 or not np.all(np.isfinite(data)):
            continue
        postjam = data[data[:, 9] > 0.5]
        if len(postjam) < 5:
            continue
        fit_window = postjam[-25:]
        slope = float(np.polyfit(fit_window[:, 5], fit_window[:, 3], 1)[0])
        if not np.isfinite(slope):
            continue
        phi_final = float(fit_window[-1, 5])
        rows.append(
            {
                **match.groupdict(),
                "phi_local": phi_final,
                "P_local": float(fit_window[-1, 3]),
                "B_grow": phi_final * slope,
                "source_file": path.name,
            }
        )
    if not rows:
        raise SystemExit("No usable STEPLOG files found for local B_grow estimation")
    frame = pd.DataFrame(rows)
    cast_metadata(frame)
    for column in ("phi_local", "P_local", "B_grow"):
        frame[column] = frame[column].astype(float)
    return frame


def filter_main(frame, main_L):
    return frame[(frame["L"] == main_L) & (frame["p0"].isin(DEFAULT_MAIN_P0S))].copy()


def summarize_curves(frame, value_column):
    grouped = frame.groupby(["p0", "dphi"], as_index=False)[value_column].agg(
        median="median",
        q25=lambda s: s.quantile(0.25),
        q75=lambda s: s.quantile(0.75),
    )
    return grouped


def build_curve_summary(growth_final, bext_main, shear_main):
    pressure = summarize_curves(growth_final, "P").rename(
        columns={"median": "P_median", "q25": "P_q25", "q75": "P_q75"}
    )
    phi = summarize_curves(growth_final, "phi").rename(
        columns={"median": "phi_median", "q25": "phi_q25", "q75": "phi_q75"}
    )
    free_bud = summarize_curves(growth_final, "u").rename(
        columns={"median": "u_median", "q25": "u_q25", "q75": "u_q75"}
    )
    delta_z = summarize_curves(growth_final, "DeltaZ").rename(
        columns={"median": "DeltaZ_median", "q25": "DeltaZ_q25", "q75": "DeltaZ_q75"}
    )
    bext = summarize_curves(bext_main, "B_ext").rename(
        columns={"median": "Bext_median", "q25": "Bext_q25", "q75": "Bext_q75"}
    )
    shear = summarize_curves(shear_main, "G").rename(
        columns={"median": "G_median", "q25": "G_q25", "q75": "G_q75"}
    )
    curves = pressure.merge(phi, on=["p0", "dphi"])
    curves = curves.merge(free_bud, on=["p0", "dphi"])
    curves = curves.merge(delta_z, on=["p0", "dphi"])
    curves = curves.merge(bext, on=["p0", "dphi"])
    curves = curves.merge(shear, on=["p0", "dphi"])
    return curves.sort_values(["p0", "dphi"]).reset_index(drop=True)


def build_structural_branch(growth_main):
    jam = (
        growth_main[growth_main["stage"] == "jamm"][["seed", "L", "dphi", "p0", "u"]]
        .rename(columns={"u": "u_jamm"})
        .copy()
    )
    final = growth_main[growth_main["stage"] == "final"].copy()
    merged = final.merge(jam, on=["seed", "L", "dphi", "p0"], how="inner")

    no_feedback = merged[merged["p0"] == -1.0]
    baseline = (
        no_feedback.groupby("dphi", as_index=False)[["P", "DeltaZ"]]
        .median()
        .sort_values("P")
        .reset_index(drop=True)
    )
    work = merged[merged["p0"] != -1.0].copy()
    work["DeltaZ0_interp"] = np.interp(
        work["P"].to_numpy(),
        baseline["P"].to_numpy(),
        baseline["DeltaZ"].to_numpy(),
        left=baseline["DeltaZ"].iloc[0],
        right=baseline["DeltaZ"].iloc[-1],
    )
    work["reservoir_depletion"] = work["u_jamm"] - work["u"]
    work["structural_excess"] = work["DeltaZ"] - work["DeltaZ0_interp"]

    median_points = (
        work.groupby(["p0", "dphi"], as_index=False)[["reservoir_depletion", "structural_excess"]]
        .median()
        .sort_values(["p0", "dphi"])
        .reset_index(drop=True)
    )
    fit_slope = float(
        np.dot(work["reservoir_depletion"], work["structural_excess"])
        / np.dot(work["reservoir_depletion"], work["reservoir_depletion"])
    )
    return work, median_points, fit_slope


def linear_fit(x_values, y_values):
    slope, intercept = np.polyfit(x_values, y_values, 1)
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "corr": float(np.corrcoef(x_values, y_values)[0, 1]),
    }


def build_moduli_coefficients(bext_main, shear_main, growth_final):
    bext_work = bext_main.copy()
    bext_work["B_over_Z"] = bext_work["B_ext"] / bext_work["Z"]
    bext_work["body_from_B"] = 8.0 * bext_work["B_over_Z"]

    shear_work = shear_main[["seed", "L", "Ly", "dphi", "p0", "G", "source_file"]].rename(
        columns={"source_file": "shear_source_file"}
    )
    shear_work = shear_work.merge(
        growth_final[
            [
                "seed",
                "L",
                "Ly",
                "dphi",
                "p0",
                "N",
                "Ziso_total",
                "Nc",
                "Nf",
                "Nu",
                "phi",
                "u",
                "Z",
                "DeltaZ",
                "source_file",
            ]
        ].rename(columns={"source_file": "packing_source_file"}),
        on=["seed", "L", "Ly", "dphi", "p0"],
        how="inner",
    )
    shear_work = shear_work.merge(
        estimate_contact_m4(shear_work),
        on=["seed", "L", "Ly", "dphi", "p0"],
        how="inner",
    )
    shear_work = shear_work[(shear_work["DeltaZ"] > MIN_RATIO_DELTA_Z) & (shear_work["M4"] > 0)].copy()
    shear_work["G_over_DeltaZ"] = shear_work["G"] / shear_work["DeltaZ"]
    shear_work["M4DeltaZ"] = shear_work["M4"] * shear_work["DeltaZ"]
    shear_work["G_over_M4DeltaZ"] = shear_work["G"] / shear_work["M4DeltaZ"]
    shear_work["body_from_G"] = 2.0 * shear_work["G_over_M4DeltaZ"]

    b_fit = linear_fit(bext_work["Z"].to_numpy(), bext_work["B_ext"].to_numpy())
    g_fit = linear_fit(shear_work["DeltaZ"].to_numpy(), shear_work["G"].to_numpy())
    g_m4_fit = linear_fit(shear_work["M4DeltaZ"].to_numpy(), shear_work["G"].to_numpy())
    b_fit["body_coefficient"] = 8.0 * b_fit["slope"]
    g_fit["body_coefficient"] = np.nan
    g_m4_fit["body_coefficient"] = 2.0 * g_m4_fit["slope"]

    coeff_curves = summarize_curves(bext_work, "body_from_B").rename(
        columns={"median": "body_B_median", "q25": "body_B_q25", "q75": "body_B_q75"}
    )
    g_over_m4_curves = summarize_curves(shear_work, "G_over_M4DeltaZ").rename(
        columns={"median": "G_over_M4DeltaZ_median", "q25": "G_over_M4DeltaZ_q25", "q75": "G_over_M4DeltaZ_q75"}
    )
    coeff_curves = coeff_curves.merge(
        summarize_curves(shear_work, "body_from_G").rename(
            columns={"median": "body_G_median", "q25": "body_G_q25", "q75": "body_G_q75"}
        ),
        on=["p0", "dphi"],
        how="inner",
    )
    pooled_body = float(
        np.nanmedian(np.concatenate([bext_work["body_from_B"].to_numpy(), shear_work["body_from_G"].to_numpy()]))
    )
    pooled_g_over_m4 = float(np.nanmedian(shear_work["G_over_M4DeltaZ"].to_numpy()))
    pooled_m4 = float(np.nanmedian(shear_work["M4"].to_numpy()))
    return (
        bext_work,
        shear_work,
        coeff_curves,
        g_over_m4_curves,
        b_fit,
        g_fit,
        g_m4_fit,
        pooled_body,
        pooled_g_over_m4,
        pooled_m4,
    )


def build_bgrow_partition(growth_final, bext_main, growth_moments_main, bgrow_local):
    merged = bgrow_local.merge(
        growth_final[["seed", "L", "Ly", "dphi", "p0", "u"]],
        on=["seed", "L", "Ly", "dphi", "p0"],
        how="inner",
    ).merge(
        bext_main[["seed", "L", "Ly", "dphi", "p0", "B_ext"]],
        on=["seed", "L", "Ly", "dphi", "p0"],
        how="inner",
    ).merge(
        growth_moments_main[["seed", "L", "Ly", "dphi", "p0", "chi_c"]],
        on=["seed", "L", "Ly", "dphi", "p0"],
        how="inner",
    )
    merged["Bgrow_over_Bext"] = merged["B_grow"] / merged["B_ext"]
    merged["lambda_c"] = ((1.0 - merged["u"]) * merged["chi_c"]) / (
        merged["u"] + (1.0 - merged["u"]) * merged["chi_c"]
    )

    ratio_curves = summarize_curves(merged, "chi_c").rename(
        columns={"median": "chi_c_median", "q25": "chi_c_q25", "q75": "chi_c_q75"}
    )
    ratio_curves = ratio_curves.merge(
        summarize_curves(merged, "lambda_c").rename(
            columns={"median": "lambda_c_median", "q25": "lambda_c_q25", "q75": "lambda_c_q75"}
        ),
        on=["p0", "dphi"],
        how="inner",
    )
    median_curve_rows = []
    for p0, subset in growth_final.groupby("p0"):
        subset = subset.sort_values("phi").copy()
        if len(subset) < 2:
            continue
        curve = (
            subset.groupby("dphi", as_index=False)[["phi", "P"]]
            .median()
            .sort_values("phi")
            .reset_index(drop=True)
        )
        curve["B_grow_curve"] = curve["phi"] * np.gradient(curve["P"].to_numpy(), curve["phi"].to_numpy())
        curve["p0"] = p0
        median_curve_rows.append(curve[["p0", "dphi", "B_grow_curve"]])
    if not median_curve_rows:
        raise SystemExit("Unable to construct median B_grow curves from the growth endpoints")
    median_curve = pd.concat(median_curve_rows, ignore_index=True)
    median_curve = median_curve.merge(
        summarize_curves(bext_main, "B_ext").rename(
            columns={"median": "Bext_median", "q25": "Bext_q25", "q75": "Bext_q75"}
        ),
        on=["p0", "dphi"],
        how="inner",
    )
    median_curve["ratio_curve"] = median_curve["B_grow_curve"] / median_curve["Bext_median"]
    ratio_curves = ratio_curves.merge(
        median_curve[["p0", "dphi", "ratio_curve"]],
        on=["p0", "dphi"],
        how="inner",
    )
    ratio_curves = ratio_curves.merge(
        summarize_curves(merged, "chi_c").rename(
            columns={"median": "chi_c_point_median", "q25": "chi_c_point_q25", "q75": "chi_c_point_q75"}
        ),
        on=["p0", "dphi"],
        how="inner",
    )
    ratio_curves = ratio_curves.merge(
        summarize_curves(merged, "Bgrow_over_Bext").rename(
            columns={"median": "ratio_point_median", "q25": "ratio_point_q25", "q75": "ratio_point_q75"}
        ),
        on=["p0", "dphi"],
        how="inner",
    )

    stats = {
        "corr_chi_c": float(np.corrcoef(ratio_curves["ratio_curve"], ratio_curves["chi_c_median"])[0, 1]),
        "corr_lambda_c": float(np.corrcoef(ratio_curves["ratio_curve"], ratio_curves["lambda_c_median"])[0, 1]),
        "rmse_chi_c": float(np.sqrt(np.mean((ratio_curves["ratio_curve"] - ratio_curves["chi_c_median"]) ** 2))),
        "rmse_lambda_c": float(
            np.sqrt(np.mean((ratio_curves["ratio_curve"] - ratio_curves["lambda_c_median"]) ** 2))
        ),
    }
    return merged, ratio_curves, stats


def build_growth_histogram_slice(histogram, growth_moments, main_L, hist_dphi):
    hist_slice = histogram[
        (histogram["L"] == main_L)
        & (np.isclose(histogram["dphi"], hist_dphi))
        & (histogram["p0"].isin(DEFAULT_HIST_P0S))
    ].copy()
    moment_slice = growth_moments[
        (growth_moments["L"] == main_L)
        & (growth_moments["p0"].isin(DEFAULT_MAIN_P0S))
    ].copy()
    moment_curves = summarize_curves(moment_slice, "h_gamma").rename(
        columns={"median": "h_gamma_median", "q25": "h_gamma_q25", "q75": "h_gamma_q75"}
    )
    return hist_slice, moment_curves


def add_line_with_band(axis, subset, value, low, high, color, label=None, logx=True):
    subset = subset.sort_values("dphi")
    x = subset["dphi"].to_numpy()
    y = subset[value].to_numpy()
    axis.plot(x, y, color=color, linewidth=2.3, marker="o", markersize=4, label=label)
    axis.fill_between(
        x,
        subset[low].to_numpy(),
        subset[high].to_numpy(),
        color=color,
        alpha=0.16,
        linewidth=0,
    )
    if logx:
        axis.set_xscale("log")
    axis.grid(True, alpha=0.25)


def plot_structure_branch(curves, structural_points, structural_medians, fit_slope, output_dir, formats):
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8), constrained_layout=True)

    axis = axes[0]
    for p0 in p0_order(curves["p0"].unique()):
        subset = curves[curves["p0"] == p0]
        add_line_with_band(
            axis,
            subset,
            "u_median",
            "u_q25",
            "u_q75",
            color=P0_COLORS[p0],
            label=p0_label(p0),
        )
    axis.set_xlabel(r"$\Delta \phi$")
    axis.set_ylabel(r"Free-bud fraction $u$")
    axis.set_title("Reservoir Depletion")
    axis.legend(frameon=False, fontsize=9)

    axis = axes[1]
    max_x = max(structural_points["reservoir_depletion"].max(), structural_medians["reservoir_depletion"].max())
    xline = np.linspace(0.0, max_x * 1.05, 100)
    axis.plot(xline, 2.0 * xline, color="#111827", linestyle="--", linewidth=1.5, label=r"Theory: $2(u_J-u)$")
    for p0 in p0_order(structural_points["p0"].unique()):
        raw = structural_points[structural_points["p0"] == p0]
        med = structural_medians[structural_medians["p0"] == p0]
        axis.scatter(
            raw["reservoir_depletion"],
            raw["structural_excess"],
            s=18,
            color=P0_COLORS[p0],
            alpha=0.18,
            linewidths=0,
        )
        axis.scatter(
            med["reservoir_depletion"],
            med["structural_excess"],
            s=54,
            color=P0_COLORS[p0],
            edgecolor="white",
            linewidth=0.8,
            label=p0_label(p0),
        )
    axis.text(
        0.04,
        0.95,
        rf"Through-origin fit slope $\approx {fit_slope:.2f}$",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    axis.set_xlabel(r"Reservoir depletion $u_J - u$")
    axis.set_ylabel(r"Structural excess $\Delta Z - \Delta Z_0(P)$")
    axis.set_title("Counting-Law Collapse")
    axis.grid(True, alpha=0.25)
    axis.legend(frameon=False, fontsize=8, loc="lower right")

    fig.suptitle("Structure Above Jamming", fontsize=14)
    save_figure(fig, output_dir, "fig_structure_branch", formats)


def plot_moduli_coefficients(
    bext_work,
    shear_work,
    coeff_curves,
    g_over_m4_curves,
    b_fit,
    g_fit,
    g_m4_fit,
    pooled_body,
    pooled_g_over_m4,
    pooled_m4,
    output_dir,
    formats,
):
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 9.4), constrained_layout=True)

    axis = axes[0, 0]
    xline = np.linspace(bext_work["Z"].min() * 0.95, bext_work["Z"].max() * 1.03, 100)
    axis.plot(
        xline,
        b_fit["slope"] * xline + b_fit["intercept"],
        color="#111827",
        linestyle="--",
        linewidth=1.6,
        label="Linear fit",
    )
    for p0 in p0_order(bext_work["p0"].unique()):
        subset = bext_work[bext_work["p0"] == p0]
        axis.scatter(
            subset["Z"],
            subset["B_ext"],
            s=22,
            color=P0_COLORS[p0],
            alpha=0.35,
            linewidths=0,
            label=p0_label(p0),
        )
    axis.text(
        0.03,
        0.96,
        "\n".join(
            [
                rf"$B_{{\mathrm{{ext}}}} \approx {b_fit['slope']:.3f} Z {b_fit['intercept']:+.3f}$",
                rf"$8a_B \approx {b_fit['body_coefficient']:.3f}$",
                rf"$r \approx {b_fit['corr']:.3f}$",
            ]
        ),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    axis.set_xlabel(r"Coordination $Z$")
    axis.set_ylabel(r"External bulk modulus $B_{\mathrm{ext}}$")
    axis.set_title(r"Fit: $B_{\mathrm{ext}}$ vs $Z$")
    axis.grid(True, alpha=0.25)
    axis.legend(frameon=False, fontsize=8, loc="upper left")

    axis = axes[0, 1]
    xline = np.linspace(0.0, shear_work["DeltaZ"].max() * 1.03, 100)
    axis.plot(
        xline,
        g_fit["slope"] * xline + g_fit["intercept"],
        color="#111827",
        linestyle="--",
        linewidth=1.6,
        label="Linear fit",
    )
    for p0 in p0_order(shear_work["p0"].unique()):
        subset = shear_work[shear_work["p0"] == p0]
        axis.scatter(
            subset["DeltaZ"],
            subset["G"],
            s=22,
            color=P0_COLORS[p0],
            alpha=0.35,
            linewidths=0,
            label=p0_label(p0),
        )
    axis.text(
        0.03,
        0.96,
        "\n".join(
            [
                rf"$G \approx {g_fit['slope']:.3f}\,\Delta Z {g_fit['intercept']:+.3f}$",
                rf"Fit to $M_4\Delta Z$: $2a_G \approx {g_m4_fit['body_coefficient']:.3f}$",
                rf"$r \approx {g_fit['corr']:.3f}$",
            ]
        ),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    axis.set_xlabel(r"Excess coordination $\Delta Z$")
    axis.set_ylabel(r"Shear modulus $G$")
    axis.set_title(r"Fit: $G$ vs $\Delta Z$")
    axis.grid(True, alpha=0.25)
    axis.legend(frameon=False, fontsize=8, loc="upper left")

    axis = axes[1, 0]
    for p0 in p0_order(g_over_m4_curves["p0"].unique()):
        subset = g_over_m4_curves[g_over_m4_curves["p0"] == p0]
        add_line_with_band(
            axis,
            subset,
            "G_over_M4DeltaZ_median",
            "G_over_M4DeltaZ_q25",
            "G_over_M4DeltaZ_q75",
            color=P0_COLORS[p0],
            label=p0_label(p0),
        )
    axis.axhline(pooled_g_over_m4, color="#111827", linestyle="--", linewidth=1.4)
    axis.text(
        0.03,
        0.95,
        "\n".join(
            [
                rf"Pooled median $\approx {pooled_g_over_m4:.3f}$",
                rf"Measured median $M_4 \approx {pooled_m4:.3f}$",
            ]
        ),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    axis.set_xlabel(r"$\Delta \phi$")
    axis.set_ylabel(r"$G / (M_4 \Delta Z)$")
    axis.set_title(r"Constancy of $G / (M_4 \Delta Z)$")
    axis.legend(frameon=False, fontsize=8, loc="upper left")

    axis = axes[1, 1]
    for p0 in p0_order(coeff_curves["p0"].unique()):
        subset = coeff_curves[coeff_curves["p0"] == p0].sort_values("dphi")
        axis.plot(
            subset["dphi"],
            subset["body_B_median"],
            color=P0_COLORS[p0],
            linewidth=2.1,
            marker="o",
            markersize=4,
        )
        axis.plot(
            subset["dphi"],
            subset["body_G_median"],
            color=P0_COLORS[p0],
            linewidth=1.9,
            linestyle="--",
            marker="s",
            markersize=3.5,
        )
    axis.axhline(pooled_body, color="#111827", linestyle=":", linewidth=1.4)
    axis.set_xscale("log")
    axis.set_xlabel(r"$\Delta \phi$")
    axis.set_ylabel(r"Estimated common factor $n_J k \ell_c^2$")
    axis.set_title("Common Coefficient-Body Comparison")
    axis.grid(True, alpha=0.25)
    axis.text(
        0.03,
        0.95,
        "\n".join(
            [
                r"Solid: $8B_{\mathrm{ext}}/Z$",
                r"Dashed: $2G/(M_4\Delta Z)$",
                rf"Pooled median $\approx {pooled_body:.3f}$",
            ]
        ),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    color_handles = [
        Line2D([0], [0], color=P0_COLORS[p0], linewidth=2.2, label=p0_label(p0))
        for p0 in p0_order(coeff_curves["p0"].unique())
    ]
    if color_handles:
        axis.legend(handles=color_handles, frameon=False, fontsize=8, loc="lower right")

    fig.suptitle("Elastic Moduli and Coefficient Consistency", fontsize=14)
    save_figure(fig, output_dir, "fig_moduli_coefficients", formats)


def plot_decoupling(curves, output_dir, formats):
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.6), constrained_layout=True)
    specs = [
        ("P_median", "P_q25", "P_q75", r"Pressure $P$", "Pressure Cost", True),
        ("Bext_median", "Bext_q25", "Bext_q75", r"Bulk modulus $B_{\mathrm{ext}}$", "External Compression", False),
        ("G_median", "G_q25", "G_q75", r"Shear modulus $G$", "Shear Stiffening", False),
    ]
    for axis, (value, low, high, ylabel, title, logy) in zip(axes, specs):
        for p0 in p0_order(curves["p0"].unique()):
            subset = curves[curves["p0"] == p0]
            add_line_with_band(
                axis,
                subset,
                value,
                low,
                high,
                color=P0_COLORS[p0],
                label=p0_label(p0),
            )
        axis.set_xlabel(r"$\Delta \phi$")
        axis.set_ylabel(ylabel)
        axis.set_title(title)
        if logy:
            axis.set_yscale("log")
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.05))
    fig.suptitle("Feedback Decouples Pressure Build-Up from Stiffness Growth", fontsize=14)
    save_figure(fig, output_dir, "fig_pressure_moduli_decoupling", formats)


def plot_bgrow_partition(partition_curves, partition_raw, partition_stats, output_dir, formats):
    feedback_curves = partition_curves[partition_curves["p0"].isin(DEFAULT_MAIN_P0S)].copy()
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2), constrained_layout=True)

    axis = axes[0]
    for p0 in p0_order(feedback_curves["p0"].unique()):
        subset = feedback_curves[feedback_curves["p0"] == p0]
        axis.plot(
            subset["dphi"],
            subset["ratio_curve"],
            color=P0_COLORS[p0],
            linewidth=2.3,
            marker="o",
            markersize=4,
            label=p0_label(p0),
        )
        axis.plot(
            subset["dphi"],
            subset["lambda_c_median"],
            color=P0_COLORS[p0],
            linewidth=1.7,
            linestyle="--",
        )
        axis.plot(
            subset["dphi"],
            subset["chi_c_median"],
            color=P0_COLORS[p0],
            linewidth=1.5,
            linestyle=":",
        )
    axis.set_xlabel(r"$\Delta \phi$")
    axis.set_ylabel(r"Ratio / partition factor")
    axis.set_title(r"$B_{\mathrm{grow}} / B_{\mathrm{ext}}$ vs $\lambda_c$ and $\chi_c$")
    axis.set_xscale("log")
    axis.grid(True, alpha=0.25)
    color_handles = [
        Line2D([0], [0], color=P0_COLORS[p0], linewidth=2.2, label=p0_label(p0))
        for p0 in p0_order(feedback_curves["p0"].unique())
    ]
    style_handles = [
        Line2D([0], [0], color="#111827", linewidth=2.2, linestyle="-", label=r"$B_{\mathrm{grow}}/B_{\mathrm{ext}}$"),
        Line2D([0], [0], color="#111827", linewidth=1.8, linestyle="--", label=r"$\lambda_c$"),
        Line2D([0], [0], color="#111827", linewidth=1.6, linestyle=":", label=r"$\chi_c$"),
    ]
    legend1 = axis.legend(handles=color_handles, frameon=False, fontsize=8, loc="upper left")
    axis.add_artist(legend1)
    axis.legend(handles=style_handles, frameon=False, fontsize=8, loc="lower right")

    axis = axes[1]
    axis.scatter(
        feedback_curves["lambda_c_median"],
        feedback_curves["ratio_curve"],
        s=22,
        c=feedback_curves["p0"].map(P0_COLORS),
        alpha=0.75,
        linewidths=0.0,
        label=r"$\lambda_c$",
    )
    axis.scatter(
        feedback_curves["chi_c_median"],
        feedback_curves["ratio_curve"],
        s=26,
        facecolors="none",
        edgecolors=feedback_curves["p0"].map(P0_COLORS),
        alpha=0.85,
        linewidths=0.9,
        marker="^",
        label=r"$\chi_c$",
    )
    max_value = max(
        feedback_curves["ratio_curve"].max(),
        feedback_curves["chi_c_median"].max(),
        feedback_curves["lambda_c_median"].max(),
    )
    xline = np.linspace(0.0, max_value * 1.02, 100)
    axis.plot(xline, xline, color="#111827", linestyle="--", linewidth=1.4)
    axis.text(
        0.03,
        0.96,
        "\n".join(
            [
                rf"corr with $\lambda_c$: {partition_stats['corr_lambda_c']:.3f}",
                rf"corr with $\chi_c$: {partition_stats['corr_chi_c']:.3f}",
                rf"RMSE$_{{\lambda_c}}$: {partition_stats['rmse_lambda_c']:.3f}",
                rf"RMSE$_{{\chi_c}}$: {partition_stats['rmse_chi_c']:.3f}",
            ]
        ),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    axis.set_xlabel(r"Proxy factor ($\lambda_c$ filled, $\chi_c$ open)")
    axis.set_ylabel(r"$B_{\mathrm{grow}} / B_{\mathrm{ext}}$")
    axis.set_title("Direct Partition Test")
    axis.grid(True, alpha=0.25)
    axis.legend(frameon=False, fontsize=8, loc="lower right")

    fig.suptitle("Flux Partition Between Growth and External Compression", fontsize=14)
    save_figure(fig, output_dir, "fig_bgrow_partition", formats)


def plot_growth_crossover(hist_slice, moment_curves, hist_dphi, output_dir, formats):
    hist_p0s = [value for value in p0_order(hist_slice["p0"].unique()) if value > 0.0]
    fig = plt.figure(figsize=(max(13.2, 3.2 * len(hist_p0s)), 7.0), constrained_layout=True)
    grid = gridspec.GridSpec(2, len(hist_p0s), figure=fig, height_ratios=[1.0, 0.9])
    bins = np.linspace(0.0, 1.0, 33)

    for index, p0 in enumerate(hist_p0s):
        axis = fig.add_subplot(grid[0, index])
        subset = hist_slice[hist_slice["p0"] == p0]["g"].to_numpy()
        axis.hist(
            subset,
            bins=bins,
            density=True,
            histtype="stepfilled",
            alpha=0.55,
            linewidth=1.2,
            color=P0_COLORS[p0],
            edgecolor=P0_COLORS[p0],
        )
        axis.axvline(1.0, color="#111827", linestyle="--", linewidth=1.2)
        axis.set_xlim(0.0, 1.02)
        axis.set_xlabel(r"Normalized growth $g = \gamma / \gamma_0$")
        axis.set_ylabel("Density")
        axis.set_title(f"{p0_label(p0)} at " + rf"$\Delta \phi = {hist_dphi:.2f}$")
        axis.grid(True, alpha=0.20)
        axis.text(
            0.03,
            0.94,
            f"pooled buds = {len(subset)}",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8,
        )

    axis = fig.add_subplot(grid[1, :])
    for p0 in p0_order(moment_curves["p0"].unique()):
        subset = moment_curves[moment_curves["p0"] == p0]
        add_line_with_band(
            axis,
            subset,
            "h_gamma_median",
            "h_gamma_q25",
            "h_gamma_q75",
            color=P0_COLORS[p0],
            label=p0_label(p0),
        )
    axis.set_xlabel(r"$\Delta \phi$")
    axis.set_ylabel(r"Heterogeneity $h_\gamma$")
    axis.set_title("Derived Growth-Rate Heterogeneity")
    axis.legend(frameon=False, ncol=4, fontsize=9, loc="upper left")

    fig.suptitle("Growth-Rate Distribution Crossover", fontsize=14)
    save_figure(fig, output_dir, "fig_growth_crossover", formats)


def write_claim_map(output_dir, main_L, hist_dphi, skipped_nc_exists):
    lines = [
        "# Paper-support plot map",
        "",
        "Supported now with existing growth/shear/B_ext outputs:",
        "- `fig_structure_branch`: depletion of the free-bud reservoir and the counting-law relation `DeltaZ - DeltaZ0(P) ~ 2(u_J - u)`.",
        "- `fig_moduli_coefficients`: fitted `B_ext ~ Z` and `G ~ DeltaZ` relations, plus coefficient-body consistency checks.",
        "- `fig_pressure_moduli_decoupling`: suppressed pressure build-up under strong feedback alongside continued growth of `B_ext` and `G`.",
        "- `fig_bgrow_partition`: `B_grow / B_ext` compared against both `chi_c` and the full partition factor `lambda_c` extracted from the saved bud-level rates.",
        f"- `fig_growth_crossover`: bimodal-to-unimodal growth-rate crossover at `L = {main_L}` and `dphi = {hist_dphi:.2f}`, plus the heterogeneity diagnostic `h_gamma`.",
        "",
        "Not yet supported by the current output folder:",
        "- Threshold-distribution objects such as `S(a*)` or direct completion-threshold histograms.",
        "- Secondary-arrest quantities that need lineage/depletion runs (`phi_2(P0)`, `P_2(P0)`, `B_2(P0)`, `G_2(P0)`).",
        "- Any direct test that needs the dedicated lineage outputs (`POSTJAMM_SUMMARY_*`, `TRANSITIONS_*`, division-injection terms), because `output/lineage_growth/` is still empty in this workspace.",
    ]
    if skipped_nc_exists:
        lines.extend(
            [
                "",
                "Data note:",
                "- Some `NC_*` files were incomplete and were skipped. The skipped list is in `skipped_nc_files.txt` in this folder.",
            ]
        )
    (output_dir / "claim_map.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = parse_args()
    formats = plot_formats(args.formats)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.style.use("seaborn-v0_8-whitegrid")
    matplotlib.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "#fbfbf8",
            "axes.edgecolor": "#2f2f2f",
            "axes.labelcolor": "#1f1f1f",
            "xtick.color": "#1f1f1f",
            "ytick.color": "#1f1f1f",
            "grid.color": "#9aa5b1",
            "font.size": 11,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
        }
    )

    growth, skipped_nc_files = load_growth_endpoints()
    bext = load_bext_endpoints()
    shear = load_shear_endpoints()
    growth_moments, histogram = load_growth_rate_distributions()
    bgrow_local = load_bgrow_local_slopes()

    growth_main = filter_main(growth, args.main_L)
    bext_main = filter_main(bext, args.main_L)
    shear_main = filter_main(shear, args.main_L)
    growth_moments_main = filter_main(growth_moments, args.main_L)

    growth_final = growth_main[growth_main["stage"] == "final"].copy()
    curves = build_curve_summary(growth_final, bext_main, shear_main)
    structural_points, structural_medians, fit_slope = build_structural_branch(growth_main)
    (
        bext_coeff,
        shear_coeff,
        coeff_curves,
        g_over_m4_curves,
        b_fit,
        g_fit,
        g_m4_fit,
        pooled_body,
        pooled_g_over_m4,
        pooled_m4,
    ) = build_moduli_coefficients(bext_main, shear_main, growth_final)
    partition_raw, partition_curves, partition_stats = build_bgrow_partition(
        growth_final, bext_main, growth_moments_main, filter_main(bgrow_local, args.main_L)
    )
    hist_slice, moment_curves = build_growth_histogram_slice(
        histogram, growth_moments, args.main_L, args.hist_dphi
    )

    if skipped_nc_files:
        (output_dir / "skipped_nc_files.txt").write_text(
            "\n".join(skipped_nc_files) + "\n", encoding="utf-8"
        )
    if hist_slice.empty:
        raise SystemExit(
            f"No bud-growth histogram slice found for L={args.main_L} and dphi={args.hist_dphi:.6g}"
        )

    save_dataframe(growth_main, output_dir / f"growth_endpoints_L{args.main_L}_main.csv")
    save_dataframe(bext_main, output_dir / f"bext_endpoints_L{args.main_L}_main.csv")
    save_dataframe(shear_main, output_dir / f"shear_endpoints_L{args.main_L}_main.csv")
    save_dataframe(curves, output_dir / f"curve_summary_L{args.main_L}_main.csv")
    save_dataframe(structural_points, output_dir / f"structural_branch_points_L{args.main_L}_main.csv")
    save_dataframe(bext_coeff, output_dir / f"bext_coefficients_L{args.main_L}_main.csv")
    save_dataframe(shear_coeff, output_dir / f"shear_coefficients_L{args.main_L}_main.csv")
    save_dataframe(coeff_curves, output_dir / f"moduli_coefficient_curves_L{args.main_L}_main.csv")
    save_dataframe(partition_raw, output_dir / f"bgrow_partition_points_L{args.main_L}_main.csv")
    save_dataframe(partition_curves, output_dir / f"bgrow_partition_curves_L{args.main_L}_main.csv")
    save_dataframe(hist_slice, output_dir / f"growth_histogram_slice_L{args.main_L}_dphi{args.hist_dphi:.2f}.csv")
    save_dataframe(moment_curves, output_dir / f"growth_moments_L{args.main_L}_main.csv")
    save_dataframe(
        pd.DataFrame(
            [
                {
                    "quantity": "B_ext_vs_Z",
                    "slope": b_fit["slope"],
                    "intercept": b_fit["intercept"],
                    "corr": b_fit["corr"],
                    "body_coefficient": b_fit["body_coefficient"],
                    "aux_value": np.nan,
                    "aux_label": "",
                },
                {
                    "quantity": "G_vs_DeltaZ",
                    "slope": g_fit["slope"],
                    "intercept": g_fit["intercept"],
                    "corr": g_fit["corr"],
                    "body_coefficient": g_fit["body_coefficient"],
                    "aux_value": np.nan,
                    "aux_label": "",
                },
                {
                    "quantity": "G_vs_M4DeltaZ",
                    "slope": g_m4_fit["slope"],
                    "intercept": g_m4_fit["intercept"],
                    "corr": g_m4_fit["corr"],
                    "body_coefficient": g_m4_fit["body_coefficient"],
                    "aux_value": pooled_m4,
                    "aux_label": "median_M4",
                },
                {
                    "quantity": "B_grow_partition",
                    "slope": np.nan,
                    "intercept": np.nan,
                    "corr": partition_stats["corr_lambda_c"],
                    "body_coefficient": np.nan,
                    "aux_value": partition_stats["corr_chi_c"],
                    "aux_label": "corr_chi_c",
                },
            ]
        ),
        output_dir / "fit_summary.csv",
    )

    plot_structure_branch(
        curves,
        structural_points,
        structural_medians,
        fit_slope,
        output_dir,
        formats,
    )
    plot_moduli_coefficients(
        bext_coeff,
        shear_coeff,
        coeff_curves,
        g_over_m4_curves,
        b_fit,
        g_fit,
        g_m4_fit,
        pooled_body,
        pooled_g_over_m4,
        pooled_m4,
        output_dir,
        formats,
    )
    plot_decoupling(curves, output_dir, formats)
    plot_bgrow_partition(partition_curves, partition_raw, partition_stats, output_dir, formats)
    plot_growth_crossover(hist_slice, moment_curves, args.hist_dphi, output_dir, formats)

    write_claim_map(
        output_dir,
        args.main_L,
        args.hist_dphi,
        skipped_nc_exists=bool(skipped_nc_files),
    )

    print(
        "Generated manuscript-support plots "
        f"for L={args.main_L} with {len(curves)} curve rows, "
        f"{len(structural_points)} structural points, "
        f"{len(partition_raw)} B_grow/B_ext points, "
        f"and {len(hist_slice)} pooled bud-growth samples.",
        flush=True,
    )


if __name__ == "__main__":
    main()
