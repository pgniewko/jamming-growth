#!/usr/bin/env python3

import argparse
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
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from run_lineage_growth import (
    FIXED as LINEAGE_FIXED,
    parse_divlog,
    parse_lineage,
    parse_lineage_snapshot,
    parse_log,
    parse_postjamm_summary,
    parse_transitions,
)

DEFAULT_INPUT_DIR = REPO_ROOT / "output" / "lineage_growth" / "growth"
DEFAULT_LOG_DIR = REPO_ROOT / "output" / "logs" / "lineage_growth"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "output" / "plots" / "lineage_support"
RUN_KEYS = ["seed", "L", "Ly", "dphi", "p0"]
PRESETS = {
    "pilot-large": {
        "sizes": ["15", "20"],
        "p0s": ["2e-3", "1e-2"],
        "dphis": ["1.2e-1"],
        "seeds": ["1201", "1202"],
    },
    "full-large": {
        "sizes": ["15", "20"],
        "p0s": ["-1", "2e-3", "5e-3", "1e-2"],
        "dphis": ["1.2e-1"],
        "seeds": ["1201", "1202", "1203", "1204", "1205", "1206"],
    },
}
P0_LABELS = {
    -1.0: "No feedback",
    2e-3: r"$P_0 = 2 \times 10^{-3}$",
    5e-3: r"$P_0 = 5 \times 10^{-3}$",
    1e-2: r"$P_0 = 10^{-2}$",
}
P0_COLORS = {
    -1.0: "#4c566a",
    2e-3: "#0f766e",
    5e-3: "#c56a00",
    1e-2: "#b42318",
}
SIZE_LINESTYLES = {
    15: "-",
    20: "--",
}
SUFFIX_RE = re.compile(
    r"^v(?P<version>[^_]+)_ar(?P<ar>[^_]+)_div_(?P<divtype>[^_]+)"
    r"_desync(?P<desync>[^_]+)_seed_(?P<seed>[^_]+)"
    r"_Lx(?P<L>[^_]+)_Ly(?P<Ly>[^_]+)_att(?P<att>[^_]+)"
    r"_dphi(?P<dphi>[^_]+)_P(?P<p0>.+)\.dat$"
)


def parse_csv_list(value):
    items = [item.strip() for item in value.split(",") if item.strip()]
    if not items:
        raise argparse.ArgumentTypeError("expected a non-empty comma-separated list")
    return items


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize large-system lineage outputs into pilot-review tables and figures."
    )
    parser.add_argument("--preset", choices=sorted(PRESETS), default="pilot-large")
    parser.add_argument("--sizes", type=parse_csv_list, help="Comma-separated box sizes.")
    parser.add_argument("--p0s", type=parse_csv_list, help="Comma-separated feedback strengths.")
    parser.add_argument("--dphis", type=parse_csv_list, help="Comma-separated overcompressions.")
    parser.add_argument("--seeds", type=parse_csv_list, help="Comma-separated seeds.")
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--log-dir", type=Path, default=DEFAULT_LOG_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--formats", default="png,pdf")
    args = parser.parse_args()
    preset = PRESETS[args.preset]
    for field in ("sizes", "p0s", "dphis", "seeds"):
        if getattr(args, field) is None:
            setattr(args, field, list(preset[field]))
    return args


def plot_formats(value):
    formats = [item.strip().lower() for item in value.split(",") if item.strip()]
    if not formats:
        raise ValueError("at least one output format is required")
    return formats


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


def save_dataframe(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def save_figure(fig, output_dir, stem, formats):
    for ext in formats:
        fig.savefig(output_dir / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def record_rejection(rejections, dataset, source, reason):
    rejections.append(
        {
            "dataset": dataset,
            "source_file": source,
            "reason": str(reason),
        }
    )


def save_rejections(rejections, output_dir):
    save_dataframe(pd.DataFrame(rejections, columns=["dataset", "source_file", "reason"]), output_dir / "rejected_inputs.csv")


def metadata_from_path(path, prefix):
    match = SUFFIX_RE.match(path.name.removeprefix(prefix))
    if match is None:
        return None
    return match.groupdict()


def lineage_basename(seed, size, dphi, p0):
    return (
        f"v{LINEAGE_FIXED['version']}_ar{LINEAGE_FIXED['ar']}_div_{LINEAGE_FIXED['divtype']}"
        f"_desync{LINEAGE_FIXED['desync']}_seed_{seed}_Lx{size}_Ly{size}"
        f"_att{LINEAGE_FIXED['att']}_dphi{dphi}_P{p0}.dat"
    )


def lineage_paths(input_dir, log_dir, name):
    return {
        "lineage_jamm": input_dir / f"LINEAGE_LF_JAMM_{name}",
        "divlog": input_dir / f"DIVLOG_{name}",
        "transitions": input_dir / f"TRANSITIONS_{name}",
        "postjamm_summary": input_dir / f"POSTJAMM_SUMMARY_{name}",
        "stdout_log": log_dir / f"stdout_{name[:-4]}.log",
    }


def lineage_run_complete(paths):
    try:
        if not parse_log(paths["stdout_log"]):
            return False
        jamm_count = parse_lineage_snapshot(paths["lineage_jamm"])
        if jamm_count <= 0:
            return False
        if not parse_lineage(paths["lineage_jamm"], jamm_count):
            return False
        if not parse_divlog(paths["divlog"]):
            return False
        if not parse_transitions(paths["transitions"]):
            return False
        if not parse_postjamm_summary(paths["postjamm_summary"]):
            return False
    except (FileNotFoundError, OSError, ValueError):
        return False
    return True


def completed_run_names(args, input_dir, log_dir):
    completed = set()
    rejections = []
    for size in args.sizes:
        for p0 in args.p0s:
            for dphi in args.dphis:
                for seed in args.seeds:
                    name = lineage_basename(seed, size, dphi, p0)
                    if lineage_run_complete(lineage_paths(input_dir, log_dir, name)):
                        completed.add(name)
                    else:
                        record_rejection(rejections, "lineage_run", name, "run is missing required files or failed validation")
    return completed, rejections


def p0_order(values):
    ordered = [value for value in (-1.0, 2e-3, 5e-3, 1e-2) if value in values]
    leftovers = sorted(set(values) - set(ordered))
    return ordered + leftovers


def p0_label(value):
    return P0_LABELS.get(value, f"$P_0 = {value:.0e}$")


def p0_color(value):
    return P0_COLORS.get(value, "#1f2937")


def load_lineage_jamm(input_dir, valid_names):
    rows = []
    for path in sorted(input_dir.glob("LINEAGE_LF_JAMM_*.dat")):
        if path.name.removeprefix("LINEAGE_LF_JAMM_") not in valid_names:
            continue
        meta = metadata_from_path(path, "LINEAGE_LF_JAMM_")
        if meta is None:
            continue
        with path.open("r", encoding="utf-8") as handle:
            header = handle.readline()
            if not header.startswith("# index cell_id parent_id"):
                raise ValueError(f"unexpected lineage snapshot header in {path.name}")
            for line in handle:
                if not line.strip():
                    continue
                fields = line.split()
                if len(fields) != 8:
                    raise ValueError(f"unexpected lineage row width in {path.name}")
                rows.append(
                    {
                        **meta,
                        "index": int(fields[0]),
                        "cell_id": int(fields[1]),
                        "parent_id": int(fields[2]),
                        "alpha": float(fields[3]),
                        "bud_diameter": float(fields[4]),
                        "bud_contacts": int(fields[5]),
                        "bud_unconstrained": int(fields[6]),
                        "bud_compressed": int(fields[7]),
                        "source_file": path.name,
                    }
                )
    if not rows:
        raise SystemExit(f"No lineage jamming snapshots found under {input_dir}")
    frame = pd.DataFrame(rows)
    cast_metadata(frame)
    frame["tracked_initial_free"] = frame["bud_unconstrained"].astype(int)
    return frame


def load_transitions(input_dir, valid_names):
    rows = []
    for path in sorted(input_dir.glob("TRANSITIONS_*.dat")):
        if path.name.removeprefix("TRANSITIONS_") not in valid_names:
            continue
        meta = metadata_from_path(path, "TRANSITIONS_")
        if meta is None:
            continue
        with path.open("r", encoding="utf-8") as handle:
            header = handle.readline()
            if not header.startswith("# step phi cell_id parent_id"):
                raise ValueError(f"unexpected transitions header in {path.name}")
            for line in handle:
                if not line.strip():
                    continue
                fields = line.split()
                if len(fields) != 12:
                    raise ValueError(f"unexpected transitions row width in {path.name}")
                rows.append(
                    {
                        **meta,
                        "step": int(fields[0]),
                        "phi": float(fields[1]),
                        "cell_id": int(fields[2]),
                        "parent_id": int(fields[3]),
                        "initial_bud_diameter": float(fields[4]),
                        "bud_diameter": float(fields[5]),
                        "delta_bud_area": float(fields[6]),
                        "bud_contacts": int(fields[7]),
                        "bud_unconstrained": int(fields[8]),
                        "bud_compressed": int(fields[9]),
                        "event_code": int(fields[10]),
                        "event_label": fields[11],
                        "source_file": path.name,
                    }
                )
    if not rows:
        raise SystemExit(f"No lineage transition logs found under {input_dir}")
    frame = pd.DataFrame(rows)
    cast_metadata(frame)
    return frame


def load_postjamm_summary(input_dir, valid_names):
    rows = []
    for path in sorted(input_dir.glob("POSTJAMM_SUMMARY_*.dat")):
        if path.name.removeprefix("POSTJAMM_SUMMARY_") not in valid_names:
            continue
        meta = metadata_from_path(path, "POSTJAMM_SUMMARY_")
        if meta is None:
            continue
        with path.open("r", encoding="utf-8") as handle:
            header = handle.readline()
            if not header.startswith("# step phi P N Nc Nf Nu Ziso total_growthrate chi_c"):
                raise ValueError(f"unexpected post-jamming summary header in {path.name}")
            for line in handle:
                if not line.strip():
                    continue
                fields = line.split()
                if len(fields) != 15:
                    raise ValueError(f"unexpected post-jamming summary row width in {path.name}")
                rows.append(
                    {
                        **meta,
                        "step": int(fields[0]),
                        "phi": float(fields[1]),
                        "P": float(fields[2]),
                        "N": int(fields[3]),
                        "Nc": int(fields[4]),
                        "Nf": int(fields[5]),
                        "Nu": int(fields[6]),
                        "Ziso": int(fields[7]),
                        "total_growthrate": float(fields[8]),
                        "chi_c": float(fields[9]),
                        "n_compressed": int(fields[10]),
                        "n_initial_free_total": int(fields[11]),
                        "n_initial_free_active": int(fields[12]),
                        "n_initial_free_completed": int(fields[13]),
                        "n_initial_free_divided": int(fields[14]),
                        "source_file": path.name,
                    }
                )
    if not rows:
        raise SystemExit(f"No post-jamming lineage summaries found under {input_dir}")
    frame = pd.DataFrame(rows)
    cast_metadata(frame)
    return frame


def load_divlog_summary(input_dir, valid_names):
    rows = []
    for path in sorted(input_dir.glob("DIVLOG_*.dat")):
        if path.name.removeprefix("DIVLOG_") not in valid_names:
            continue
        meta = metadata_from_path(path, "DIVLOG_")
        if meta is None:
            continue
        event_rows = 0
        with path.open("r", encoding="utf-8") as handle:
            header = handle.readline()
            if not header.startswith("# step phi parent_cell_id"):
                raise ValueError(f"unexpected DIVLOG header in {path.name}")
            for line in handle:
                if not line.strip():
                    continue
                fields = line.split()
                if len(fields) != 6:
                    raise ValueError(f"unexpected DIVLOG row width in {path.name}")
                int(fields[0])
                float(fields[1])
                int(fields[2])
                int(fields[3])
                int(fields[4])
                int(fields[5])
                event_rows += 1
        rows.append({**meta, "divlog_events": event_rows, "source_file": path.name})
    if not rows:
        raise SystemExit(f"No lineage division logs found under {input_dir}")
    frame = pd.DataFrame(rows)
    cast_metadata(frame)
    return frame


def filter_frame(frame, args):
    mask = (
        frame["L"].isin([int(item) for item in args.sizes])
        & frame["seed"].isin([int(item) for item in args.seeds])
        & frame["dphi"].isin([float(item) for item in args.dphis])
        & frame["p0"].isin([float(item) for item in args.p0s])
    )
    return frame[mask].copy()


def first_event(transitions, event_code, prefix):
    subset = transitions[transitions["event_code"] == event_code].sort_values("step")
    subset = subset.drop_duplicates(RUN_KEYS + ["cell_id"], keep="first")
    return subset[RUN_KEYS + ["cell_id", "step", "phi", "delta_bud_area", "bud_compressed"]].rename(
        columns={
            "step": f"{prefix}_step",
            "phi": f"{prefix}_phi",
            "delta_bud_area": f"{prefix}_delta_bud_area",
            "bud_compressed": f"{prefix}_bud_compressed",
        }
    )


def build_threshold_events(lineage_jamm, transitions):
    tracked = lineage_jamm[lineage_jamm["tracked_initial_free"] == 1].copy()
    tracked = tracked.rename(columns={"bud_diameter": "jamm_bud_diameter"})
    tracked["tracked_at_jamming"] = 1
    first_contact = first_event(transitions, 1, "first_contact")
    completion = first_event(transitions, 2, "completion")
    division = first_event(transitions, 3, "division")

    thresholds = tracked[RUN_KEYS + ["cell_id", "parent_id", "jamm_bud_diameter", "alpha"]].copy()
    thresholds = thresholds.merge(first_contact, on=RUN_KEYS + ["cell_id"], how="left")
    thresholds = thresholds.merge(completion, on=RUN_KEYS + ["cell_id"], how="left")
    thresholds = thresholds.merge(division, on=RUN_KEYS + ["cell_id"], how="left")
    thresholds["has_first_contact"] = thresholds["first_contact_step"].notna()
    thresholds["has_completion"] = thresholds["completion_step"].notna()
    thresholds["has_division"] = thresholds["division_step"].notna()
    thresholds["completion_before_division"] = thresholds["has_completion"] & (
        ~thresholds["has_division"] | (thresholds["completion_step"] < thresholds["division_step"])
    )
    thresholds["divided_before_completion"] = thresholds["has_division"] & (
        ~thresholds["has_completion"] | (thresholds["division_step"] < thresholds["completion_step"])
    )
    thresholds["usable_completion"] = thresholds["completion_before_division"]
    thresholds["a1"] = thresholds["first_contact_delta_bud_area"]
    thresholds["a2"] = thresholds["completion_delta_bud_area"]
    thresholds["delta_a_stage2"] = thresholds["a2"] - thresholds["a1"]
    return thresholds


def build_depletion_curves(postjamm, lineage_jamm):
    tracked_counts = (
        lineage_jamm.groupby(RUN_KEYS, as_index=False)["tracked_initial_free"]
        .sum()
        .rename(columns={"tracked_initial_free": "tracked_total_jamm"})
    )
    curves = postjamm.merge(tracked_counts, on=RUN_KEYS, how="inner")
    curves = curves.sort_values(RUN_KEYS + ["step"]).reset_index(drop=True)
    first_rows = (
        curves.groupby(RUN_KEYS, as_index=False)
        .first()[RUN_KEYS + ["phi", "P", "N", "n_initial_free_total"]]
        .rename(
            columns={
                "phi": "phi_j",
                "P": "P_j",
                "N": "N_j",
                "n_initial_free_total": "n_initial_free_total_j",
            }
        )
    )
    curves = curves.merge(first_rows, on=RUN_KEYS, how="inner")
    if not np.array_equal(
        curves["tracked_total_jamm"].to_numpy(),
        curves["n_initial_free_total"].to_numpy(),
    ):
        mismatch = curves[curves["tracked_total_jamm"] != curves["n_initial_free_total"]]
        if not mismatch.empty:
            example = mismatch.iloc[0]
            raise ValueError(
                "tracked reservoir mismatch between LINEAGE_LF_JAMM and POSTJAMM_SUMMARY for "
                f"L={example['L']} seed={example['seed']} p0={example['p0']}"
            )
    total = curves["n_initial_free_total"].replace(0, np.nan)
    curves["active_fraction"] = curves["n_initial_free_active"] / total
    curves["completed_fraction"] = curves["n_initial_free_completed"] / total
    curves["divided_fraction"] = curves["n_initial_free_divided"] / total
    curves["delta_phi"] = curves["phi"] - curves["phi_j"]
    curves["u_j_tracked"] = curves["n_initial_free_total"] / curves["N_j"]
    curves["accounting_ok"] = (
        curves["n_initial_free_active"]
        + curves["n_initial_free_completed"]
        + curves["n_initial_free_divided"]
        == curves["n_initial_free_total"]
    )
    if not curves["accounting_ok"].all():
        bad = curves.loc[~curves["accounting_ok"]].iloc[0]
        raise ValueError(
            "tracked-reservoir accounting failed for "
            f"L={bad['L']} seed={bad['seed']} p0={bad['p0']} step={bad['step']}"
        )
    return curves


def build_secondary_summary(curves, thresholds, divlog_summary):
    rows = []
    threshold_summary = (
        thresholds.groupby(RUN_KEYS, as_index=False)
        .agg(
            n_tracked=("cell_id", "count"),
            n_first_contact=("has_first_contact", "sum"),
            n_completion=("usable_completion", "sum"),
            n_divided=("divided_before_completion", "sum"),
            a2_median=("a2", "median"),
            a2_q25=("a2", lambda s: s.quantile(0.25)),
            a2_q75=("a2", lambda s: s.quantile(0.75)),
        )
        .fillna(np.nan)
    )
    curves = curves.merge(threshold_summary, on=RUN_KEYS, how="left")
    curves = curves.merge(divlog_summary[RUN_KEYS + ["divlog_events"]], on=RUN_KEYS, how="left")
    for key, subset in curves.groupby(RUN_KEYS, sort=True):
        subset = subset.sort_values("step").reset_index(drop=True)
        exhausted = subset[subset["n_initial_free_active"] == 0]
        exhausted_by_endpoint = not exhausted.empty
        target_row = exhausted.iloc[0] if exhausted_by_endpoint else subset.iloc[-1]
        rows.append(
            {
                "seed": key[0],
                "L": key[1],
                "Ly": key[2],
                "dphi": key[3],
                "p0": key[4],
                "phi_j": float(subset["phi_j"].iloc[0]),
                "u_j_tracked": float(subset["u_j_tracked"].iloc[0]),
                "tracked_total": int(subset["n_initial_free_total"].iloc[0]),
                "phi_2_obs": float(target_row["phi"]) if exhausted_by_endpoint else np.nan,
                "delta_phi_2_obs": float(target_row["delta_phi"]) if exhausted_by_endpoint else np.nan,
                "P_2_obs": float(target_row["P"]) if exhausted_by_endpoint else np.nan,
                "bar_chi_until_target": float(
                    subset.loc[: target_row.name, "chi_c"].mean()
                ),
                "endpoint_active_fraction": float(subset["active_fraction"].iloc[-1]),
                "endpoint_completed_fraction": float(subset["completed_fraction"].iloc[-1]),
                "endpoint_divided_fraction": float(subset["divided_fraction"].iloc[-1]),
                "exhausted_by_endpoint": int(exhausted_by_endpoint),
                "divlog_events": int(subset["divlog_events"].iloc[0]),
                "n_completion": int(subset["n_completion"].iloc[0]),
                "n_divided": int(subset["n_divided"].iloc[0]),
                "a2_median": float(subset["a2_median"].iloc[0]) if pd.notna(subset["a2_median"].iloc[0]) else np.nan,
                "a2_q25": float(subset["a2_q25"].iloc[0]) if pd.notna(subset["a2_q25"].iloc[0]) else np.nan,
                "a2_q75": float(subset["a2_q75"].iloc[0]) if pd.notna(subset["a2_q75"].iloc[0]) else np.nan,
            }
        )
    summary = pd.DataFrame(rows)
    cast_metadata(summary)
    return summary.sort_values(RUN_KEYS).reset_index(drop=True)


def add_line_with_band(axis, x, y, y_low, y_high, color, linestyle="-", label=None):
    axis.plot(x, y, color=color, linewidth=2.2, linestyle=linestyle, label=label)
    axis.fill_between(x, y_low, y_high, color=color, alpha=0.16, linewidth=0)
    axis.grid(True, alpha=0.25)


def summarize_interpolated_curves(curves, value_column, n_points=150):
    rows = []
    for (L, p0), subset in curves.groupby(["L", "p0"], sort=True):
        max_delta = float(subset["delta_phi"].max())
        grid = np.linspace(0.0, max_delta, n_points)
        series = []
        for _, run in subset.groupby(RUN_KEYS, sort=True):
            run = run.sort_values("delta_phi")
            series.append(
                np.interp(
                    grid,
                    run["delta_phi"].to_numpy(),
                    run[value_column].to_numpy(),
                    left=run[value_column].iloc[0],
                    right=np.nan,
                )
            )
        values = np.asarray(series, dtype=float)
        rows.append(
            pd.DataFrame(
                {
                    "L": L,
                    "p0": p0,
                    "delta_phi": grid,
                    f"{value_column}_median": np.nanmedian(values, axis=0),
                    f"{value_column}_q25": np.nanquantile(values, 0.25, axis=0),
                    f"{value_column}_q75": np.nanquantile(values, 0.75, axis=0),
                }
            )
        )
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def plot_thresholds(thresholds, output_dir, formats):
    usable = thresholds[thresholds["usable_completion"] & thresholds["a2"].notna()].copy()
    sizes = sorted(usable["L"].unique())
    fig, axes = plt.subplots(1, len(sizes), figsize=(6.4 * len(sizes), 4.8), constrained_layout=True)
    axes = np.atleast_1d(axes)
    for axis, L in zip(axes, sizes):
        size_subset = usable[usable["L"] == L]
        for p0 in p0_order(size_subset["p0"].unique()):
            values = np.sort(size_subset[size_subset["p0"] == p0]["a2"].to_numpy())
            if len(values) == 0:
                continue
            survival = 1.0 - np.arange(1, len(values) + 1) / len(values)
            axis.step(values, survival, where="post", color=p0_color(p0), linewidth=2.2, label=p0_label(p0))
        axis.set_xlabel(r"Completion threshold $a^\ast$")
        axis.set_ylabel(r"Survival $S(a^\ast)$")
        axis.set_title(f"L = {L}")
        axis.grid(True, alpha=0.25)
    if len(axes):
        axes[0].legend(frameon=False, fontsize=9, loc="upper right")
    fig.suptitle("Completion-Threshold Survival by System Size", fontsize=14)
    save_figure(fig, output_dir, "fig_lineage_thresholds", formats)


def plot_two_stage(thresholds, output_dir, formats):
    paired = thresholds[thresholds["has_first_contact"] & thresholds["has_completion"]].copy()
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.0), constrained_layout=True)

    axis = axes[0]
    if not paired.empty:
        all_values = np.concatenate([paired["a1"].to_numpy(), paired["a2"].to_numpy()])
        xline = np.linspace(0.0, float(np.nanmax(all_values)) * 1.02, 100)
        axis.plot(xline, xline, color="#111827", linestyle="--", linewidth=1.4)
        for L in sorted(paired["L"].unique()):
            subset = paired[paired["L"] == L]
            axis.scatter(
                subset["a1"],
                subset["a2"],
                s=32,
                c=subset["p0"].map(p0_color),
                alpha=0.75,
                linewidths=0.0,
                marker="o" if L == 15 else "^",
                label=f"L = {L}",
            )
    axis.set_xlabel(r"First-contact threshold $a_1^\ast$")
    axis.set_ylabel(r"Completion threshold $a_2^\ast$")
    axis.set_title("Two-Stage Threshold Check")
    axis.grid(True, alpha=0.25)
    size_handles = [
        Line2D([0], [0], marker="o" if L == 15 else "^", color="#111827", linestyle="", label=f"L = {L}")
        for L in sorted(paired["L"].unique())
    ]
    color_handles = [
        Line2D([0], [0], color=p0_color(p0), linewidth=2.2, label=p0_label(p0))
        for p0 in p0_order(paired["p0"].unique())
    ]
    if size_handles:
        legend1 = axis.legend(handles=size_handles, frameon=False, fontsize=8, loc="upper left")
        axis.add_artist(legend1)
    if color_handles:
        axis.legend(handles=color_handles, frameon=False, fontsize=8, loc="lower right")

    axis = axes[1]
    gaps = paired[paired["delta_a_stage2"].notna()].copy()
    for L in sorted(gaps["L"].unique()):
        subset_L = gaps[gaps["L"] == L]
        for p0 in p0_order(subset_L["p0"].unique()):
            values = np.sort(subset_L[subset_L["p0"] == p0]["delta_a_stage2"].to_numpy())
            if len(values) == 0:
                continue
            ecdf = np.arange(1, len(values) + 1) / len(values)
            axis.step(
                values,
                ecdf,
                where="post",
                color=p0_color(p0),
                linestyle=SIZE_LINESTYLES.get(L, "-"),
                linewidth=2.0,
            )
    axis.set_xlabel(r"Gap $a_2^\ast - a_1^\ast$")
    axis.set_ylabel("ECDF")
    axis.set_title("Completion Gap Distribution")
    axis.grid(True, alpha=0.25)

    fig.suptitle("First Contact Versus Mechanical Completion", fontsize=14)
    save_figure(fig, output_dir, "fig_lineage_two_stage", formats)


def plot_depletion(curves, output_dir, formats):
    summary = summarize_interpolated_curves(curves, "active_fraction")
    sizes = sorted(curves["L"].unique())
    fig, axes = plt.subplots(1, len(sizes), figsize=(6.5 * len(sizes), 4.8), constrained_layout=True)
    axes = np.atleast_1d(axes)
    for axis, L in zip(axes, sizes):
        size_curves = curves[curves["L"] == L]
        size_summary = summary[summary["L"] == L]
        for p0 in p0_order(size_curves["p0"].unique()):
            subset = size_curves[size_curves["p0"] == p0]
            for _, run in subset.groupby(RUN_KEYS, sort=True):
                axis.plot(
                    run["delta_phi"],
                    run["active_fraction"],
                    color=p0_color(p0),
                    linewidth=1.1,
                    alpha=0.30,
                )
                exhausted = run[run["n_initial_free_active"] == 0]
                if not exhausted.empty:
                    first = exhausted.iloc[0]
                    axis.scatter(first["delta_phi"], first["active_fraction"], color=p0_color(p0), s=26, zorder=3)
            med = size_summary[size_summary["p0"] == p0].sort_values("delta_phi")
            if med.empty:
                continue
            add_line_with_band(
                axis,
                med["delta_phi"].to_numpy(),
                med["active_fraction_median"].to_numpy(),
                med["active_fraction_q25"].to_numpy(),
                med["active_fraction_q75"].to_numpy(),
                p0_color(p0),
                label=p0_label(p0),
            )
        axis.set_xlabel(r"$\phi - \phi_J$")
        axis.set_ylabel(r"Tracked active fraction $u_{\mathrm{track}} / u_J$")
        axis.set_ylim(-0.03, 1.03)
        axis.set_title(f"L = {L}")
    if len(axes):
        axes[0].legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle("Depletion of the Initial Free-Bud Reservoir", fontsize=14)
    save_figure(fig, output_dir, "fig_lineage_depletion", formats)


def plot_injection_check(summary, output_dir, formats):
    summary = summary.sort_values(["L", "p0"]).reset_index(drop=True)
    labels = [f"L={int(row.L)}\n{p0_label(row.p0)}" for row in summary.itertuples()]
    x = np.arange(len(summary))

    fig, ax = plt.subplots(figsize=(max(9.5, 1.6 * len(summary)), 5.0), constrained_layout=True)
    ax.bar(x, summary["endpoint_completed_fraction"], color="#0f766e", label="Completed")
    ax.bar(
        x,
        summary["endpoint_divided_fraction"],
        bottom=summary["endpoint_completed_fraction"],
        color="#c56a00",
        label="Divided",
    )
    ax.bar(
        x,
        summary["endpoint_active_fraction"],
        bottom=summary["endpoint_completed_fraction"] + summary["endpoint_divided_fraction"],
        color="#94a3b8",
        label="Still active",
    )
    ax.set_xticks(x, labels)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("Fraction of tracked initial-free buds")
    ax.set_title("Endpoint Partition of the Initial Free-Bud Reservoir")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    save_figure(fig, output_dir, "fig_lineage_injection_check", formats)


def plot_secondary_summary(summary, output_dir, formats):
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.9), constrained_layout=True)
    exhausted = summary[summary["exhausted_by_endpoint"] == 1].copy()
    not_exhausted = summary[summary["exhausted_by_endpoint"] == 0].copy()
    for axis, value, ylabel, title in zip(
        axes,
        ["delta_phi_2_obs", "P_2_obs"],
        [r"$\phi_2 - \phi_J$", r"$P_2$"],
        ["Observed Secondary-Arrest Location", "Observed Pressure at Secondary Arrest"],
    ):
        for L in sorted(summary["L"].unique()):
            subset = exhausted[exhausted["L"] == L].sort_values("p0")
            if subset.empty:
                continue
            x = np.arange(len(subset))
            axis.plot(
                x,
                subset[value],
                color="#111827",
                linestyle=SIZE_LINESTYLES.get(L, "-"),
                linewidth=2.0,
                marker="o",
                markersize=5,
                label=f"L = {L}",
            )
            axis.set_xticks(x, [p0_label(item) for item in subset["p0"]])
        if not not_exhausted.empty:
            axis.text(
                0.03,
                0.06,
                f"Not exhausted by endpoint: {len(not_exhausted)} run(s)",
                transform=axis.transAxes,
                ha="left",
                va="bottom",
                fontsize=9,
            )
        axis.set_ylabel(ylabel)
        axis.set_title(title)
        axis.grid(True, alpha=0.25)
    if len(summary["L"].unique()):
        axes[0].legend(frameon=False, fontsize=9, loc="upper right")
    fig.suptitle("Pilot Secondary-Arrest Summary", fontsize=14)
    save_figure(fig, output_dir, "fig_lineage_secondary_summary", formats)


def write_summary_note(summary, output_dir, preset_name):
    exhausted = int(summary["exhausted_by_endpoint"].sum())
    lines = [
        "# Lineage pilot summary",
        "",
        f"- Preset: `{preset_name}`",
        f"- Runs summarized: {len(summary)}",
        f"- Exhausted by endpoint: {exhausted}/{len(summary)}",
        f"- Sizes: {', '.join(str(item) for item in sorted(summary['L'].unique()))}",
        f"- Feedbacks: {', '.join(f'{item:g}' for item in sorted(summary['p0'].unique()))}",
    ]
    (output_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = parse_args()
    formats = plot_formats(args.formats)
    input_dir = args.input_dir.resolve()
    log_dir = args.log_dir.resolve()
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

    valid_names, rejections = completed_run_names(args, input_dir, log_dir)
    if not valid_names:
        raise SystemExit("No completed lineage runs passed validation for the requested parameter set")

    lineage_jamm = filter_frame(load_lineage_jamm(input_dir, valid_names), args)
    transitions = filter_frame(load_transitions(input_dir, valid_names), args)
    postjamm = filter_frame(load_postjamm_summary(input_dir, valid_names), args)
    divlog_summary = filter_frame(load_divlog_summary(input_dir, valid_names), args)

    thresholds = build_threshold_events(lineage_jamm, transitions)
    curves = build_depletion_curves(postjamm, lineage_jamm)
    summary = build_secondary_summary(curves, thresholds, divlog_summary)

    save_dataframe(thresholds, output_dir / "threshold_events.csv")
    save_dataframe(curves, output_dir / "depletion_curves.csv")
    save_dataframe(summary, output_dir / "secondary_summary.csv")
    save_rejections(rejections, output_dir)

    plot_thresholds(thresholds, output_dir, formats)
    plot_two_stage(thresholds, output_dir, formats)
    plot_depletion(curves, output_dir, formats)
    plot_injection_check(summary, output_dir, formats)
    plot_secondary_summary(summary, output_dir, formats)
    write_summary_note(summary, output_dir, args.preset)

    print(
        "Generated lineage-support plots "
        f"for {len(summary)} run(s) with "
        f"{int(summary['exhausted_by_endpoint'].sum())} exhausted by the endpoint.",
        flush=True,
    )


if __name__ == "__main__":
    main()
