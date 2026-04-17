#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mplconfig_"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.mkdtemp(prefix="xdgcache_"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from pipeline_config import REPO_ROOT, basename, iter_job_params, seed_range
from pipeline_paths import growth_paths
from pipeline_validate import growth_done


P0_ORDER = ["-1", "1e-1", "1e-2", "1e-3", "1e-4", "1e-5"]
P0_LABELS = {
    "-1": "no_fb",
    "1e-1": r"$10^{-1}$",
    "1e-2": r"$10^{-2}$",
    "1e-3": r"$10^{-3}$",
    "1e-4": r"$10^{-4}$",
    "1e-5": r"$10^{-5}$",
}
P0_COLORS = {
    "-1": "#4c566a",
    "1e-1": "#2f5d8a",
    "1e-2": "#7b4db2",
    "1e-3": "#1f9e89",
    "1e-4": "#84b547",
    "1e-5": "#f28e2b",
}
SIZE_STYLES = {
    15: {"color": "#2f5d8a", "marker": "o"},
    20: {"color": "#c77259", "marker": "s"},
    25: {"color": "#1f9e89", "marker": "^"},
}
POSTJAM_TOL = 1e-12
BOOTSTRAP_SEED = 12345


@dataclass(frozen=True)
class PostSummaryRow:
    step: int
    phi: float
    n: int
    initial_free_active: int
    postjam_divisions_total: int


@dataclass(frozen=True)
class EventRow:
    step: int
    phi: float


@dataclass
class CandidateRecord:
    lx: int
    p0: str
    seed: str
    dphi_target: float
    post_rows: list[PostSummaryRow]
    div_events: list[EventRow]
    track_events: list[EventRow]

    @property
    def phi_j(self) -> float:
        return self.post_rows[0].phi

    @property
    def step_j(self) -> int:
        return self.post_rows[0].step

    @property
    def phi_final(self) -> float:
        return self.post_rows[-1].phi

    @property
    def delta_phi_final(self) -> float:
        return self.phi_final - self.phi_j


@dataclass
class RunSummary:
    lx: int
    p0: str
    seed: str
    dphi_target: float
    n_j: int
    phi_j: float
    delta_phi_final: float
    phi_2_obs: float | None
    delta_phi_2_obs: float | None
    exhausted_by_endpoint: int
    n_div_all_at_phi2: int | None
    n_div_all_final: int
    n_div_all_after_phi2: int | None
    frac_before_phi2: float | None
    n_div_track_at_phi2: int | None
    n_div_track_final: int
    divlog_postjam_count: int
    postjamm_summary_final_count: int
    postjamm_summary_count_gap: int
    all_event_deltas: np.ndarray
    track_event_deltas: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate supplementary post-jamming division-count diagnostics.")
    parser.add_argument("--panel-a-size", type=int, default=20)
    parser.add_argument("--panel-b-sizes", nargs="+", type=int, default=[20])
    parser.add_argument("--output-dir", default=str(REPO_ROOT / "paper" / "revision-figures"))
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--bin-width", type=float, default=0.0025)
    parser.add_argument("--grid-step", type=float, default=0.0025)
    parser.add_argument("--bootstrap-reps", type=int, default=2000)
    parser.add_argument("--seed-start", type=int, help="Optional inclusive seed start.")
    parser.add_argument("--seed-stop", type=int, help="Optional inclusive seed stop.")
    args = parser.parse_args()
    if (args.seed_start is None) ^ (args.seed_stop is None):
        raise SystemExit("--seed-start and --seed-stop must be given together")
    if args.seed_start is not None and args.seed_stop < args.seed_start:
        raise SystemExit("--seed-stop must be greater than or equal to --seed-start")
    if args.bin_width <= 0.0:
        raise SystemExit("--bin-width must be positive")
    if args.grid_step <= 0.0:
        raise SystemExit("--grid-step must be positive")
    if args.bootstrap_reps <= 0:
        raise SystemExit("--bootstrap-reps must be positive")
    return args


def configure_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 4,
            "ytick.major.size": 4,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def parse_post_rows(path: Path) -> list[PostSummaryRow]:
    rows: list[PostSummaryRow] = []
    with path.open("r", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            if not line.strip():
                continue
            ff = line.split()
            rows.append(
                PostSummaryRow(
                    step=int(ff[0]),
                    phi=float(ff[1]),
                    n=int(ff[3]),
                    initial_free_active=int(ff[12]),
                    postjam_divisions_total=int(ff[18]),
                )
            )
    return rows


def parse_div_events(path: Path) -> list[EventRow]:
    events: list[EventRow] = []
    with path.open("r", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            if not line.strip():
                continue
            ff = line.split()
            events.append(EventRow(step=int(ff[0]), phi=float(ff[1])))
    return events


def parse_track_events(path: Path) -> list[EventRow]:
    events: list[EventRow] = []
    with path.open("r", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            if not line.strip():
                continue
            ff = line.split()
            if int(ff[10]) != 3:
                continue
            events.append(EventRow(step=int(ff[0]), phi=float(ff[1])))
    return events


def select_postjam_event_deltas(events: list[EventRow], step_j: int, phi_j: float) -> np.ndarray:
    deltas = [
        max(0.0, event.phi - phi_j)
        for event in events
        if event.step >= step_j and event.phi >= phi_j - POSTJAM_TOL
    ]
    if not deltas:
        return np.asarray([], dtype=float)
    deltas.sort()
    return np.asarray(deltas, dtype=float)


def load_candidate_records(sizes: list[str], seeds: list[str] | None) -> list[CandidateRecord]:
    records: list[CandidateRecord] = []
    for params in iter_job_params(sizes=sizes, seeds=seeds):
        name = basename(params)
        growth = growth_paths(name)
        if not growth_done(growth):
            continue
        post_rows = parse_post_rows(growth["postjamm_summary"])
        if not post_rows:
            continue
        records.append(
            CandidateRecord(
                lx=int(params["lx"]),
                p0=params["p0"],
                seed=params["seed"],
                dphi_target=float(params["dphi"]),
                post_rows=post_rows,
                div_events=parse_div_events(growth["divlog"]),
                track_events=parse_track_events(growth["transitions"]),
            )
        )
    return records


def choose_representative_records(records: list[CandidateRecord]) -> list[CandidateRecord]:
    selected: dict[tuple[int, str, str], CandidateRecord] = {}
    for record in records:
        key = (record.lx, record.p0, record.seed)
        current = selected.get(key)
        if current is None:
            selected[key] = record
            continue
        if record.delta_phi_final > current.delta_phi_final + POSTJAM_TOL:
            selected[key] = record
            continue
        if math.isclose(record.delta_phi_final, current.delta_phi_final, abs_tol=POSTJAM_TOL) and record.dphi_target > current.dphi_target:
            selected[key] = record
    return sorted(selected.values(), key=lambda row: (row.lx, float(row.p0), int(row.seed)))


def count_leq(values: np.ndarray, cutoff: float) -> int:
    return int(np.searchsorted(values, cutoff + POSTJAM_TOL, side="right"))


def summarize_record(record: CandidateRecord) -> RunSummary:
    phi_j = record.phi_j
    step_j = record.step_j
    phi_2_obs = next((row.phi for row in record.post_rows if row.initial_free_active == 0), None)
    delta_phi_2_obs = None if phi_2_obs is None else phi_2_obs - phi_j
    all_event_deltas = select_postjam_event_deltas(record.div_events, step_j, phi_j)
    track_event_deltas = select_postjam_event_deltas(record.track_events, step_j, phi_j)
    n_div_all_final = int(all_event_deltas.size)
    n_div_track_final = int(track_event_deltas.size)
    summary_final = record.post_rows[-1].postjam_divisions_total

    n_div_all_at_phi2 = None
    n_div_all_after_phi2 = None
    frac_before_phi2 = None
    n_div_track_at_phi2 = None
    exhausted_by_endpoint = 0
    if delta_phi_2_obs is not None:
        exhausted_by_endpoint = 1
        n_div_all_at_phi2 = count_leq(all_event_deltas, delta_phi_2_obs)
        n_div_all_after_phi2 = n_div_all_final - n_div_all_at_phi2
        n_div_track_at_phi2 = count_leq(track_event_deltas, delta_phi_2_obs)
        if n_div_all_final > 0:
            frac_before_phi2 = n_div_all_at_phi2 / n_div_all_final

    return RunSummary(
        lx=record.lx,
        p0=record.p0,
        seed=record.seed,
        dphi_target=record.dphi_target,
        n_j=record.post_rows[0].n,
        phi_j=phi_j,
        delta_phi_final=record.delta_phi_final,
        phi_2_obs=phi_2_obs,
        delta_phi_2_obs=delta_phi_2_obs,
        exhausted_by_endpoint=exhausted_by_endpoint,
        n_div_all_at_phi2=n_div_all_at_phi2,
        n_div_all_final=n_div_all_final,
        n_div_all_after_phi2=n_div_all_after_phi2,
        frac_before_phi2=frac_before_phi2,
        n_div_track_at_phi2=n_div_track_at_phi2,
        n_div_track_final=n_div_track_final,
        divlog_postjam_count=n_div_all_final,
        postjamm_summary_final_count=summary_final,
        postjamm_summary_count_gap=n_div_all_final - summary_final,
        all_event_deltas=all_event_deltas,
        track_event_deltas=track_event_deltas,
    )


def quantile_triplet(values: list[float]) -> tuple[float, float, float]:
    array = np.asarray(values, dtype=float)
    return (
        float(np.median(array)),
        float(np.percentile(array, 25.0)),
        float(np.percentile(array, 75.0)),
    )


def bootstrap_mean_ci(values: list[float], reps: int, seed: int = BOOTSTRAP_SEED) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    if len(array) <= 1:
        value = float(array[0]) if len(array) == 1 else float("nan")
        return value, value
    rng = np.random.default_rng(seed)
    means = np.empty(reps, dtype=float)
    for index in range(reps):
        sample = rng.choice(array, size=len(array), replace=True)
        means[index] = float(np.mean(sample))
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def bootstrap_median_ci(values: list[float], reps: int, seed: int = BOOTSTRAP_SEED) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    if len(array) <= 1:
        value = float(array[0]) if len(array) == 1 else float("nan")
        return value, value
    rng = np.random.default_rng(seed)
    medians = np.empty(reps, dtype=float)
    for index in range(reps):
        sample = rng.choice(array, size=len(array), replace=True)
        medians[index] = float(np.median(sample))
    return float(np.percentile(medians, 2.5)), float(np.percentile(medians, 97.5))


def build_cumulative_rows(group: list[RunSummary], grid_step: float) -> list[dict]:
    max_endpoint = max(run.delta_phi_final for run in group)
    grid = np.arange(0.0, max_endpoint + 0.5 * grid_step, grid_step, dtype=float)
    delta_phi_2_values = [run.delta_phi_2_obs for run in group if run.delta_phi_2_obs is not None]
    median_delta_phi_2 = float(np.median(delta_phi_2_values)) if delta_phi_2_values else None
    rows: list[dict] = []
    for delta_phi in grid:
        covered = [run for run in group if run.delta_phi_final >= delta_phi - POSTJAM_TOL]
        if not covered:
            continue
        all_counts = [count_leq(run.all_event_deltas, float(delta_phi)) for run in covered]
        track_counts = [count_leq(run.track_event_deltas, float(delta_phi)) for run in covered]
        all_med, all_q1, all_q3 = quantile_triplet(all_counts)
        track_med, track_q1, track_q3 = quantile_triplet(track_counts)
        rows.append(
            {
                "lx": covered[0].lx,
                "p0": covered[0].p0,
                "delta_phi": float(delta_phi),
                "n_runs": len(covered),
                "median_delta_phi_2_obs": median_delta_phi_2,
                "all_median": all_med,
                "all_q1": all_q1,
                "all_q3": all_q3,
                "track_mean": float(np.mean(track_counts)),
                "track_median": track_med,
                "track_q1": track_q1,
                "track_q3": track_q3,
            }
        )
    return rows


def build_rate_rows(group: list[RunSummary], bin_width: float, bootstrap_reps: int) -> list[dict]:
    max_endpoint = max(run.delta_phi_final for run in group)
    edges = np.arange(0.0, max_endpoint + bin_width, bin_width, dtype=float)
    if edges[-1] < max_endpoint + POSTJAM_TOL:
        edges = np.append(edges, edges[-1] + bin_width)
    delta_phi_2_values = [run.delta_phi_2_obs for run in group if run.delta_phi_2_obs is not None]
    median_delta_phi_2 = float(np.median(delta_phi_2_values)) if delta_phi_2_values else None
    rows: list[dict] = []
    for left, right in zip(edges[:-1], edges[1:]):
        covered = [run for run in group if run.delta_phi_final >= right - POSTJAM_TOL]
        if not covered:
            continue
        counts = []
        for run in covered:
            lo = np.searchsorted(run.all_event_deltas, left - POSTJAM_TOL, side="left")
            hi = np.searchsorted(run.all_event_deltas, right - POSTJAM_TOL, side="left")
            counts.append((hi - lo) / bin_width)
        mean_rate = float(np.mean(counts))
        ci_low, ci_high = bootstrap_mean_ci(counts, bootstrap_reps)
        rows.append(
            {
                "lx": covered[0].lx,
                "p0": covered[0].p0,
                "bin_left": float(left),
                "bin_right": float(right),
                "n_runs": len(covered),
                "median_delta_phi_2_obs": median_delta_phi_2,
                "rate_mean": mean_rate,
                "rate_ci_low": ci_low,
                "rate_ci_high": ci_high,
            }
        )
    return rows


def build_phi2_rows(group: list[RunSummary], bootstrap_reps: int) -> list[dict]:
    exhausted = [run for run in group if run.exhausted_by_endpoint == 1 and run.n_div_all_at_phi2 is not None]
    if not exhausted:
        return []
    values = [float(run.n_div_all_at_phi2) for run in exhausted if run.n_div_all_at_phi2 is not None]
    median, q1, q3 = quantile_triplet(values)
    ci_low, ci_high = bootstrap_median_ci(values, bootstrap_reps)
    return [
        {
            "lx": exhausted[0].lx,
            "p0": exhausted[0].p0,
            "n_runs": len(exhausted),
            "median": median,
            "q1": q1,
            "q3": q3,
            "ci_low": ci_low,
            "ci_high": ci_high,
        }
    ]


def group_by_size_and_p0(summaries: list[RunSummary]) -> dict[tuple[int, str], list[RunSummary]]:
    grouped: dict[tuple[int, str], list[RunSummary]] = {}
    for summary in summaries:
        grouped.setdefault((summary.lx, summary.p0), []).append(summary)
    for value in grouped.values():
        value.sort(key=lambda row: int(row.seed))
    return grouped


def format_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, int):
        return str(value)
    if isinstance(value, str):
        return value
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.12g}"
    return str(value)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: format_value(row.get(name)) for name in fieldnames})


def add_panel_label(ax: plt.Axes, text: str) -> None:
    ax.text(-0.22, 1.05, text, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def save_figure(fig: plt.Figure, output_dir: Path, stem: str, dpi: int) -> None:
    fig.savefig(output_dir / f"{stem}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_grouped_lookup(rows: list[dict], key_name: str) -> dict[tuple[int, str], list[dict]]:
    grouped: dict[tuple[int, str], list[dict]] = {}
    for row in rows:
        grouped.setdefault((int(row["lx"]), row["p0"]), []).append(row)
    for values in grouped.values():
        values.sort(key=lambda row: row[key_name])
    return grouped


def average_nj_by_size(summaries: list[RunSummary]) -> dict[int, float]:
    grouped: dict[int, list[int]] = {}
    for summary in summaries:
        grouped.setdefault(summary.lx, []).append(summary.n_j)
    return {size: float(np.mean(values)) for size, values in grouped.items()}


def make_figure(
    output_dir: Path,
    panel_a_size: int,
    panel_b_sizes: list[int],
    cumulative_rows: list[dict],
    phi2_rows: list[dict],
    avg_nj_by_size: dict[int, float],
    dpi: int,
) -> None:
    cumulative_lookup = build_grouped_lookup(cumulative_rows, "delta_phi")
    phi2_lookup = build_grouped_lookup(phi2_rows, "median")

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(6.3, 7.2), constrained_layout=False)
    fig.subplots_adjust(hspace=0.42, top=0.9)

    present_p0s_a = [p0 for p0 in P0_ORDER if (panel_a_size, p0) in cumulative_lookup]
    x_max = 0.0
    for p0 in present_p0s_a:
        curves = cumulative_lookup[(panel_a_size, p0)]
        x = np.asarray([row["delta_phi"] for row in curves], dtype=float)
        y = np.asarray([row["all_median"] for row in curves], dtype=float)
        y_q1 = np.asarray([row["all_q1"] for row in curves], dtype=float)
        y_q3 = np.asarray([row["all_q3"] for row in curves], dtype=float)
        ax_top.plot(x, y, lw=1.8, color=P0_COLORS[p0])
        ax_top.fill_between(x, y_q1, y_q3, color=P0_COLORS[p0], alpha=0.18, linewidth=0.0)
        x_max = max(x_max, float(x[-1]))
        delta_phi_2 = curves[0]["median_delta_phi_2_obs"]
        if delta_phi_2 is not None:
            ax_top.axvline(delta_phi_2, color=P0_COLORS[p0], lw=1.0, ls="--", alpha=0.8)

    ax_top.set_xlim(0.0, max(x_max, 1e-6))
    ax_top.set_xlabel(r"$\Delta \phi = \phi - \phi_J$")
    ax_top.set_ylabel("cumulative divisions per run")
    add_panel_label(ax_top, "A")

    x_positions = np.arange(len(P0_ORDER), dtype=float)
    present_l_sizes = [size for size in panel_b_sizes if any((size, p0) in phi2_lookup for p0 in P0_ORDER)]
    if present_l_sizes:
        size = present_l_sizes[0]
        norm = avg_nj_by_size[size]
        xs = []
        median = []
        ci_low = []
        ci_high = []
        point_p0s = []
        for idx, p0 in enumerate(P0_ORDER):
            if (size, p0) not in phi2_lookup:
                continue
            row = phi2_lookup[(size, p0)][0]
            xs.append(x_positions[idx])
            median.append(row["median"] / norm)
            ci_low.append(row["ci_low"] / norm)
            ci_high.append(row["ci_high"] / norm)
            point_p0s.append(p0)
        if xs:
            xs_arr = np.asarray(xs, dtype=float)
            median_arr = np.asarray(median, dtype=float)
            ci_low_arr = np.asarray(ci_low, dtype=float)
            ci_high_arr = np.asarray(ci_high, dtype=float)
            ax_bot.plot(xs_arr, median_arr, color="#8b8174", lw=1.0, alpha=0.75, zorder=1)
            for idx, p0 in enumerate(point_p0s):
                ax_bot.errorbar(
                    [xs_arr[idx]],
                    [median_arr[idx]],
                    yerr=np.asarray([[max(0.0, median_arr[idx] - ci_low_arr[idx])], [max(0.0, ci_high_arr[idx] - median_arr[idx])]], dtype=float),
                    fmt="o",
                    color=P0_COLORS[p0],
                    mfc=P0_COLORS[p0],
                    mec="#222222",
                    mew=0.7,
                    ms=6.2,
                    capsize=2.5,
                    lw=1.0,
                    zorder=3,
                )

    ax_bot.set_xticks(x_positions, [P0_LABELS[p0] for p0 in P0_ORDER], rotation=20)
    ax_bot.set_xlabel(r"$P_0$")
    ax_bot.set_ylabel(r"$N_{\mathrm{div}}(\phi_2) / \langle N_J \rangle$")
    ax_bot.set_title(r"cumulative cell divisions between $\phi_J$ and $\phi_2$", pad=8)
    add_panel_label(ax_bot, "B")

    legend_handles = [
        Line2D([0], [0], color=P0_COLORS[p0], lw=2.0, label=P0_LABELS[p0])
        for p0 in present_p0s_a
    ]
    if legend_handles:
        fig.legend(handles=legend_handles, loc="upper center", ncols=min(len(legend_handles), 6), fontsize=9, bbox_to_anchor=(0.5, 0.985))
    save_figure(fig, output_dir, "fig_lineage_division_counts", dpi)


def main() -> None:
    args = parse_args()
    configure_style()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    seeds = None
    if args.seed_start is not None:
        seeds = seed_range(args.seed_start, args.seed_stop)

    panel_b_sizes = [int(size) for size in args.panel_b_sizes]
    load_sizes = sorted({str(args.panel_a_size), *[str(size) for size in panel_b_sizes]}, key=int)
    candidates = load_candidate_records(load_sizes, seeds)
    representatives = choose_representative_records(candidates)
    summaries = [summarize_record(record) for record in representatives]
    grouped = group_by_size_and_p0(summaries)

    cumulative_rows: list[dict] = []
    rate_rows: list[dict] = []
    phi2_rows: list[dict] = []
    for key in sorted(grouped, key=lambda item: (item[0], float(item[1]))):
        group = grouped[key]
        cumulative_rows.extend(build_cumulative_rows(group, args.grid_step))
        rate_rows.extend(build_rate_rows(group, args.bin_width, args.bootstrap_reps))
        phi2_rows.extend(build_phi2_rows(group, args.bootstrap_reps))

    summary_rows = [
        {
            "lx": summary.lx,
            "p0": summary.p0,
            "seed": summary.seed,
            "dphi_target": summary.dphi_target,
            "n_j": summary.n_j,
            "phi_j": summary.phi_j,
            "phi_2_obs": summary.phi_2_obs,
            "delta_phi_2_obs": summary.delta_phi_2_obs,
            "delta_phi_final": summary.delta_phi_final,
            "exhausted_by_endpoint": summary.exhausted_by_endpoint,
            "N_div_all_at_phi2": summary.n_div_all_at_phi2,
            "N_div_all_final": summary.n_div_all_final,
            "N_div_all_after_phi2": summary.n_div_all_after_phi2,
            "frac_before_phi2": summary.frac_before_phi2,
            "N_div_track_at_phi2": summary.n_div_track_at_phi2,
            "N_div_track_final": summary.n_div_track_final,
            "divlog_postjam_count": summary.divlog_postjam_count,
            "postjamm_summary_final_count": summary.postjamm_summary_final_count,
            "postjamm_summary_count_gap": summary.postjamm_summary_count_gap,
        }
        for summary in summaries
    ]

    write_csv(
        output_dir / "division_counts_cumulative.csv",
        [
            "lx",
            "p0",
            "delta_phi",
            "n_runs",
            "median_delta_phi_2_obs",
            "all_median",
            "all_q1",
            "all_q3",
            "track_mean",
            "track_median",
            "track_q1",
            "track_q3",
        ],
        cumulative_rows,
    )
    write_csv(
        output_dir / "division_counts_rate.csv",
        [
            "lx",
            "p0",
            "bin_left",
            "bin_right",
            "n_runs",
            "median_delta_phi_2_obs",
            "rate_mean",
            "rate_ci_low",
            "rate_ci_high",
        ],
        rate_rows,
    )
    write_csv(
        output_dir / "division_counts_summary.csv",
        [
            "lx",
            "p0",
            "seed",
            "dphi_target",
            "n_j",
            "phi_j",
            "phi_2_obs",
            "delta_phi_2_obs",
            "delta_phi_final",
            "exhausted_by_endpoint",
            "N_div_all_at_phi2",
            "N_div_all_final",
            "N_div_all_after_phi2",
            "frac_before_phi2",
            "N_div_track_at_phi2",
            "N_div_track_final",
            "divlog_postjam_count",
            "postjamm_summary_final_count",
            "postjamm_summary_count_gap",
        ],
        summary_rows,
    )

    if not summaries:
        raise SystemExit("No validated runs found for the requested sizes.")
    avg_nj_by_size = average_nj_by_size(summaries)
    make_figure(output_dir, int(args.panel_a_size), panel_b_sizes, cumulative_rows, phi2_rows, avg_nj_by_size, args.dpi)

    mismatch_count = sum(1 for summary in summaries if summary.postjamm_summary_count_gap != 0)
    print(f"Selected representative runs: {len(summaries)}")
    print(f"Runs with DIVLOG / POSTJAMM_SUMMARY count mismatch: {mismatch_count}")
    print(f"Wrote outputs to: {output_dir}")


if __name__ == "__main__":
    main()
