#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import math
import os
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mplconfig_"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.mkdtemp(prefix="xdgcache_"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from pipeline_config import DEFAULT_DPHI_PROBE, REPO_ROOT, basename, iter_job_params, seed_range
from pipeline_paths import bext_paths, growth_paths, shear_paths
from pipeline_validate import bext_done, growth_done, shear_done


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
DPHI_RANGE_STYLES = [
    ("1e-4 to 1e-3", 1e-4, 1.01e-3, "o"),
    ("3e-3 to 1e-2", 3e-3, 1.01e-2, "s"),
    ("3e-2 to 1.2e-1", 3e-2, 1.21e-1, "^"),
]
P_MAX = 0.1
SECONDARY_ARREST_ACTIVE_CUTOFF = 1
GROWTH_CROSSOVER_DPHI = 3e-3
GROWTH_HISTOGRAM_P0_ORDER = ["1e-1", "1e-2", "1e-3", "1e-4"]
GROWTH_HISTOGRAM_TEXT = {
    "1e-1": r"$P_0=10^{-1}$",
    "1e-2": r"$P_0=10^{-2}$",
    "1e-3": r"$P_0=10^{-3}$",
    "1e-4": r"$P_0=10^{-4}$",
}
@dataclass
class PostRow:
    phi: float
    pressure: float
    n: int
    nc: int
    nf: int
    nu: int
    ziso: int
    chi_c: float
    initial_free_total: int
    initial_free_active: int
    postjam_divisions_total: int

    @property
    def nonfloating(self) -> int:
        return self.n - self.nf

    @property
    def delta_z(self) -> float:
        return (self.nc - self.ziso) / self.nonfloating

    @property
    def u_total(self) -> float:
        return self.nu / self.nonfloating

    @property
    def u_j(self) -> float:
        return self.initial_free_total / self.nonfloating

    @property
    def active_fraction(self) -> float:
        if self.initial_free_total == 0:
            return 0.0
        return self.initial_free_active / self.initial_free_total


@dataclass
class CohortRow:
    status_code: int
    completion_delta_area: float
    division_delta_area: float
    final_delta_area: float


@dataclass
class BextRow:
    bext: float


@dataclass
class JobRecord:
    lx: int
    p0: str
    dphi: float
    seed: str
    post_rows: list[PostRow]
    cohort_rows: list[CohortRow]
    bext_row: BextRow | None
    shear_phi2_data_path: Path | None
    stats_path: Path

    @property
    def phi_j(self) -> float:
        return self.post_rows[0].phi

    @property
    def endpoint(self) -> PostRow:
        return self.post_rows[-1]

    @property
    def delta_phi(self) -> float:
        return self.endpoint.phi - self.phi_j

    @property
    def final_divisions(self) -> int:
        return self.endpoint.postjam_divisions_total


@dataclass
class ArrestPair:
    p0: str
    seed: str
    dphi: float
    phi2: float
    p2: float
    g2: float
    kept: bool = True


@dataclass
class ArrestSelection:
    p0: str
    seed: str
    record: JobRecord
    row_cutoff: PostRow


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate publication-quality paper figures for L=15.")
    parser.add_argument("--size", default="15")
    parser.add_argument("--output-dir", default=str(REPO_ROOT / "paper" / "figures"))
    parser.add_argument("--dphi-probe", type=float, default=DEFAULT_DPHI_PROBE)
    parser.add_argument("--max-postjam-divisions", type=int, default=10)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--seed-start", type=int, help="Optional inclusive seed start.")
    parser.add_argument("--seed-stop", type=int, help="Optional inclusive seed stop.")
    args = parser.parse_args()
    if (args.seed_start is None) ^ (args.seed_stop is None):
        raise SystemExit("--seed-start and --seed-stop must be given together")
    if args.seed_start is not None and args.seed_stop < args.seed_start:
        raise SystemExit("--seed-stop must be greater than or equal to --seed-start")
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


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def mean(values: list[float]) -> float:
    return float(sum(values) / len(values))


def stderr(values: list[float]) -> float:
    if len(values) <= 1:
        return 0.0
    return float(np.std(values, ddof=1) / math.sqrt(len(values)))


def bootstrap_median_stderr(values: list[float], n_boot: int = 2000, seed: int = 12345) -> float:
    if len(values) <= 1:
        return 0.0
    array = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    medians = np.empty(n_boot, dtype=float)
    for index in range(n_boot):
        sample = rng.choice(array, size=len(array), replace=True)
        medians[index] = float(np.median(sample))
    return float(np.std(medians, ddof=1))


def tukey_filter(values: list[tuple[str, float]], fence_scale: float = 1.5) -> list[tuple[str, float]]:
    if len(values) < 4:
        return values
    array = np.asarray([value for _seed, value in values], dtype=float)
    q1 = float(np.percentile(array, 25.0))
    q3 = float(np.percentile(array, 75.0))
    iqr = q3 - q1
    if iqr <= 0.0:
        return values
    lo = q1 - fence_scale * iqr
    hi = q3 + fence_scale * iqr
    filtered = [(seed, value) for seed, value in values if lo <= value <= hi]
    return filtered if filtered else values


def parse_post_rows(path: Path) -> list[PostRow]:
    rows: list[PostRow] = []
    with path.open("r", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            if not line.strip():
                continue
            ff = line.split()
            rows.append(
                PostRow(
                    phi=float(ff[1]),
                    pressure=float(ff[2]),
                    n=int(ff[3]),
                    nc=int(ff[4]),
                    nf=int(ff[5]),
                    nu=int(ff[6]),
                    ziso=int(ff[7]),
                    chi_c=float(ff[9]),
                    initial_free_total=int(ff[11]),
                    initial_free_active=int(ff[12]),
                    postjam_divisions_total=int(ff[18]),
                )
            )
    return rows


def parse_cohort_rows(path: Path) -> list[CohortRow]:
    rows: list[CohortRow] = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            ff = line.split()
            rows.append(
                CohortRow(
                    status_code=int(ff[3]),
                    completion_delta_area=float(ff[10]),
                    division_delta_area=float(ff[13]),
                    final_delta_area=float(ff[17]),
                )
            )
    return rows


def parse_bext(path: Path) -> BextRow:
    with path.open("r", encoding="utf-8") as handle:
        next(handle)
        ff = next(handle).split()
    return BextRow(bext=float(ff[5]))


def estimate_growth_bulk(post_rows: list[PostRow], n_tail: int = 5) -> float:
    tail = post_rows[-min(len(post_rows), n_tail) :]
    phi = np.asarray([row.phi for row in tail], dtype=float)
    pressure = np.asarray([row.pressure for row in tail], dtype=float)
    if len(phi) < 2:
        return float("nan")
    slope = np.polyfit(phi, pressure, deg=1)[0]
    return float(phi[-1] * slope)


def estimate_shear_modulus(path: Path, n_head: int = 20) -> float:
    strain = []
    stress = []
    topology_ref = None
    with path.open("r", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            if not line.strip():
                continue
            ff = line.split()
            topology = (
                int(ff[4]),
                int(ff[5]),
                int(ff[6]),
                int(ff[7]),
                int(ff[9]),
                int(ff[10]),
                int(ff[11]),
            )
            if topology_ref is None:
                topology_ref = topology
            elif topology != topology_ref:
                break
            strain.append(float(ff[0]))
            stress.append(float(ff[2]))
            if len(strain) >= n_head:
                break
    if len(strain) < 2:
        return float("nan")
    return float(np.polyfit(np.asarray(strain), np.asarray(stress), deg=1)[0])


def parse_growth_rates(path: Path) -> list[float]:
    growth_rates: list[float] = []
    with path.open("r", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            if not line.strip():
                continue
            ff = line.split()
            if int(ff[3]) != 1:
                continue
            innate_rate = float(ff[6])
            if innate_rate <= 0.0:
                continue
            growth_rates.append(float(ff[5]) / innate_rate)
    return growth_rates


def load_records(size: str, dphi_probe: float, seeds: list[str] | None = None) -> list[JobRecord]:
    records: list[JobRecord] = []
    for params in iter_job_params(sizes=[size], seeds=seeds):
        name = basename(params)
        growth = growth_paths(name)
        if not growth_done(growth):
            continue
        shear_phi2 = shear_paths(name, source_tag="PHI2")
        bext = bext_paths(name, dphi_probe)
        cohort_path = growth["cohort_gz"] if growth["cohort_gz"].is_file() else growth["cohort"]
        records.append(
            JobRecord(
                lx=int(params["lx"]),
                p0=params["p0"],
                dphi=float(params["dphi"]),
                seed=params["seed"],
                post_rows=parse_post_rows(growth["postjamm_summary"]),
                cohort_rows=parse_cohort_rows(cohort_path),
                bext_row=parse_bext(bext["data"]) if bext_done(bext) else None,
                shear_phi2_data_path=shear_phi2["g_data"] if shear_done(shear_phi2) else None,
                stats_path=growth["stats_frame"],
            )
        )
    return records


def dphi_range_label(value: float) -> str:
    for label, lo, hi, _marker in DPHI_RANGE_STYLES:
        if lo <= value < hi:
            return label
    return "other"


def dphi_range_marker(value: float) -> str:
    for _label, lo, hi, marker in DPHI_RANGE_STYLES:
        if lo <= value < hi:
            return marker
    return "D"


def build_depletion_curves(records: list[JobRecord], max_postjam_divisions: int) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    per_seed: dict[tuple[str, str], list[JobRecord]] = defaultdict(list)
    for record in records:
        per_seed[(record.p0, record.seed)].append(record)

    pooled: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for (p0, _seed), group in per_seed.items():
        eligible = [record for record in group if record.final_divisions <= max_postjam_divisions]
        if not eligible:
            continue
        chosen = max(eligible, key=lambda record: record.delta_phi)
        for row in chosen.post_rows:
            pooled[p0].append((row.phi - chosen.phi_j, row.active_fraction))

    curves = {}
    for p0, points in pooled.items():
        points.sort(key=lambda item: item[0])
        x = np.asarray([item[0] for item in points], dtype=float)
        y = np.asarray([item[1] for item in points], dtype=float)
        if len(x) == 0:
            continue
        bins = np.linspace(0.0, float(x.max()), 26)
        centers = 0.5 * (bins[:-1] + bins[1:])
        curve_y = []
        curve_x = []
        for left, right, center in zip(bins[:-1], bins[1:], centers):
            mask = (x >= left) & (x < right if right < bins[-1] else x <= right)
            if np.any(mask):
                curve_x.append(center)
                curve_y.append(float(np.mean(y[mask])))
        curves[p0] = (np.asarray(curve_x), np.asarray(curve_y))
    return curves


def first_active_cutoff_row(record: JobRecord, active_cutoff: int) -> PostRow | None:
    return next((row for row in record.post_rows if row.initial_free_active <= active_cutoff), None)


def choose_arrest_selection(
    group: list[JobRecord],
    max_postjam_divisions: int,
    active_cutoff: int = SECONDARY_ARREST_ACTIVE_CUTOFF,
    require_phi2_shear: bool = False,
) -> ArrestSelection | None:
    candidates: list[ArrestSelection] = []
    for record in group:
        if require_phi2_shear and record.shear_phi2_data_path is None:
            continue
        if record.final_divisions > max_postjam_divisions:
            continue
        row_cutoff = first_active_cutoff_row(record, active_cutoff)
        if row_cutoff is None:
            continue
        candidates.append(
            ArrestSelection(
                p0=record.p0,
                seed=record.seed,
                record=record,
                row_cutoff=row_cutoff,
            )
        )
    if not candidates:
        return None
    key = (
        (lambda selection: selection.record.dphi)
        if require_phi2_shear
        else (lambda selection: selection.record.delta_phi)
    )
    return min(candidates, key=key)


def build_arrest_selections(
    records: list[JobRecord],
    max_postjam_divisions: int,
    active_cutoff: int = SECONDARY_ARREST_ACTIVE_CUTOFF,
    require_phi2_shear: bool = False,
) -> list[ArrestSelection]:
    per_seed: dict[tuple[str, str], list[JobRecord]] = defaultdict(list)
    for record in records:
        per_seed[(record.p0, record.seed)].append(record)

    selections: list[ArrestSelection] = []
    for group in per_seed.values():
        chosen = choose_arrest_selection(
            group,
            max_postjam_divisions,
            active_cutoff=active_cutoff,
            require_phi2_shear=require_phi2_shear,
        )
        if chosen is not None:
            selections.append(chosen)

    selections.sort(key=lambda selection: (float(selection.p0), selection.seed))
    return selections


def baseline_delta_z_by_pressure(records: list[JobRecord]) -> tuple[np.ndarray, np.ndarray]:
    grouped: dict[float, list[JobRecord]] = defaultdict(list)
    for record in records:
        if record.p0 != "-1":
            continue
        grouped[record.dphi].append(record)

    pressure_points = []
    delta_z_points = []
    for dphi in sorted(grouped):
        group = grouped[dphi]
        pressure_points.append(mean([record.endpoint.pressure for record in group]))
        delta_z_points.append(mean([record.endpoint.delta_z for record in group]))

    if not pressure_points:
        return np.asarray([], dtype=float), np.asarray([], dtype=float)

    order = np.argsort(np.asarray(pressure_points, dtype=float))
    return (
        np.asarray(pressure_points, dtype=float)[order],
        np.asarray(delta_z_points, dtype=float)[order],
    )


def build_structural_points(records: list[JobRecord]) -> list[dict]:
    baseline_pressure, baseline_delta_z = baseline_delta_z_by_pressure(records)
    if len(baseline_pressure) == 0:
        return []

    grouped: dict[tuple[str, float], list[JobRecord]] = defaultdict(list)
    for record in records:
        if record.p0 == "-1":
            continue
        grouped[(record.p0, record.dphi)].append(record)

    points = []
    for (p0, dphi), group in grouped.items():
        predicted = []
        observed = []
        for record in group:
            dz0 = float(np.interp(record.endpoint.pressure, baseline_pressure, baseline_delta_z))
            predicted.append(dz0 + 2.0 * (record.endpoint.u_j - record.endpoint.u_total))
            observed.append(record.endpoint.delta_z)
        points.append(
            {
                "p0": p0,
                "dphi": dphi,
                "range_label": dphi_range_label(dphi),
                "marker": dphi_range_marker(dphi),
                "pred_mean": mean(predicted),
                "pred_err": stderr(predicted),
                "obs_mean": mean(observed),
                "obs_err": stderr(observed),
                "count": len(group),
            }
        )
    points.sort(key=lambda row: (float(row["p0"]), row["dphi"]))
    return points


def cohort_time_and_event(row: CohortRow) -> tuple[float, bool] | None:
    if row.status_code == 1 and row.completion_delta_area > 0.0:
        return row.completion_delta_area, True
    if row.status_code == 2 and row.division_delta_area > 0.0:
        return row.division_delta_area, False
    if row.status_code == 0 and row.final_delta_area > 0.0:
        return row.final_delta_area, False
    return None


def kaplan_meier_a_epsilon_and_integral(record: JobRecord, active_cutoff: int) -> tuple[float, float] | None:
    observations = []
    for row in record.cohort_rows:
        parsed = cohort_time_and_event(row)
        if parsed is not None:
            observations.append(parsed)
    if not observations:
        return None

    initial_free_total = record.post_rows[0].initial_free_total
    if initial_free_total <= 0:
        return None
    if initial_free_total <= active_cutoff:
        return 0.0, 0.0

    target_survival = active_cutoff / initial_free_total
    times = sorted({time for time, _event in observations})
    survival = 1.0
    prev_time = 0.0
    integral = 0.0
    at_risk = len(observations)

    for time in times:
        integral += survival * (time - prev_time)
        events = sum(1 for obs_time, event in observations if obs_time == time and event)
        censored = sum(1 for obs_time, event in observations if obs_time == time and not event)
        if at_risk <= 0:
            break
        if events > 0:
            survival_after = survival * (1.0 - events / at_risk)
            if survival_after <= target_survival:
                return time, integral
            survival = survival_after
        at_risk -= events + censored
        prev_time = time
    return None


def predicted_phi2(
    record: JobRecord,
    row_cutoff: PostRow,
    active_cutoff: int = SECONDARY_ARREST_ACTIVE_CUTOFF,
) -> float | None:
    km = kaplan_meier_a_epsilon_and_integral(record, active_cutoff)
    if km is None:
        return None
    a_eps, integral_s = km
    nonfloating = record.post_rows[0].nonfloating
    n_j = nonfloating / float(record.lx * record.lx)
    u_j = record.post_rows[0].u_j
    chi_bar = mean([row.chi_c for row in record.post_rows if row.phi <= row_cutoff.phi])
    return record.phi_j + n_j * (chi_bar * a_eps + (1.0 - chi_bar) * u_j * integral_s)


def build_phi2_summary(selections: list[ArrestSelection], active_cutoff: int = SECONDARY_ARREST_ACTIVE_CUTOFF) -> list[dict]:
    raw = []
    for selection in selections:
        pred = predicted_phi2(selection.record, selection.row_cutoff, active_cutoff=active_cutoff)
        if pred is None:
            continue
        raw.append(
            {
                "p0": selection.p0,
                "seed": selection.seed,
                "phi2_meas": selection.row_cutoff.phi,
                "p2_meas": selection.row_cutoff.pressure,
                "phi2_pred": pred,
            }
        )

    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in raw:
        grouped[row["p0"]].append(row)

    summary = []
    for p0 in P0_ORDER:
        values = grouped.get(p0, [])
        if not values:
            continue
        meas = [row["phi2_meas"] for row in values]
        pred = [row["phi2_pred"] for row in values]
        p2 = [row["p2_meas"] for row in values]
        summary.append(
            {
                "p0": p0,
                "count": len(values),
                "phi2_meas": mean(meas),
                "phi2_meas_err": stderr(meas),
                "phi2_pred": mean(pred),
                "phi2_pred_err": stderr(pred),
                "p2": mean(p2),
                "p2_err": stderr(p2),
                "p2_norm": mean(p2) / P_MAX,
                "p2_norm_err": stderr(p2) / P_MAX,
            }
        )
    return summary


def build_mechanics_points(records: list[JobRecord]) -> list[dict]:
    grouped: dict[tuple[str, float], list[JobRecord]] = defaultdict(list)
    for record in records:
        if record.bext_row is None:
            continue
        grouped[(record.p0, record.dphi)].append(record)

    rows = []
    for (p0, dphi), group in grouped.items():
        predicted = []
        measured = []
        for record in group:
            endpoint = record.endpoint
            denominator = endpoint.u_total + (1.0 - endpoint.u_total) * endpoint.chi_c
            if denominator <= 0.0:
                continue
            lambda_c = ((1.0 - endpoint.u_total) * endpoint.chi_c) / denominator
            predicted.append(lambda_c * record.bext_row.bext)
            measured.append(estimate_growth_bulk(record.post_rows))
        if not predicted:
            continue
        rows.append(
            {
                "p0": p0,
                "dphi": dphi,
                "marker": dphi_range_marker(dphi),
                "pred_mean": mean(predicted),
                "pred_err": stderr(predicted),
                "meas_mean": mean(measured),
                "meas_err": stderr(measured),
                "count": len(predicted),
            }
        )
    rows.sort(key=lambda row: (float(row["p0"]), row["dphi"]))
    return rows


def build_paired_arrest_points(
    selections: list[ArrestSelection],
) -> list[ArrestPair]:
    pairs: list[ArrestPair] = []
    for selection in selections:
        if selection.record.shear_phi2_data_path is None:
            continue
        value = estimate_shear_modulus(selection.record.shear_phi2_data_path)
        if not math.isfinite(value):
            continue
        pairs.append(
            ArrestPair(
                p0=selection.p0,
                seed=selection.seed,
                dphi=selection.record.dphi,
                phi2=selection.row_cutoff.phi,
                p2=selection.row_cutoff.pressure,
                g2=value,
            )
        )

    pairs.sort(key=lambda pair: (float(pair.p0), pair.seed))
    return pairs


def apply_g2_tukey_filter(pairs: list[ArrestPair]) -> None:
    grouped: dict[str, list[ArrestPair]] = defaultdict(list)
    for pair in pairs:
        pair.kept = pair.g2 > 0.0 and math.isfinite(pair.g2)
        if pair.kept:
            grouped[pair.p0].append(pair)

    for group in grouped.values():
        filtered = tukey_filter([(pair.seed, pair.g2) for pair in group])
        keep = {seed for seed, _value in filtered}
        for pair in group:
            pair.kept = pair.seed in keep


def build_g2_summary(pairs: list[ArrestPair]) -> list[dict]:
    raw: dict[str, list[float]] = defaultdict(list)
    for pair in pairs:
        if pair.kept:
            raw[pair.p0].append(pair.g2)

    summary = []
    for p0 in P0_ORDER:
        values = raw.get(p0, [])
        if not values:
            continue
        summary.append(
            {
                "p0": p0,
                "g2": float(np.median(np.asarray(values, dtype=float))),
                "g2_err": bootstrap_median_stderr(values),
                "count": len(values),
            }
        )
    return summary


def build_paired_g2_over_p2_summary(pairs: list[ArrestPair]) -> list[dict]:
    grouped: dict[str, list[float]] = defaultdict(list)
    for pair in pairs:
        if not pair.kept or pair.p2 <= 0.0:
            continue
        value = pair.g2 / pair.p2
        if not math.isfinite(value) or value <= 0.0:
            continue
        grouped[pair.p0].append(value)

    summary = []
    for p0 in P0_ORDER:
        values = grouped.get(p0, [])
        if not values:
            continue
        summary.append(
            {
                "p0": p0,
                "ratio": float(np.median(np.asarray(values, dtype=float))),
                "ratio_err": bootstrap_median_stderr(values),
                "count": len(values),
            }
        )
    return summary


def build_paired_g2_vs_p2_summary(pairs: list[ArrestPair]) -> list[dict]:
    grouped: dict[str, list[ArrestPair]] = defaultdict(list)
    for pair in pairs:
        if pair.kept and pair.p2 > 0.0 and pair.g2 > 0.0:
            grouped[pair.p0].append(pair)

    summary = []
    for p0 in P0_ORDER:
        values = grouped.get(p0, [])
        if not values:
            continue
        p2_values = [pair.p2 for pair in values]
        g2_values = [pair.g2 for pair in values]
        summary.append(
            {
                "p0": p0,
                "count": len(values),
                "p2": float(np.median(np.asarray(p2_values, dtype=float))),
                "p2_err": bootstrap_median_stderr(p2_values),
                "g2": float(np.median(np.asarray(g2_values, dtype=float))),
                "g2_err": bootstrap_median_stderr(g2_values),
            }
        )
    return summary


def build_lambda_c_summary(selections: list[ArrestSelection]) -> list[dict]:
    raw: dict[str, list[float]] = defaultdict(list)
    for selection in selections:
        lambda_values = []
        for row in selection.record.post_rows:
            if row.phi > selection.row_cutoff.phi:
                break
            denominator = row.u_total + (1.0 - row.u_total) * row.chi_c
            if denominator <= 0.0:
                continue
            lambda_values.append(((1.0 - row.u_total) * row.chi_c) / denominator)
        if lambda_values:
            raw[selection.p0].append(mean(lambda_values))

    summary = []
    for p0 in P0_ORDER:
        values = raw.get(p0, [])
        if not values:
            continue
        summary.append(
            {
                "p0": p0,
                "lambda_c_bar": mean(values),
                "lambda_c_bar_err": stderr(values),
                "count": len(values),
            }
        )
    return summary


def build_growth_hgamma_summary(records: list[JobRecord], target_dphi: float) -> list[dict]:
    raw: dict[str, list[dict]] = defaultdict(list)
    for record in records:
        if not math.isclose(record.dphi, target_dphi, rel_tol=0.0, abs_tol=1e-12):
            continue
        growth_rates = parse_growth_rates(record.stats_path)
        if not growth_rates:
            continue
        mean_g = mean(growth_rates)
        variance_g = float(np.var(np.asarray(growth_rates, dtype=float)))
        h_gamma = variance_g / (mean_g * mean_g) if mean_g > 0.0 else float("nan")
        spike_fraction = sum(1 for value in growth_rates if math.isclose(value, 1.0, rel_tol=0.0, abs_tol=1e-12)) / len(growth_rates)
        raw[record.p0].append({"h_gamma": h_gamma, "spike_fraction": spike_fraction})

    summary = []
    for p0 in P0_ORDER:
        rows = raw.get(p0, [])
        if not rows:
            continue
        h_values = [row["h_gamma"] for row in rows]
        spike_values = [row["spike_fraction"] for row in rows]
        summary.append(
            {
                "p0": p0,
                "dphi": target_dphi,
                "h_gamma": mean(h_values),
                "h_gamma_err": stderr(h_values),
                "spike_fraction": mean(spike_values),
                "spike_fraction_err": stderr(spike_values),
                "count": len(rows),
            }
        )
    return summary


def build_growth_histograms(records: list[JobRecord], target_dphi: float) -> dict[str, dict]:
    pooled: dict[str, list[float]] = defaultdict(list)
    for record in records:
        if not math.isclose(record.dphi, target_dphi, rel_tol=0.0, abs_tol=1e-12):
            continue
        if record.p0 not in GROWTH_HISTOGRAM_P0_ORDER:
            continue
        pooled[record.p0].extend(parse_growth_rates(record.stats_path))

    edges = np.asarray([0.0, 0.03, 0.06, 0.10, 0.15, 0.22, 0.32, 0.45, 0.60, 0.78, 0.90, 0.97, 0.995], dtype=float)
    histograms: dict[str, dict] = {}
    for p0 in GROWTH_HISTOGRAM_P0_ORDER:
        values = pooled.get(p0, [])
        if not values:
            continue
        spike_count = sum(1 for value in values if math.isclose(value, 1.0, rel_tol=0.0, abs_tol=1e-12))
        continuous = [value for value in values if value < 1.0 - 1e-12]
        counts, _ = np.histogram(np.asarray(continuous, dtype=float), bins=edges)
        total = len(values)
        widths = np.diff(edges)
        density = counts.astype(float) / (total * widths)
        histograms[p0] = {
            "dphi": target_dphi,
            "edges": edges,
            "density": density,
            "spike_fraction": spike_count / total,
            "count": total,
        }
    return histograms


def add_panel_label(ax, text: str) -> None:
    ax.text(-0.18, 1.05, text, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def save_figure(fig: plt.Figure, output_dir: Path, stem: str, dpi: int) -> None:
    fig.savefig(output_dir / f"{stem}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_structure_figure(
    output_dir: Path,
    structural_points: list[dict],
    depletion_curves: dict[str, tuple[np.ndarray, np.ndarray]],
    phi2_summary: list[dict],
    dpi: int,
    stem: str = "figure2_structure_l15",
) -> None:
    fig = plt.figure(figsize=(7.1, 4.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1.0], height_ratios=[1.0, 0.95], wspace=0.42, hspace=0.45)
    ax_struct = fig.add_subplot(gs[:, 0])
    ax_deplete = fig.add_subplot(gs[0, 1])
    ax_phi2 = fig.add_subplot(gs[1, 1])

    add_panel_label(ax_struct, "A")
    add_panel_label(ax_deplete, "B")
    add_panel_label(ax_phi2, "C")

    pred = np.asarray([row["pred_mean"] for row in structural_points])
    obs = np.asarray([row["obs_mean"] for row in structural_points])
    lo = float(min(pred.min(), obs.min()) * 0.92)
    hi = float(max(pred.max(), obs.max()) * 1.03)
    ax_struct.plot([lo, hi], [lo, hi], ls="--", lw=1.0, color="0.55", zorder=1)
    for row in structural_points:
        ax_struct.errorbar(
            row["pred_mean"],
            row["obs_mean"],
            fmt=row["marker"],
            ms=5.8,
            mec="white",
            mew=0.5,
            color=P0_COLORS[row["p0"]],
            alpha=0.95,
            zorder=3,
        )
    ax_struct.set_xlim(lo, hi)
    ax_struct.set_ylim(lo, hi)
    ax_struct.set_xlabel(r"Eq. (6) prediction for $\Delta Z$")
    ax_struct.set_ylabel(r"measured $\Delta Z$")
    ax_struct.set_title(r"Structural branch from $\Delta Z = \Delta Z_0 + 2(u_J-u)$", pad=6)

    color_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=P0_COLORS[p0], markeredgecolor="white", markeredgewidth=0.5, markersize=6, label=P0_LABELS[p0])
        for p0 in P0_ORDER
        if any(row["p0"] == p0 for row in structural_points)
    ]
    marker_handles = [
        Line2D(
            [0],
            [0],
            marker=marker,
            color="0.25",
            lw=0,
            markersize=6,
            markerfacecolor="0.25",
            markeredgecolor="0.25",
            markeredgewidth=0.0,
            label=label,
        )
        for label, _lo, _hi, marker in DPHI_RANGE_STYLES
    ]
    leg1 = ax_struct.legend(handles=color_handles, title=r"$P_0$", loc="upper left", bbox_to_anchor=(0.0, 1.02), ncol=2, handletextpad=0.35, columnspacing=0.8)
    ax_struct.add_artist(leg1)
    ax_struct.legend(handles=marker_handles, title=r"$\Delta\phi$ range", loc="lower right", handletextpad=0.35)

    for p0 in P0_ORDER:
        if p0 not in depletion_curves:
            continue
        x, y = depletion_curves[p0]
        ax_deplete.plot(x, y, lw=2.0, color=P0_COLORS[p0], label=P0_LABELS[p0])
    ax_deplete.set_xlabel(r"$\phi-\phi_J$")
    ax_deplete.set_ylabel(r"active initial-free fraction")
    ax_deplete.set_ylim(-0.03, 1.03)
    ax_deplete.set_title("Reservoir depletion", pad=6)
    ax_deplete.legend(loc="upper right", fontsize=8, ncol=2, handlelength=2.2, columnspacing=0.9)

    x = np.arange(len(phi2_summary))
    meas = [row["phi2_meas"] for row in phi2_summary]
    pred_y = [row["phi2_pred"] for row in phi2_summary]
    ax_phi2.errorbar(x, meas, yerr=[row["phi2_meas_err"] for row in phi2_summary], fmt="o", color="#222222", ms=5, capsize=2, label="simulation")
    ax_phi2.errorbar(x, pred_y, yerr=[row["phi2_pred_err"] for row in phi2_summary], fmt="s", mfc="white", mec="#b23a48", color="#b23a48", ms=5, capsize=2, label="theory")
    ax_phi2.set_xticks(x, [P0_LABELS[row["p0"]] for row in phi2_summary])
    ax_phi2.set_ylabel(r"$\phi_2$")
    ax_phi2.set_xlabel(r"$P_0$")
    ax_phi2.set_title("Secondary arrest", pad=6)
    ax_phi2.legend(loc="best", fontsize=8)

    save_figure(fig, output_dir, stem, dpi)


def make_mechanics_figure(
    output_dir: Path,
    mechanics_points: list[dict],
    phi2_summary: list[dict],
    g2_summary: list[dict],
    paired_arrest_points: list[ArrestPair],
    lambda_c_summary: list[dict],
    dpi: int,
    stem: str = "figure3_mechanics_l15",
) -> None:
    fig = plt.figure(figsize=(7.1, 4.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 0.9], height_ratios=[1.0, 1.0], wspace=0.42, hspace=0.45)
    ax_id = fig.add_subplot(gs[:, 0])
    ax_p2 = fig.add_subplot(gs[0, 1])
    ax_g2 = fig.add_subplot(gs[1, 1])

    add_panel_label(ax_id, "A")
    add_panel_label(ax_p2, "B")
    add_panel_label(ax_g2, "C")

    positive = [row for row in mechanics_points if row["pred_mean"] > 0.0 and row["meas_mean"] > 0.0]
    xs = np.asarray([row["pred_mean"] for row in positive])
    ys = np.asarray([row["meas_mean"] for row in positive])
    lo = 0.0
    hi = float(max(xs.max(), ys.max()) * 1.15)
    ax_id.plot([lo, hi], [lo, hi], ls="--", lw=1.0, color="0.55", zorder=1)
    for row in positive:
        ax_id.errorbar(
            row["pred_mean"],
            row["meas_mean"],
            xerr=row["pred_err"],
            yerr=row["meas_err"],
            fmt=row["marker"],
            ms=6.2,
            mec="white",
            mew=0.75,
            elinewidth=1.0,
            capsize=2.0,
            color=P0_COLORS[row["p0"]],
            alpha=0.95,
            zorder=3,
        )
    ax_id.set_xlim(lo, hi)
    ax_id.set_ylim(lo, hi)
    ax_id.set_xlabel(r"$\lambda_c B_{\mathrm{ext}}$")
    ax_id.set_ylabel(r"$B_{\mathrm{grow}}$")
    ax_id.set_title(r"Flux-partition identity $B_{\mathrm{grow}}=\lambda_c B_{\mathrm{ext}}$", pad=6)

    color_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=P0_COLORS[p0], markeredgecolor="white", markeredgewidth=0.5, markersize=6, label=P0_LABELS[p0])
        for p0 in P0_ORDER
        if any(row["p0"] == p0 for row in positive)
    ]
    marker_handles = [
        Line2D([0], [0], marker=marker, color="0.25", lw=0, markersize=6, label=label)
        for label, _lo, _hi, marker in DPHI_RANGE_STYLES
    ]
    leg1 = ax_id.legend(handles=color_handles, title=r"$P_0$", loc="upper left", bbox_to_anchor=(0.0, 1.02), ncol=2, handletextpad=0.35, columnspacing=0.8)
    ax_id.add_artist(leg1)
    ax_id.legend(handles=marker_handles, title=r"$\Delta\phi$ range", loc="lower right", handletextpad=0.35)

    x = np.arange(len(phi2_summary))
    ax_p2.errorbar(
        x,
        [row["p2_norm"] for row in phi2_summary],
        yerr=[row["p2_norm_err"] for row in phi2_summary],
        fmt="o-",
        color="#222222",
        lw=1.3,
        ms=5,
        capsize=2,
    )
    ax_p2.set_yscale("log")
    ax_p2.set_xticks(x, [P0_LABELS[row["p0"]] for row in phi2_summary])
    ax_p2.set_ylabel(r"$P_2/P_{\max}$")
    ax_p2.set_xlabel(r"$P_0$")
    ax_p2.set_title("Pressure at arrest", pad=6)

    inset = ax_p2.inset_axes([0.10, 0.14, 0.42, 0.42])
    xin = np.arange(len(lambda_c_summary))
    inset.errorbar(
        xin,
        [row["lambda_c_bar"] for row in lambda_c_summary],
        yerr=[row["lambda_c_bar_err"] for row in lambda_c_summary],
        fmt="o-",
        color="#444444",
        lw=1.1,
        ms=3.6,
        capsize=2,
    )
    inset.set_xticks(xin, [P0_LABELS[row["p0"]] for row in lambda_c_summary])
    inset.set_ylim(0.0, 1.02)
    inset.yaxis.set_label_position("left")
    inset.yaxis.tick_left()
    inset.set_ylabel(r"$\bar{\lambda}_c$", fontsize=7, labelpad=1)
    inset.yaxis.set_label_coords(0.18, 0.5)
    inset.set_yticks([0.0, 0.5, 1.0])
    inset.set_yticklabels(["0", "0.5", "1.0"])
    inset.set_title("mean compressed-\nsector flux", fontsize=6.5, pad=1.5)
    inset.tick_params(axis="both", labelsize=6, pad=1)
    inset.spines["top"].set_visible(False)
    inset.spines["left"].set_visible(True)
    inset.spines["right"].set_visible(False)

    xg = np.arange(len(g2_summary))
    xg_by_p0 = {row["p0"]: index for index, row in enumerate(g2_summary)}
    jitter_rng = np.random.default_rng(12345)
    for pair in paired_arrest_points:
        if not pair.kept or pair.p0 not in xg_by_p0:
            continue
        x_value = xg_by_p0[pair.p0] + float(jitter_rng.uniform(-0.085, 0.085))
        ax_g2.plot(
            x_value,
            pair.g2,
            "o",
            ms=3.5,
            color="#9aa1a8",
            alpha=0.5,
            mec="none",
            zorder=1,
        )
    ax_g2.errorbar(
        xg,
        [row["g2"] for row in g2_summary],
        yerr=[row["g2_err"] for row in g2_summary],
        fmt="o-",
        color="#222222",
        lw=1.3,
        ms=5,
        capsize=2,
    )
    ax_g2.set_xticks(xg, [P0_LABELS[row["p0"]] for row in g2_summary])
    ax_g2.set_ylabel(r"$G_2$")
    ax_g2.set_xlabel(r"$P_0$")
    ax_g2.set_title("Shear stiffness at arrest", pad=6)

    save_figure(fig, output_dir, stem, dpi)


def make_growth_crossover_figure(output_dir: Path, hgamma_summary: list[dict], histograms: dict[str, dict], dpi: int) -> None:
    fig = plt.figure(figsize=(7.1, 3.2))
    ax_h = fig.add_subplot(1, 1, 1)

    plot_rows = [row for row in hgamma_summary if row["p0"] in ["-1", *GROWTH_HISTOGRAM_P0_ORDER, "1e-5"]]
    x = np.arange(len(plot_rows))
    x_by_p0 = {}
    ax_h.plot(x, [row["h_gamma"] for row in plot_rows], color="0.35", lw=1.1, zorder=1)
    for index, row in enumerate(plot_rows):
        x_by_p0[row["p0"]] = index
        ax_h.errorbar(
            index,
            row["h_gamma"],
            yerr=row["h_gamma_err"],
            fmt="o",
            ms=5.2,
            mec="white",
            mew=0.5,
            elinewidth=0.9,
            capsize=2.0,
            color=P0_COLORS[row["p0"]],
            zorder=3,
        )
    ax_h.set_xticks(x, [P0_LABELS[row["p0"]] for row in plot_rows])
    ax_h.set_ylabel(r"$h_\gamma$", fontsize=12)
    ax_h.set_xlabel(r"$P_0$", fontsize=12)
    ax_h.set_ylim(-0.08, 6.25)
    ax_h.set_title(r"Growth-rate crossover at $\Delta\phi=3\times10^{-3}$", pad=6)

    inset_positions = {
        "1e-1": [0.10, 0.60, 0.15, 0.22],
        "1e-2": [0.31, 0.60, 0.15, 0.22],
        "1e-3": [0.54, 0.60, 0.15, 0.22],
        "1e-4": [0.83, 0.12, 0.16, 0.24],
    }

    for p0 in GROWTH_HISTOGRAM_P0_ORDER:
        histogram = histograms[p0]
        axis = ax_h.inset_axes(inset_positions[p0])
        edges = histogram["edges"]
        density = histogram["density"]
        widths = np.diff(edges)
        axis.bar(
            edges[:-1],
            density,
            width=widths,
            align="edge",
            color=P0_COLORS[p0],
            edgecolor="white",
            linewidth=0.35,
        )
        spike_height = max(density.max() * 0.18, histogram["spike_fraction"] * 10.0)
        axis.plot(
            [1.0, 1.0],
            [0.0, spike_height],
            color=P0_COLORS[p0],
            lw=1.7,
            solid_capstyle="round",
        )
        axis.text(0.02, 0.82, GROWTH_HISTOGRAM_TEXT[p0], transform=axis.transAxes, fontsize=8)
        axis.set_facecolor("white")
        axis.patch.set_alpha(1.0)
        axis.set_ylim(0.0, max(density.max(), spike_height) * 1.08)
        axis.locator_params(axis="y", nbins=2)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.tick_params(axis="both", labelsize=6, pad=1)
        axis.set_xlim(0.0, 1.02)
        axis.set_xticks([0.0, 0.5, 1.0])
        axis.set_yticks([])
        axis.set_xlabel("")
        axis.set_ylabel("")
        x_frac, y_frac, width, _height = inset_positions[p0]
        point_y = next(row["h_gamma"] for row in plot_rows if row["p0"] == p0)
        ax_h.annotate(
            "",
            xy=(x_by_p0[p0], point_y),
            xycoords="data",
            xytext=(x_frac + 0.5 * width, y_frac),
            textcoords=ax_h.transAxes,
            arrowprops={"arrowstyle": "->", "color": "0.45", "lw": 0.65, "linestyle": "--", "mutation_scale": 8},
            zorder=2,
        )

    ax_h.text(0.52, 0.86, r"Representative pooled $\Pi(g)$", transform=ax_h.transAxes, ha="center", fontsize=9)

    save_figure(fig, output_dir, "figureS2_growth_crossover_l15", dpi)


def make_paired_arrest_figure(
    output_dir: Path,
    paired_arrest_points: list[ArrestPair],
    g2_vs_p2_summary: list[dict],
    g2_over_p2_summary: list[dict],
    dpi: int,
    stem: str = "figureS1_paired_arrest_l15",
) -> None:
    fig = plt.figure(figsize=(3.45, 4.75))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.06, 0.74], hspace=0.34)
    ax_scatter = fig.add_subplot(gs[0, 0])
    ax_ratio = fig.add_subplot(gs[1, 0])

    add_panel_label(ax_scatter, "A")
    add_panel_label(ax_ratio, "B")

    for pair in paired_arrest_points:
        if not pair.kept or pair.p2 <= 0.0 or pair.g2 <= 0.0:
            continue
        ax_scatter.scatter(
            pair.p2,
            pair.g2,
            s=28,
            color=P0_COLORS[pair.p0],
            alpha=0.45,
            edgecolors="white",
            linewidths=0.4,
            zorder=2,
        )

    for row in g2_vs_p2_summary:
        ax_scatter.errorbar(
            row["p2"],
            row["g2"],
            xerr=row["p2_err"],
            yerr=row["g2_err"],
            fmt="o",
            ms=6.2,
            mec="white",
            mew=0.6,
            elinewidth=1.0,
            capsize=2.0,
            color=P0_COLORS[row["p0"]],
            zorder=4,
        )
    ax_scatter.set_xscale("log")
    ax_scatter.set_xlabel(r"$P_2$")
    ax_scatter.set_ylabel(r"$G_2$")
    legend_handles = [
        Line2D(
            [0],
            [0],
            lw=0,
            color="none",
            label=r"$P_0$",
        )
    ]
    legend_handles.extend(
        [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=P0_COLORS[p0],
            markeredgecolor="white",
            markeredgewidth=0.5,
            markersize=6,
            label=P0_LABELS[p0],
        )
        for p0 in P0_ORDER
        if any(pair.kept and pair.p0 == p0 and pair.p2 > 0.0 and pair.g2 > 0.0 for pair in paired_arrest_points)
        ]
    )
    ax_scatter.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.14),
        borderaxespad=0.0,
        fontsize=8,
        ncol=4,
        columnspacing=0.8,
        handletextpad=0.35,
        labelspacing=0.45,
        alignment="left",
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        facecolor="white",
        edgecolor="0.4",
        borderpad=0.3,
    )

    xr = np.arange(len(g2_over_p2_summary))
    ax_ratio.errorbar(
        xr,
        [row["ratio"] for row in g2_over_p2_summary],
        yerr=[row["ratio_err"] for row in g2_over_p2_summary],
        fmt="o-",
        color="#222222",
        lw=1.3,
        ms=5,
        capsize=2,
        zorder=2,
    )
    for index, row in enumerate(g2_over_p2_summary):
        ax_ratio.plot(index, row["ratio"], "o", ms=5.5, color=P0_COLORS[row["p0"]], mec="white", mew=0.5, zorder=3)
    ax_ratio.set_yscale("log")
    ax_ratio.set_xticks(xr, [P0_LABELS[row["p0"]] for row in g2_over_p2_summary])
    ax_ratio.set_ylabel(r"$G_2/P_2$")
    ax_ratio.set_xlabel(r"$P_0$")

    save_figure(fig, output_dir, stem, dpi)


def main() -> None:
    args = parse_args()
    configure_style()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    seeds = None
    if args.seed_start is not None:
        seeds = seed_range(args.seed_start, args.seed_stop)

    records = load_records(args.size, args.dphi_probe, seeds=seeds)
    if not records:
        if seeds is None:
            raise SystemExit(f"No completed records found for L={args.size}")
        raise SystemExit(
            f"No completed records found for L={args.size} and seeds {args.seed_start}-{args.seed_stop}"
        )

    structural_points = build_structural_points(records)
    depletion_curves = build_depletion_curves(records, args.max_postjam_divisions)
    measured_arrest_selections = build_arrest_selections(
        records,
        args.max_postjam_divisions,
        active_cutoff=SECONDARY_ARREST_ACTIVE_CUTOFF,
    )
    paired_arrest_selections = build_arrest_selections(
        records,
        args.max_postjam_divisions,
        active_cutoff=SECONDARY_ARREST_ACTIVE_CUTOFF,
        require_phi2_shear=True,
    )
    phi2_summary = build_phi2_summary(measured_arrest_selections, active_cutoff=SECONDARY_ARREST_ACTIVE_CUTOFF)
    mechanics_points = build_mechanics_points(records)
    paired_arrest_points = build_paired_arrest_points(paired_arrest_selections)
    apply_g2_tukey_filter(paired_arrest_points)
    g2_summary = build_g2_summary(paired_arrest_points)
    g2_over_p2_summary = build_paired_g2_over_p2_summary(paired_arrest_points)
    g2_vs_p2_summary = build_paired_g2_vs_p2_summary(paired_arrest_points)
    lambda_c_summary = build_lambda_c_summary(measured_arrest_selections)
    hgamma_summary = build_growth_hgamma_summary(records, GROWTH_CROSSOVER_DPHI)
    growth_histograms = build_growth_histograms(records, GROWTH_CROSSOVER_DPHI)

    make_structure_figure(output_dir, structural_points, depletion_curves, phi2_summary, args.dpi)
    make_mechanics_figure(output_dir, mechanics_points, phi2_summary, g2_summary, paired_arrest_points, lambda_c_summary, args.dpi)
    make_growth_crossover_figure(output_dir, hgamma_summary, growth_histograms, args.dpi)
    make_paired_arrest_figure(output_dir, paired_arrest_points, g2_vs_p2_summary, g2_over_p2_summary, args.dpi)

    if seeds is None:
        print(f"Loaded records: {len(records)}")
    else:
        print(f"Loaded records: {len(records)} for seeds {args.seed_start}-{args.seed_stop}")
    print(f"Wrote figures to: {output_dir}")


if __name__ == "__main__":
    main()
