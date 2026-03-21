#!/usr/bin/env python3

import argparse
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

from run_growth_shear import basename, build_paths, growth_done, job_params, shear_done


DIMENSIONS = (
    ("lx", "L (size)"),
    ("p0", "P0 (feedback pressure)"),
    ("dphi", "dphi"),
)

PAIRWISE_DIMENSIONS = (
    ("lx", "p0", "L x P0"),
    ("lx", "dphi", "L x dphi"),
    ("p0", "dphi", "P0 x dphi"),
)


def parse_args():
    parser = argparse.ArgumentParser(description="Monitor completed growth+shear jobs.")
    parser.add_argument(
        "--interval-seconds",
        type=int,
        default=300,
        help="Polling interval in seconds. Default: 300",
    )
    parser.add_argument(
        "--write-reports",
        action="store_true",
        help="Also write rendered reports to disk.",
    )
    parser.add_argument(
        "--latest-report",
        type=Path,
        help="Path for the latest rendered report.",
    )
    parser.add_argument(
        "--history-dir",
        type=Path,
        help="Directory for timestamped report snapshots.",
    )
    parser.add_argument(
        "--once",
        action="store_true",
        help="Render one report and exit.",
    )
    return parser.parse_args()


def format_percent(completed, total):
    if total == 0:
        return "0.0%"
    return f"{(completed / total) * 100:5.1f}%"


def render_table(headers, rows):
    widths = [len(str(header)) for header in headers]
    for row in rows:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], len(str(cell)))

    def render_row(row):
        return " | ".join(str(cell).ljust(widths[index]) for index, cell in enumerate(row))

    separator = "-+-".join("-" * width for width in widths)
    lines = [render_row(headers), separator]
    lines.extend(render_row(row) for row in rows)
    return "\n".join(lines)


def collect_status():
    jobs = list(job_params())
    expected_single = {key: Counter() for key, _label in DIMENSIONS}
    completed_single = {key: Counter() for key, _label in DIMENSIONS}
    expected_pair = {(left, right): Counter() for left, right, _label in PAIRWISE_DIMENSIONS}
    completed_pair = {(left, right): Counter() for left, right, _label in PAIRWISE_DIMENSIONS}

    completed_jobs = 0
    for params in jobs:
        for key, _label in DIMENSIONS:
            expected_single[key][params[key]] += 1
        for left, right, _label in PAIRWISE_DIMENSIONS:
            expected_pair[(left, right)][(params[left], params[right])] += 1

        paths = build_paths(basename(params))
        if growth_done(paths) and shear_done(paths):
            completed_jobs += 1
            for key, _label in DIMENSIONS:
                completed_single[key][params[key]] += 1
            for left, right, _label in PAIRWISE_DIMENSIONS:
                completed_pair[(left, right)][(params[left], params[right])] += 1

    return {
        "jobs": jobs,
        "completed_jobs": completed_jobs,
        "expected_single": expected_single,
        "completed_single": completed_single,
        "expected_pair": expected_pair,
        "completed_pair": completed_pair,
    }

def render_single_dimension_tables(status):
    sections = []
    for key, label in DIMENSIONS:
        rows = []
        for value in sorted(status["expected_single"][key], key=_sort_key):
            expected = status["expected_single"][key][value]
            completed = status["completed_single"][key][value]
            rows.append((value, completed, expected, format_percent(completed, expected)))
        sections.append(
            f"{label}\n{render_table((label, 'Completed', 'Total', 'Percent'), rows)}"
        )
    return "\n\n".join(sections)


def render_pairwise_table(left, right, label, status):
    left_values = sorted(status["expected_single"][left], key=_sort_key)
    right_values = sorted(status["expected_single"][right], key=_sort_key)
    headers = [f"{label}"] + list(right_values) + ["Row total"]

    rows = []
    for left_value in left_values:
        row = [left_value]
        row_completed = 0
        row_expected = 0
        for right_value in right_values:
            expected = status["expected_pair"][(left, right)][(left_value, right_value)]
            completed = status["completed_pair"][(left, right)][(left_value, right_value)]
            row.append(f"{completed}/{expected}")
            row_completed += completed
            row_expected += expected
        row.append(f"{row_completed}/{row_expected}")
        rows.append(tuple(row))

    total_row = ["Column total"]
    total_completed = 0
    total_expected = 0
    for right_value in right_values:
        column_completed = sum(
            status["completed_pair"][(left, right)][(left_value, right_value)]
            for left_value in left_values
        )
        column_expected = sum(
            status["expected_pair"][(left, right)][(left_value, right_value)]
            for left_value in left_values
        )
        total_row.append(f"{column_completed}/{column_expected}")
        total_completed += column_completed
        total_expected += column_expected
    total_row.append(f"{total_completed}/{total_expected}")
    rows.append(tuple(total_row))

    return f"{label}\n{render_table(headers, rows)}"


def _sort_key(value):
    try:
        return (0, float(value))
    except ValueError:
        return (1, value)


def render_report(status):
    timestamp = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    total_jobs = len(status["jobs"])
    completed_jobs = status["completed_jobs"]

    sections = [
        f"Completed job monitor snapshot\nTimestamp: {timestamp}",
        (
            "Summary\n"
            f"Completed growth+shear jobs: {completed_jobs}/{total_jobs} "
            f"({format_percent(completed_jobs, total_jobs).strip()})"
        ),
        render_single_dimension_tables(status),
    ]

    for left, right, label in PAIRWISE_DIMENSIONS:
        sections.append(render_pairwise_table(left, right, label, status))

    return "\n\n".join(sections) + "\n"


def write_report(report_text, latest_report, history_dir):
    if latest_report is None or history_dir is None:
        raise SystemExit("--latest-report and --history-dir are required with --write-reports")
    latest_report.parent.mkdir(parents=True, exist_ok=True)
    history_dir.mkdir(parents=True, exist_ok=True)
    latest_report.write_text(report_text)
    snapshot_name = datetime.now().astimezone().strftime("%Y%m%dT%H%M%S%z.txt")
    (history_dir / snapshot_name).write_text(report_text)


def main():
    args = parse_args()

    while True:
        status = collect_status()
        report_text = render_report(status)
        print(report_text, flush=True)
        if args.write_reports:
            write_report(report_text, args.latest_report, args.history_dir)
        if args.once:
            return
        time.sleep(args.interval_seconds)


if __name__ == "__main__":
    main()
