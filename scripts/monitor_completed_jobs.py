#!/usr/bin/env python3

from collections import Counter
from datetime import datetime

from run_growth_shear import (
    all_job_params,
    basename,
    build_paths,
    dphi_allowed,
    growth_done,
    job_params,
    shear_done,
)


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


def pair_is_allowed(left, right, left_value, right_value):
    if left == "p0" and right == "dphi":
        return dphi_allowed(left_value, right_value)
    if left == "dphi" and right == "p0":
        return dphi_allowed(right_value, left_value)
    return True


def collect_status():
    all_jobs = list(all_job_params())
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
        "all_jobs": all_jobs,
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
            if not pair_is_allowed(left, right, left_value, right_value):
                row.append("n/a")
                continue
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
            if pair_is_allowed(left, right, left_value, right_value)
        )
        column_expected = sum(
            status["expected_pair"][(left, right)][(left_value, right_value)]
            for left_value in left_values
            if pair_is_allowed(left, right, left_value, right_value)
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
    total_grid_jobs = len(status["all_jobs"])
    total_jobs = len(status["jobs"])
    completed_jobs = status["completed_jobs"]
    total_grid_points = len(
        {(params["lx"], params["p0"], params["dphi"]) for params in status["all_jobs"]}
    )
    active_grid_points = len(
        {(params["lx"], params["p0"], params["dphi"]) for params in status["jobs"]}
    )

    sections = [
        f"Completed job monitor snapshot\nTimestamp: {timestamp}",
        (
            "Summary\n"
            f"Completed growth+shear jobs: {completed_jobs}/{total_jobs} "
            f"({format_percent(completed_jobs, total_jobs).strip()})\n"
            f"Scheduled jobs after dphi filter: {total_jobs}/{total_grid_jobs}\n"
            f"Active grid points (L, P0, dphi): {active_grid_points}/{total_grid_points}"
        ),
        render_single_dimension_tables(status),
    ]

    for left, right, label in PAIRWISE_DIMENSIONS:
        sections.append(render_pairwise_table(left, right, label, status))

    return "\n\n".join(sections) + "\n"

def main():
    status = collect_status()
    report_text = render_report(status)
    print(report_text, flush=True)


if __name__ == "__main__":
    main()
