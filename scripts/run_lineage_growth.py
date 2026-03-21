#!/usr/bin/env python3

import argparse
import gzip
import os
import subprocess
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = REPO_ROOT / "output"
GROWTH_ROOT = OUTPUT_ROOT / "lineage_growth"
GROWTH_DIR = GROWTH_ROOT / "growth"
LOG_DIR = OUTPUT_ROOT / "logs" / "lineage_growth"
GROWTH_SCRIPT = REPO_ROOT / "examples" / "run_jamming.sh"
LINEAGE_EXE = REPO_ROOT / "bin" / "jamming_by_growth_lineage"

FIXED = {
    "att": "0.0",
    "rate0": "0.002",
    "skip": "25",
    "desync": "0.4",
    "ar": "1.01",
    "divtype": "4",
    "version": "1.0",
}

PRESETS = {
    "pilot-large": {
        "sizes": ["15", "20"],
        "p0s": ["2e-3", "1e-2"],
        "dphis": ["1.2e-1"],
        "seeds": ["1201", "1202"],
        "n_cpus": 4,
    },
    "full-large": {
        "sizes": ["15", "20"],
        "p0s": ["-1", "2e-3", "5e-3", "1e-2"],
        "dphis": ["1.2e-1"],
        "seeds": ["1201", "1202", "1203", "1204", "1205", "1206"],
        "n_cpus": 6,
    },
}

EVENT_CODES = {
    "first_bud_contact": 1,
    "free_to_constrained": 2,
    "initial_free_division": 3,
}
STEPLOG_HEADER = (
    "# step N fret_per_particle P dt phi total_growthrate "
    "before_jamming at_jamming above_jamming"
)


def parse_csv_list(value):
    items = [item.strip() for item in value.split(",")]
    items = [item for item in items if item]
    if not items:
        raise argparse.ArgumentTypeError("expected a non-empty comma-separated list")
    return items


def parse_args():
    default_cpus = max(1, (os.cpu_count() or 1) - 1)
    parser = argparse.ArgumentParser(
        description="Run a compact growth-only lineage stage for an explicit parameter set."
    )
    parser.add_argument(
        "--preset",
        choices=sorted(PRESETS),
        help="Named default parameter set. Explicit --sizes/--p0s/--dphis/--seeds override preset values.",
    )
    parser.add_argument("--sizes", type=parse_csv_list, help="Comma-separated box sizes.")
    parser.add_argument("--p0s", type=parse_csv_list, help="Comma-separated feedback strengths.")
    parser.add_argument("--dphis", type=parse_csv_list, help="Comma-separated overcompressions.")
    parser.add_argument("--seeds", type=parse_csv_list, help="Comma-separated positive seeds.")
    parser.add_argument(
        "--n-cpus",
        type=int,
        default=None,
        help=f"Number of concurrent jobs. Default: preset value or {default_cpus}",
    )
    parser.add_argument("--force", action="store_true", help="Rerun completed jobs.")
    parser.add_argument(
        "--keep-all-output",
        action="store_true",
        help="Keep transient growth-style outputs and the pre-jamming lineage snapshot after successful jobs.",
    )
    args = parser.parse_args()
    preset = PRESETS.get(args.preset, {})
    for field in ("sizes", "p0s", "dphis", "seeds"):
        if getattr(args, field) is None:
            value = preset.get(field)
            if value is None:
                parser.error(
                    "either provide --preset or explicitly pass all of "
                    "--sizes/--p0s/--dphis/--seeds"
                )
            setattr(args, field, list(value))
    if args.n_cpus is None:
        args.n_cpus = int(preset.get("n_cpus", default_cpus))
    return args


def job_params(args):
    for size in args.sizes:
        for p0 in args.p0s:
            for dphi in args.dphis:
                for seed in args.seeds:
                    yield {"lx": size, "p0": p0, "dphi": dphi, "seed": seed}


def basename(params):
    return (
        f"v{FIXED['version']}_ar{FIXED['ar']}_div_{FIXED['divtype']}_desync{FIXED['desync']}"
        f"_seed_{params['seed']}_Lx{params['lx']}_Ly{params['lx']}"
        f"_att{FIXED['att']}_dphi{params['dphi']}_P{params['p0']}.dat"
    )


def file_with_gz(path):
    gz_path = path.with_name(path.name + ".gz")
    if path.is_file() and path.stat().st_size > 0:
        return path
    if gz_path.is_file() and gz_path.stat().st_size > 0:
        return gz_path
    return None


def remove_if_exists(path):
    if path.exists():
        path.unlink()


def open_text(path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_packing(path):
    with open_text(path) as handle:
        header = handle.readline().split()
        if len(header) != 2:
            raise ValueError("missing or invalid packing header")
        count = int(header[0])
        float(header[1])
        if count <= 0:
            raise ValueError("packing must contain at least one cell")
        rows = 0
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 5:
                raise ValueError(f"expected 5 columns per packing row, found {len(fields)}")
            for field in fields:
                float(field)
            rows += 1
    if rows != count:
        raise ValueError(f"header count {count} does not match {rows} packing rows")
    return rows


def parse_trajectory(path):
    with open_text(path) as handle:
        frames = 0
        while True:
            header_line = handle.readline()
            if header_line == "":
                break
            if not header_line.strip():
                continue
            header = header_line.split()
            if len(header) != 2:
                raise ValueError("missing or invalid trajectory header")
            count = int(header[0])
            float(header[1])
            if count <= 0:
                raise ValueError("trajectory frame must contain at least one lobe")
            for _ in range(count):
                row = handle.readline()
                if row == "":
                    raise ValueError("truncated trajectory frame")
                fields = row.split()
                if len(fields) != 4:
                    raise ValueError(f"expected 4 columns per trajectory row, found {len(fields)}")
                float(fields[0])
                float(fields[1])
                float(fields[2])
                int(fields[3])
            frames += 1
    if frames == 0:
        raise ValueError("trajectory contains no frames")
    return frames


def parse_stats(path, expected_rows):
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().split()
        if len(header) != 3:
            raise ValueError("missing or invalid stats header")
        row_count = int(header[0])
        float(header[1])
        float(header[2])
        if row_count != expected_rows:
            raise ValueError(f"stats header count {row_count} does not match expected {expected_rows}")
        rows = 0
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 8:
                raise ValueError(f"expected 8 columns per stats row, found {len(fields)}")
            float(fields[0])
            float(fields[1])
            float(fields[2])
            int(fields[3])
            float(fields[4])
            float(fields[5])
            float(fields[6])
            float(fields[7])
            rows += 1
    if rows != expected_rows:
        raise ValueError(f"stats file has {rows} rows, expected {expected_rows}")


def parse_nc(path, expected_jamm_count, expected_frame_count):
    fields = path.read_text(encoding="utf-8").split()
    if len(fields) != 26:
        raise ValueError(f"expected 26 NC values, found {len(fields)}")
    values = list(map(float, fields))
    for offset, expected_count in ((0, expected_jamm_count), (13, expected_frame_count)):
        count = int(values[offset])
        if count != expected_count:
            raise ValueError(f"NC count {count} does not match expected {expected_count}")
        for index in range(offset, offset + 8):
            int(values[index])
        for index in range(offset + 8, offset + 13):
            float(values[index])


def parse_steplog(path, expected_final_count):
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip()
        if header != STEPLOG_HEADER:
            raise ValueError("unexpected STEPLOG header")
        rows = 0
        last_count = None
        saw_postjam = False
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 10:
                raise ValueError(f"expected 10 columns per steplog row, found {len(fields)}")
            int(fields[0])
            last_count = int(fields[1])
            float(fields[2])
            float(fields[3])
            float(fields[4])
            float(fields[5])
            float(fields[6])
            before_jamming = int(fields[7])
            at_jamming = int(fields[8])
            above_jamming = int(fields[9])
            if before_jamming not in (0, 1) or at_jamming not in (0, 1) or above_jamming not in (0, 1):
                raise ValueError("invalid jamming-state flags in steplog")
            saw_postjam = saw_postjam or above_jamming == 1
            rows += 1
    if rows == 0:
        raise ValueError("STEPLOG contains no data rows")
    if last_count != expected_final_count:
        raise ValueError(f"STEPLOG final N {last_count} does not match expected {expected_final_count}")
    if not saw_postjam:
        raise ValueError("STEPLOG never reaches the post-jamming regime")


def parse_lineage_snapshot(path):
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header.startswith("# index cell_id parent_id"):
            raise ValueError("unexpected lineage snapshot header")
        rows = 0
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 8:
                raise ValueError(f"expected 8 columns per lineage row, found {len(fields)}")
            int(fields[0])
            int(fields[1])
            int(fields[2])
            float(fields[3])
            float(fields[4])
            int(fields[5])
            int(fields[6])
            int(fields[7])
            rows += 1
    if rows == 0:
        raise ValueError("lineage snapshot contains no rows")
    return rows


def parse_lineage(path, expected_count):
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header.startswith("# index cell_id parent_id"):
            return False
        rows = 0
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 8:
                return False
            int(fields[0])
            int(fields[1])
            int(fields[2])
            float(fields[3])
            float(fields[4])
            int(fields[5])
            int(fields[6])
            int(fields[7])
            rows += 1
    return rows == expected_count


def parse_divlog(path):
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header.startswith("# step phi parent_cell_id"):
            return False
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 6:
                return False
            int(fields[0])
            float(fields[1])
            int(fields[2])
            int(fields[3])
            int(fields[4])
            int(fields[5])
    return True


def parse_transitions(path):
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header.startswith("# step phi cell_id parent_id"):
            return False
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 12:
                return False
            int(fields[0])
            float(fields[1])
            int(fields[2])
            int(fields[3])
            float(fields[4])
            float(fields[5])
            float(fields[6])
            int(fields[7])
            int(fields[8])
            int(fields[9])
            code = int(fields[10])
            label = fields[11]
            if EVENT_CODES.get(label) != code:
                return False
    return True


def parse_postjamm_summary(path):
    rows = 0
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
        if not header.startswith("# step phi P N Nc Nf Nu Ziso total_growthrate chi_c"):
            return False
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 15:
                return False
            int(fields[0])
            float(fields[1])
            float(fields[2])
            int(fields[3])
            int(fields[4])
            int(fields[5])
            int(fields[6])
            int(fields[7])
            float(fields[8])
            float(fields[9])
            int(fields[10])
            int(fields[11])
            int(fields[12])
            int(fields[13])
            int(fields[14])
            rows += 1
    return rows > 0


def parse_log(path):
    return path.is_file() and path.stat().st_size > 0


def build_paths(name):
    return {
        "frame": GROWTH_DIR / f"LF_DPHI_{name}",
        "jamm": GROWTH_DIR / f"LF_JAMM_{name}",
        "stats_frame": GROWTH_DIR / f"STATS_LF_DPHI_{name}",
        "stats_jamm": GROWTH_DIR / f"STATS_LF_JAMM_{name}",
        "lineage_frame": GROWTH_DIR / f"LINEAGE_LF_DPHI_{name}",
        "lineage_jamm": GROWTH_DIR / f"LINEAGE_LF_JAMM_{name}",
        "divlog": GROWTH_DIR / f"DIVLOG_{name}",
        "transitions": GROWTH_DIR / f"TRANSITIONS_{name}",
        "postjamm_summary": GROWTH_DIR / f"POSTJAMM_SUMMARY_{name}",
        "nc": GROWTH_DIR / f"NC_{name}",
        "steplog": GROWTH_DIR / f"STEPLOG_{name}",
        "traj": GROWTH_DIR / name,
        "stdout_log": LOG_DIR / f"stdout_{name[:-4]}.log",
    }


def clean_job(paths):
    for path in (
        paths["frame"],
        paths["frame"].with_name(paths["frame"].name + ".gz"),
        paths["jamm"],
        paths["jamm"].with_name(paths["jamm"].name + ".gz"),
        paths["stats_frame"],
        paths["stats_jamm"],
        paths["lineage_frame"],
        paths["lineage_jamm"],
        paths["divlog"],
        paths["transitions"],
        paths["postjamm_summary"],
        paths["nc"],
        paths["steplog"],
        paths["traj"],
        paths["traj"].with_name(paths["traj"].name + ".gz"),
        paths["stdout_log"],
    ):
        remove_if_exists(path)


def clean_finished_job(paths, keep_all_output):
    if keep_all_output:
        return
    for path in (
        paths["frame"],
        paths["frame"].with_name(paths["frame"].name + ".gz"),
        paths["jamm"],
        paths["jamm"].with_name(paths["jamm"].name + ".gz"),
        paths["stats_frame"],
        paths["stats_jamm"],
        paths["lineage_frame"],
        paths["nc"],
        paths["steplog"],
        paths["traj"],
        paths["traj"].with_name(paths["traj"].name + ".gz"),
    ):
        remove_if_exists(path)


def lineage_done(paths):
    if not parse_log(paths["stdout_log"]):
        return False
    try:
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
        if paths["lineage_frame"].is_file():
            frame_count = parse_lineage_snapshot(paths["lineage_frame"])
            if not parse_lineage(paths["lineage_frame"], frame_count):
                return False

        for path, parser in (
            (paths["frame"], parse_packing),
            (paths["jamm"], parse_packing),
            (paths["traj"], parse_trajectory),
        ):
            gz_path = path.with_name(path.name + ".gz")
            if path.is_file():
                return False
            if gz_path.is_file():
                parser(gz_path)

        if paths["stats_frame"].is_file():
            if not paths["lineage_frame"].is_file():
                return False
            parse_stats(paths["stats_frame"], frame_count * 2)
        if paths["stats_jamm"].is_file():
            parse_stats(paths["stats_jamm"], jamm_count * 2)
        if paths["nc"].is_file():
            if not paths["lineage_frame"].is_file():
                return False
            parse_nc(paths["nc"], jamm_count, frame_count)
        if paths["steplog"].is_file():
            if not paths["lineage_frame"].is_file():
                return False
            parse_steplog(paths["steplog"], frame_count)
    except (FileNotFoundError, OSError, ValueError):
        return False
    return True


def run_growth(params):
    cmd = [
        "bash",
        str(GROWTH_SCRIPT),
        "--exe",
        str(LINEAGE_EXE),
        "--growth-dir",
        str(GROWTH_DIR),
        "--log-dir",
        str(LOG_DIR),
        "--p0",
        params["p0"],
        "--lx",
        params["lx"],
        "--dphi",
        params["dphi"],
        "--seed",
        params["seed"],
        "--att",
        FIXED["att"],
        "--rate0",
        FIXED["rate0"],
        "--skip",
        FIXED["skip"],
        "--desync",
        FIXED["desync"],
        "--ar",
        FIXED["ar"],
        "--divtype",
        FIXED["divtype"],
        "--version",
        FIXED["version"],
    ]
    subprocess.run(cmd, cwd=REPO_ROOT, check=True)


def run_job(params, force, keep_all_output):
    name = basename(params)
    paths = build_paths(name)
    if force:
        clean_job(paths)
    elif lineage_done(paths):
        clean_finished_job(paths, keep_all_output)
        return name
    else:
        clean_job(paths)

    try:
        run_growth(params)
        clean_finished_job(paths, keep_all_output)
        if not lineage_done(paths):
            raise RuntimeError(f"incomplete lineage-growth outputs for {name}")
    except Exception:
        clean_job(paths)
        raise
    return name


def main():
    args = parse_args()
    if args.n_cpus < 1:
        raise SystemExit("--n-cpus must be at least 1")
    if not GROWTH_SCRIPT.is_file():
        raise SystemExit("Missing growth wrapper script.")
    if not LINEAGE_EXE.is_file():
        raise SystemExit(f"Missing executable: {LINEAGE_EXE}. Run 'make all' first.")

    GROWTH_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    futures = {}
    with ThreadPoolExecutor(max_workers=args.n_cpus) as executor:
        for params in job_params(args):
            future = executor.submit(run_job, params, args.force, args.keep_all_output)
            futures[future] = params
        while futures:
            done, _ = wait(futures, return_when=FIRST_COMPLETED)
            for future in done:
                params = futures.pop(future)
                try:
                    future.result()
                except Exception as exc:
                    for pending in futures:
                        pending.cancel()
                    raise SystemExit(
                        "Failed job: "
                        f"lx={params['lx']} p0={params['p0']} "
                        f"dphi={params['dphi']} seed={params['seed']}\n{exc}"
                    ) from exc


if __name__ == "__main__":
    main()
