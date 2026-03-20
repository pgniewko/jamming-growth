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
    parser.add_argument("--sizes", required=True, type=parse_csv_list, help="Comma-separated box sizes.")
    parser.add_argument("--p0s", required=True, type=parse_csv_list, help="Comma-separated feedback strengths.")
    parser.add_argument("--dphis", required=True, type=parse_csv_list, help="Comma-separated overcompressions.")
    parser.add_argument("--seeds", required=True, type=parse_csv_list, help="Comma-separated positive seeds.")
    parser.add_argument(
        "--n-cpus",
        type=int,
        default=default_cpus,
        help=f"Number of concurrent jobs. Default: {default_cpus}",
    )
    parser.add_argument("--force", action="store_true", help="Rerun completed jobs.")
    return parser.parse_args()


def job_params(args):
    for seed in args.seeds:
        for size in args.sizes:
            for p0 in args.p0s:
                for dphi in args.dphis:
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
            return None
        count = int(header[0])
        float(header[1])
        rows = 0
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 5:
                return None
            for field in fields:
                float(field)
            rows += 1
    return rows if rows == count else None


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
        paths["nc"],
        paths["steplog"],
        paths["traj"],
        paths["traj"].with_name(paths["traj"].name + ".gz"),
        paths["stdout_log"],
    ):
        remove_if_exists(path)


def job_complete(paths):
    frame_path = file_with_gz(paths["frame"])
    jamm_path = file_with_gz(paths["jamm"])
    if frame_path is None or jamm_path is None:
        return False
    try:
        frame_count = parse_packing(frame_path)
        jamm_count = parse_packing(jamm_path)
        if frame_count is None or jamm_count is None:
            return False
        if not parse_lineage(paths["lineage_frame"], frame_count):
            return False
        if not parse_lineage(paths["lineage_jamm"], jamm_count):
            return False
        if not parse_divlog(paths["divlog"]):
            return False
    except (FileNotFoundError, OSError, ValueError):
        return False
    if not parse_log(paths["stdout_log"]):
        return False
    if not parse_log(paths["nc"]) or not parse_log(paths["steplog"]):
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


def run_job(params, force):
    name = basename(params)
    paths = build_paths(name)
    if force:
        clean_job(paths)
    if not job_complete(paths):
        clean_job(paths)
        run_growth(params)
        if not job_complete(paths):
            raise RuntimeError(f"incomplete lineage-growth outputs for {name}")
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
            future = executor.submit(run_job, params, args.force)
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
