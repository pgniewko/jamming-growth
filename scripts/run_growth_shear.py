#!/usr/bin/env python3

import argparse
import os
import subprocess
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = REPO_ROOT / "output"
GROWTH_SCRIPT = REPO_ROOT / "examples" / "run_jamming.sh"
SHEAR_SCRIPT = REPO_ROOT / "examples" / "run_shear.sh"

SIZES = ["8", "15", "20", "25"]
P0S = ["-1", "1e-2", "5e-3", "2e-3", "1e-3", "5e-4", "1e-4", "1e-5"]
DPHIS = ["1e-4", "3e-4", "1e-3", "3e-3", "1e-2", "2e-2", "4e-2", "6e-2", "8e-2", "1e-1", "1.2e-1"]
SEEDS = [str(seed) for seed in range(1201, 1221)]

FIXED = {
    "att": "0.0",
    "rate0": "0.002",
    "skip": "25",
    "desync": "0.4",
    "ar": "1.01",
    "divtype": "4",
    "version": "1.0",
    "strain_step": "1e-6",
    "shear_steps": "5000",
}


def parse_args():
    default_cpus = max(1, (os.cpu_count() or 1) - 1)
    parser = argparse.ArgumentParser(description="Run the full growth + shear numerical sweep.")
    parser.add_argument(
        "--n-cpus",
        type=int,
        default=default_cpus,
        help=f"Number of concurrent jobs. Default: {default_cpus}",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rerun completed jobs.",
    )
    return parser.parse_args()


def job_params():
    for size in SIZES:
        for p0 in P0S:
            for dphi in DPHIS:
                for seed in SEEDS:
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


def clean_growth(paths):
    for path in (
        paths["growth_frame"],
        paths["growth_frame_gz"],
        paths["growth_jamm"],
        paths["growth_jamm_gz"],
        paths["growth_stats"],
        paths["growth_stats_jamm"],
        paths["growth_nc"],
        paths["growth_steplog"],
        paths["growth_traj"],
        paths["growth_traj_gz"],
        paths["growth_log"],
    ):
        remove_if_exists(path)


def clean_shear(paths):
    for path in (
        paths["shear_input_local"],
        paths["shear_input_local_gz"],
        paths["shear_g_data"],
        paths["shear_traj"],
        paths["shear_traj_gz"],
        paths["shear_log"],
    ):
        remove_if_exists(path)


def build_paths(name):
    growth_dir = OUTPUT_ROOT / "growth"
    shear_dir = OUTPUT_ROOT / "shear"
    growth_log_dir = OUTPUT_ROOT / "logs" / "growth"
    shear_log_dir = OUTPUT_ROOT / "logs" / "shear"
    return {
        "growth_frame": growth_dir / f"LF_DPHI_{name}",
        "growth_frame_gz": growth_dir / f"LF_DPHI_{name}.gz",
        "growth_jamm": growth_dir / f"LF_JAMM_{name}",
        "growth_jamm_gz": growth_dir / f"LF_JAMM_{name}.gz",
        "growth_stats": growth_dir / f"STATS_LF_DPHI_{name}",
        "growth_stats_jamm": growth_dir / f"STATS_LF_JAMM_{name}",
        "growth_nc": growth_dir / f"NC_{name}",
        "growth_steplog": growth_dir / f"STEPLOG_{name}",
        "growth_traj": growth_dir / name,
        "growth_traj_gz": growth_dir / f"{name}.gz",
        "growth_log": growth_log_dir / f"stdout_{name[:-4]}.log",
        "shear_input_local": shear_dir / f"LF_DPHI_{name}",
        "shear_input_local_gz": shear_dir / f"LF_DPHI_{name}.gz",
        "shear_g_data": shear_dir / f"G_data_LF_DPHI_{name}",
        "shear_traj": shear_dir / f"SHEAR_TRAJ_LF_DPHI_{name}",
        "shear_traj_gz": shear_dir / f"SHEAR_TRAJ_LF_DPHI_{name}.gz",
        "shear_log": shear_log_dir / f"stdout_{name[:-4]}.log",
    }


def growth_done(paths):
    return file_with_gz(paths["growth_frame"]) is not None


def shear_done(paths):
    return paths["shear_g_data"].is_file() and paths["shear_g_data"].stat().st_size > 0 and file_with_gz(paths["shear_traj"]) is not None


def run_script(script, extra_args):
    cmd = ["bash", str(script), "--results-root", str(OUTPUT_ROOT), *extra_args]
    subprocess.run(cmd, cwd=REPO_ROOT, check=True)


def run_growth(params, force):
    name = basename(params)
    paths = build_paths(name)
    if force:
        clean_growth(paths)
    if not growth_done(paths):
        clean_growth(paths)
        run_script(
            GROWTH_SCRIPT,
            [
                "--p0", params["p0"],
                "--lx", params["lx"],
                "--dphi", params["dphi"],
                "--seed", params["seed"],
                "--att", FIXED["att"],
                "--rate0", FIXED["rate0"],
                "--skip", FIXED["skip"],
                "--desync", FIXED["desync"],
                "--ar", FIXED["ar"],
                "--divtype", FIXED["divtype"],
                "--version", FIXED["version"],
            ],
        )
    return name, paths


def run_shear(params, paths, force):
    if force:
        clean_shear(paths)
    if not shear_done(paths):
        clean_shear(paths)
        run_script(
            SHEAR_SCRIPT,
            [
                "--p0", params["p0"],
                "--lx", params["lx"],
                "--dphi", params["dphi"],
                "--seed", params["seed"],
                "--att", FIXED["att"],
                "--desync", FIXED["desync"],
                "--ar", FIXED["ar"],
                "--divtype", FIXED["divtype"],
                "--version", FIXED["version"],
                "--strain-step", FIXED["strain_step"],
                "--shear-steps", FIXED["shear_steps"],
            ],
        )


def run_job(params, force):
    name, paths = run_growth(params, force)
    run_shear(params, paths, force)
    return name


def main():
    args = parse_args()
    if args.n_cpus < 1:
        raise SystemExit("--n-cpus must be at least 1")
    if not GROWTH_SCRIPT.is_file() or not SHEAR_SCRIPT.is_file():
        raise SystemExit("Missing example wrapper scripts.")

    jobs = list(job_params())
    futures = {}
    with ThreadPoolExecutor(max_workers=args.n_cpus) as executor:
        for params in jobs:
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
                        f"lx={params['lx']} p0={params['p0']} dphi={params['dphi']} seed={params['seed']}\n{exc}"
                    ) from exc


if __name__ == "__main__":
    main()
