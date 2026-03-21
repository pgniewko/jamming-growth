#!/usr/bin/env python3

import argparse
import gzip
import os
import signal
import subprocess
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = REPO_ROOT / "output"
GROWTH_SCRIPT = REPO_ROOT / "examples" / "run_jamming.sh"
SHEAR_SCRIPT = REPO_ROOT / "examples" / "run_shear.sh"

SIZES = ["8", "15", "20", "25"]
P0S = ["-1", "1e-1", "1e-2", "1e-3", "1e-4", "1e-5"]
DPHIS = ["1e-4", "3e-4", "1e-3", "3e-3", "1e-2", "3e-2", "6e-2", "1e-1", "1.2e-1"]
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

GROWTH_STEPLOG_HEADER = (
    "# step N fret_per_particle P dt phi total_growthrate "
    "before_jamming at_jamming above_jamming"
)
SHEAR_G_DATA_HEADER = (
    "# strain shear_stress delta_shear_stress N Nc Nf Nu Ziso phi Nmm Nbb Nmb"
)
DEFAULT_JOB_TIMEOUT_SECONDS = 10800
EXIT_MIN_DT = 10
EXIT_MAX_POSTJAM_STEPS = 11


class JobTimeoutError(Exception):
    pass


class OutputValidationError(Exception):
    def __init__(self, stage, message):
        super().__init__(message)
        self.stage = stage


def dphi_allowed(p0, dphi):
    p0_value = float(p0)
    dphi_value = float(dphi)
    if p0_value > 1e-3 or p0_value <= 0.0:
        return True
    if p0_value == 1e-3:
        return dphi_value <= 6e-2
    if p0_value == 1e-4:
        return dphi_value <= 6e-2
    if p0_value == 1e-5:
        return dphi_value <= 3e-2
    return True


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
    parser.add_argument(
        "--keep-all-output",
        action="store_true",
        help="Keep shear trajectory files after successful jobs.",
    )
    parser.add_argument(
        "--job-timeout-seconds",
        type=int,
        default=DEFAULT_JOB_TIMEOUT_SECONDS,
        help=(
            "Kill a growth job if it runs longer than this many seconds. "
            "Use 0 to disable. Default: "
            f"{DEFAULT_JOB_TIMEOUT_SECONDS}"
        ),
    )
    args = parser.parse_args()
    if args.job_timeout_seconds < 0:
        raise SystemExit("--job-timeout-seconds must be non-negative")
    return args


def all_job_params():
    for size in SIZES:
        for p0 in P0S:
            for dphi in DPHIS:
                for seed in SEEDS:
                    yield {"lx": size, "p0": p0, "dphi": dphi, "seed": seed}


def job_params():
    for params in all_job_params():
        if dphi_allowed(params["p0"], params["dphi"]):
            yield params


def basename(params):
    return (
        f"v{FIXED['version']}_ar{FIXED['ar']}_div_{FIXED['divtype']}_desync{FIXED['desync']}"
        f"_seed_{params['seed']}_Lx{params['lx']}_Ly{params['lx']}"
        f"_att{FIXED['att']}_dphi{params['dphi']}_P{params['p0']}.dat"
    )


def open_text(path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def valid_file_with_gz(path, validator):
    candidates = (path, path.with_name(path.name + ".gz"))
    for candidate in candidates:
        if not candidate.is_file() or candidate.stat().st_size == 0:
            continue
        try:
            validator(candidate)
        except (OSError, ValueError):
            continue
        return candidate
    return None


def valid_gzip_file(path, validator):
    candidate = path.with_name(path.name + ".gz")
    if not candidate.is_file() or candidate.stat().st_size == 0:
        return None
    try:
        validator(candidate)
    except (OSError, ValueError):
        return None
    return candidate


def parse_cell_packing(path):
    with open_text(path) as handle:
        header = handle.readline().split()
        if len(header) != 2:
            raise ValueError("missing or invalid packing header")
        count = int(header[0])
        phi = float(header[1])
        if count <= 0:
            raise ValueError("packing must contain at least one cell")
        rows = 0
        for line in handle:
            text = line.strip()
            if not text:
                continue
            fields = text.split()
            if len(fields) != 5:
                raise ValueError(f"expected 5 columns per cell row, found {len(fields)}")
            for field in fields:
                float(field)
            rows += 1
    if rows != count:
        raise ValueError(f"header count {count} does not match {rows} cell rows")
    return {"count": count, "phi": phi}


def parse_lobe_packing(path):
    with open_text(path) as handle:
        frame_count = 0
        while True:
            header_line = handle.readline()
            if header_line == "":
                break
            if not header_line.strip():
                continue
            header = header_line.split()
            if len(header) != 2:
                raise ValueError("missing or invalid lobe-packing header")
            count = int(header[0])
            float(header[1])
            if count <= 0:
                raise ValueError("lobe packing must contain at least one particle")
            for _ in range(count):
                row = handle.readline()
                if row == "":
                    raise ValueError("truncated lobe-packing frame")
                fields = row.split()
                if len(fields) != 4:
                    raise ValueError(f"expected 4 columns per lobe row, found {len(fields)}")
                float(fields[0])
                float(fields[1])
                float(fields[2])
                int(fields[3])
            frame_count += 1
    if frame_count == 0:
        raise ValueError("lobe packing contains no frames")
    return frame_count


def parse_growth_stats(path, expected_rows):
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
            text = line.strip()
            if not text:
                continue
            fields = text.split()
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
        if header != GROWTH_STEPLOG_HEADER:
            raise ValueError("unexpected STEPLOG header")
        rows = 0
        last_count = None
        saw_postjam = False
        for line in handle:
            text = line.strip()
            if not text:
                continue
            if text.startswith("#"):
                continue
            fields = text.split()
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


def parse_stdout_log(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty log file: {path.name}")


def parse_shear_g_data(path):
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip()
        if header != SHEAR_G_DATA_HEADER:
            raise ValueError("unexpected G_data header")
        rows = 0
        first_n = None
        previous_strain = None
        for line in handle:
            text = line.strip()
            if not text:
                continue
            fields = text.split()
            if len(fields) != 12:
                raise ValueError(f"expected 12 columns per G_data row, found {len(fields)}")
            strain = float(fields[0])
            float(fields[1])
            float(fields[2])
            count = int(fields[3])
            int(fields[4])
            int(fields[5])
            int(fields[6])
            int(fields[7])
            float(fields[8])
            int(fields[9])
            int(fields[10])
            int(fields[11])
            if count <= 0:
                raise ValueError("invalid particle count in G_data")
            if previous_strain is not None and strain < previous_strain:
                raise ValueError("non-monotonic strain sequence in G_data")
            previous_strain = strain
            if first_n is None:
                first_n = count
            elif count != first_n:
                raise ValueError("particle count changes across G_data rows")
            rows += 1
    expected_rows = int(FIXED["shear_steps"])
    if rows != expected_rows:
        raise ValueError(f"G_data has {rows} rows, expected {expected_rows}")
    return first_n


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


def clean_failed_growth(paths):
    for path in (
        paths["growth_frame"],
        paths["growth_frame_gz"],
        paths["growth_jamm"],
        paths["growth_jamm_gz"],
        paths["growth_stats"],
        paths["growth_stats_jamm"],
        paths["growth_nc"],
        paths["growth_traj"],
        paths["growth_traj_gz"],
    ):
        remove_if_exists(path)


def clean_failed_shear(paths):
    for path in (
        paths["shear_input_local"],
        paths["shear_input_local_gz"],
        paths["shear_g_data"],
        paths["shear_traj"],
        paths["shear_traj_gz"],
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


def clean_finished_shear(paths, keep_all_output):
    if keep_all_output:
        return
    remove_if_exists(paths["shear_input_local"])
    remove_if_exists(paths["shear_input_local_gz"])
    remove_if_exists(paths["shear_traj"])
    remove_if_exists(paths["shear_traj_gz"])


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
    try:
        frame_path = valid_gzip_file(paths["growth_frame"], parse_cell_packing)
        jamm_path = valid_gzip_file(paths["growth_jamm"], parse_cell_packing)
        traj_path = valid_gzip_file(paths["growth_traj"], parse_lobe_packing)
        if frame_path is None or jamm_path is None or traj_path is None:
            return False
        frame_info = parse_cell_packing(frame_path)
        jamm_info = parse_cell_packing(jamm_path)
        parse_lobe_packing(traj_path)
        parse_growth_stats(paths["growth_stats"], frame_info["count"] * 2)
        parse_growth_stats(paths["growth_stats_jamm"], jamm_info["count"] * 2)
        parse_nc(paths["growth_nc"], jamm_info["count"], frame_info["count"])
        parse_steplog(paths["growth_steplog"], frame_info["count"])
        parse_stdout_log(paths["growth_log"])
        return True
    except (FileNotFoundError, OSError, ValueError):
        return False


def shear_done(paths):
    try:
        expected_count = parse_shear_g_data(paths["shear_g_data"])
        parse_stdout_log(paths["shear_log"])
        for candidate in (paths["shear_traj"], paths["shear_traj_gz"]):
            if candidate.is_file():
                frame_count = parse_lobe_packing(candidate)
                if frame_count != int(FIXED["shear_steps"]):
                    raise ValueError(
                        f"shear trajectory has {frame_count} frames, expected {FIXED['shear_steps']}"
                    )
                break
        return expected_count > 0
    except (FileNotFoundError, OSError, ValueError):
        return False


def run_script(script, extra_args, timeout_seconds=None):
    cmd = ["bash", str(script), "--results-root", str(OUTPUT_ROOT), *extra_args]
    process = subprocess.Popen(cmd, cwd=REPO_ROOT, start_new_session=True)
    try:
        returncode = process.wait(timeout=timeout_seconds)
    except subprocess.TimeoutExpired as exc:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
        raise JobTimeoutError(f"timed out after {timeout_seconds} seconds") from exc
    if returncode != 0:
        raise subprocess.CalledProcessError(returncode, cmd)


def classify_growth_failure(exc):
    if isinstance(exc, JobTimeoutError):
        return "TIMEOUT"
    if isinstance(exc, OutputValidationError) and exc.stage == "growth":
        return "INVALID_OUTPUTS"
    if isinstance(exc, subprocess.CalledProcessError):
        if exc.returncode == EXIT_MIN_DT:
            return "MIN_DT"
        if exc.returncode == EXIT_MAX_POSTJAM_STEPS:
            return "MAX_POSTJAM_STEPS"
    return None


def classify_shear_failure(exc):
    if isinstance(exc, OutputValidationError) and exc.stage == "shear":
        return "INVALID_OUTPUTS"
    return None


def run_growth(params, force, timeout_seconds):
    name = basename(params)
    paths = build_paths(name)
    if force:
        clean_growth(paths)
    elif growth_done(paths):
        return {"status": "completed", "name": name, "paths": paths}
    else:
        clean_growth(paths)

    try:
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
            timeout_seconds=timeout_seconds or None,
        )
        if not growth_done(paths):
            raise OutputValidationError("growth", f"incomplete growth outputs for {name}")
    except Exception as exc:
        reason = classify_growth_failure(exc)
        if reason is not None:
            clean_failed_growth(paths)
            return {
                "status": "failed",
                "stage": "growth",
                "reason": reason,
                "name": name,
                "paths": paths,
            }
        clean_growth(paths)
        raise
    return {"status": "completed", "name": name, "paths": paths}


def run_shear(params, paths, force, keep_all_output):
    if force:
        clean_shear(paths)
    elif shear_done(paths):
        clean_finished_shear(paths, keep_all_output)
        return {"status": "completed", "name": basename(params), "paths": paths}
    else:
        clean_shear(paths)

    try:
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
        if not shear_done(paths):
            raise OutputValidationError("shear", f"incomplete shear outputs for {basename(params)}")
    except Exception as exc:
        reason = classify_shear_failure(exc)
        if reason is not None:
            clean_failed_shear(paths)
            return {
                "status": "failed",
                "stage": "shear",
                "reason": reason,
                "name": basename(params),
                "paths": paths,
            }
        clean_shear(paths)
        raise
    clean_finished_shear(paths, keep_all_output)
    return {"status": "completed", "name": basename(params), "paths": paths}


def run_job(params, force, keep_all_output, timeout_seconds):
    growth_result = run_growth(params, force, timeout_seconds)
    if growth_result["status"] == "failed":
        return growth_result
    shear_result = run_shear(params, growth_result["paths"], force, keep_all_output)
    if shear_result["status"] == "failed":
        return shear_result
    return {"status": "completed", "name": growth_result["name"]}


def main():
    args = parse_args()
    if args.n_cpus < 1:
        raise SystemExit("--n-cpus must be at least 1")
    if not GROWTH_SCRIPT.is_file() or not SHEAR_SCRIPT.is_file():
        raise SystemExit("Missing example wrapper scripts.")

    jobs = list(job_params())
    futures = {}
    completed_jobs = 0
    failed_jobs = []
    with ThreadPoolExecutor(max_workers=args.n_cpus) as executor:
        for params in jobs:
            future = executor.submit(
                run_job,
                params,
                args.force,
                args.keep_all_output,
                args.job_timeout_seconds,
            )
            futures[future] = params
        while futures:
            done, _ = wait(futures, return_when=FIRST_COMPLETED)
            for future in done:
                params = futures.pop(future)
                try:
                    result = future.result()
                    if result["status"] == "failed":
                        failed_jobs.append((params, result["stage"], result["reason"]))
                        print(
                            f"Failed {result['stage']} job: "
                            f"lx={params['lx']} p0={params['p0']} "
                            f"dphi={params['dphi']} seed={params['seed']} "
                            f"reason={result['reason']}"
                        )
                    else:
                        completed_jobs += 1
                except Exception as exc:
                    for pending in futures:
                        pending.cancel()
                    raise SystemExit(
                        "Failed job: "
                        f"lx={params['lx']} p0={params['p0']} dphi={params['dphi']} seed={params['seed']}\n{exc}"
                    ) from exc
    print(f"Completed jobs: {completed_jobs}")
    print(f"Failed jobs: {len(failed_jobs)}")


if __name__ == "__main__":
    main()
