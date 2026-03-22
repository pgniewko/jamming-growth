import gzip
import os
import shutil
import signal
import subprocess
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from threading import Lock

from pipeline_cleanup import clean_bext, clean_growth, clean_shear, normalize_bext_success, normalize_growth_success, normalize_shear_success
from pipeline_config import BEXT_DIR, BEXT_EXE, DEFAULT_DPHI_PROBE, DEFAULT_JOB_TIMEOUT_SECONDS, EXIT_MAX_POSTJAM_STEPS, EXIT_MIN_DT, FIXED, GROWTH_SCRIPT, OUTPUT_ROOT, REPO_ROOT, SHEAR_SCRIPT
from pipeline_paths import bext_paths, ensure_stage_dirs, growth_paths, shear_paths
from pipeline_validate import bext_done, growth_done, shear_done
from pipeline_config import basename


class JobTimeoutError(Exception):
    pass


class OutputValidationError(Exception):
    def __init__(self, stage, message):
        super().__init__(message)
        self.stage = stage


ACTIVE_PGIDS = set()
ACTIVE_PGIDS_LOCK = Lock()


def register_process(process):
    with ACTIVE_PGIDS_LOCK:
        ACTIVE_PGIDS.add(process.pid)


def unregister_process(process):
    with ACTIVE_PGIDS_LOCK:
        ACTIVE_PGIDS.discard(process.pid)


def terminate_process_group(pid, sig):
    try:
        os.killpg(pid, sig)
    except (PermissionError, ProcessLookupError):
        pass


def terminate_active_processes():
    with ACTIVE_PGIDS_LOCK:
        pids = tuple(ACTIVE_PGIDS)
    for pid in pids:
        terminate_process_group(pid, signal.SIGTERM)
    for pid in pids:
        terminate_process_group(pid, signal.SIGKILL)


def run_command(cmd, timeout_seconds=None, cwd=None, stdout_handle=None, input_bytes=None):
    process = subprocess.Popen(
        cmd,
        cwd=cwd or REPO_ROOT,
        start_new_session=True,
        stdin=subprocess.PIPE if input_bytes is not None else None,
        stdout=stdout_handle,
        stderr=subprocess.STDOUT if stdout_handle is not None else None,
    )
    register_process(process)
    try:
        if input_bytes is None:
            returncode = process.wait(timeout=timeout_seconds)
        else:
            process.communicate(input=input_bytes, timeout=timeout_seconds)
            returncode = process.returncode
    except subprocess.TimeoutExpired as exc:
        terminate_process_group(process.pid, signal.SIGTERM)
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            terminate_process_group(process.pid, signal.SIGKILL)
            process.wait()
        raise JobTimeoutError(f"timed out after {timeout_seconds} seconds") from exc
    finally:
        unregister_process(process)
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


def classify_derived_failure(exc, stage):
    if isinstance(exc, JobTimeoutError):
        return "TIMEOUT"
    if isinstance(exc, OutputValidationError) and exc.stage == stage:
        return "INVALID_OUTPUTS"
    return None


def run_growth_stage(params, force=False, save_all_data=False, timeout_seconds=DEFAULT_JOB_TIMEOUT_SECONDS):
    name = basename(params)
    paths = growth_paths(name)
    if force:
        clean_growth(paths)
    elif growth_done(paths):
        normalize_growth_success(paths, save_all_data)
        return {"status": "skipped", "name": name, "paths": paths}
    else:
        clean_growth(paths)

    cmd = [
        "bash",
        str(GROWTH_SCRIPT),
        "--results-root",
        str(OUTPUT_ROOT),
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
    if save_all_data:
        cmd.append("--save-all-data")
    try:
        run_command(cmd, timeout_seconds=timeout_seconds or None)
        normalize_growth_success(paths, save_all_data)
        if not growth_done(paths):
            raise OutputValidationError("growth", f"incomplete growth outputs for {name}")
    except Exception as exc:
        reason = classify_growth_failure(exc)
        if reason is not None:
            clean_growth(paths)
            return {"status": "failed", "stage": "growth", "reason": reason, "name": name, "paths": paths}
        clean_growth(paths)
        raise
    return {"status": "completed", "name": name, "paths": paths}


def run_shear_stage(params, force=False, save_all_data=False, timeout_seconds=None):
    name = basename(params)
    paths = shear_paths(name)
    if force:
        clean_shear(paths)
    elif shear_done(paths):
        normalize_shear_success(paths, save_all_data)
        return {"status": "skipped", "name": name, "paths": paths}
    else:
        clean_shear(paths)

    cmd = [
        "bash",
        str(SHEAR_SCRIPT),
        "--results-root",
        str(OUTPUT_ROOT),
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
        "--desync",
        FIXED["desync"],
        "--ar",
        FIXED["ar"],
        "--divtype",
        FIXED["divtype"],
        "--version",
        FIXED["version"],
        "--strain-step",
        FIXED["strain_step"],
        "--shear-steps",
        FIXED["shear_steps"],
    ]
    if save_all_data:
        cmd.append("--save-all-data")
    try:
        run_command(cmd, timeout_seconds=timeout_seconds or None)
        normalize_shear_success(paths, save_all_data)
        if not shear_done(paths):
            raise OutputValidationError("shear", f"incomplete shear outputs for {name}")
    except Exception as exc:
        reason = classify_derived_failure(exc, "shear")
        if reason is not None:
            clean_shear(paths)
            return {"status": "failed", "stage": "shear", "reason": reason, "name": name, "paths": paths}
        clean_shear(paths)
        raise
    return {"status": "completed", "name": name, "paths": paths}


def stage_bext_input(source, destination):
    if source.suffix == ".gz":
        with gzip.open(source, "rb") as src, destination.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copyfile(source, destination)


def growth_input_source(name):
    growth = growth_paths(name)
    if growth["frame_gz"].is_file():
        return growth["frame_gz"]
    if growth["frame"].is_file():
        return growth["frame"]
    raise FileNotFoundError(f"missing growth input for {name}")


def run_bext_stage(params, dphi_probe, force=False, save_all_data=False, timeout_seconds=None):
    name = basename(params)
    paths = bext_paths(name, dphi_probe)
    if force:
        clean_bext(paths)
    elif bext_done(paths):
        normalize_bext_success(paths, save_all_data)
        return {"status": "skipped", "name": name, "paths": paths}
    else:
        clean_bext(paths)

    source = growth_input_source(name)
    stage_bext_input(source, paths["input_local"])
    try:
        with paths["log"].open("wb") as log_handle:
            run_command(
                [str(BEXT_EXE)],
                cwd=BEXT_DIR,
                timeout_seconds=timeout_seconds or None,
                stdout_handle=log_handle,
                input_bytes="\n".join(
                    [params["lx"], params["lx"], FIXED["att"], paths["input_local"].name, f"{dphi_probe:.16g}", ""]
                ).encode(),
            )
        normalize_bext_success(paths, save_all_data)
        if not bext_done(paths):
            raise OutputValidationError("bext", f"incomplete B_ext outputs for {name}")
    except Exception as exc:
        reason = classify_derived_failure(exc, "bext")
        clean_bext(paths)
        if reason is not None:
            return {"status": "failed", "stage": "bext", "reason": reason, "name": name, "paths": paths}
        raise
    return {"status": "completed", "name": name, "paths": paths}


def run_pipeline_job(params, force=False, save_all_data=False, timeout_seconds=DEFAULT_JOB_TIMEOUT_SECONDS, dphi_probe=DEFAULT_DPHI_PROBE):
    ensure_stage_dirs()
    name = basename(params)
    if not force and not growth_done(growth_paths(name)):
        clean_shear(shear_paths(name))
        clean_bext(bext_paths(name, dphi_probe))
    growth_result = run_growth_stage(
        params,
        force=force,
        save_all_data=save_all_data,
        timeout_seconds=timeout_seconds,
    )
    if growth_result["status"] == "failed":
        clean_shear(shear_paths(growth_result["name"]))
        clean_bext(bext_paths(growth_result["name"], dphi_probe))
        return growth_result

    shear_result = run_shear_stage(params, force=force, save_all_data=save_all_data)
    if shear_result["status"] == "failed":
        return shear_result

    bext_result = run_bext_stage(params, dphi_probe=dphi_probe, force=force, save_all_data=save_all_data)
    if bext_result["status"] == "failed":
        return bext_result

    if (
        growth_result["status"] == "skipped"
        and shear_result["status"] == "skipped"
        and bext_result["status"] == "skipped"
    ):
        return {"status": "skipped", "name": growth_result["name"]}
    return {"status": "completed", "name": growth_result["name"]}


def execute_param_jobs(job_iterable, runner, max_workers, failure_formatter):
    jobs = list(job_iterable)
    futures = {}
    completed_jobs = 0
    skipped_jobs = 0
    failed_jobs = []
    executor = ThreadPoolExecutor(max_workers=max_workers)
    try:
        for params in jobs:
            future = executor.submit(runner, params)
            futures[future] = params
        while futures:
            done, _ = wait(futures, return_when=FIRST_COMPLETED)
            for future in done:
                params = futures.pop(future)
                try:
                    result = future.result()
                    if result["status"] == "failed":
                        failed_jobs.append((params, result["stage"], result["reason"]))
                        print(failure_formatter(params, result), flush=True)
                    elif result["status"] == "skipped":
                        skipped_jobs += 1
                    else:
                        completed_jobs += 1
                except Exception as exc:
                    for pending in futures:
                        pending.cancel()
                    raise SystemExit(
                        "Failed job: "
                        f"lx={params['lx']} p0={params['p0']} dphi={params['dphi']} seed={params['seed']}\n{exc}"
                    ) from exc
    except KeyboardInterrupt as exc:
        for pending in futures:
            pending.cancel()
        terminate_active_processes()
        executor.shutdown(wait=False, cancel_futures=True)
        raise SystemExit("Interrupted; terminated active pipeline jobs.") from exc
    finally:
        executor.shutdown(wait=True, cancel_futures=True)
    return completed_jobs, skipped_jobs, failed_jobs
