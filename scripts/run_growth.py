#!/usr/bin/env python3

import argparse

from pipeline_cli import add_common_job_arguments, finalize_common_args
from pipeline_config import GROWTH_SCRIPT, job_params
from pipeline_paths import ensure_stage_dirs
from pipeline_run import execute_param_jobs, run_growth_stage


def parse_args():
    parser = argparse.ArgumentParser(description="Run the growth sweep.")
    add_common_job_arguments(parser, include_timeout=True)
    return finalize_common_args(parser.parse_args())


def main():
    args = parse_args()
    if not GROWTH_SCRIPT.is_file():
        raise SystemExit("Missing growth wrapper script.")
    ensure_stage_dirs()

    completed, skipped, failed = execute_param_jobs(
        job_params(),
        lambda params: run_growth_stage(
            params,
            force=args.force,
            save_all_data=args.save_all_data,
            timeout_seconds=args.job_timeout_seconds,
        ),
        args.n_cpus,
        lambda params, result: (
            "Failed growth job: "
            f"lx={params['lx']} p0={params['p0']} dphi={params['dphi']} seed={params['seed']} "
            f"reason={result['reason']}"
        ),
    )
    print(f"Completed growth jobs: {completed}")
    print(f"Skipped growth jobs: {skipped}")
    print(f"Failed growth jobs: {len(failed)}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
