#!/usr/bin/env python3

import argparse

from pipeline_cli import add_common_job_arguments, completed_growth_params, finalize_common_args
from pipeline_config import SHEAR_SCRIPT
from pipeline_paths import ensure_stage_dirs
from pipeline_run import execute_param_jobs, run_shear_stage


def parse_args():
    parser = argparse.ArgumentParser(description="Run the shear pass over completed growth jobs.")
    add_common_job_arguments(parser)
    return finalize_common_args(parser.parse_args())


def main():
    args = parse_args()
    if not SHEAR_SCRIPT.is_file():
        raise SystemExit("Missing shear wrapper script.")
    ensure_stage_dirs()

    growth_ready = completed_growth_params()
    completed, skipped, failed = execute_param_jobs(
        growth_ready,
        lambda params: run_shear_stage(
            params,
            force=args.force,
            save_all_data=args.save_all_data,
        ),
        args.n_cpus,
        lambda params, result: (
            "Failed shear job: "
            f"lx={params['lx']} p0={params['p0']} dphi={params['dphi']} seed={params['seed']} "
            f"reason={result['reason']}"
        ),
    )
    print(f"Queued shear jobs: {len(growth_ready)}")
    print(f"Completed shear jobs: {completed}")
    print(f"Skipped shear jobs: {skipped}")
    print(f"Failed shear jobs: {len(failed)}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
