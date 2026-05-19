#!/usr/bin/env python3

import argparse

from pipeline_cli import add_common_job_arguments, completed_growth_params, finalize_common_args
from pipeline_paths import ensure_stage_dirs
from pipeline_run import execute_param_jobs, run_bext_stage


def parse_args():
    parser = argparse.ArgumentParser(description="Run the B_ext probe over completed growth jobs.")
    add_common_job_arguments(parser, include_probe=True)
    return finalize_common_args(parser.parse_args())


def main():
    args = parse_args()
    ensure_stage_dirs()

    growth_ready = completed_growth_params()
    completed, skipped, failed = execute_param_jobs(
        growth_ready,
        lambda params: run_bext_stage(
            params,
            dphi_probe=args.dphi_probe,
            force=args.force,
            save_all_data=args.save_all_data,
        ),
        args.n_cpus,
        lambda params, result: (
            "Failed B_ext job: "
            f"lx={params['lx']} p0={params['p0']} dphi={params['dphi']} seed={params['seed']} "
            f"reason={result['reason']}"
        ),
    )
    print(f"Queued B_ext jobs: {len(growth_ready)}")
    print(f"Completed B_ext jobs: {completed}")
    print(f"Skipped B_ext jobs: {skipped}")
    print(f"Failed B_ext jobs: {len(failed)}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
