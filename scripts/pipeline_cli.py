import argparse
import os

from pipeline_config import DEFAULT_DPHI_PROBE, DEFAULT_JOB_TIMEOUT_SECONDS, add_save_all_output_arguments, basename, job_params, resolve_save_all_data
from pipeline_paths import growth_paths
from pipeline_validate import growth_done


def default_cpus():
    return max(1, (os.cpu_count() or 1) - 1)


def add_common_job_arguments(parser, include_timeout=False, include_probe=False):
    cpus = default_cpus()
    parser.add_argument(
        "--n-cpus",
        type=int,
        default=cpus,
        help=f"Number of concurrent jobs. Default: {cpus}",
    )
    parser.add_argument("--force", action="store_true", help="Rerun completed jobs.")
    add_save_all_output_arguments(parser)
    if include_timeout:
        parser.add_argument(
            "--job-timeout-seconds",
            type=int,
            default=DEFAULT_JOB_TIMEOUT_SECONDS,
            help=(
                "Kill a growth job if it runs longer than this many seconds. "
                f"Default: {DEFAULT_JOB_TIMEOUT_SECONDS}"
            ),
        )
    if include_probe:
        parser.add_argument(
            "--dphi-probe",
            type=float,
            default=DEFAULT_DPHI_PROBE,
            help=f"Fixed probe compression used for all B_ext jobs. Default: {DEFAULT_DPHI_PROBE:g}",
        )


def finalize_common_args(args):
    if args.n_cpus < 1:
        raise SystemExit("--n-cpus must be at least 1")
    if hasattr(args, "job_timeout_seconds") and args.job_timeout_seconds < 0:
        raise SystemExit("--job-timeout-seconds must be non-negative")
    if hasattr(args, "dphi_probe") and args.dphi_probe <= 0.0:
        raise SystemExit("--dphi-probe must be positive")
    args.save_all_data = resolve_save_all_data(args)
    return args


def completed_growth_params():
    return [params for params in job_params() if growth_done(growth_paths(basename(params)))]
