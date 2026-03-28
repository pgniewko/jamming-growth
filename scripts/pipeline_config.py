from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = REPO_ROOT / "output"

GROWTH_SCRIPT = REPO_ROOT / "examples" / "run_jamming.sh"
SHEAR_SCRIPT = REPO_ROOT / "examples" / "run_shear.sh"
BEXT_EXE = REPO_ROOT / "bin" / "box_compress_bext"

GROWTH_DIR = OUTPUT_ROOT / "growth"
SHEAR_DIR = OUTPUT_ROOT / "shear"
BEXT_DIR = OUTPUT_ROOT / "bext"

GROWTH_LOG_DIR = OUTPUT_ROOT / "logs" / "growth"
SHEAR_LOG_DIR = OUTPUT_ROOT / "logs" / "shear"
BEXT_LOG_DIR = OUTPUT_ROOT / "logs" / "bext"

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
BEXT_HEADER = (
    "# phi0 P0 phi1 P1 dphi_probe B_ext Lx0 Ly0 Lx1 Ly1 N Nc Nf Nu Ziso Nmm Nbb Nmb fret0 fret1"
)

DEFAULT_JOB_TIMEOUT_SECONDS = 10800
DEFAULT_DPHI_PROBE = 1e-8
EXIT_MIN_DT = 10
EXIT_MAX_POSTJAM_STEPS = 11


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


def seed_range(start, stop):
    if stop < start:
        raise ValueError("seed stop must be greater than or equal to seed start")
    return [str(seed) for seed in range(start, stop + 1)]


def iter_all_job_params(sizes=None, p0s=None, dphis=None, seeds=None):
    sizes = SIZES if sizes is None else [str(size) for size in sizes]
    p0s = P0S if p0s is None else [str(p0) for p0 in p0s]
    dphis = DPHIS if dphis is None else [str(dphi) for dphi in dphis]
    seeds = SEEDS if seeds is None else [str(seed) for seed in seeds]
    for size in sizes:
        for p0 in p0s:
            for dphi in dphis:
                for seed in seeds:
                    yield {"lx": size, "p0": p0, "dphi": dphi, "seed": seed}


def iter_job_params(sizes=None, p0s=None, dphis=None, seeds=None):
    for params in iter_all_job_params(sizes=sizes, p0s=p0s, dphis=dphis, seeds=seeds):
        if dphi_allowed(params["p0"], params["dphi"]):
            yield params


def all_job_params():
    yield from iter_all_job_params()


def job_params():
    yield from iter_job_params()


def basename(params):
    return (
        f"v{FIXED['version']}_ar{FIXED['ar']}_div_{FIXED['divtype']}_desync{FIXED['desync']}"
        f"_seed_{params['seed']}_Lx{params['lx']}_Ly{params['lx']}"
        f"_att{FIXED['att']}_dphi{params['dphi']}_P{params['p0']}.dat"
    )


def probe_tag(value):
    return f"{value:.4E}"


def add_save_all_output_arguments(parser):
    parser.add_argument(
        "--save-all-data",
        action="store_true",
        help="Retain large optional outputs such as trajectories and retained probe frames.",
    )


def resolve_save_all_data(args):
    return bool(getattr(args, "save_all_data", False))
