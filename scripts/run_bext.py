#!/usr/bin/env python3

import argparse
import gzip
import os
import re
import shutil
import subprocess
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = REPO_ROOT / "output"
GROWTH_DIR = OUTPUT_ROOT / "growth"
BEXT_DIR = OUTPUT_ROOT / "bext"
LOG_DIR = OUTPUT_ROOT / "logs" / "bext"
EXE = REPO_ROOT / "bin" / "box_compress_bext"
HEADER = "# phi0 P0 phi1 P1 dphi_probe B_ext Lx0 Ly0 Lx1 Ly1 N Nc Nf Nu Ziso Nmm Nbb Nmb fret0 fret1"
NAME_RE = re.compile(
    r"^(?P<prefix>.+)_Lx(?P<lx>[^_]+)_Ly(?P<ly>[^_]+)_att(?P<att>[^_]+)_dphi(?P<dphi>[^_]+)_P(?P<p0>.+)\.dat$"
)


def parse_args():
    default_cpus = max(1, (os.cpu_count() or 1) - 1)
    parser = argparse.ArgumentParser(description="Run the B_ext probe over existing growth packings.")
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
        "--dphi-probe",
        type=float,
        default=1e-6,
        help="Fixed probe compression used for all B_ext jobs. Default: 1e-6",
    )
    return parser.parse_args()


def probe_tag(value):
    return f"{value:.4E}"


def iter_growth_packings():
    seen = set()
    for path in sorted(GROWTH_DIR.glob("LF_DPHI_*.dat")):
        seen.add(path.name)
        yield path
    for path in sorted(GROWTH_DIR.glob("LF_DPHI_*.dat.gz")):
        plain_name = path.name[:-3]
        if plain_name not in seen:
            yield path


def parse_growth_metadata(path):
    name = path.name
    if name.endswith(".gz"):
        name = name[:-3]
    if not name.startswith("LF_DPHI_"):
        return None
    base = name[len("LF_DPHI_") :]
    match = NAME_RE.match(base)
    if match is None:
        return None
    return {
        "input_name": name,
        "base_name": base,
        "lx": match.group("lx"),
        "ly": match.group("ly"),
        "att": match.group("att"),
    }


def build_paths(meta, dphi_probe):
    tag = probe_tag(dphi_probe)
    stem = meta["base_name"][:-4]
    input_local = BEXT_DIR / meta["input_name"]
    return {
        "input_local": input_local,
        "input_local_gz": input_local.with_name(input_local.name + ".gz"),
        "data": BEXT_DIR / f"B_ext_data_dphiprobe{tag}_{meta['input_name']}",
        "frame": BEXT_DIR / f"LF_BEXT_dphiprobe{tag}_{meta['input_name']}",
        "log": LOG_DIR / f"stdout_{stem}_dphiprobe{tag}.log",
    }


def parse_data_line(line):
    parts = line.split()
    if len(parts) != 20:
        raise ValueError(f"expected 20 columns, found {len(parts)}")
    for index in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 18, 19):
        float(parts[index])
    for index in (10, 11, 12, 13, 14, 15, 16, 17):
        int(parts[index])


def bext_done(path):
    if not path.is_file() or path.stat().st_size == 0:
        return False
    with path.open() as handle:
        header = handle.readline().strip()
        if header != HEADER:
            return False
        for line in handle:
            text = line.strip()
            if not text:
                continue
            try:
                parse_data_line(text)
            except ValueError:
                return False
            return True
    return False


def remove_if_exists(path):
    if path.exists():
        path.unlink()


def clean_job(paths):
    for path in paths.values():
        remove_if_exists(path)


def stage_input(source, destination):
    if source.suffix == ".gz":
        with gzip.open(source, "rb") as src, destination.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copyfile(source, destination)


def run_job(source, meta, dphi_probe, force):
    paths = build_paths(meta, dphi_probe)
    if force:
        clean_job(paths)
    elif bext_done(paths["data"]):
        return paths["data"].name

    clean_job(paths)
    stage_input(source, paths["input_local"])

    with paths["log"].open("wb") as log_handle:
        subprocess.run(
            [str(EXE)],
            cwd=BEXT_DIR,
            check=True,
            input="\n".join(
                [
                    meta["lx"],
                    meta["ly"],
                    meta["att"],
                    meta["input_name"],
                    f"{dphi_probe:.16g}",
                    "",
                ]
            ).encode(),
            stdout=log_handle,
            stderr=subprocess.STDOUT,
        )

    if not bext_done(paths["data"]):
        raise RuntimeError(f"Invalid B_ext output: {paths['data']}")
    return paths["data"].name


def main():
    args = parse_args()
    if args.n_cpus < 1:
        raise SystemExit("--n-cpus must be at least 1")
    if args.dphi_probe <= 0.0:
        raise SystemExit("--dphi-probe must be positive")
    if not EXE.is_file():
        raise SystemExit(f"Missing executable: {EXE}. Run 'make bin/box_compress_bext' first.")
    if not GROWTH_DIR.is_dir():
        raise SystemExit(f"Missing growth output directory: {GROWTH_DIR}")

    BEXT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    jobs = []
    for source in iter_growth_packings():
        meta = parse_growth_metadata(source)
        if meta is not None:
            jobs.append((source, meta))
    if not jobs:
        raise SystemExit(f"No valid LF_DPHI growth packings found under {GROWTH_DIR}")

    futures = {}
    with ThreadPoolExecutor(max_workers=args.n_cpus) as executor:
        for source, meta in jobs:
            future = executor.submit(run_job, source, meta, args.dphi_probe, args.force)
            futures[future] = meta["input_name"]
        while futures:
            done, _ = wait(futures, return_when=FIRST_COMPLETED)
            for future in done:
                name = futures.pop(future)
                try:
                    future.result()
                except Exception as exc:
                    for pending in futures:
                        pending.cancel()
                    raise SystemExit(f"Failed B_ext job: {name}\n{exc}") from exc


if __name__ == "__main__":
    main()
