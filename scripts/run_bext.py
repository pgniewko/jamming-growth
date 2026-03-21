#!/usr/bin/env python3

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
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
        help="Fixed probe compression used for all B_ext jobs. Default: 1e-6 (was 1e-8 in Gniewek, 2019)",
    )
    parser.add_argument(
        "--keep-all-output",
        action="store_true",
        help="Keep staged input packings and compressed LF_BEXT frames after successful jobs.",
    )
    return parser.parse_args()


def probe_tag(value):
    return f"{value:.4E}"


def log(message):
    print(f"[run_bext] {message}", flush=True)


def iter_growth_packings():
    for path in sorted(GROWTH_DIR.glob("LF_DPHI_*.dat.gz")):
        if path.stat().st_size > 0:
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
    frame = BEXT_DIR / f"LF_BEXT_dphiprobe{tag}_{meta['input_name']}"
    return {
        "input_local": input_local,
        "input_local_gz": input_local.with_name(input_local.name + ".gz"),
        "data": BEXT_DIR / f"B_ext_data_dphiprobe{tag}_{meta['input_name']}",
        "frame": frame,
        "frame_gz": frame.with_name(frame.name + ".gz"),
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


def open_text(path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_packing(path):
    with open_text(path) as handle:
        header = handle.readline().split()
        if len(header) != 2:
            raise ValueError("missing or invalid packing header")
        count = int(header[0])
        float(header[1])
        rows = 0
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 5:
                raise ValueError(f"expected 5 columns per particle row, found {len(fields)}")
            for field in fields:
                float(field)
            rows += 1
    if rows != count:
        raise ValueError(f"header count {count} does not match {rows} particle rows")


def validate_packing(path):
    if not path.is_file():
        return f"missing file: {path.name}"
    if path.stat().st_size == 0:
        return f"empty file: {path.name}"
    try:
        parse_packing(path)
    except (OSError, ValueError) as exc:
        return f"invalid packing {path.name}: {exc}"
    return None


def validate_bext_data(path):
    if not path.is_file():
        raise ValueError(f"missing file: {path.name}")
    if path.stat().st_size == 0:
        raise ValueError(f"empty file: {path.name}")
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip()
        if header != HEADER:
            raise ValueError("unexpected B_ext header")
        rows = 0
        for line in handle:
            text = line.strip()
            if not text:
                continue
            parse_data_line(text)
            rows += 1
    if rows != 1:
        raise ValueError(f"expected exactly 1 B_ext data row, found {rows}")


def validate_log(path):
    if not path.is_file():
        raise ValueError(f"missing log file: {path.name}")
    if path.stat().st_size == 0:
        raise ValueError(f"empty log file: {path.name}")


def bext_done(paths):
    try:
        validate_bext_data(paths["data"])
        validate_log(paths["log"])
        if paths["frame"].is_file():
            raise ValueError(f"uncompressed retained B_ext frame: {paths['frame'].name}")
        if paths["frame_gz"].is_file():
            reason = validate_packing(paths["frame_gz"])
            if reason is not None:
                raise ValueError(reason)
        return True
    except (OSError, ValueError):
        return False


def remove_if_exists(path):
    if path.exists():
        path.unlink()


def clean_job(paths):
    for path in paths.values():
        remove_if_exists(path)


def clean_finished_job(paths, keep_all_output):
    remove_if_exists(paths["input_local"])
    remove_if_exists(paths["input_local_gz"])
    if keep_all_output:
        return
    remove_if_exists(paths["frame"])
    remove_if_exists(paths["frame_gz"])


def clean_failed_job(paths):
    remove_if_exists(paths["input_local"])
    remove_if_exists(paths["input_local_gz"])
    remove_if_exists(paths["data"])
    remove_if_exists(paths["frame"])
    remove_if_exists(paths["frame_gz"])
    remove_if_exists(paths["log"])


def stage_input(source, destination):
    if source.suffix == ".gz":
        with gzip.open(source, "rb") as src, destination.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copyfile(source, destination)


def gzip_output(path):
    if not path.is_file():
        return
    gz_path = path.with_name(path.name + ".gz")
    with path.open("rb") as src, gzip.open(gz_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    path.unlink()


def report_invalid_sources(invalid_sources):
    if not invalid_sources:
        return
    print(
        f"Skipping {len(invalid_sources)} invalid growth packing(s):",
        file=sys.stderr,
    )
    preview = invalid_sources[:10]
    for source, reason in preview:
        print(f"  - {source.name}: {reason}", file=sys.stderr)
    remaining = len(invalid_sources) - len(preview)
    if remaining > 0:
        print(f"  ... and {remaining} more", file=sys.stderr)


def run_job(source, meta, dphi_probe, force, keep_all_output):
    paths = build_paths(meta, dphi_probe)
    if force:
        clean_job(paths)
    else:
        if bext_done(paths):
            clean_finished_job(paths, keep_all_output)
            return ("skipped", paths["data"].name)
        clean_job(paths)
    stage_input(source, paths["input_local"])

    try:
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

    except Exception:
        clean_failed_job(paths)
        raise
    if keep_all_output:
        gzip_output(paths["frame"])
    clean_finished_job(paths, keep_all_output)
    if not bext_done(paths):
        clean_failed_job(paths)
        raise RuntimeError(f"Invalid B_ext output: {paths['data']}")
    return ("completed", paths["data"].name)


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
    invalid_sources = []
    skipped_existing = 0
    for source in iter_growth_packings():
        meta = parse_growth_metadata(source)
        if meta is not None:
            reason = validate_packing(source)
            if reason is None:
                if not args.force and bext_done(build_paths(meta, args.dphi_probe)):
                    skipped_existing += 1
                else:
                    jobs.append((source, meta))
            else:
                invalid_sources.append((source, reason))
    report_invalid_sources(invalid_sources)
    log(
        "Discovered "
        f"{len(jobs) + skipped_existing} completed growth packings, "
        f"{skipped_existing} with existing valid B_ext results, "
        f"{len(jobs)} queued to run."
    )
    if invalid_sources:
        log(f"Ignored {len(invalid_sources)} invalid completed-growth candidate(s).")
    if not jobs:
        if invalid_sources:
            raise SystemExit(
                f"No valid LF_DPHI growth packings found under {GROWTH_DIR}; "
                "all discovered candidates were empty or invalid."
            )
        log("Nothing to do.")
        return

    futures = {}
    completed = 0
    total = len(jobs)
    with ThreadPoolExecutor(max_workers=args.n_cpus) as executor:
        for source, meta in jobs:
            future = executor.submit(
                run_job,
                source,
                meta,
                args.dphi_probe,
                args.force,
                args.keep_all_output,
            )
            futures[future] = meta["input_name"]
        while futures:
            done, _ = wait(futures, return_when=FIRST_COMPLETED)
            for future in done:
                name = futures.pop(future)
                try:
                    status, _ = future.result()
                    if status == "completed":
                        completed += 1
                        log(f"Completed {completed}/{total}: {name}")
                except Exception as exc:
                    for pending in futures:
                        pending.cancel()
                    raise SystemExit(f"Failed B_ext job: {name}\n{exc}") from exc
    log(
        f"Finished B_ext pass. Completed {completed} new job(s), "
        f"reused {skipped_existing} existing result(s)."
    )


if __name__ == "__main__":
    main()
