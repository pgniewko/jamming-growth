from pathlib import Path

from pipeline_config import BEXT_DIR, BEXT_LOG_DIR, GROWTH_DIR, GROWTH_LOG_DIR, SHEAR_DIR, SHEAR_LOG_DIR, probe_tag


def normalize_source_tag(source_tag):
    tag = str(source_tag).upper()
    if tag not in {"DPHI", "JAMM", "PHI2"}:
        raise ValueError(f"unsupported source tag: {source_tag}")
    return tag


def growth_paths(name):
    return {
        "frame": GROWTH_DIR / f"LF_DPHI_{name}",
        "frame_gz": GROWTH_DIR / f"LF_DPHI_{name}.gz",
        "jamm": GROWTH_DIR / f"LF_JAMM_{name}",
        "jamm_gz": GROWTH_DIR / f"LF_JAMM_{name}.gz",
        "phi2_frame": GROWTH_DIR / f"LF_PHI2_{name}",
        "phi2_frame_gz": GROWTH_DIR / f"LF_PHI2_{name}.gz",
        "lineage_frame": GROWTH_DIR / f"LINEAGE_LF_DPHI_{name}",
        "lineage_frame_gz": GROWTH_DIR / f"LINEAGE_LF_DPHI_{name}.gz",
        "lineage_jamm": GROWTH_DIR / f"LINEAGE_LF_JAMM_{name}",
        "lineage_jamm_gz": GROWTH_DIR / f"LINEAGE_LF_JAMM_{name}.gz",
        "stats_frame": GROWTH_DIR / f"STATS_LF_DPHI_{name}",
        "stats_jamm": GROWTH_DIR / f"STATS_LF_JAMM_{name}",
        "divlog": GROWTH_DIR / f"DIVLOG_{name}",
        "transitions": GROWTH_DIR / f"TRANSITIONS_{name}",
        "postjamm_summary": GROWTH_DIR / f"POSTJAMM_SUMMARY_{name}",
        "cohort": GROWTH_DIR / f"COHORT_INITIAL_FREE_{name}",
        "cohort_gz": GROWTH_DIR / f"COHORT_INITIAL_FREE_{name}.gz",
        "nc": GROWTH_DIR / f"NC_{name}",
        "steplog": GROWTH_DIR / f"STEPLOG_{name}",
        "traj": GROWTH_DIR / name,
        "traj_gz": GROWTH_DIR / f"{name}.gz",
        "log": GROWTH_LOG_DIR / f"stdout_{name[:-4]}.log",
    }


def shear_paths(name, source_tag="DPHI"):
    tag = normalize_source_tag(source_tag)
    stem = name[:-4]
    if tag == "DPHI":
        log_name = f"stdout_{stem}.log"
    else:
        log_name = f"stdout_{tag.lower()}_{stem}.log"
    return {
        "input_local": SHEAR_DIR / f"LF_{tag}_{name}",
        "input_local_gz": SHEAR_DIR / f"LF_{tag}_{name}.gz",
        "g_data": SHEAR_DIR / f"G_data_LF_{tag}_{name}",
        "traj": SHEAR_DIR / f"SHEAR_TRAJ_LF_{tag}_{name}",
        "traj_gz": SHEAR_DIR / f"SHEAR_TRAJ_LF_{tag}_{name}.gz",
        "log": SHEAR_LOG_DIR / log_name,
    }


def bext_paths(name, dphi_probe):
    tag = probe_tag(dphi_probe)
    stem = name[:-4]
    input_local = BEXT_DIR / f"LF_DPHI_{name}"
    frame = BEXT_DIR / f"LF_BEXT_dphiprobe{tag}_LF_DPHI_{name}"
    return {
        "input_local": input_local,
        "input_local_gz": input_local.with_name(input_local.name + ".gz"),
        "data": BEXT_DIR / f"B_ext_data_dphiprobe{tag}_LF_DPHI_{name}",
        "frame": frame,
        "frame_gz": frame.with_name(frame.name + ".gz"),
        "log": BEXT_LOG_DIR / f"stdout_{stem}_dphiprobe{tag}.log",
    }


def ensure_stage_dirs():
    for directory in (GROWTH_DIR, SHEAR_DIR, BEXT_DIR, GROWTH_LOG_DIR, SHEAR_LOG_DIR, BEXT_LOG_DIR):
        directory.mkdir(parents=True, exist_ok=True)
