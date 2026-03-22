from pathlib import Path

from pipeline_config import BEXT_DIR, BEXT_LOG_DIR, GROWTH_DIR, GROWTH_LOG_DIR, SHEAR_DIR, SHEAR_LOG_DIR, probe_tag


def growth_paths(name):
    return {
        "frame": GROWTH_DIR / f"LF_DPHI_{name}",
        "frame_gz": GROWTH_DIR / f"LF_DPHI_{name}.gz",
        "jamm": GROWTH_DIR / f"LF_JAMM_{name}",
        "jamm_gz": GROWTH_DIR / f"LF_JAMM_{name}.gz",
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


def shear_paths(name):
    return {
        "input_local": SHEAR_DIR / f"LF_DPHI_{name}",
        "input_local_gz": SHEAR_DIR / f"LF_DPHI_{name}.gz",
        "g_data": SHEAR_DIR / f"G_data_LF_DPHI_{name}",
        "traj": SHEAR_DIR / f"SHEAR_TRAJ_LF_DPHI_{name}",
        "traj_gz": SHEAR_DIR / f"SHEAR_TRAJ_LF_DPHI_{name}.gz",
        "log": SHEAR_LOG_DIR / f"stdout_{name[:-4]}.log",
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
