import gzip
import shutil


def remove_if_exists(path):
    if path.exists():
        path.unlink()


def gzip_output(path):
    if not path.is_file():
        return
    gz_path = path.with_name(path.name + ".gz")
    with path.open("rb") as src, gzip.open(gz_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    path.unlink()


def clean_growth(paths):
    for key in (
        "frame",
        "frame_gz",
        "jamm",
        "jamm_gz",
        "lineage_frame",
        "lineage_frame_gz",
        "lineage_jamm",
        "lineage_jamm_gz",
        "stats_frame",
        "stats_jamm",
        "divlog",
        "transitions",
        "postjamm_summary",
        "nc",
        "steplog",
        "traj",
        "traj_gz",
        "log",
    ):
        remove_if_exists(paths[key])


def clean_shear(paths):
    for key in ("input_local", "input_local_gz", "g_data", "traj", "traj_gz", "log"):
        remove_if_exists(paths[key])


def clean_bext(paths):
    for key in ("input_local", "input_local_gz", "data", "frame", "frame_gz", "log"):
        remove_if_exists(paths[key])


def normalize_growth_success(paths, save_all_data):
    for key in ("frame", "jamm", "lineage_frame", "lineage_jamm"):
        gzip_output(paths[key])
    if save_all_data:
        gzip_output(paths["traj"])
    else:
        remove_if_exists(paths["traj"])
        remove_if_exists(paths["traj_gz"])


def normalize_shear_success(paths, save_all_data):
    remove_if_exists(paths["input_local"])
    remove_if_exists(paths["input_local_gz"])
    if save_all_data:
        gzip_output(paths["traj"])
    else:
        remove_if_exists(paths["traj"])
        remove_if_exists(paths["traj_gz"])


def normalize_bext_success(paths, save_all_data):
    remove_if_exists(paths["input_local"])
    remove_if_exists(paths["input_local_gz"])
    if save_all_data:
        gzip_output(paths["frame"])
    else:
        remove_if_exists(paths["frame"])
        remove_if_exists(paths["frame_gz"])
