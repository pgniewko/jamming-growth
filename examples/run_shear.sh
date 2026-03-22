#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
EXE_FILE="${ROOT_DIR}/bin/shear_yeast_linearshear"

att="0.0"
desync="0.4"
ar="1.01"
divtype="4"
P0="1e-3"
Lx="8"
dphi="1e-3"
seed="1234"
version="1.0"
strain_step="1e-6"
shear_steps="5000"
results_root="${ROOT_DIR}/output/debug"
force=0
save_all_data=0

compress_if_present() {
    local path="$1"
    if [[ -f "${path}" ]]; then
        gzip -f "${path}"
    fi
}

usage() {
    cat <<'EOF'
Usage: run_shear.sh [options]

Runs one shear job from a saved growth state.

Options:
  --results-root PATH   Root directory for saved outputs.
  --p0 VALUE            Feedback strength. Default: 1e-3
  --lx VALUE            Square box size. Default: 8
  --dphi VALUE          Overcompression, e.g. 1e-3.
  --seed VALUE          Positive seed used in filenames.
  --att VALUE           Attraction range. Default: 0.0
  --desync VALUE        Growth desynchronization. Default: 0.4
  --ar VALUE            Initial bud aspect ratio. Default: 1.01
  --divtype VALUE       Division type. Default: 4
  --version VALUE       Output version. Default: 1.0
  --strain-step VALUE   Shear strain increment passed to the Fortran code. Default: 1e-6
  --shear-steps VALUE   Number of shear steps. Default: 5000
  --force               Overwrite existing outputs.
  --save-all-data       Retain the shear trajectory.
  -h, --help            Show this help.
EOF
}

cleanup_outputs() {
    rm -f "${local_input}" "${local_input}.gz" "${g_data_file}" "${shear_traj_file}" "${shear_traj_file}.gz" "${stdout_log}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --results-root) results_root="$2"; shift 2 ;;
        --p0) P0="$2"; shift 2 ;;
        --lx) Lx="$2"; shift 2 ;;
        --dphi) dphi="$2"; shift 2 ;;
        --seed) seed="$2"; shift 2 ;;
        --att) att="$2"; shift 2 ;;
        --desync) desync="$2"; shift 2 ;;
        --ar) ar="$2"; shift 2 ;;
        --divtype) divtype="$2"; shift 2 ;;
        --version) version="$2"; shift 2 ;;
        --strain-step|--ddelrx) strain_step="$2"; shift 2 ;;
        --shear-steps) shear_steps="$2"; shift 2 ;;
        --force) force=1; shift ;;
        --save-all-data) save_all_data=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ ! -x "${EXE_FILE}" ]]; then
    echo "Missing executable: ${EXE_FILE}. Run 'make all' first." >&2
    exit 1
fi

if [[ "${results_root}" != /* ]]; then
    results_root="${PWD}/${results_root}"
fi

growth_dir="${results_root}/growth"
shear_dir="${results_root}/shear"
log_dir="${results_root}/logs/shear"
mkdir -p "${shear_dir}" "${log_dir}"

Ly="${Lx}"
basename="v${version}_ar${ar}_div_${divtype}_desync${desync}_seed_${seed}"
basename="${basename}_Lx${Lx}_Ly${Ly}_att${att}_dphi${dphi}_P${P0}.dat"
input_file="${growth_dir}/LF_DPHI_${basename}"
local_input="${shear_dir}/LF_DPHI_${basename}"
stdout_log="${log_dir}/stdout_${basename%.dat}.log"
g_data_file="${shear_dir}/G_data_LF_DPHI_${basename}"
shear_traj_file="${shear_dir}/SHEAR_TRAJ_LF_DPHI_${basename}"

if [[ ! -f "${input_file}" && ! -f "${input_file}.gz" ]]; then
    echo "Missing growth configuration: ${input_file}" >&2
    exit 1
fi

if [[ ${force} -eq 0 && -s "${g_data_file}" && -s "${stdout_log}" ]]; then
    if [[ ${save_all_data} -eq 0 || -s "${shear_traj_file}" || -s "${shear_traj_file}.gz" ]]; then
    echo "Skipping existing shear run: ${g_data_file}"
    exit 0
    fi
fi

cleanup_outputs

if [[ -f "${input_file}" ]]; then
    cp -f "${input_file}" "${local_input}"
else
    gzip -dc "${input_file}.gz" > "${local_input}"
fi

if [[ $(wc -l < "${local_input}") -lt 2 ]]; then
    echo "Input file is too short: ${local_input}" >&2
    exit 1
fi

(
    cd "${shear_dir}"
    "${EXE_FILE}" <<EOF >"${stdout_log}" 2>&1
${Lx}
${Ly}
${att}
LF_DPHI_${basename}
${strain_step}
${shear_steps}
EOF
)

compress_if_present "${local_input}"
if [[ ${save_all_data} -eq 1 ]]; then
    compress_if_present "${shear_traj_file}"
else
    rm -f "${shear_traj_file}" "${shear_traj_file}.gz"
fi
rm -f "${local_input}" "${local_input}.gz"

echo "Saved shear outputs in ${shear_dir}"
echo "Saved shear stdout in ${stdout_log}"
