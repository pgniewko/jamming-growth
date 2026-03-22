#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
EXE_FILE="${ROOT_DIR}/bin/jamming_by_growth"

att="0.0"
rate0="0.002"
skip="25"
desync="0.4"
ar="1.01"
divtype="4"
P0="1e-3"
Lx="8"
dphi="1e-3"
seed="1234"
version="1.0"
results_root="${ROOT_DIR}/output/debug"
growth_dir=""
log_dir=""
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
Usage: run_jamming.sh [options]

Runs one growth/jamming job.

Options:
  --results-root PATH   Root directory for saved outputs.
  --growth-dir PATH     Growth output directory. Overrides --results-root.
  --log-dir PATH        Stdout log directory. Overrides --results-root.
  --p0 VALUE            Feedback strength. Default: 1e-3
  --lx VALUE            Square box size. Default: 8
  --dphi VALUE          Overcompression, e.g. 1e-3.
  --seed VALUE          Positive seed used in filenames.
  --att VALUE           Attraction range. Default: 0.0
  --rate0 VALUE         Growth rate. Default: 0.002
  --skip VALUE          Trajectory save interval. Default: 25
  --desync VALUE        Growth desynchronization. Default: 0.4
  --ar VALUE            Initial bud aspect ratio. Default: 1.01
  --divtype VALUE       Division type. Default: 4
  --version VALUE       Output version. Default: 1.0
  --force               Overwrite existing outputs.
  --save-all-data       Retain the raw growth trajectory.
  -h, --help            Show this help.
EOF
}

cleanup_outputs() {
    rm -f \
        "${growth_dir}/${basename}" \
        "${growth_dir}/${basename}.gz" \
        "${growth_dir}/LF_JAMM_${basename}" \
        "${growth_dir}/LF_JAMM_${basename}.gz" \
        "${growth_dir}/LF_DPHI_${basename}" \
        "${growth_dir}/LF_DPHI_${basename}.gz" \
        "${growth_dir}/LINEAGE_LF_JAMM_${basename}" \
        "${growth_dir}/LINEAGE_LF_JAMM_${basename}.gz" \
        "${growth_dir}/LINEAGE_LF_DPHI_${basename}" \
        "${growth_dir}/LINEAGE_LF_DPHI_${basename}.gz" \
        "${stats_jamm}" \
        "${stats_frame}" \
        "${nc_file}" \
        "${steplog_file}" \
        "${divlog_file}" \
        "${transitions_file}" \
        "${postjamm_file}" \
        "${stdout_log}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --results-root) results_root="$2"; shift 2 ;;
        --growth-dir) growth_dir="$2"; shift 2 ;;
        --log-dir) log_dir="$2"; shift 2 ;;
        --p0) P0="$2"; shift 2 ;;
        --lx) Lx="$2"; shift 2 ;;
        --dphi) dphi="$2"; shift 2 ;;
        --seed) seed="$2"; shift 2 ;;
        --att) att="$2"; shift 2 ;;
        --rate0) rate0="$2"; shift 2 ;;
        --skip) skip="$2"; shift 2 ;;
        --desync) desync="$2"; shift 2 ;;
        --ar) ar="$2"; shift 2 ;;
        --divtype) divtype="$2"; shift 2 ;;
        --version) version="$2"; shift 2 ;;
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

if [[ -z "${growth_dir}" ]]; then
    growth_dir="${results_root}/growth"
elif [[ "${growth_dir}" != /* ]]; then
    growth_dir="${PWD}/${growth_dir}"
fi

if [[ -z "${log_dir}" ]]; then
    log_dir="${results_root}/logs/growth"
elif [[ "${log_dir}" != /* ]]; then
    log_dir="${PWD}/${log_dir}"
fi

mkdir -p "${growth_dir}" "${log_dir}"

Ly="${Lx}"
seed_fortran=$((-1 * seed))
dphi_fortran=$(printf '%s' "${dphi}" | tr '[:upper:]' '[:lower:]' | sed 's/e/d/g')
basename="v${version}_ar${ar}_div_${divtype}_desync${desync}_seed_${seed}"
basename="${basename}_Lx${Lx}_Ly${Ly}_att${att}_dphi${dphi}_P${P0}.dat"
stdout_log="${log_dir}/stdout_${basename%.dat}.log"
final_frame="${growth_dir}/LF_DPHI_${basename}.gz"
jam_frame="${growth_dir}/LF_JAMM_${basename}.gz"
lineage_frame="${growth_dir}/LINEAGE_LF_DPHI_${basename}.gz"
lineage_jamm="${growth_dir}/LINEAGE_LF_JAMM_${basename}.gz"
stats_frame="${growth_dir}/STATS_LF_DPHI_${basename}"
stats_jamm="${growth_dir}/STATS_LF_JAMM_${basename}"
nc_file="${growth_dir}/NC_${basename}"
steplog_file="${growth_dir}/STEPLOG_${basename}"
divlog_file="${growth_dir}/DIVLOG_${basename}"
transitions_file="${growth_dir}/TRANSITIONS_${basename}"
postjamm_file="${growth_dir}/POSTJAMM_SUMMARY_${basename}"

if [[ ${force} -eq 0 && -s "${final_frame}" && -s "${jam_frame}" && -s "${lineage_frame}" && -s "${lineage_jamm}" && -s "${stats_frame}" && -s "${stats_jamm}" && -s "${nc_file}" && -s "${steplog_file}" && -s "${divlog_file}" && -s "${transitions_file}" && -s "${postjamm_file}" && -s "${stdout_log}" ]]; then
    echo "Skipping existing growth run: ${final_frame}"
    exit 0
fi

cleanup_outputs

(
    cd "${growth_dir}"
    "${EXE_FILE}" <<EOF >"${stdout_log}" 2>&1
${ar}
${Lx}
${Ly}
${divtype}
${P0}
${att}
${rate0}
${desync}
${seed_fortran}
${skip}
${dphi_fortran}
${basename}
EOF
)

compress_if_present "${growth_dir}/${basename}"
compress_if_present "${growth_dir}/LF_JAMM_${basename}"
compress_if_present "${growth_dir}/LF_DPHI_${basename}"
compress_if_present "${growth_dir}/LINEAGE_LF_JAMM_${basename}"
compress_if_present "${growth_dir}/LINEAGE_LF_DPHI_${basename}"

if [[ ${save_all_data} -eq 0 ]]; then
    rm -f "${growth_dir}/${basename}" "${growth_dir}/${basename}.gz"
fi

echo "Saved growth outputs in ${growth_dir}"
echo "Saved growth stdout in ${stdout_log}"
