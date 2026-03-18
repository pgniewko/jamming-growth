#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
VENV_DIR="${SCRIPT_DIR}/.venv"
PIP_CACHE_DIR="${SCRIPT_DIR}/.pip-cache"
PYTHON_BIN=""
RUN_BUILD=1
RUN_VERIFY=1
UPGRADE_PIP=1
INSTALL_SYSTEM_DEPS=0
SYSTEM_PACKAGE_MANAGER=""
MISSING_SYSTEM_DEPS=()

usage() {
    cat <<'EOF'
Usage:
  bash ./install.sh [options]

Options:
  --python PATH        Python interpreter to use for venv creation.
  --install-system-deps
                       Install missing system tools with brew or apt-get.
  --no-build           Skip `make all`.
  --no-verify          Skip verification checks.
  --no-upgrade-pip     Do not upgrade pip inside the venv.
  -h, --help           Show this help.

Notes:
  - The installer checks for Python 3, the `venv` module, `make`, and
    `gfortran` before it creates the project environment.
  - `bash ./install.sh` creates the venv, installs dependencies, builds the
    binaries, and verifies Python imports from `requirements.txt`.
  - The installer never activates the environment for you. Activate it
    separately with `source .venv/bin/activate` when needed.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --python) PYTHON_BIN="$2"; shift 2 ;;
        --install-system-deps) INSTALL_SYSTEM_DEPS=1; shift ;;
        --no-build) RUN_BUILD=0; shift ;;
        --no-verify) RUN_VERIFY=0; shift ;;
        --no-upgrade-pip) UPGRADE_PIP=0; shift ;;
        -h|--help) usage; return 0 2>/dev/null || exit 0 ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            return 1 2>/dev/null || exit 1
            ;;
    esac
done

detect_package_manager() {
    if command -v brew >/dev/null 2>&1; then
        SYSTEM_PACKAGE_MANAGER="brew"
    elif command -v apt-get >/dev/null 2>&1; then
        SYSTEM_PACKAGE_MANAGER="apt-get"
    else
        SYSTEM_PACKAGE_MANAGER=""
    fi
}

add_missing_dep() {
    local dep="$1"
    for existing in "${MISSING_SYSTEM_DEPS[@]:-}"; do
        if [[ "${existing}" == "${dep}" ]]; then
            return 0
        fi
    done
    MISSING_SYSTEM_DEPS+=("${dep}")
}

check_python_venv() {
    local candidate="$1"
    [[ -n "${candidate}" ]] || return 1
    "${candidate}" -c "import sys; assert sys.version_info >= (3, 10); import venv" >/dev/null 2>&1
}

install_system_deps() {
    detect_package_manager
    if [[ ${#MISSING_SYSTEM_DEPS[@]} -eq 0 ]]; then
        return 0
    fi

    if [[ -z "${SYSTEM_PACKAGE_MANAGER}" ]]; then
        echo "Missing system dependencies: ${MISSING_SYSTEM_DEPS[*]}" >&2
        echo "No supported package manager was found. Install them manually and rerun." >&2
        return 1
    fi

    if [[ ${INSTALL_SYSTEM_DEPS} -eq 0 ]]; then
        echo "Missing system dependencies: ${MISSING_SYSTEM_DEPS[*]}" >&2
        if [[ "${SYSTEM_PACKAGE_MANAGER}" == "brew" ]]; then
            echo "Install them with: brew install python gcc make" >&2
        else
            echo "Install them with: sudo apt-get update && sudo apt-get install -y python3 python3-venv make gfortran" >&2
        fi
        echo "Or rerun this script with --install-system-deps." >&2
        return 1
    fi

    if [[ "${SYSTEM_PACKAGE_MANAGER}" == "brew" ]]; then
        brew install python gcc make
    else
        sudo apt-get update
        sudo apt-get install -y python3 python3-venv make gfortran
    fi
}

verify_python_imports() {
    local req_file="${SCRIPT_DIR}/requirements.txt"

    "${VENV_PYTHON}" - "${req_file}" <<'PY'
import importlib
import re
import sys
from pathlib import Path

req_file = Path(sys.argv[1])
if not req_file.exists():
    print(f"Missing requirements file: {req_file}", file=sys.stderr)
    sys.exit(1)

def normalize_import_target(raw: str) -> str:
    name = re.split(r'[<>=!~\[]', raw, 1)[0].strip()
    name = name.replace("-", "_")
    return name

targets = []
for line in req_file.read_text().splitlines():
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        continue
    if stripped.startswith(("-e ", "--")):
        print(f"Unsupported requirement entry for import verification: {stripped}", file=sys.stderr)
        sys.exit(1)
    targets.append((stripped, normalize_import_target(stripped)))

if not targets:
    print("No Python dependencies listed in requirements.txt; skipping import verification.")
    sys.exit(0)

failures = []
for requirement, module_name in targets:
    try:
        importlib.import_module(module_name)
    except Exception as exc:
        failures.append((requirement, module_name, exc))

if failures:
    print("Python dependency verification failed:", file=sys.stderr)
    for requirement, module_name, exc in failures:
        print(f"  - {requirement} -> import {module_name!r} failed: {exc}", file=sys.stderr)
    sys.exit(1)

print("Verified Python imports:")
for requirement, module_name in targets:
    print(f"  - {requirement} -> {module_name}")
PY
}

if [[ -z "${PYTHON_BIN}" ]]; then
    if command -v python3 >/dev/null 2>&1 && check_python_venv "$(command -v python3)"; then
        PYTHON_BIN=$(command -v python3)
    elif command -v python >/dev/null 2>&1 && check_python_venv "$(command -v python)"; then
        PYTHON_BIN=$(command -v python)
    fi
fi

if [[ -z "${PYTHON_BIN}" ]]; then
    add_missing_dep "python3+venv"
fi

if ! command -v make >/dev/null 2>&1; then
    add_missing_dep "make"
fi

if ! command -v gfortran >/dev/null 2>&1; then
    add_missing_dep "gfortran"
fi

install_system_deps || return 1 2>/dev/null || exit 1

if [[ -z "${PYTHON_BIN}" ]]; then
    if command -v python3 >/dev/null 2>&1 && check_python_venv "$(command -v python3)"; then
        PYTHON_BIN=$(command -v python3)
    elif command -v python >/dev/null 2>&1 && check_python_venv "$(command -v python)"; then
        PYTHON_BIN=$(command -v python)
    else
        echo "No usable Python 3 interpreter with the venv module was found after system dependency installation." >&2
        return 1 2>/dev/null || exit 1
    fi
fi

echo "Using Python interpreter: ${PYTHON_BIN}"
mkdir -p "${PIP_CACHE_DIR}"
export PIP_CACHE_DIR
"${PYTHON_BIN}" -m venv "${VENV_DIR}"

VENV_PYTHON="${VENV_DIR}/bin/python"
VENV_PIP="${VENV_DIR}/bin/pip"

if [[ ${UPGRADE_PIP} -eq 1 ]]; then
    "${VENV_PYTHON}" -m pip install --upgrade pip
fi

"${VENV_PIP}" install -r "${SCRIPT_DIR}/requirements.txt"

if [[ ${RUN_BUILD} -eq 1 ]]; then
    make -C "${SCRIPT_DIR}" all
fi

if [[ ${RUN_VERIFY} -eq 1 ]]; then
    verify_python_imports
fi

echo
echo "Setup complete."
echo "To activate the environment later, run:"
echo "  source \"${VENV_DIR}/bin/activate\""
