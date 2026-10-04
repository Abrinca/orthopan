#!/usr/bin/env bash
set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${DATA_DIR:-${REPO_ROOT}/resources/orthofinder-demo}"
# OrthoFinder needs its own Python 3.12 env, separate from the project's venv.
VENV_DIR="${ORTHOFINDER_VENV:-${REPO_ROOT}/scratch/.venv}"
CORE_DIR="${DATA_DIR}/faa-core"
ADD_DIR="${DATA_DIR}/faa-add"
OUTPUT_DIR="${DATA_DIR}/results"
THREADS="$(nproc 2>/dev/null || echo 4)"

# Check input directories
if [ ! -d "${CORE_DIR}" ]; then
    echo "ERROR: Core FASTA directory '${CORE_DIR}' does not exist." >&2
    exit 1
fi

if [ ! -d "${ADD_DIR}" ]; then
    echo "ERROR: Additional FASTA directory '${ADD_DIR}' does not exist." >&2
    exit 1
fi

# Verify and activate virtual environment
if [ ! -d "${VENV_DIR}" ]; then
    echo "ERROR: Virtual environment not found at '${VENV_DIR}'." >&2
    echo "Please set up the environment first:" >&2
    echo "  uv venv --python python3.12" >&2
    echo "  source .venv/bin/activate" >&2
    echo "  uv pip install git+https://github.com/OrthoFinder/OrthoFinder.git" >&2
    exit 1
fi

# shellcheck disable=SC1091
source "${VENV_DIR}/bin/activate"

if ! command -v orthofinder >/dev/null 2>&1; then
    echo "ERROR: 'orthofinder' command not found in virtual environment." >&2
    exit 1
fi

# Python compatibility fix for OrthoFinder with NumPy 2.x:
# OrthoFinder upstream still calls `np.chararray` in sample_genes.py, which was
# moved to `np.char.chararray` in NumPy 2.0+. We provide sitecustomize.py so the
# alias is defined automatically when python starts.
COMPAT_DIR="${VENV_DIR}/../.compat"
mkdir -p "${COMPAT_DIR}"
cat << 'EOF' > "${COMPAT_DIR}/sitecustomize.py"
import numpy as np
if not hasattr(np, 'chararray') and hasattr(np, 'char') and hasattr(np.char, 'chararray'):
    np.chararray = np.char.chararray
EOF
export PYTHONPATH="${COMPAT_DIR}${PYTHONPATH:+:${PYTHONPATH}}"

VENV_SP="$(python3 -c 'import site; print(site.getsitepackages()[0])' 2>/dev/null || true)"
if [ -n "${VENV_SP}" ] && [ -d "${VENV_SP}" ]; then
    cp "${COMPAT_DIR}/sitecustomize.py" "${VENV_SP}/sitecustomize.py" 2>/dev/null || true
fi

echo "=========================================="
echo "Running OrthoFinder v2 Workflow (Local venv)"
echo "OrthoFinder:      $(command -v orthofinder)"
echo "Version:          $(orthofinder -v 2>&1 | tr '\n' ' ')"
echo "Core FASTA:       ${CORE_DIR}"
echo "Additional FASTA: ${ADD_DIR}"
echo "Output dir:       ${OUTPUT_DIR}"
echo "Threads:          ${THREADS}"
echo "=========================================="

# Check if completed core analysis already exists
EXISTING_CORE=""
if [ -d "${OUTPUT_DIR}" ]; then
    for candidate in "${OUTPUT_DIR}"/Results_*; do
        if [ -d "${candidate}" ] && grep -q "Started OrthoFinder" "${candidate}/Log.txt" 2>/dev/null && grep -q "OrthoFinder run completed" "${candidate}/Log.txt" 2>/dev/null && ! grep -q -- "--assign" "${candidate}/Log.txt" 2>/dev/null; then
            EXISTING_CORE="${candidate}"
            break
        fi
    done
fi

if [ -n "${EXISTING_CORE}" ]; then
    echo ""
    echo "=== Step 1: Found completed core analysis ==="
    echo "Reusing core results: ${EXISTING_CORE}"
    echo "(To re-run Step 1 from scratch, remove '${OUTPUT_DIR}' before running)"
    RESULTS_DIR="${EXISTING_CORE}"
else
    # OrthoFinder requires that the non-default output directory does not already exist
    if [ -d "${OUTPUT_DIR}" ]; then
        echo "ERROR: Output directory '${OUTPUT_DIR}' already exists but contains no completed core results." >&2
        echo "OrthoFinder requires that the non-default output directory does not exist yet." >&2
        echo "Please remove or rename the existing directory before running." >&2
        exit 1
    fi

    echo ""
    echo "=== Step 1: Run core OrthoFinder analysis ==="
    orthofinder -f "${CORE_DIR}" -o "${OUTPUT_DIR}" -t "${THREADS}" "$@"

    # Identify the generated Results directory inside OUTPUT_DIR
    RESULTS_DIR="$(find "${OUTPUT_DIR}" -maxdepth 1 -type d -name "Results_*" | grep -v "_[0-9]\+$" | sort | tail -n 1)"
    if [ -z "${RESULTS_DIR}" ]; then
        echo "ERROR: Could not find Results_* directory inside '${OUTPUT_DIR}'." >&2
        exit 1
    fi
fi

echo ""
echo "=== Step 2: Assign additional species to existing core results ==="
echo "Using core results: ${RESULTS_DIR}"
orthofinder --assign "${ADD_DIR}" --core "${RESULTS_DIR}" -t "${THREADS}" "$@"
