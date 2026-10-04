#!/usr/bin/env bash
set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${DATA_DIR:-${REPO_ROOT}/resources/orthofinder-demo}"
IMAGE="docker.io/arnstrm2/orthofinder:3.1.5"
OUTPUT_DIR="${DATA_DIR}/results"
THREADS="$(nproc 2>/dev/null || echo 4)"

# OrthoFinder requires that the non-default output directory does not already exist
if [ -d "${OUTPUT_DIR}" ]; then
    echo "ERROR: Output directory '${OUTPUT_DIR}' already exists." >&2
    echo "OrthoFinder requires that the non-default output directory does not exist yet." >&2
    echo "Please remove or rename the existing directory before running." >&2
    exit 1
fi

echo "=========================================="
echo "Running OrthoFinder with Podman"
echo "Image:       ${IMAGE}"
echo "Input FASTA: ${DATA_DIR}/faa-core"
echo "Output dir:  ${OUTPUT_DIR}"
echo "Threads:     ${THREADS}"
echo "=========================================="

# Note: In rootless Podman, omitting --user runs as root inside the container's
# user namespace, which maps 1:1 to your user on the host. This gives write access
# to /data while ensuring generated files are owned by your host user.
podman run --rm \
    -v "${DATA_DIR}:/data:Z" \
    -w /data \
    "${IMAGE}" \
    orthofinder -f /data/faa-core -o /data/results -t "${THREADS}" "$@"
