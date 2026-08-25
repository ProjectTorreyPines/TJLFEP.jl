#!/bin/bash -l
#SBATCH -A m3739
#SBATCH -q debug
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -C gpu
#SBATCH -J TJLFEP_smoke
#SBATCH -o smoke_test_%j.out
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32

set -euo pipefail

module load cudatoolkit/12.9
module load julia/1.11.7

# Shared m3739 tree: staged sources (src/), full depot, and published sysimages.
CFS_DIR="${TJLFEP_CFS_DIR:-/global/cfs/cdirs/m3739/TJLFEP}"

# Your writable depot MUST come first; the CFS depot provides read-only packages/artifacts.
export JULIA_DEPOT_PATH="${SCRATCH}/.julia:${CFS_DIR}/depot"
export JULIA_CUDA_USE_COMPAT=false

TJLFEP_ROOT="${TJLFEP_ROOT:-${CFS_DIR}/src/TJLFEP}"
# GACODE/IMAS/TurbulentTransport are TJLFEP weak deps; only the FUSE project resolves them
# (incl. the patched IMASdd/IMASggd), so JIT mode must run under that project.
SMOKE_PROJECT="${SMOKE_PROJECT:-${CFS_DIR}/src/FUSE}"
# Optional GPU sysimage (build once with batch_build_gpu_sysimage_generic.sh). Missing -> JIT.
SYSIMAGE="${TJLFEP_GPU_SYSIMAGE:-${CFS_DIR}/TJLFEP_gpu_generic_sysimage.so}"

cd "${TJLFEP_ROOT}/build"
echo "=== TJLFEP smoke test ==="
echo "host: $(hostname)"
echo "date: $(date)"

if [[ -f "${SYSIMAGE}" ]]; then
    ls -lh "${SYSIMAGE}"
    SYSIMG_ARGS=(--sysimage="${SYSIMAGE}")
else
    echo "sysimage not found at ${SYSIMAGE} -> running with JIT"
    SYSIMG_ARGS=()
fi

julia --project="${SMOKE_PROJECT}" \
    "${SYSIMG_ARGS[@]}" \
    -t "${SLURM_CPUS_PER_TASK:-32}" \
    run/smoke_test.jl

echo "=== smoke test finished ==="
