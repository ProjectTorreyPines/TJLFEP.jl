#!/bin/bash -l
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -C cpu
#SBATCH -J TJLFEP_sysimage
#SBATCH -o build_sysimage_%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G

set -euo pipefail

module load julia/1.11.7

export JULIA_DEPOT_PATH="${PSCRATCH}/.julia"
mkdir -p "${JULIA_DEPOT_PATH}/compiled"

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"

echo "=== TJLFEP sysimage build ==="
echo "host: $(hostname)"
echo "date: $(date)"
echo "julia: $(which julia)"
julia --version
echo "JULIA_DEPOT_PATH=${JULIA_DEPOT_PATH}"
echo "TJLF branch:"
git -C "${TJLFEP_ROOT}/../TJLF" rev-parse --abbrev-ref HEAD 2>/dev/null || true
git -C "${TJLFEP_ROOT}/../TJLF" log -1 --oneline 2>/dev/null || true

# Ensure env matches gpu_new TJLF before compiling
julia --project="${TJLFEP_ROOT}" -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

julia --project="${TJLFEP_ROOT}" -t "${SLURM_CPUS_PER_TASK:-32}" build_sysimage.jl

SO="${TJLFEP_ROOT}/build/noTJLF_TJLFEP_sysimage.so"
if [[ ! -f "${SO}" ]]; then
    echo "ERROR: sysimage not found at ${SO}"
    exit 1
fi
ls -lh "${SO}"
echo "=== sysimage build OK: ${SO} ==="
