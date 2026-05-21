#!/bin/bash -l
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -C gpu
#SBATCH -J TJLFEP_smoke_ns
#SBATCH -o smoke_nosysimage_%j.out
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32

set -euo pipefail

module load cudatoolkit/12.9
module load julia/1.11.7

export JULIA_DEPOT_PATH="${PSCRATCH}/.julia"
mkdir -p "${JULIA_DEPOT_PATH}/compiled"

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"

echo "=== TJLFEP smoke (no sysimage) ==="
echo "host: $(hostname)"
echo "date: $(date)"
julia --version
echo "JULIA_DEPOT_PATH=${JULIA_DEPOT_PATH}"

echo "--- Pkg.instantiate ---"
julia --project="${TJLFEP_ROOT}" -e 'using Pkg; Pkg.instantiate()'

echo "--- smoke_test.jl ---"
julia --project="${TJLFEP_ROOT}" \
    -t "${SLURM_CPUS_PER_TASK:-32}" \
    smoke_test.jl

echo "=== smoke finished OK ==="
