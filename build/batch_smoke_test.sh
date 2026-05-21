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

export JULIA_DEPOT_PATH="${PSCRATCH}/.julia"

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
SYSIMAGE="${TJLFEP_ROOT}/build/noTJLF_TJLFEP_sysimage.so"

if [[ ! -f "${SYSIMAGE}" ]]; then
    echo "ERROR: sysimage missing: ${SYSIMAGE}"
    echo "Run batch_build_sysimage.sh first."
    exit 1
fi

cd "${TJLFEP_ROOT}/build"
echo "=== TJLFEP smoke test ==="
echo "host: $(hostname)"
echo "date: $(date)"
ls -lh "${SYSIMAGE}"

julia --project="${TJLFEP_ROOT}" \
    --sysimage="${SYSIMAGE}" \
    -t "${SLURM_CPUS_PER_TASK:-32}" \
    smoke_test.jl

echo "=== smoke test finished ==="
