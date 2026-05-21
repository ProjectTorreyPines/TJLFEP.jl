#!/bin/bash -l
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 02:00:00
#SBATCH -C cpu
#SBATCH -J TJLFEP_smoke_cpu
#SBATCH -o smoke_nosysimage_cpu_%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G

set -euo pipefail

module load julia/1.11.7

export JULIA_DEPOT_PATH="${PSCRATCH}/.julia"
mkdir -p "${JULIA_DEPOT_PATH}/compiled"

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"

echo "=== TJLFEP smoke (no sysimage, CPU) ==="
echo "host: $(hostname)"
echo "date: $(date)"
julia --version
echo "JULIA_DEPOT_PATH=${JULIA_DEPOT_PATH}"

echo "--- Pkg.instantiate + TJLF path check ---"
julia --project="${TJLFEP_ROOT}" -e '
using Pkg
Pkg.instantiate()
using TJLF
tjlf = pathof(TJLF)
if !occursin("/dev/TJLF", tjlf)
    error("Expected dev TJLF (gpu_new); got: $tjlf")
end
if !isdefined(TJLF, :pick_device)
    error("TJLF missing pick_device — use gpu_new at ../TJLF")
end
println("TJLF: ", tjlf)
println("pick_device(:auto) = ", TJLF.pick_device(:auto))
'

echo "--- smoke_test.jl ---"
julia --project="${TJLFEP_ROOT}" \
    -t "${SLURM_CPUS_PER_TASK:-32}" \
    smoke_test.jl

echo "=== smoke finished OK ==="
