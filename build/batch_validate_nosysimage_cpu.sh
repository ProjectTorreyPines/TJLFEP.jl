#!/bin/bash -l
# Full DIII-D validation (20 radii) without sysimage on CPU — allow long walltime.
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH -C cpu
#SBATCH -J TJLFEP_val_cpu
#SBATCH -o validate_nosysimage_cpu_%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G

set -euo pipefail

module load julia/1.11.7

export JULIA_DEPOT_PATH="${PSCRATCH}/.julia"
TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
export CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
# Generated Julia inputs only — not under CASE_DIR (ref outputs stay untouched).
export FILE_DIR="${FILE_DIR:-${TJLFEP_ROOT}/build/fileInput_${SLURM_JOB_ID:-local}}"
cd "${TJLFEP_ROOT}/build"

run_julia() {
  stdbuf -oL -eL julia "$@"
}

echo "=== TJLFEP validate (file-based, CPU) CASE_DIR=$CASE_DIR ==="
echo "host: $(hostname)"
julia --version

run_julia --project="${TJLFEP_ROOT}" -e '
using Pkg; Pkg.instantiate(); using TJLF
@assert occursin("/dev/TJLF", pathof(TJLF))
@assert isdefined(TJLF, :pick_device)
println("TJLF: ", pathof(TJLF))
flush(stdout)
'

run_julia --project="${TJLFEP_ROOT}" -t "${SLURM_CPUS_PER_TASK:-32}" validate_nosysimage.jl

echo "=== validate finished ==="
