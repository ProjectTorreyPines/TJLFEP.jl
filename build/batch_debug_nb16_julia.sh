#!/bin/bash -l
# Julia TJLFEP: N_BASIS=16, SCAN_N=1, single node (ir=2).
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 02:00:00
#SBATCH -C cpu
#SBATCH -J TJLFEP_nb16
#SBATCH -o debug_nb16_julia_%j.out
#SBATCH -e debug_nb16_julia_%j.err
#SBATCH --cpus-per-task=64

set -euo pipefail

module load julia/1.11.7
export JULIA_DEPOT_PATH="${PSCRATCH}/.julia"
export TJLFEP_DEBUG=0
export TJLFEP_FILE_ONLY=1

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
export CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
export TGLFEP_FILE="${TJLFEP_ROOT}/build/debug_nb16/input.TGLFEP"
export FILE_DIR="${TJLFEP_ROOT}/build/debug_nb16/fileInput_${SLURM_JOB_ID:-local}"

cd "${TJLFEP_ROOT}/build"

run_julia() { stdbuf -oL -eL julia "$@"; }

echo "=== TJLFEP nb16 single radius (file-only, SCAN_N=1) ==="
run_julia --project="${TJLFEP_ROOT}" -t "${SLURM_CPUS_PER_TASK:-64}" debug_compare_nb16.jl
echo "=== done ==="
