#!/bin/bash -l
# Julia TJLFEP: N_BASIS=32, SCAN_N=1, single node (ir=2).
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 02:00:00
#SBATCH -C cpu
#SBATCH -J TJLFEP_nb32_cpu
#SBATCH -o debug_nb32_julia_%j.out
#SBATCH -e debug_nb32_julia_%j.err
#SBATCH --cpus-per-task=64

set -euo pipefail

module load julia/1.11.7
export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"
export TJLFEP_DEBUG=0
export TJLFEP_FILE_ONLY=1

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
export CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
export GACODE_DUMP="${GACODE_DUMP:-${CASE_DIR}/input.gacode}"
export TGLFEP_FILE="${TJLFEP_ROOT}/build/debug_nb32/input.TGLFEP"
export FILE_DIR="${TJLFEP_ROOT}/build/debug_nb32/fileInput_${SLURM_JOB_ID:-local}"

cd "${TJLFEP_ROOT}/build"

run_julia() { stdbuf -oL -eL julia "$@"; }

echo "=== TJLFEP nb32 single radius (file-only, SCAN_N=1) ==="
run_julia --project="${TJLFEP_ROOT}" -t "${SLURM_CPUS_PER_TASK:-64}" debug_compare_nb32.jl
echo "=== done ==="
