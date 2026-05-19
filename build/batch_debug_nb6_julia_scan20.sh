#!/bin/bash -l
# Julia file-based debug: N_BASIS=6, SCAN_N=20, single node threads.
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 02:00:00
#SBATCH -C cpu
#SBATCH -J TJLFEP_nb6_s20
#SBATCH -o debug_nb6_julia20_%j.out
#SBATCH -e debug_nb6_julia20_%j.err
#SBATCH --cpus-per-task=32

set -euo pipefail

module load julia/1.11.7
export JULIA_DEPOT_PATH="${PSCRATCH}/.julia"
export TJLFEP_DEBUG=0

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
export CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
export TGLFEP_FILE="${TJLFEP_ROOT}/build/debug_nb6/input_scan20.TGLFEP"
export FILE_DIR="${TJLFEP_ROOT}/build/debug_nb6/fileInput_scan20_${SLURM_JOB_ID:-local}"

cd "${TJLFEP_ROOT}/build"
run_julia() { stdbuf -oL -eL julia "$@"; }

run_julia --project="${TJLFEP_ROOT}" -t "${SLURM_CPUS_PER_TASK:-32}" debug_compare_nb6_scan20.jl
