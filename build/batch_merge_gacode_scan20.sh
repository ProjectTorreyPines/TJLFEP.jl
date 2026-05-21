#!/bin/bash -l
# Merge array task outputs -> alpha_dndr_crit.input, alpha_dpdr_crit.input
# (Only needed for batch_run_gacode_scan20_gpu_array.sh; 5-node job merges in-job.)
#SBATCH -A m3739_g
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -C cpu
#SBATCH -J TJLFEP_s20_merge
#SBATCH -o gacode_scan20_merge_%j.out
#SBATCH -e gacode_scan20_merge_%j.err
#SBATCH --cpus-per-task=8

set -euo pipefail

module load julia/1.11.7
export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"
export TJLFEP_FILE_ONLY=1

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
export CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
export GACODE_FILE="${GACODE_FILE:-${CASE_DIR}/input.gacode}"
export TGLFEP_FILE="${TGLFEP_FILE:-${TJLFEP_ROOT}/build/debug_nb6/input_scan20.TGLFEP}"

# Set by submit script or manually: OUT_DIR=.../gacode_scan20_<ARRAY_JOB_ID>_tasks
: "${OUT_DIR:?set OUT_DIR to gacode_scan20_<array_jobid>_tasks}"

cd "${TJLFEP_ROOT}/build"
echo "=== merge OUT_DIR=${OUT_DIR} ==="

stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" merge_gacode_scan20_array.jl

echo "=== merge done ==="
