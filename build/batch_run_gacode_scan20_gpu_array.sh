#!/bin/bash -l
# SCAN_N=20 via Slurm array: one GPU + CPU threads per radius (gacode + TGLFEP only).
# Prefer batch_run_gacode_scan20_gpu_5nodes.sh (5 nodes, 20 tasks) for fewer nodes.
#
# Strategy: 20 independent tasks (--array=0-19), each with:
#   - 1 GPU (EV solve via TJLFEPCUDAExt; rest on CPU threads)
#   - 32 CPU threads for kwscale_scan / TJLF
#
# After array completes, merge:
#   ARRAY_JOB=<id> sbatch batch_merge_gacode_scan20.sh
# or: ./submit_gacode_scan20_gpu.sh  (submits array + merge with dependency)
#
#SBATCH -A m3739_g
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 02:00:00
#SBATCH -C gpu
#SBATCH -J TJLFEP_s20_gpu
#SBATCH -o gacode_scan20_gpu_%A_%a.out
#SBATCH -e gacode_scan20_gpu_%A_%a.err
#SBATCH --array=0-19
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-task=1

set -euo pipefail

module load cudatoolkit/12.9
module load julia/1.11.7

export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"
export TJLFEP_FILE_ONLY=1
export USE_GPU=1
export TJLFEP_DEBUG=0
export TJLFEP_PRINTOUT=0

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
export CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
export GACODE_FILE="${GACODE_FILE:-${CASE_DIR}/input.gacode}"
export TGLFEP_FILE="${TGLFEP_FILE:-${TJLFEP_ROOT}/build/debug_nb6/input_scan20.TGLFEP}"
export OUT_DIR="${OUT_DIR:-${TJLFEP_ROOT}/build/gacode_scan20_${SLURM_ARRAY_JOB_ID}_tasks}"

cd "${TJLFEP_ROOT}/build"

echo "=== array task ${SLURM_ARRAY_TASK_ID:-?} / job ${SLURM_ARRAY_JOB_ID:-?} ==="
echo "host=$(hostname) OUT_DIR=${OUT_DIR}"
nvidia-smi -L 2>/dev/null | head -1 || true

stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    -t "${SLURM_CPUS_PER_TASK:-32}" run_gacode_scan20_array_task.jl

echo "=== array task done ==="
