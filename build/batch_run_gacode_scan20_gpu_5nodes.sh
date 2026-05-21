#!/bin/bash -l
# SCAN_N=20 on 5 GPU nodes: 20 tasks, 4 per node, 1 A100 per task (gacode + TGLFEP only).
# Merge runs on the batch head node CPU after all srun tasks finish (no separate Slurm job).
#
#   sbatch batch_run_gacode_scan20_gpu_5nodes.sh
# or: ./submit_gacode_scan20_gpu.sh
#
#SBATCH -A m3739_g
#SBATCH -q premium
#SBATCH -N 5
#SBATCH -n 20
#SBATCH -t 02:00:00
#SBATCH -C gpu
#SBATCH -J TJLFEP_s20_5gpu
#SBATCH -o gacode_scan20_gpu5_%j.out
#SBATCH -e gacode_scan20_gpu5_%j.err
#SBATCH --ntasks-per-node=4
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
# Always use this job id (do not inherit OUT_DIR from a prior manual merge).
export OUT_DIR="${TJLFEP_ROOT}/build/gacode_scan20_${SLURM_JOB_ID}_tasks"

cd "${TJLFEP_ROOT}/build"

echo "=== SCAN_N=20 on ${SLURM_NNODES:-?} GPU nodes, ${SLURM_NTASKS:-?} tasks ==="
echo "host=$(hostname) OUT_DIR=${OUT_DIR}"
nvidia-smi -L 2>/dev/null | head -4 || true

srun --export=ALL --label -n "${SLURM_NTASKS:-20}" --cpu-bind=cores \
    stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    -t "${SLURM_CPUS_PER_TASK:-32}" run_gacode_scan20_array_task.jl

echo "=== all ${SLURM_NTASKS:-20} tasks done; merging on CPU ==="
export USE_GPU=0
# shellcheck source=julia_sysimage.inc.sh
source "${TJLFEP_ROOT}/build/julia_sysimage.inc.sh"
stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    "${JULIA_SYSIMAGE_ARGS[@]}" \
    -t 8 merge_gacode_scan20_array.jl

echo "=== scan + merge done; outputs in ${OUT_DIR} ==="
