#!/bin/bash -l
# N_BASIS=32, SCAN_N=20 on 5 GPU nodes: 20 tasks, 4 per node, 1 GPU per task.
#SBATCH -A m3739_g
#SBATCH -q premium
#SBATCH -N 5
#SBATCH -n 20
#SBATCH -t 06:00:00
#SBATCH -C gpu
#SBATCH -J TJLFEP_nb32_s20gpu
#SBATCH -o gacode_nb32_scan20_gpu5_%j.out
#SBATCH -e gacode_nb32_scan20_gpu5_%j.err
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
export TGLFEP_FILE="${TGLFEP_FILE:-${TJLFEP_ROOT}/build/debug_nb32/input_scan20.TGLFEP}"
export OUT_DIR="${TJLFEP_ROOT}/build/gacode_nb32_scan20_${SLURM_JOB_ID}_tasks"

cd "${TJLFEP_ROOT}/build"

echo "=== nb32 SCAN_N=20 on ${SLURM_NNODES:-5} GPU nodes, ${SLURM_NTASKS:-20} tasks ==="
echo "TGLFEP_FILE=${TGLFEP_FILE} OUT_DIR=${OUT_DIR}"
nvidia-smi -L 2>/dev/null | head -4 || true

srun --export=ALL --label -n "${SLURM_NTASKS:-20}" --cpu-bind=cores \
    stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    -t "${SLURM_CPUS_PER_TASK:-32}" run_gacode_scan20_array_task.jl

echo "=== merging on CPU ==="
export USE_GPU=0
# shellcheck source=julia_sysimage.inc.sh
source "${TJLFEP_ROOT}/build/julia_sysimage.inc.sh"
stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    "${JULIA_SYSIMAGE_ARGS[@]}" \
    -t 8 merge_gacode_scan20_array.jl

echo "=== nb32 scan20 GPU done; outputs in ${OUT_DIR} ==="
