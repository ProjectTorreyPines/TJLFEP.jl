#!/bin/bash -l
# N_BASIS=32, SCAN_N=1 from input.gacode + debug_nb32/input.TGLFEP on 1 GPU node.
#SBATCH -A m3739_g
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 02:00:00
#SBATCH -C gpu
#SBATCH -J TJLFEP_nb32_gpu
#SBATCH -o run_gacode_nb32_gpu_%j.out
#SBATCH -e run_gacode_nb32_gpu_%j.err
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32

set -euo pipefail

module load cudatoolkit/12.9
module load julia/1.11.7

export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"
export TJLFEP_FILE_ONLY=1
export USE_GPU=1
export TJLFEP_DEBUG=0

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
export CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
export GACODE_FILE="${GACODE_FILE:-${CASE_DIR}/input.gacode}"
export TGLFEP_FILE="${TGLFEP_FILE:-${TJLFEP_ROOT}/build/debug_nb32/input.TGLFEP}"

cd "${TJLFEP_ROOT}/build"

echo "=== runTHD_from_gacode N_BASIS=32 GPU ==="
echo "host=$(hostname) date=$(date)"
echo "TGLFEP_FILE=${TGLFEP_FILE}"
nvidia-smi -L 2>/dev/null || true

stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    -t "${SLURM_CPUS_PER_TASK:-32}" run_gacode_nb32.jl

echo "=== done ==="
