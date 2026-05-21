#!/bin/bash -l
# nb6 SCAN_N=1 from input.gacode + input.TGLFEP on 1 GPU node.
#SBATCH -A m3739_g
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -C gpu
#SBATCH -J TJLFEP_gacode_nb6_gpu
#SBATCH -o run_gacode_nb6_gpu_%j.out
#SBATCH -e run_gacode_nb6_gpu_%j.err
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32

set -euo pipefail

module load cudatoolkit/12.9
module load julia/1.11.7

export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"
export TJLFEP_FILE_ONLY=1
export USE_GPU=1

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"

echo "=== runTHD_from_gacode nb6 GPU ==="
echo "host=$(hostname) date=$(date)"
nvidia-smi -L 2>/dev/null || true

stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    -t "${SLURM_CPUS_PER_TASK:-32}" run_gacode_nb6.jl

echo "=== done ==="
