#!/bin/bash -l
# Timing: Julia TJLFEP GPU, N_BASIS=6, SCAN_N=1 (gacode path).
#SBATCH -A m3739_g
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -C gpu
#SBATCH -J time_nb6_jgpu
#SBATCH -o time_nb6_julia_gpu_%j.out
#SBATCH -e time_nb6_julia_gpu_%j.err
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32

set -euo pipefail

module load cudatoolkit/12.9
module load julia/1.11.7
export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"
export TJLFEP_FILE_ONLY=1
export TJLFEP_DEBUG=0
export USE_GPU=1

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"

echo "host=$(hostname) gpu=$(nvidia-smi -L 2>/dev/null | head -1 || echo none)"
stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    -t "${SLURM_CPUS_PER_TASK:-32}" time_nb6.jl
