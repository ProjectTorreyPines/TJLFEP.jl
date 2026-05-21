#!/bin/bash -l
# Timing: Julia TJLFEP CPU, N_BASIS=6, SCAN_N=1 (gacode path).
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -C cpu
#SBATCH -J time_nb6_jcpu
#SBATCH -o time_nb6_julia_cpu_%j.out
#SBATCH -e time_nb6_julia_cpu_%j.err
#SBATCH --cpus-per-task=8

set -euo pipefail

module load julia/1.11.7
export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"
export TJLFEP_FILE_ONLY=1
export TJLFEP_DEBUG=0
export USE_GPU=0

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"

stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    -t "${SLURM_CPUS_PER_TASK:-8}" time_nb6.jl
