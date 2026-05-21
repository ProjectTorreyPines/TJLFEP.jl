#!/bin/bash -l
# Timing: Julia GPU SCAN_N=20, N_BASIS=8, 5 nodes, 20 tasks + CPU merge.
#SBATCH -A m3739_g
#SBATCH -q premium
#SBATCH -N 5
#SBATCH -n 20
#SBATCH -t 03:00:00
#SBATCH -C gpu
#SBATCH -J time_nb8_jgpu
#SBATCH -o time_scan20_nb8_julia_gpu_%j.out
#SBATCH -e time_scan20_nb8_julia_gpu_%j.err
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
export TGLFEP_FILE="${TGLFEP_FILE:-${TJLFEP_ROOT}/build/debug_nb8/input_scan20.TGLFEP}"
export OUT_DIR="${TJLFEP_ROOT}/build/gacode_nb8_scan20_time_${SLURM_JOB_ID}_tasks"

cd "${TJLFEP_ROOT}/build"

JOB_T0=$(date +%s.%N)
echo "TIMING_START backend=julia device=gpu path=gacode nodes=${SLURM_NNODES:-5} tasks=${SLURM_NTASKS:-20} SCAN_N=20 N_BASIS=8"
echo "OUT_DIR=${OUT_DIR}"
nvidia-smi -L 2>/dev/null | head -4 || true

T0=$(date +%s.%N)
srun --export=ALL --label -n "${SLURM_NTASKS:-20}" --cpu-bind=cores \
    stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    -t "${SLURM_CPUS_PER_TASK:-32}" run_gacode_scan20_array_task.jl
T1=$(date +%s.%N)
SCAN_S=$(python3 -c "print(f'{float(\"${T1}\") - float(\"${T0}\"): .3f}')")
echo "TIMING_RESULT backend=julia device=gpu phase=scan seconds=${SCAN_S} SCAN_N=20 N_BASIS=8 nodes=${SLURM_NNODES:-5} tasks=${SLURM_NTASKS:-20}"

export USE_GPU=0
# shellcheck source=julia_sysimage.inc.sh
source "${TJLFEP_ROOT}/build/julia_sysimage.inc.sh"
T0=$(date +%s.%N)
stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    "${JULIA_SYSIMAGE_ARGS[@]}" \
    -t 8 merge_gacode_scan20_array.jl
T1=$(date +%s.%N)
MERGE_S=$(python3 -c "print(f'{float(\"${T1}\") - float(\"${T0}\"): .3f}')")
echo "TIMING_RESULT backend=julia device=gpu phase=merge seconds=${MERGE_S} SCAN_N=20 N_BASIS=8"

JOB_T1=$(date +%s.%N)
TOTAL_S=$(python3 -c "print(f'{float(\"${JOB_T1}\") - float(\"${JOB_T0}\"): .3f}')")
echo "TIMING_RESULT backend=julia device=gpu phase=total_job seconds=${TOTAL_S} SCAN_N=20 N_BASIS=8 nodes=${SLURM_NNODES:-5} tasks=${SLURM_NTASKS:-20}"
echo "=== done; outputs in ${OUT_DIR} ==="
