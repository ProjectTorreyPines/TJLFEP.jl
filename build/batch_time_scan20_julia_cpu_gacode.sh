#!/bin/bash -l
# Timing: Julia CPU SCAN_N=20, 10 nodes, 2 tasks/node (20 srun tasks) + merge (gacode path).
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 10
#SBATCH -n 20
#SBATCH -t 02:00:00
#SBATCH -C cpu
#SBATCH -J time_s20_jcpu_g
#SBATCH -o time_scan20_julia_cpu_gacode_%j.out
#SBATCH -e time_scan20_julia_cpu_gacode_%j.err
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64

set -euo pipefail

module load julia/1.11.7

export JULIA_DEPOT_PATH="${PSCRATCH}/.julia${JULIA_DEPOT_PATH:+:${JULIA_DEPOT_PATH}}"
export TJLFEP_FILE_ONLY=1
export USE_GPU=0
export TJLFEP_DEBUG=0
export TJLFEP_PRINTOUT=0

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
export CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
export GACODE_FILE="${GACODE_FILE:-${CASE_DIR}/input.gacode}"
export TGLFEP_FILE="${TGLFEP_FILE:-${TJLFEP_ROOT}/build/debug_nb6/input_scan20.TGLFEP}"
export OUT_DIR="${TJLFEP_ROOT}/build/gacode_scan20_time_cpu_${SLURM_JOB_ID}_tasks"

cd "${TJLFEP_ROOT}/build"

JOB_T0=$(date +%s.%N)
echo "TIMING_START backend=julia device=cpu path=gacode nodes=${SLURM_NNODES:-10} tasks=${SLURM_NTASKS:-20} tasks_per_node=2 SCAN_N=20 N_BASIS=6"
echo "OUT_DIR=${OUT_DIR}"

T0=$(date +%s.%N)
srun --export=ALL --label -n "${SLURM_NTASKS:-20}" --cpu-bind=cores \
    stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    -t "${SLURM_CPUS_PER_TASK:-64}" run_gacode_scan20_array_task.jl
T1=$(date +%s.%N)
SCAN_S=$(python3 -c "print(f'{float(\"${T1}\") - float(\"${T0}\"): .3f}')")
echo "TIMING_RESULT backend=julia device=cpu path=gacode phase=scan seconds=${SCAN_S} SCAN_N=20 N_BASIS=6 nodes=${SLURM_NNODES:-10} tasks=${SLURM_NTASKS:-20}"

T0=$(date +%s.%N)
stdbuf -oL -eL julia --startup-file=no --project="${TJLFEP_ROOT}" \
    -t 8 merge_gacode_scan20_array.jl
T1=$(date +%s.%N)
MERGE_S=$(python3 -c "print(f'{float(\"${T1}\") - float(\"${T0}\"): .3f}')")
echo "TIMING_RESULT backend=julia device=cpu path=gacode phase=merge seconds=${MERGE_S} SCAN_N=20 N_BASIS=6"

JOB_T1=$(date +%s.%N)
TOTAL_S=$(python3 -c "print(f'{float(\"${JOB_T1}\") - float(\"${JOB_T0}\"): .3f}')")
echo "TIMING_RESULT backend=julia device=cpu path=gacode phase=total_job seconds=${TOTAL_S} SCAN_N=20 N_BASIS=6 nodes=${SLURM_NNODES:-10} tasks=${SLURM_NTASKS:-20}"
echo "=== done; outputs in ${OUT_DIR} ==="
