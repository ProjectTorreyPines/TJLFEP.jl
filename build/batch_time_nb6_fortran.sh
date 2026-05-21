#!/bin/bash -l
# Timing: Fortran TGLF-EP, N_BASIS=6, SCAN_N=1, 1 MPI rank.
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -C cpu
#SBATCH -J time_nb6_frt
#SBATCH -o time_nb6_fortran_%j.out
#SBATCH -e time_nb6_fortran_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

set -euo pipefail

export GACODE_ROOT=/pscratch/sd/t/tneiser/gacode_cpu/gacode
export GACODE_ADD_ROOT=/pscratch/sd/t/tneiser/gacode_cpu/gacode_add
export GACODE_PLATFORM=PERLMUTTER_CPU
export TGLFEP_DIR="${GACODE_ADD_ROOT}/TGLF-EP"
export TGLFEP_DEBUG=0

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
DEBUG_DIR="${TJLFEP_ROOT}/build/debug_nb6"
DRIVER="${TGLFEP_DIR}/TGLFEP_driver"
RUN_DIR="${TJLFEP_ROOT}/build/timing_runs/fortran_nb6_${SLURM_JOB_ID}"

set +u
# shellcheck source=/dev/null
. "${GACODE_ROOT}/shared/bin/gacode_setup"
# shellcheck source=/dev/null
. "${GACODE_ROOT}/platform/env/env.${GACODE_PLATFORM}"
set -u

mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"

ln -sf "${DEBUG_DIR}/input.TGLFEP" input.TGLFEP
ln -sf "${CASE_DIR}/input.gacode" input.gacode

echo "TIMING_START backend=fortran device=cpu SCAN_N=1 N_BASIS=6 ntasks=1"
T0=$(date +%s.%N)
srun -n 1 "${DRIVER}"
T1=$(date +%s.%N)
ELAPSED=$(python3 -c "print(f'{float(\"${T1}\") - float(\"${T0}\"): .3f}')")
echo "TIMING_RESULT backend=fortran device=cpu seconds=${ELAPSED} SCAN_N=1 N_BASIS=6 ntasks=1"
