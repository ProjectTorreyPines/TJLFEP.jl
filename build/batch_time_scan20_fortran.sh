#!/bin/bash -l
# Timing: Fortran SCAN_N=20, 10 nodes, 20 MPI ranks.
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 10
#SBATCH -t 02:00:00
#SBATCH -C cpu
#SBATCH -J time_s20_frt
#SBATCH -o time_scan20_fortran_%j.out
#SBATCH -e time_scan20_fortran_%j.err
#SBATCH --ntasks=20
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
RUN_DIR="${TJLFEP_ROOT}/build/timing_runs/fortran_scan20_${SLURM_JOB_ID}"

set +u
# shellcheck source=/dev/null
. "${GACODE_ROOT}/shared/bin/gacode_setup"
# shellcheck source=/dev/null
. "${GACODE_ROOT}/platform/env/env.${GACODE_PLATFORM}"
set -u

mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"
ln -sf "${DEBUG_DIR}/input_scan20.TGLFEP" input.TGLFEP
ln -sf "${CASE_DIR}/input.gacode" input.gacode

JOB_T0=$(date +%s.%N)
echo "TIMING_START backend=fortran device=cpu path=gacode nodes=${SLURM_NNODES:-10} tasks=${SLURM_NTASKS:-20} SCAN_N=20 N_BASIS=6"

T0=$(date +%s.%N)
set +e
srun -n "${SLURM_NTASKS:-20}" "${DRIVER}"
SRUN_RC=$?
set -e
T1=$(date +%s.%N)
SCAN_S=$(python3 -c "print(f'{float(\"${T1}\") - float(\"${T0}\"): .3f}')")
JOB_T1=$(date +%s.%N)
TOTAL_S=$(python3 -c "print(f'{float(\"${JOB_T1}\") - float(\"${JOB_T0}\"): .3f}')")

echo "TIMING_RESULT backend=fortran device=cpu phase=scan seconds=${SCAN_S} SCAN_N=20 N_BASIS=6 nodes=${SLURM_NNODES:-10} tasks=${SLURM_NTASKS:-20} srun_exit=${SRUN_RC}"
echo "TIMING_RESULT backend=fortran device=cpu phase=total_job seconds=${TOTAL_S} SCAN_N=20 N_BASIS=6 nodes=${SLURM_NNODES:-10} tasks=${SLURM_NTASKS:-20}"
[ "${SRUN_RC}" -ne 0 ] && echo "NOTE: srun exited ${SRUN_RC} (often MPI teardown); check outputs in ${RUN_DIR}"
echo "srun_exit_code=${SRUN_RC} (MPI teardown may segfault after successful scan)"
