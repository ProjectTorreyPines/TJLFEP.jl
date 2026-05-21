#!/bin/bash -l
# Fortran: N_BASIS=16, SCAN_N=20, 10 nodes, 20 MPI ranks.
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 10
#SBATCH -t 04:00:00
#SBATCH -C cpu
#SBATCH -J TGLFEP_nb16_s20
#SBATCH -o debug_nb16_fortran20_10n_%j.out
#SBATCH -e debug_nb16_fortran20_10n_%j.err
#SBATCH --ntasks=20
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8

set -euo pipefail

export GACODE_ROOT=/pscratch/sd/t/tneiser/gacode_cpu/gacode
export GACODE_ADD_ROOT=/pscratch/sd/t/tneiser/gacode_cpu/gacode_add
export GACODE_PLATFORM=PERLMUTTER_CPU
export TGLFEP_DIR="${GACODE_ADD_ROOT}/TGLF-EP"
export TGLFEP_DEBUG=0

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
DEBUG_DIR="${TJLFEP_ROOT}/build/debug_nb16"
DRIVER="${TGLFEP_DIR}/TGLFEP_driver"
RUN_DIR="${TJLFEP_ROOT}/build/fortran_runs/debug_nb16_scan20_10n_${SLURM_JOB_ID:-local}"

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

echo "=== Fortran nb16 SCAN_N=20 on ${SLURM_NNODES:-10} nodes, ntasks=${SLURM_NTASKS:-20} ==="
srun -n "${SLURM_NTASKS:-20}" "${DRIVER}"
echo "=== done -> ${RUN_DIR} ==="
