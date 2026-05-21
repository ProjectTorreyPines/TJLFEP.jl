#!/bin/bash -l
# Fortran TGLF-EP reference run (user build under gacode_cpu/gacode_add).
# Matches ~/.bashrc.ext: GACODE_ROOT, GACODE_ADD_ROOT, TGLFEP_DIR.
#SBATCH -A m3739
#SBATCH -q premium
#SBATCH -N 10
#SBATCH -t 01:00:00
#SBATCH -C cpu
#SBATCH -J TGLFEP_fortran
#SBATCH -o fortran_tglfep_%j.out
#SBATCH -e fortran_tglfep_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=128

set -euo pipefail

# Non-interactive batch: set explicitly (bashrc.ext only exports these for -i shells).
export GACODE_ROOT=/pscratch/sd/t/tneiser/gacode_cpu/gacode
export GACODE_ADD_ROOT=/pscratch/sd/t/tneiser/gacode_cpu/gacode_add
export GACODE_PLATFORM=PERLMUTTER_CPU
export TGLFEP_DIR="${GACODE_ADD_ROOT}/TGLF-EP"

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
CASE_DIR="${CASE_DIR:-${TJLFEP_ROOT}/src/DIIIDfiles/202017C42_500ms_v3.1}"
DRIVER="${TGLFEP_DIR}/TGLFEP_driver"
# All writes go here — never into CASE_DIR (preserves archived out.scalefactor_r*, etc.).
RUN_DIR="${TJLFEP_ROOT}/build/fortran_runs/${SLURM_JOB_ID:-local}"

echo "=== Fortran TGLF-EP (gacode_cpu) ==="
echo "host: $(hostname)"
echo "date: $(date)"
echo "GACODE_ROOT=${GACODE_ROOT}"
echo "GACODE_ADD_ROOT=${GACODE_ADD_ROOT}"
echo "TGLFEP_DIR=${TGLFEP_DIR}"
echo "DRIVER=${DRIVER}"
echo "CASE_DIR=${CASE_DIR}"

if [[ ! -x "${DRIVER}" ]]; then
    echo "ERROR: TGLFEP_driver not found or not executable: ${DRIVER}"
    exit 1
fi

echo "RUN_DIR=${RUN_DIR} (outputs only; CASE_DIR is read-only)"

# GACODE runtime (same as interactive: source gacode_setup + platform env).
set +u
# shellcheck source=/dev/null
. "${GACODE_ROOT}/shared/bin/gacode_setup"
# shellcheck source=/dev/null
. "${GACODE_ROOT}/platform/env/env.${GACODE_PLATFORM}"
set -u

mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"

for f in input.TGLFEP input.gacode; do
    if [[ ! -f "${CASE_DIR}/${f}" ]]; then
        echo "ERROR: missing ${CASE_DIR}/${f}"
        exit 1
    fi
    ln -sf "${CASE_DIR}/${f}" "${f}"
done

NTASKS="${SLURM_NTASKS:-1280}"
echo "Running in ${RUN_DIR} on ${SLURM_NNODES:-10} nodes, ${NTASKS} MPI tasks (SCAN_N=20)"
srun -n "${NTASKS}" "${DRIVER}"

echo "=== Fortran run finished ==="
ls -la out.TGLFEP out.scalefactor_r* alpha_*_crit.input 2>/dev/null | head -30 || true
