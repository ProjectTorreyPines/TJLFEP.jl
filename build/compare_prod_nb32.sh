#!/bin/bash
# Compare N_BASIS=32 / SCAN_N=20 Fortran vs Julia outputs after Slurm jobs finish.
# Usage:
#   FORTRAN_DIR=.../fortran_runs/prod_nb32_<FJOB> \
#   JULIA_DIR=.../validate_out_<JJOB>_files_dist \
#     ./compare_prod_nb32.sh

set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
F_DIR="${FORTRAN_DIR:?set FORTRAN_DIR to fortran_runs/prod_nb32_<jobid>}"
J_DIR="${JULIA_DIR:?set JULIA_DIR to validate_out_<jobid>_files_dist}"

julia --project="${TJLFEP_ROOT}" "${TJLFEP_ROOT}/src/DIIIDfiles/compare_fortran_julia.jl" "${J_DIR}"
