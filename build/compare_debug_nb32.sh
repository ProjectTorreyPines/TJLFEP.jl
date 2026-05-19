#!/bin/bash
# Compare nb32 single-radius Fortran vs Julia outputs.
# Usage:
#   FORTRAN_DIR=.../fortran_runs/debug_nb32_<FJOB> \
#   JULIA_DIR=.../debug_out_nb32_<JJOB> \
#     ./compare_debug_nb32.sh

set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
export FORTRAN_DIR="${FORTRAN_DIR:?set FORTRAN_DIR}"
J_DIR="${JULIA_DIR:?set JULIA_DIR}"

julia --project="${TJLFEP_ROOT}" "${TJLFEP_ROOT}/src/DIIIDfiles/compare_fortran_julia.jl" "${J_DIR}"
