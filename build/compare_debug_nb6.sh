#!/bin/bash
# Compare nb6 single-radius Fortran vs Julia α profiles.
# Usage:
#   FORTRAN_DIR=.../fortran_runs/debug_nb6_<FJOB> \
#   JULIA_DIR=.../debug_out_nb6_<JJOB> \
#   PLOT_OUTDIR=.../compare_nb6_plots \
#     ./compare_debug_nb6.sh

set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
export FORTRAN_DIR="${FORTRAN_DIR:?set FORTRAN_DIR}"
export JULIA_DIR="${JULIA_DIR:?set JULIA_DIR}"
export PLOT_OUTDIR="${PLOT_OUTDIR:-${TJLFEP_ROOT}/build/compare_nb6_plots}"
export FORTRAN_LABEL="${FORTRAN_LABEL:-Fortran}"
export JULIA_LABEL="${JULIA_LABEL:-Julia}"
export PLOT_TITLE="${PLOT_TITLE:-N_BASIS=6, SCAN_N=1, ir=2}"

julia --project="${TJLFEP_ROOT}" "${TJLFEP_ROOT}/build/plot_nb32_scan_compare.jl"
