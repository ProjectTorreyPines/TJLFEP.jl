#!/bin/bash -l
# Submit Fortran + Julia nb6 debug (SCAN_N=1) in parallel; record job IDs for nb32 chaining.
set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"
mkdir -p timing_runs

FRT=$(sbatch --parsable batch_debug_nb6_fortran.sh)
JUL=$(sbatch --parsable batch_debug_nb6_julia.sh)

echo "Submitted nb6 debug (parallel):"
echo "  Fortran: ${FRT}  -> debug_nb6_fortran_${FRT}.out"
echo "  Julia:   ${JUL}  -> debug_nb6_julia_${JUL}.out"
echo "${FRT} ${JUL}" > timing_runs/last_nb6_debug_jobs.txt
