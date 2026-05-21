#!/bin/bash -l
# Submit N_BASIS=32 Fortran + Julia debug after nb6 jobs finish.
#
# Usage:
#   ./submit_nb32_after_nb6.sh                    # uses timing_runs/last_nb6_debug_jobs.txt
#   ./submit_nb32_after_nb6.sh JOB1 JOB2          # explicit nb6 job IDs
#   NB6_JOBS="123 456" ./submit_nb32_after_nb6.sh
set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"
mkdir -p timing_runs

DEP_ARGS=()
if [[ $# -gt 0 ]]; then
    DEP="afterok:$(IFS=:; echo "$*")"
    DEP_ARGS=(--dependency="${DEP}")
elif [[ -n "${NB6_JOBS:-}" ]]; then
    DEP="afterok:$(echo "${NB6_JOBS}" | tr ' ' ':')"
    DEP_ARGS=(--dependency="${DEP}")
elif [[ -f timing_runs/last_nb6_debug_jobs.txt ]]; then
    read -r FRT_JOB JUL_JOB < timing_runs/last_nb6_debug_jobs.txt
    DEP_ARGS=(--dependency=afterok:"${FRT_JOB}:${JUL_JOB}")
    echo "Chaining on nb6 jobs: ${FRT_JOB} ${JUL_JOB}"
else
    echo "WARNING: no nb6 job IDs; submitting nb32 without dependency."
fi

F32=$(sbatch --parsable "${DEP_ARGS[@]}" batch_debug_nb32_fortran.sh)
J32=$(sbatch --parsable "${DEP_ARGS[@]}" batch_debug_nb32_julia.sh)

echo "Submitted nb32 (N_BASIS=32, SCAN_N=1):"
echo "  Fortran: ${F32}  -> debug_nb32_fortran_${F32}.out"
echo "  Julia:   ${J32}  -> debug_nb32_julia_${J32}.out"
echo "${F32} ${J32}" > timing_runs/last_nb32_debug_jobs.txt
