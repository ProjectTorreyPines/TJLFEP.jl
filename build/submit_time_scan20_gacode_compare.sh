#!/bin/bash -l
# SCAN_N=20 timing: Fortran + Julia CPU gacode (10n) + Julia GPU gacode (5n).
set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"
mkdir -p timing_runs

FRT=$(sbatch --parsable batch_time_scan20_fortran.sh)
CPU=$(sbatch --parsable batch_time_scan20_julia_cpu_gacode.sh)
GPU=$(sbatch --parsable batch_time_scan20_julia_gpu.sh)

echo "Submitted SCAN_N=20 gacode-path timing:"
echo "  Fortran 10n/20t:     ${FRT}  -> time_scan20_fortran_${FRT}.out"
echo "  Julia CPU gacode 10n: ${CPU}  -> time_scan20_julia_cpu_gacode_${CPU}.out"
echo "  Julia GPU gacode 5n:  ${GPU}  -> time_scan20_julia_gpu_${GPU}.out"
echo "${FRT} ${CPU} ${GPU}" > timing_runs/last_time_scan20_gacode_jobs.txt
