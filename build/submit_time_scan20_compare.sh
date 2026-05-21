#!/bin/bash -l
set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"
mkdir -p timing_runs

FRT=$(sbatch --parsable batch_time_scan20_fortran.sh)
CPU=$(sbatch --parsable batch_time_scan20_julia_cpu.sh)
GPU=$(sbatch --parsable batch_time_scan20_julia_gpu.sh)

echo "Submitted SCAN_N=20 timing:"
echo "  Fortran 10n/20t:        ${FRT}  -> time_scan20_fortran_${FRT}.out"
echo "  Julia CPU 10n/2t dist:  ${CPU}  -> time_scan20_julia_cpu_${CPU}.out"
echo "  Julia GPU 5n/4t gacode: ${GPU}  -> time_scan20_julia_gpu_${GPU}.out"
echo "${FRT} ${CPU} ${GPU}" > timing_runs/last_time_scan20_jobs.txt
echo ""
echo "Comparable layout: Fortran 10n/20t | Julia CPU 10n/2t dist | Julia GPU 5n/4t"
echo "If Julia CPU failed, rerun: ./submit_time_scan20_missing.sh"
