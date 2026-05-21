#!/bin/bash -l
set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"
mkdir -p timing_runs

FRT=$(sbatch --parsable batch_time_nb6_fortran.sh)
CPU=$(sbatch --parsable batch_time_nb6_julia_cpu.sh)
GPU=$(sbatch --parsable batch_time_nb6_julia_gpu.sh)

echo "Submitted nb6 timing (SCAN_N=1, N_BASIS=6, input.gacode path):"
echo "  Fortran:    ${FRT}  -> time_nb6_fortran_${FRT}.out"
echo "  Julia CPU:  ${CPU}  -> time_nb6_julia_cpu_${CPU}.out"
echo "  Julia GPU:  ${GPU}  -> time_nb6_julia_gpu_${GPU}.out"
echo "${FRT} ${CPU} ${GPU}" > timing_runs/last_time_nb6_jobs.txt
