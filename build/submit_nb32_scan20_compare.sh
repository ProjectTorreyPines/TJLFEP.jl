#!/bin/bash -l
# N_BASIS=32, SCAN_N=20: Fortran 10n | Julia CPU 10n dist | Julia GPU 5n
set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"
mkdir -p timing_runs

FRT=$(sbatch --parsable batch_debug_nb32_fortran_scan20_10n.sh)
CPU=$(sbatch --parsable batch_debug_nb32_julia_scan20_10n.sh)
GPU=$(sbatch --parsable batch_run_gacode_nb32_scan20_gpu_5nodes.sh)

echo "Submitted N_BASIS=32 SCAN_N=20:"
echo "  Fortran 10n/20t:  ${FRT}  -> debug_nb32_fortran20_10n_${FRT}.out"
echo "  Julia CPU 10n:    ${CPU}  -> debug_nb32_julia20_10n_${CPU}.out"
echo "  Julia GPU 5n/20t: ${GPU}  -> gacode_nb32_scan20_gpu5_${GPU}.out"
echo "${FRT} ${CPU} ${GPU}" > timing_runs/last_nb32_scan20_jobs.txt
