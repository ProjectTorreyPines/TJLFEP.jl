#!/bin/bash -l
# N_BASIS=16, SCAN_N=20: Fortran 10n premium | Julia CPU 10n premium | Julia GPU 5n
set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"
mkdir -p timing_runs

FRT=$(sbatch --parsable batch_debug_nb16_fortran_scan20_10n.sh)
CPU=$(sbatch --parsable batch_debug_nb16_julia_scan20_10n.sh)
GPU=$(sbatch --parsable batch_run_gacode_nb16_scan20_gpu_5nodes.sh)

echo "Submitted N_BASIS=16 SCAN_N=20:"
echo "  Fortran 10n/20t (premium):  ${FRT}  -> debug_nb16_fortran20_10n_${FRT}.out"
echo "  Julia CPU 10n (premium):    ${CPU}  -> debug_nb16_julia20_10n_${CPU}.out"
echo "  Julia GPU 5n/20t:           ${GPU}  -> gacode_nb16_scan20_gpu5_${GPU}.out"
echo "${FRT} ${CPU} ${GPU}" > timing_runs/last_nb16_scan20_jobs.txt
