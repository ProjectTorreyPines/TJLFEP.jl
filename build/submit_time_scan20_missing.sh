#!/bin/bash -l
# Submit only missing SCAN_N=20 timing: Julia CPU 10 nodes (SlurmClusterManager).
# Fortran (10n) and Julia GPU (5n) are in timing_runs/last_time_scan20_jobs.txt when complete.
set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"
mkdir -p timing_runs

CPU=$(sbatch --parsable batch_time_scan20_julia_cpu.sh)
echo "Submitted Julia CPU distributed timing (10 nodes, 20 workers): ${CPU}"
echo "  log: time_scan20_julia_cpu_${CPU}.out"
echo "${CPU}" >> timing_runs/last_time_scan20_cpu_resubmit.txt
