#!/bin/bash -l
# Parse TIMING_RESULT from nb6 timing job logs.
set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"

if [[ $# -ge 3 ]]; then
  FRT_JOB=$1
  CPU_JOB=$2
  GPU_JOB=$3
elif [[ -f timing_runs/last_time_nb6_jobs.txt ]]; then
  read -r FRT_JOB CPU_JOB GPU_JOB < timing_runs/last_time_nb6_jobs.txt
else
  echo "usage: $0 [fortran_job cpu_job gpu_job]"
  exit 1
fi

parse() {
  local f=$1
  if [[ ! -f $f ]]; then
    echo "MISSING $f"
    return
  fi
  grep -m1 '^TIMING_RESULT' "$f" || echo "NO TIMING_RESULT in $f"
}

echo "=== N_BASIS=6 SCAN_N=1 wall time (runTHD_from_gacode / TGLFEP_driver) ==="
echo ""
parse "time_nb6_fortran_${FRT_JOB}.out"
parse "time_nb6_julia_cpu_${CPU_JOB}.out"
parse "time_nb6_julia_gpu_${GPU_JOB}.out"
echo ""
echo "Slurm elapsed (s):"
sacct -j "${FRT_JOB},${CPU_JOB},${GPU_JOB}" --format=JobID,JobName,Elapsed,State -n 2>/dev/null || true
