#!/bin/bash -l
set -euo pipefail

TJLFEP_ROOT="${TJLFEP_ROOT:-/pscratch/sd/t/tneiser/.julia/dev/TJLFEP}"
cd "${TJLFEP_ROOT}/build"

if [[ $# -ge 3 ]]; then
  FRT_JOB=$1
  CPU_JOB=$2
  GPU_JOB=$3
elif [[ -f timing_runs/last_time_scan20_jobs.txt ]]; then
  read -r FRT_JOB CPU_JOB GPU_JOB < timing_runs/last_time_scan20_jobs.txt
  CPU_LOG="time_scan20_julia_cpu_${CPU_JOB}.out"
else
  echo "usage: $0 [fortran_job cpu_job gpu_job]"
  exit 1
fi
CPU_LOG="${CPU_LOG:-time_scan20_julia_cpu_${CPU_JOB}.out}"

show() {
  local label=$1
  local f=$2
  echo "--- $label ($f) ---"
  if [[ ! -f $f ]]; then
    echo "  MISSING"
    return
  fi
  grep '^TIMING_RESULT' "$f" | sed 's/^/  /' || echo "  (no TIMING_RESULT lines)"
}

echo "=== SCAN_N=20 N_BASIS=6 timing summary ==="
echo "CPU: Fortran MPI 10 nodes (20 ranks, 2/node)"
echo "Julia CPU: SlurmClusterManager 10 nodes (20 workers, 2/node, 64 threads/worker)"
echo "Julia GPU: gacode srun 5 nodes (20 tasks, 4/node, 1 GPU/task)"
echo ""
show "Fortran" "time_scan20_fortran_${FRT_JOB}.out"
echo ""
show "Julia CPU (distributed)" "${CPU_LOG}"
echo ""
show "Julia GPU" "time_scan20_julia_gpu_${GPU_JOB}.out"
echo ""
echo "Slurm elapsed:"
sacct -j "${FRT_JOB},${CPU_JOB},${GPU_JOB}" --format=JobID,JobName,Elapsed,State -n 2>/dev/null || true

# Quick comparison on phase=scan or phase=compute
CPU_LOG="${CPU_LOG}" FRT_JOB="${FRT_JOB}" CPU_JOB="${CPU_JOB}" GPU_JOB="${GPU_JOB}" python3 <<'PY'
import os, re
from pathlib import Path

def parse(path, phase):
    if not Path(path).exists():
        return None
    for line in Path(path).read_text().splitlines():
        if not line.startswith("TIMING_RESULT"):
            continue
        if f"phase={phase}" not in line:
            continue
        m = re.search(r"seconds=([0-9.]+)", line.replace(" ", ""))
        if m:
            return float(m.group(1))
    return None

frt, cpu, gpu = os.environ["FRT_JOB"], os.environ["CPU_JOB"], os.environ["GPU_JOB"]
cpu_log = os.environ.get("CPU_LOG", f"time_scan20_julia_cpu_{cpu}.out")
jobs = {
    "Fortran scan": (f"time_scan20_fortran_{frt}.out", "scan"),
    "Julia CPU compute": (cpu_log, "compute"),
    "Julia CPU worker_setup": (cpu_log, "worker_setup"),
    "Julia CPU precompile": (cpu_log, "precompile"),
    "Julia CPU total": (cpu_log, "total_job"),
    "Julia GPU scan": (f"time_scan20_julia_gpu_{gpu}.out", "scan"),
    "Julia GPU merge": (f"time_scan20_julia_gpu_{gpu}.out", "merge"),
    "Julia GPU total": (f"time_scan20_julia_gpu_{gpu}.out", "total_job"),
}
print("=== Comparable seconds ===")
vals = {}
for name, (path, phase) in jobs.items():
    s = parse(path, phase)
    vals[name] = s
    print(f"  {name}: {s if s is not None else 'n/a'} s")

frt_s = vals.get("Fortran scan")
jcpu = vals.get("Julia CPU compute")
jgscan = vals.get("Julia GPU scan")
if frt_s and jcpu:
    print(f"\n  Julia CPU compute / Fortran scan = {jcpu/frt_s:.2f}x")
if frt_s and jgscan:
    print(f"  Julia GPU scan / Fortran scan = {jgscan/frt_s:.2f}x")
if jcpu and jgscan:
    print(f"  Julia GPU scan / Julia CPU compute = {jgscan/jcpu:.2f}x")
PY
